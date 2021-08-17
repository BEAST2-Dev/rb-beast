package beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.DefaultMixtureSiteModel;
import beast.evolution.sitemodel.MixtureSiteModel;
import beast.evolution.sitemodel.MixtureSubstModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;



@Description("Tree likelihood that can handle a general mixture model efficiently")
public class MixtureTreeLikelihood extends TreeLikelihood {
	
	
	DefaultMixtureSiteModel siteModel;
	MixtureLikelihoodCore mixtureLikelihoodCore;
	// list of categories that need (re)calculation
	List<Integer> classes;
	public static List<Integer> ALL;
	
	
	@Override
	public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getNrTaxa() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = (SubstitutionModel.Base) m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int nStateCount = dataInput.get().getMaxStateCount();
        int nPatterns = dataInput.get().getPatternCount();
        if (nStateCount == 4) {
            likelihoodCore = new MixtureBeerLikelihoodCore4();
        } else {
            likelihoodCore = new MixtureBeerLikelihoodCore(nStateCount);
        }
        System.out.println("TreeLikelihood uses " + likelihoodCore.getClass().getName());

//        m_fProportionInvariant = m_siteModel.getProportianInvariant();
//        m_siteModel.setPropInvariantIsCategory(false);
//        if (m_fProportionInvariant > 0) {
//            calcConstantPatternIndices(nPatterns, nStateCount);
//        }

        initCore();

        patternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        matrixSize = (nStateCount + 1) * (nStateCount + 1);
        probabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get() instanceof AscertainedAlignment) {
            useAscertainedSitePatterns = true;
        }

        siteModel = (DefaultMixtureSiteModel) siteModelInput.get();
        mixtureLikelihoodCore = (MixtureLikelihoodCore) likelihoodCore;
        ALL = new ArrayList<Integer>();
        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
        	ALL.add(i);
        }
        classes = ALL;
	}

    /* Assumes there IS a branch rate model as opposed to traverse() */
	@Override
	protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int iNode = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;
        m_branchLengths[iNode] = branchTime;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != storedBranchLengths[iNode])) {
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(iNode);
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                siteModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities, i);
                likelihoodCore.setNodeMatrix(iNode, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(iNode);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(iNode);
                }

                if (siteModel.integrateAcrossCategories()) {
                	mixtureLikelihoodCore.calculatePartials(childNum1, childNum2, iNode, classes);
                } else {
                    throw new IllegalArgumentException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    //final double[] frequencies = //m_pFreqs.get().
                    final double[][] frequencies = new double[siteModel.getCategoryCount()][];
                    for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                    	frequencies[i] = siteModel.getFrequencies(i);
                    }

                    final double[] proportions = siteModel.getCategoryProportions(node);
                    ((MixtureLikelihoodCore)likelihoodCore).integratePartialsMixture(node.getNr(), proportions, m_fRootPartials, frequencies, patternLogLikelihoods);

//                    if (m_iConstantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        m_fProportionInvariant = siteModel.getProportianInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : m_iConstantPattern) {
//                            m_fRootPartials[i] += m_fProportionInvariant;
//                        }
//                    }

//                    IMPLEMENTATIONS NOT CORRECT:
//                    ((FrequenciesMixture)m_likelihoodCore).calculateLogLikelihoods(m_fRootPartials, frequencies, m_fPatternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverseWithBRM

	
	@Override
	protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
    	classes = new ArrayList<Integer>();

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
        	classes = ALL;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
        	classes = ALL;
            return true;
        }
        if (treeInput.get().somethingIsDirty()) {
        	classes = ALL;
        	return true;
        }
        if (siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
    		for (int k = 0; k < siteModel.getCategoryCount(); k++) {
    			if (siteModel.hasDirtySubstModel(k)) {
    				classes.add(k);
    			}
    		}
            return true;
        }
        return false;
	}
	
}
