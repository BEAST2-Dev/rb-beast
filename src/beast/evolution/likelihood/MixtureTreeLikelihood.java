package beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.StrictClockModel;
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
	public void initAndValidate() throws Exception {
        // sanity check: alignment should have same #taxa as tree
        if (m_data.get().getNrTaxa() != m_tree.get().getLeafNodeCount()) {
            throw new Exception("The number of nodes in the tree does not match the number of sequences");
        }
        // No Beagle instance was found, so we use the good old java likelihood core
        m_beagle = null;

        int nodeCount = m_tree.get().getNodeCount();
        if (!(m_pSiteModel.get() instanceof SiteModel.Base)) {
        	throw new Exception ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) m_pSiteModel.get();
        m_siteModel.setDataType(m_data.get().getDataType());
        m_substitutionModel = (SubstitutionModel.Base) m_siteModel.m_pSubstModel.get();

        if (m_pBranchRateModel.get() != null) {
            m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        m_StoredBranchLengths = new double[nodeCount];

        int nStateCount = m_data.get().getMaxStateCount();
        int nPatterns = m_data.get().getPatternCount();
        if (nStateCount == 4) {
            m_likelihoodCore = new MixtureBeerLikelihoodCore4();
        } else {
            m_likelihoodCore = new MixtureBeerLikelihoodCore(nStateCount);
        }
        System.out.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

//        m_fProportionInvariant = m_siteModel.getProportianInvariant();
//        m_siteModel.setPropInvariantIsCategory(false);
//        if (m_fProportionInvariant > 0) {
//            calcConstantPatternIndices(nPatterns, nStateCount);
//        }

        initCore();

        m_fPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        m_nMatrixSize = (nStateCount + 1) * (nStateCount + 1);
        m_fProbabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        Arrays.fill(m_fProbabilities, 1.0);

        if (m_data.get() instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }

        siteModel = (DefaultMixtureSiteModel) m_pSiteModel.get();
        mixtureLikelihoodCore = (MixtureLikelihoodCore) m_likelihoodCore;
        ALL = new ArrayList<Integer>();
        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
        	ALL.add(i);
        }
        classes = ALL;
	}

    /* Assumes there IS a branch rate model as opposed to traverse() */
	@Override
    int traverse(final Node node) throws Exception {

        int update = (node.isDirty() | m_nHasDirt);

        final int iNode = node.getNr();

        final double branchRate = m_branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;
        m_branchLengths[iNode] = branchTime;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[iNode])) {
            final Node parent = node.getParent();
            m_likelihoodCore.setNodeMatrixForUpdate(iNode);
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                siteModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, m_fProbabilities, i);
                m_likelihoodCore.setNodeMatrix(iNode, i, m_fProbabilities);
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

                m_likelihoodCore.setNodePartialsForUpdate(iNode);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    m_likelihoodCore.setNodeStatesForUpdate(iNode);
                }

                if (siteModel.integrateAcrossCategories()) {
                	mixtureLikelihoodCore.calculatePartials(childNum1, childNum2, iNode, classes);
                } else {
                    throw new Exception("Error TreeLikelihood 201: Site categories not supported");
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
                    ((MixtureLikelihoodCore)m_likelihoodCore).integratePartialsMixture(node.getNr(), proportions, m_fRootPartials, frequencies, m_fPatternLogLikelihoods);

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
        m_nHasDirt = Tree.IS_CLEAN;
    	classes = new ArrayList<Integer>();

        if (m_data.get().isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_FILTHY;
        	classes = ALL;
            return true;
        }
        if (m_branchRateModel != null && m_branchRateModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
        	classes = ALL;
            return true;
        }
        if (m_tree.get().somethingIsDirty()) {
        	classes = ALL;
        	return true;
        }
        if (siteModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
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
