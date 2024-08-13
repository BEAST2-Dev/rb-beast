package rbbeast.evolution.likelihood;


import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import rbbeast.evolution.sitemodel.SingleCategorySiteModel;
import rbbeast.evolution.sitemodel.VariableCategorySiteModel;

@Description("Tree likelihood that allows the number of rate categories to change. "
		+ "It uses a TreeLikelihood for each category, then integrates over the resulting root probabilites.")
public class VariableCategoryTreeLikelihood2 extends TreeLikelihood {
	
    /** private list of likelihoods, to notify framework of TreeLikelihoods being created **/
    final private Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("*","",new ArrayList<>());

	// private List<TreeLikelihood> categoryLikelihoods;
	private int nrOfStates;
	// private int nrOfPatterns, nrOfStates, categoryCount;

	private List<Integer> classes;
	
	@Override
	public void initAndValidate() {
		m_siteModel = (SiteModel.Base) siteModelInput.get();
		// revert to TreeLikelihood default behaviour when site model is not a VariableCategorySiteModel
		if (!(m_siteModel instanceof VariableCategorySiteModel)) {
			super.initAndValidate();
			return;
		}		
		
		// Site model must be a VariableCategorySiteModel once we get here
		// initialise, but make sure no beagle instance is created
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        System.setProperty("java.only", "true");
		super.initAndValidate();
        System.setProperty("java.only", forceJava + "");

		// clean up memory
		likelihoodCore = null;
		
//		nrOfPatterns = dataInput.get().getPatternCount();
		nrOfStates = substitutionModel.getStateCount();
//		categoryCount = m_siteModel.getCategoryCount();
		
		// create TreeLikelihoods, one for each potential category
//		categoryLikelihoods = new ArrayList<>();
//		for (int i = 0; i < ((VariableCategorySiteModel) m_siteModel).maxCategoryCountInput.get(); i++) {
//			createNewTreeLikelihood();
//		}

	
		classes = new ArrayList<>();
        for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
        	classes.add(i);
        }
		
        if (nrOfStates == 4) {
            likelihoodCore = new MixtureBeerLikelihoodCore4();
        } else {
            likelihoodCore = new MixtureBeerLikelihoodCore(nrOfStates);
        }
        System.out.println("TreeLikelihood uses " + likelihoodCore.getClass().getName());

//        m_fProportionInvariant = m_siteModel.getProportianInvariant();
//        m_siteModel.setPropInvariantIsCategory(false);
//        if (m_fProportionInvariant > 0) {
//            calcConstantPatternIndices(nPatterns, nStateCount);
//        }

        initCore();

	}
	
	
    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                ((VariableCategorySiteModel) m_siteModel).maxCategoryCountInput.get(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

//	private void createNewTreeLikelihood() {
//		int i = categoryLikelihoods.size();		
//		TreeLikelihood treelikelihood = new TreeLikelihood();
//		treelikelihood.setID(getID() + i);
//		treelikelihood.getOutputs().add(this);
//		SingleCategorySiteModel siteModel = new SingleCategorySiteModel();
//		siteModel.initByName("category", i, "variableCategorySiteModel", m_siteModel, "substModel", m_siteModel.substModelInput.get());
//		
//		treelikelihood.initByName("data", dataInput.get(), 
//				"tree", treeInput.get(), 
//				"siteModel", siteModel, 
//				"branchRateModel", branchRateModelInput.get(), 
//				"useAmbiguities", m_useAmbiguities.get(),
//				"rootFrequencies", rootFrequenciesInput.get(),
//				"scaling" , Scaling.none // TODO: deal with scaling in calculatLogP
//	            // note the associated TODO in calculateLogLikelihoods()
//				);
//		treelikelihood.getOutputs().add(this);
//		likelihoodsInput.get().add(treelikelihood);
//		categoryLikelihoods.add(treelikelihood);
//	}

//	@Override
//	public double calculateLogP() {
//		if (!(m_siteModel instanceof VariableCategorySiteModel)) {
//			return super.calculateLogP();
//		}
//		
//		categoryCount = m_siteModel.getCategoryCount();
//		double [][] rootPartials = new double[categoryCount][];
//		for (int i = 0; i < categoryCount; i++) {
//			TreeLikelihood tl = categoryLikelihoods.get(i);
//			tl.calculateLogP();
//			rootPartials[i] = tl.getRootPartials();
//		}
//		
//		// integrate over category weights
//        final double[] proportions = m_siteModel.getCategoryProportions(null);
//		calculateIntegratePartials(rootPartials, proportions, m_fRootPartials);
//
//        double[] rootFrequencies = substitutionModel.getFrequencies();
//        if (rootFrequenciesInput.get() != null) {
//            rootFrequencies = rootFrequenciesInput.get().getFreqs();
//        }
//        calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
//
//        calcLogP();
//        return logP;
//	}
//	
	
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
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
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

                if (m_siteModel.integrateAcrossCategories()) {
                	((MixtureLikelihoodCore)likelihoodCore).calculatePartials(childNum1, childNum2, iNode, classes.size());
                } else {
                    throw new IllegalArgumentException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    //final double[] frequencies = //m_pFreqs.get().
                    final double[][] frequencies = new double[m_siteModel.getCategoryCount()][];
                    frequencies[0] = substitutionModel.getFrequencies();
                    for (int i = 1; i < m_siteModel.getCategoryCount(); i++) {
                    	frequencies[i] = frequencies[0];
                    }

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
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


//	@Override
//    protected void calcLogP() {
//        logP = 0.0;
//        if (useAscertainedSitePatterns) {
//            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
//            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
//                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
//            }
//        } else {
//            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
//                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
//            }
//        }
//    }
//
//	private void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
//        int v = 0;
//        for (int k = 0; k < nrOfPatterns; k++) {
//
//            double sum = 0.0;
//            for (int i = 0; i < nrOfStates; i++) {
//
//                sum += frequencies[i] * partials[v];
//                v++;
//            }
//            // TODO: deal with scaling in calculatLogP
//            // note the associated TODO in createNewTreeLikelihood()
//            outLogLikelihoods[k] = Math.log(sum); // + getLogScalingFactor(k);
//        }
//    }
//	
//	private void calculateIntegratePartials(double[][] rootPartials, double[] proportions, double[] outPartials) {
//        int u = 0;
//        int v = 0;
//        double [] inPartials = rootPartials[0];
//        for (int k = 0; k < nrOfPatterns; k++) {
//            for (int i = 0; i < nrOfStates; i++) {
//                outPartials[u] = inPartials[v] * proportions[0];
//                u++;
//                v++;
//            }
//        }
//
//        for (int l = 1; l < categoryCount; l++) {
//            u = 0;
//            v = 0;
//            inPartials = rootPartials[l];
//            for (int k = 0; k < nrOfPatterns; k++) {
//                for (int i = 0; i < nrOfStates; i++) {
//                    outPartials[u] += inPartials[v] * proportions[l];
//                    u++;
//                    v++;
//                }
//            }
//        }
//    }

	
	@Override
    public List<Input<?>> listInputs() {
    	List<Input<?>> list =  super.listInputs();
    	if (!ProgramStatus.name.equals("BEAUti") && System.getProperty("beast.is.junit.testing") == null) {
    		// do not expose internal likelihoods to BEAUti or junit tests
    		list.add(likelihoodsInput);
    	}
    	return list;
    }


	@Override
	protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
        if (m_siteModel.getCategoryCount() != classes.size()) {
        	setUpClasses();
            hasDirt = Tree.IS_FILTHY;
            return true;
        }

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
        	setUpClasses();
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
        	setUpClasses();
            return true;
        }
        if (treeInput.get().somethingIsDirty()) {
        	setUpClasses();
        	return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        return false;
	}


	private void setUpClasses() {
    	classes.clear();
		for (int k = 0; k < m_siteModel.getCategoryCount(); k++) {
			classes.add(k);
		}
	}
	

}
