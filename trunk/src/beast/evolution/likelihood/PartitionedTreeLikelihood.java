/*
* File TreeLikelihood.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/


package beast.evolution.likelihood;


import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.RejectedExecutionException;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Calculates the likelihood of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
//public class PartitionedTreeLikelihood extends Distribution {
public class PartitionedTreeLikelihood extends TreeLikelihood {

//    public Input<Alignment> m_data = new Input<Alignment>("data", "sequence data for the beast.tree", Validate.REQUIRED);
//    public Input<Tree> m_tree = new Input<Tree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
//    public Input<BranchRateModel.Base> m_pBranchRateModel = new Input<BranchRateModel.Base>("branchRateModel",
//            "A model describing the rates on the branches of the beast.tree.");
//    public Input<Boolean> m_useAmbiguities = new Input<Boolean>("useAmbiguities", "flag to indicate leafs that sites containing ambigue states should be handled instead of ignored (the default)", false);

    public Input<PartitionProvider> partitionProviderInput = new Input<PartitionProvider>("partitionProvider", "provides information on how to split the alignment into partitions", Validate.REQUIRED);
    
    public PartitionedTreeLikelihood() {
    	m_pSiteModel.setRule(Validate.OPTIONAL);
    }
    
    /**
     * calculation engine *
     */
    protected PartitionedLikelihoodCore [] m_likelihoodCore;
    BeagleTreeLikelihood m_beagle;

    /**
     * Plugin associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    SubstitutionModel.Base [] m_substitutionModel;
    protected SiteModel.Base [] m_siteModel;
    BranchRateModel.Base m_branchRateModel;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    int m_nHasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
//    double[] m_branchLengths;
//    double[] m_StoredBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    double[][] m_fPatternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    double[][] m_fRootPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    double[][] m_fProbabilities;

    int m_nMatrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    boolean m_bAscertainedSitePatterns = false;
    
    /** number of partitions used for splitting up the alignment **/
    int partitionCount;
    PartitionProvider partitionProvider;
    
    
    
    @Override
    public void initAndValidate() throws Exception {
        // sanity check: alignment should have same #taxa as tree
        if (m_data.get().getNrTaxa() != m_tree.get().getLeafNodeCount()) {
            throw new Exception("The number of nodes in the tree does not match the number of sequences");
        }
        m_beagle = null;
//        m_beagle = new BeagleTreeLikelihood();
//        m_beagle.initByName("data", m_data.get(), "tree", m_tree.get(), "siteModel", m_pSiteModel.get(), "branchRateModel", m_pBranchRateModel.get(), "useAmbiguities", m_useAmbiguities.get());
//        if (m_beagle.beagle != null) {
//            //a Beagle instance was found, so we use it
//            return;
//        }
        // No Beagle instance was found, so we use the good old java likelihood core
        m_beagle = null;

        if (m_pBranchRateModel.get() != null) {
            m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }

        if (m_pSiteModel.get() != null) {
        	m_pSiteModel.get().setDataType(m_data.get().getDataType());
        }
        
        partitionProvider = partitionProviderInput.get();
        if (partitionProvider.delayInitialisation()) {
            m_siteModel = new SiteModel[0];
            m_substitutionModel = new SubstitutionModel.Base[0];
        	return;
        }
        
        int nodeCount = m_tree.get().getNodeCount();
        partitionCount = partitionProvider.getPartitionCount();
        m_siteModel = new SiteModel[partitionCount];
        m_substitutionModel = new SubstitutionModel.Base[partitionCount];
        for (int i = 0; i < partitionCount; i++) {
            m_siteModel[i] = partitionProvider.m_pSiteModel.get().get(i);        	
            m_siteModel[i].setDataType(m_data.get().getDataType());
            m_substitutionModel[i] = m_siteModel[i].m_pSubstModel.get();
        }

        m_branchLengths = new double[nodeCount];
        m_StoredBranchLengths = new double[nodeCount];

        int nStateCount = m_data.get().getMaxStateCount();
        int nPatterns = m_data.get().getPatternCount();
        m_likelihoodCore = new PartitionedBeerLikelihoodCore[partitionCount];
        for (int i = 0; i < partitionCount; i++) {
	        if (nStateCount == 4) {
	            m_likelihoodCore[i] = new PartitionedBeerLikelihoodCore4();
	        } else {
	            m_likelihoodCore[i] = new PartitionedBeerLikelihoodCore(nStateCount);
	        }
        }
        System.out.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        for (int i = 0; i < partitionCount; i++) {
        	initCore(i);
        }

        m_fPatternLogLikelihoods = new double[partitionCount][nPatterns];
        m_fRootPartials = new double[partitionCount][nPatterns * nStateCount];
        m_nMatrixSize = (nStateCount + 1) * (nStateCount + 1);
        m_fProbabilities = new double[partitionCount][(nStateCount + 1) * (nStateCount + 1)];
        for (int i = 0; i < partitionCount; i++) {
        	Arrays.fill(m_fProbabilities[i], 1.0);
        }
        
        if (m_data.get() instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }
    }



    void initCore(int partition) {
        int nodeCount = m_tree.get().getNodeCount();
        int extNodeCount = nodeCount / 2 + 1;
        int intNodeCount = nodeCount / 2;
        
        m_likelihoodCore[partition].initialize(
                nodeCount,
                m_data.get().getPatternCount(),
                m_siteModel[partition].getCategoryCount(),
                true, m_useAmbiguities.get()
        );


        if (m_useAmbiguities.get()) {
            setPartials(partition, m_tree.get().getRoot(), m_data.get().getPatternCount());
        } else {
            setStates(partition, m_tree.get().getRoot(), m_data.get().getPatternCount());
        }
        m_nHasDirt = Tree.IS_FILTHY;
        for (int j = 0; j < intNodeCount; j++) {
            m_likelihoodCore[partition].createNodePartials(extNodeCount + j);
        }
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
//    public void sample(State state, Random random) {
//        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
//    }

    /**
     * set leaf states in likelihood core *
     */
    void setStates(int partition, Node node, int patternCount) {
        if (node.isLeaf()) {
            int i;
            int[] states = new int[patternCount];
            int iTaxon = m_data.get().getTaxonIndex(node.getID());
            for (i = 0; i < patternCount; i++) {
                states[i] = m_data.get().getPattern(iTaxon, i);
            }
            m_likelihoodCore[partition].setNodeStates(node.getNr(), states);

        } else {
            setStates(partition, node.getLeft(), patternCount);
            setStates(partition, node.getRight(), patternCount);
        }
    }

    /**
     * set leaf partials in likelihood core *
     */
    void setPartials(int partition, Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = m_data.get();
            int nStates = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * nStates];

            int k = 0;
            int iTaxon = m_data.get().getTaxonIndex(node.getID());
            for (int iPattern = 0; iPattern < patternCount; iPattern++) {
                int nState = data.getPattern(iTaxon, iPattern);
                boolean[] stateSet = data.getStateSet(nState);
                for (int iState = 0; iState < nStates; iState++) {
                    partials[k++] = (stateSet[iState] ? 1.0 : 0.0);
                }
            }
            m_likelihoodCore[partition].setNodePartials(node.getNr(), partials);

        } else {
            setPartials(partition, node.getLeft(), patternCount);
            setPartials(partition, node.getRight(), patternCount);
        }
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() throws Exception {
        if (m_beagle != null) {
            logP = m_beagle.calculateLogP();
            return logP;
        }
        Tree tree = m_tree.get();
        
//        for (int i = 0; i < partitionCount; i++) {
//        	List<Integer> partitionIndicator = partitionProvider.getPatternIndicators(i);
//        	if (partitionIndicator.size() > 0) {
//        		traverse(i, partitionIndicator, tree.getRoot());
//        	}
//        }
        for (Node node : tree.getNodesAsArray()) {
	        double branchRate = m_branchRateModel.getRateForBranch(node);
	        double branchTime = node.getLength() * branchRate;
	        m_branchLengths[node.getNr()] = branchTime;
        }

        
        threadedTraverse(tree.getRoot());
        calcLogP();
        
// TODO make scaline work
// note that the likelihood core needs updating: to know which patterns are scaled per partition
//        m_nScale++;
//        if (logP > 0 || (m_likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
//        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
//            m_nScale = 0;
//            m_fScale *= 1.01;
//            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
//            m_likelihoodCore.setUseScaling(m_fScale);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
//        }
        return logP;
    }

    void calcLogP() throws Exception {
        logP = 0.0;
        if (m_bAscertainedSitePatterns) {
        	for (int k = 0; k < partitionCount; k++) {
        		int [] weights = partitionProvider.getPatternWeights(k);
        		List<Integer> patternIndicators = partitionProvider.getPatternIndicators(k);
	            double ascertainmentCorrection = ((AscertainedAlignment) m_data.get()).getAscertainmentCorrection(m_fPatternLogLikelihoods[k]);
	            for (int i : patternIndicators) {
	                logP += (m_fPatternLogLikelihoods[k][i] - ascertainmentCorrection) * weights[i];
	            }
        	}
        } else {
        	for (int k = 0; k < partitionCount; k++) {
        		int [] weights = partitionProvider.getPatternWeights(k);
        		List<Integer> patternIndicators = partitionProvider.getPatternIndicators(k);
	            for (int i : patternIndicators) {
	                logP += m_fPatternLogLikelihoods[k][i] * weights[i];
	            }
        	}
        }
    }

    
	class CoreRunnable implements Runnable {
		int partition;
		private Node m_root;
		List<Integer> patternIndicator;
		
		CoreRunnable(int partition, Node root, List<Integer> patternIndicator) {
			    this.partition = partition;
				m_root = root;
				this.patternIndicator = patternIndicator;
		}

        public void run() {
  		  	try {
  		  		traverse(partition, patternIndicator, m_root);
  		  	} catch (Exception e) {
  		  		System.err.println("Something went wrong in a traversal thread for partition" + partition);
				e.printStackTrace();
				System.exit(0);
			}
  		    m_nCountDown.countDown();
        }

	} // CoreRunnable
    
	CountDownLatch m_nCountDown;
	
	void threadedTraverse(Node root) throws Exception {
		try {
			int nPatterns = m_fPatternLogLikelihoods.length;
			if (partitionCount >= 1) {
				m_nCountDown = new CountDownLatch(partitionCount);
				// kick off the threads
		    	for (int partition = 0; partition < partitionCount; partition++) {
		    		CoreRunnable coreRunnable = new CoreRunnable(partition, root, partitionProvider.getPatternIndicators(partition));
					BeastMCMC.g_exec.execute(coreRunnable);
		    	}
				m_nCountDown.await();
			} else {
		    	for (int iThread = 0; iThread < partitionCount; iThread++) {
		    		traverse(iThread, partitionProvider.getPatternIndicators(iThread), root);
		    	}
			}
		} catch (RejectedExecutionException e) {
			System.err.println("Reducing nr of threads to 1 because " + e.getMessage());
			// refresh thread pool
	    	for (int iThread = 0; iThread < partitionCount; iThread++) {
	    		traverse(iThread, partitionProvider.getPatternIndicators(iThread), root);
	    	}
		}
    }
    
    
    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(int partition, List<Integer> patternIndicator, Node node) throws Exception {

        int update = (node.isDirty() | m_nHasDirt);

        int iNode = node.getNr();

//        double branchRate = m_branchRateModel.getRateForBranch(node);
//        double branchTime = node.getLength() * branchRate;
//        m_branchLengths[iNode] = branchTime;
        double branchTime = m_branchLengths[iNode]; 
        double branchRate = branchTime / node.getLength();
        
        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[iNode])) {
            Node parent = node.getParent();
            m_likelihoodCore[partition].setNodeMatrixForUpdate(iNode);
            for (int i = 0; i < m_siteModel[partition].getCategoryCount(); i++) {
                double jointBranchRate = m_siteModel[partition].getRateForCategory(i, node) * branchRate;
                m_substitutionModel[partition].getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, m_fProbabilities[partition]);
                m_likelihoodCore[partition].setNodeMatrix(iNode, i, m_fProbabilities[partition]);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            Node child1 = node.getLeft(); //Two children
            int update1 = traverse(partition, patternIndicator, child1);

            Node child2 = node.getRight();
            int update2 = traverse(partition, patternIndicator, child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int childNum1 = child1.getNr();
                int childNum2 = child2.getNr();

                m_likelihoodCore[partition].setNodePartialsForUpdate(iNode);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    m_likelihoodCore[partition].setNodeStatesForUpdate(iNode);
                }

                if (m_siteModel[partition].integrateAcrossCategories()) {
                    m_likelihoodCore[partition].calculatePartials(childNum1, childNum2, iNode, patternIndicator);
                } else {
                    throw new Exception("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = //m_pFreqs.get().
                            m_substitutionModel[partition].getFrequencies();

                    double[] proportions = m_siteModel[partition].getCategoryProportions(node);
                    m_likelihoodCore[partition].integratePartials(node.getNr(), proportions, m_fRootPartials[partition], patternIndicator);

                    m_likelihoodCore[partition].calculateLogLikelihoods(m_fRootPartials[partition], frequencies, m_fPatternLogLikelihoods[partition], patternIndicator);
                }

            }
        }
        return update;
    } // traverseWithBRM

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        if (m_beagle != null) {
            return m_beagle.requiresRecalculation();
        }
        m_nHasDirt = Tree.IS_CLEAN;

        if (m_data.get().isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_FILTHY;
            return true;
        }
        for (SiteModel.Base siteModel : m_siteModel) {
	        if (siteModel.isDirtyCalculation()) {
	            m_nHasDirt = Tree.IS_DIRTY;
	            return true;
	        }
        }
        if (partitionProvider.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (m_branchRateModel != null && m_branchRateModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return m_tree.get().somethingIsDirty();
    }

    @Override
    public void store() {
        if (m_beagle != null) {
            m_beagle.store();
            super.store();
            return;
        }
        if (m_likelihoodCore != null) {
        	for (int i = 0; i < partitionCount; i++) {
        		m_likelihoodCore[i].store();
        	}
        }
        super.store();
 //       System.arraycopy(m_branchLengths, 0, m_StoredBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
        if (m_beagle != null) {
            m_beagle.restore();
            super.restore();
            return;
        }
        if (m_likelihoodCore != null) {
        	for (int i = 0; i < partitionCount; i++) {
        		m_likelihoodCore[i].restore();
        	}
        }
        super.restore();
//        double[] tmp = m_branchLengths;
//        m_branchLengths = m_StoredBranchLengths;
//        m_StoredBranchLengths = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    public List<String> getArguments() {
        return Collections.singletonList(m_data.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    public List<String> getConditions() {
        return m_siteModel[0].getConditions();
    }

} // class TreeLikelihood
