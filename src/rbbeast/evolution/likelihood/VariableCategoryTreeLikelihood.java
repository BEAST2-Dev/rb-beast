package rbbeast.evolution.likelihood;


import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.BeagleTreeLikelihood;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.likelihood.BeerLikelihoodCore4;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import rbbeast.evolution.sitemodel.VariableCategorySiteModel;

@Description("Tree likelihood that allows the number of rate categories to change. "
		+ "It uses a TreeLikelihood for each category, then integrates over the resulting root probabilites.")
public class VariableCategoryTreeLikelihood extends TreeLikelihood {
	
	private int nrOfStates;
	private int currentCategoryCount;
	
	List<LikelihoodCore> cores;
	List<MyBeagleTreeLikelihood> beagles;
	class MyBeagleTreeLikelihood extends BeagleTreeLikelihood {
		public boolean check() {
			return requiresRecalculation();
		}
	}
	
	@Override
	public void initAndValidate() {
		m_siteModel = (SiteModel.Base) siteModelInput.get();

		beagle = new MyBeagleTreeLikelihood();
        try {
	        beagle.initByName(
                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString(),
                    "rootFrequencies", rootFrequenciesInput.get());
	        if (beagle.getBeagle() != null) {
	        	
	        	beagles = new ArrayList<>();
	        	for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
	        		beagles.add(null);
	        	}
	        	beagles.add((MyBeagleTreeLikelihood)beagle);
	        	
	            //a Beagle instance was found, so we use it
	            return;
	        }
        } catch (Exception e) {
			// ignore
		}
        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

		
		
		// revert to TreeLikelihood default behaviour when site model is not a VariableCategorySiteModel
		if (!(m_siteModel instanceof VariableCategorySiteModel)) {
			super.initAndValidate();
			return;
		}

		
		nrOfStates = dataInput.get().getMaxStateCount();
		cores = new ArrayList<>();
		cores.add(null);

		// Site model must be a VariableCategorySiteModel once we get here
		// initialise, but make sure no beagle instance is created
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        System.setProperty("java.only", "true");
		super.initAndValidate();
        System.setProperty("java.only", forceJava + "");

		// clean up memory
		likelihoodCore = null;
		
	
		currentCategoryCount = m_siteModel.getCategoryCount();
		
        initCore(currentCategoryCount);

        likelihoodCore = cores.get(currentCategoryCount);
        System.out.println("TreeLikelihood uses " + likelihoodCore.getClass().getName());
        
	}
	
		
	@Override
    protected void initCore() {
		// do nothing
    }
	
    private void initBeagle(int i) {
        while (beagles.size() <= i) {
        	beagles.add(null);
        }
        beagle = new MyBeagleTreeLikelihood();
        try {
	        beagle.initByName(
                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString(),
                    "rootFrequencies", rootFrequenciesInput.get());
        } catch (Exception e) {
			// ignore
		}        
        beagles.set(i, (MyBeagleTreeLikelihood)beagle);
    }
	
    private void initCore(int i) {
    	System.err.println("Crearting likelihood core with " + i + " categories");
        final int nodeCount = treeInput.get().getNodeCount();
        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        while (cores.size() <= i) {
        	cores.add(null);
        }
        if (nrOfStates == 4) {
        	cores.set(i, new BeerLikelihoodCore4());
        } else {
        	cores.set(i, new BeerLikelihoodCore(nrOfStates));
        }
        
        
        cores.get(i).initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                i,
                true, m_useAmbiguities.get()
        );
        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount(), cores.get(i));
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount(), cores.get(i));
        }
        for (int j = 0; j < intNodeCount; j++) {
        	cores.get(i).createNodePartials(extNodeCount + j);
        }

        hasDirt = Tree.IS_FILTHY;
    }

    
    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount, LikelihoodCore likelihoodCore) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount, likelihoodCore);
            setStates(node.getRight(), patternCount, likelihoodCore);
        }
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}

	/**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(Node node, int patternCount, LikelihoodCore likelihoodCore) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {                
                double[] tipLikelihoods = data.getTipLikelihoods(taxonIndex,patternIndex_);
                if (tipLikelihoods != null) {
                	for (int state = 0; state < states; state++) {
                		partials[k++] = tipLikelihoods[state];
                	}
                }
                else {
                	int stateCount = data.getPattern(taxonIndex, patternIndex_);
	                boolean[] stateSet = data.getStateSet(stateCount);
	                for (int state = 0; state < states; state++) {
	                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
	                }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }
	
	@Override
	public void store() {
		super.store();
		setUpClasses();
	}
	
	@Override
	public void restore() {
		setUpClasses();
		super.restore();
	}


	@Override
	protected boolean requiresRecalculation() {
        if (beagle != null) {
            return ((MyBeagleTreeLikelihood)beagle).check();
        }
		
        hasDirt = Tree.IS_CLEAN;
        if (m_siteModel.getCategoryCount() != currentCategoryCount) {
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
		currentCategoryCount = m_siteModel.getCategoryCount();
		if (beagle != null) {
			if (currentCategoryCount >= beagles.size() || beagles.get(currentCategoryCount) == null) {
				initBeagle(currentCategoryCount);
			}
			beagle = beagles.get(currentCategoryCount);
		} else {
			if (currentCategoryCount >= cores.size() || cores.get(currentCategoryCount) == null) {
				initCore(currentCategoryCount);
			}
			likelihoodCore = cores.get(currentCategoryCount);
		}
	}
	

}
