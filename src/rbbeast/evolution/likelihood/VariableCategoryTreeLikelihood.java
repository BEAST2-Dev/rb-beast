package rbbeast.evolution.likelihood;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import rbbeast.evolution.sitemodel.SingleCategorySiteModel;
import rbbeast.evolution.sitemodel.VariableCategorySiteModel;

@Description("Tree likelihood that allows the number of rate categories to change. "
		+ "It uses a TreeLikelihood for each category, then integrates over the resulting root probabilites.")
public class VariableCategoryTreeLikelihood extends TreeLikelihood {
	
    /** private list of likelihoods, to notify framework of TreeLikelihoods being created **/
    final private Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("*","",new ArrayList<>());

	private List<TreeLikelihood> categoryLikelihoods;
	private int nrOfPatterns, nrOfStates, categoryCount;

	
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
		
		nrOfPatterns = dataInput.get().getPatternCount();
		nrOfStates = substitutionModel.getStateCount();
		categoryCount = m_siteModel.getCategoryCount();
		
		// create TreeLikelihoods, one for each potential category
		categoryLikelihoods = new ArrayList<>();
		for (int i = 0; i < ((VariableCategorySiteModel) m_siteModel).maxCategoryCountInput.get(); i++) {
			createNewTreeLikelihood();
		}
	}
	
	private void createNewTreeLikelihood() {
		int i = categoryLikelihoods.size();		
		TreeLikelihood treelikelihood = new TreeLikelihood();
		treelikelihood.setID(getID() + i);
		treelikelihood.getOutputs().add(this);
		SingleCategorySiteModel siteModel = new SingleCategorySiteModel(i, siteModelInput.get());
		
		treelikelihood.initByName("data", dataInput.get(), 
				"tree", treeInput.get(), 
				"siteModel", siteModel, 
				"branchRateModel", branchRateModelInput.get(), 
				"useAmbiguities", m_useAmbiguities.get(),
				"rootFrequencies", rootFrequenciesInput.get(),
				"scaling" , Scaling.none // TODO: deal with scaling in calculatLogP
	            // note the associated TODO in calculateLogLikelihoods()
				);
		treelikelihood.getOutputs().add(this);
		likelihoodsInput.get().add(treelikelihood);
		categoryLikelihoods.add(treelikelihood);
	}

	@Override
	public double calculateLogP() {
		if (!(m_siteModel instanceof VariableCategorySiteModel)) {
			return super.calculateLogP();
		}
		
		categoryCount = m_siteModel.getCategoryCount();
		double [][] rootPartials = new double[categoryCount][];
		for (int i = 0; i < categoryCount; i++) {
			TreeLikelihood tl = categoryLikelihoods.get(i);
			tl.calculateLogP();
			tl.getRootPartials();
		}
		
		// integrate over category weights
        final double[] proportions = m_siteModel.getCategoryProportions(null);
		calculateIntegratePartials(rootPartials, proportions, m_fRootPartials);

        double[] rootFrequencies = substitutionModel.getFrequencies();
        if (rootFrequenciesInput.get() != null) {
            rootFrequencies = rootFrequenciesInput.get().getFreqs();
        }
        calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);

        calcLogP();
        return logP;
	}
	
	@Override
    protected void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }
    }

	private void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {

                sum += frequencies[i] * partials[v];
                v++;
            }
            // TODO: deal with scaling in calculatLogP
            // note the associated TODO in createNewTreeLikelihood()
            outLogLikelihoods[k] = Math.log(sum); // + getLogScalingFactor(k);
        }
    }
	
	private void calculateIntegratePartials(double[][] rootPartials, double[] proportions, double[] outPartials) {
        int u = 0;
        int v = 0;
        double [] inPartials = rootPartials[0];
        for (int k = 0; k < nrOfPatterns; k++) {
            for (int i = 0; i < nrOfStates; i++) {
                outPartials[u] = inPartials[v] * proportions[0];
                u++;
                v++;
            }
        }

        for (int l = 1; l < categoryCount; l++) {
            u = 0;
            inPartials = rootPartials[l];
            for (int k = 0; k < nrOfPatterns; k++) {
                for (int i = 0; i < nrOfStates; i++) {
                    outPartials[u] += inPartials[v] * proportions[l];
                    u++;
                    v++;
                }
            }
        }
    }

	
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
		return super.requiresRecalculation();
	}
}
