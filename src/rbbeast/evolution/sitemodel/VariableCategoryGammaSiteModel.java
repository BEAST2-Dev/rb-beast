package rbbeast.evolution.sitemodel;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.inference.parameter.IntegerParameter;

@Description("As the gamma site heterogeity model, but with the ability to change catagory counts")
public class VariableCategoryGammaSiteModel extends SiteModel {
    final public Input<IntegerParameter> categoryCountInput =
            new Input<>("categoryCount", "gamma category count -- can vary during MCMC and overrides gammaCategoryCount input", new IntegerParameter("1"));

    private int propInvarCount;

    @Override
    protected void refresh() {
        propInvarCount = 0;
        if (shapeParameter != null) {
            categoryCount = categoryCountInput.get().getValue();
            if (categoryCount < 1) {
            	if (categoryCount < 0) {
            		Log.warning.println("SiteModel: Invalid category count (" + categoryCount + ") Setting category count to 1");
            	}
                categoryCount = 1;
                propInvarCount = 1;
            }
        } else {
            categoryCount = 1;
            propInvarCount = 1;
        }

        if (/*invarParameter != null && */invarParameter.getValue() > 0) {
            if (invarParameter.getValue() >= 1.0) {
            	throw new RuntimeException("Wrong value for parameter " + invarParameter.getID() +
            			". Proportion invariant should be in bewteen 0 and 1 (exclusive)");
            }
            if (hasPropInvariantCategory) {
                categoryCount += 1;
            }
        }

        categoryRates = new double[categoryCount];
        categoryProportions = new double[categoryCount];
        calculateCategoryRates(null);
        //ratesKnown = false;
    }

    @Override
    protected boolean requiresRecalculation() {
    	
    	if (categoryCount - propInvarCount != categoryCountInput.get().getValue()) {
    		refresh();
            ratesKnown = false;
    	}
    	
        // do explicit check whether any of the non-substitution model parameters changed
        if (categoryCount > 1) {
            if (shapeParameter != null && shapeParameter.somethingIsDirty() ||
                    muParameter.somethingIsDirty() ||
                    invarParameter.somethingIsDirty()) {
                ratesKnown = false;
            }
        } else {
            if (muParameter.somethingIsDirty() || !hasPropInvariantCategory && invarParameter.somethingIsDirty()) {
                ratesKnown = false;
            }
        }
//    	ratesKnown = false;
        // we only get here if something is dirty in its inputs, so always return true
        return true;
    }
    
    @Override
    public int getCategoryCount() {
    	return categoryCountInput.get().getValue() + propInvarCount;
    }

    @Override
    public void restore() {
    	if (categoryCount - propInvarCount != categoryCountInput.get().getValue()) {
    		refresh();
    	}
    	super.restore();
    }

}
