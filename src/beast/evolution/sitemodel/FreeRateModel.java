package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

@Description("Lars Jermijn et al's free rate model")
public class FreeRateModel extends SiteModel.Base {

    public Input<RealParameter> muParameterInput = new Input<RealParameter>("mutationRate", "mutation rate (defaults to 1.0)");
    public Input<RealParameter> weightInput =
            new Input<RealParameter>("weights", "weights of the various categories, should sum to 1", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);
    public Input<RealParameter> invarParameterInput =
            new Input<RealParameter>("proportionInvariant", "proportion of sites that is invariant: should be between 0 (default) and 1");

    RealParameter muParameter;
    RealParameter rateParameter;
    RealParameter weightParameter;
    RealParameter invarParameter;
    
    @Override
    public void initAndValidate() {
        muParameter = muParameterInput.get();
        if (muParameter == null) {
            muParameter = new RealParameter("1.0");
        }
        rateParameter = rateParameterInput.get();
        weightParameter = weightInput.get();
    	int dim = Math.max(rateParameter.getDimension(), weightParameter.getDimension());
        if (rateParameter.getDimension() != dim) {
        	Log.warning.println("Warning: setting dimesnio of " + rateParameter.getID() + " to " + dim);
        	rateParameter.setDimension(dim);
        }
        if (weightParameter.getDimension() != dim) {
        	Log.warning.println("Warning: setting dimesnio of " + weightParameter.getID() + " to " + dim);
        	weightParameter.setDimension(dim);
        }
        
        invarParameter = invarParameterInput.get();
        if (invarParameter == null) {
            invarParameter = new RealParameter("0.0");
            invarParameter.setBounds(Math.max(0.0, invarParameter.getLower()), Math.min(1.0, invarParameter.getUpper()));
        }

        //if (muParameter != null) {
        muParameter.setBounds(Math.max(muParameter.getLower(), 0.0), Math.min(muParameter.getUpper(), Double.POSITIVE_INFINITY));
        //}
        if (rateParameter != null) {
            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            rateParameter.setBounds(Math.max(rateParameter.getLower(), 1.0E-3), Math.min(rateParameter.getUpper(), 1.0E3));
        }


        if (/*invarParameter != null && */(invarParameter.getValue() < 0 || invarParameter.getValue() > 1)) {
            throw new IllegalArgumentException("proportion invariant should be between 0 and 1");
        }
        refresh();

        addCondition(muParameterInput);
        addCondition(invarParameterInput);
        addCondition(rateParameterInput);
    }

    @Override
    protected void refresh() {
        categoryCount = rateParameter.getDimension();

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


    // *****************************************************************
    // Interface SiteModel
    // *****************************************************************

    @Override
    public boolean integrateAcrossCategories() {
        return true;
    }

    @Override
    public int getCategoryCount() {
        return categoryCount;
    }

    @Override
    public int getCategoryOfSite(final int site, final Node node) {
        throw new IllegalArgumentException("Integrating across categories");
    }

    @Override
    public double getRateForCategory(final int category, final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        //final double mu = (muParameter != null) ? muParameter.getValue() : 1.0;

        return categoryRates[category] * muParameter.getValue();
    }


    /**
     * return category rates
     *
     * @param node rates to which the rates apply. Typically, the rates will be uniform
     *             throughout the tree and the node argument is ignored.
     */
    @Override
    public double[] getCategoryRates(final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        final double mu = muParameter.getValue();//(muParameter != null) ? muParameter.getValue() : 1.0;

        final double[] rates = new double[categoryRates.length];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = categoryRates[i] * mu;
        }

        return rates;
    }

    /**
     * @return substitution model *
     */
    @Override
    public SubstitutionModel getSubstitutionModel() {
        return substModelInput.get();
    }

    /**
     * Get the expected proportion of sites in this category.
     *
     * @param category the category number
     * @param node     node to which the proportions apply. Typically, proportions
     *                 will be uniform throughout the tree and this argument is ignored.
     * @return the proportion.
     */
    @Override
    public double getProportionForCategory(final int category, final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        return categoryProportions[category];
    }

    /**
     * Get an array of the expected proportion of sites in this category.
     *
     * @return an array of the proportion.
     */
    @Override
    public double[] getCategoryProportions(final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        return categoryProportions;
    }

    /**
     * discretisation of gamma distribution with equal proportions in each
     * category
     * @param node
     */
    protected void calculateCategoryRates(final Node node) {
        double propVariable = 1.0;
        int cat = 0;

        if (/*invarParameter != null && */invarParameter.getValue() > 0) {
            if (hasPropInvariantCategory) {
                categoryRates[0] = 0.0;
                categoryProportions[0] = invarParameter.getValue();
            }
            propVariable = 1.0 - invarParameter.getValue();
            if (hasPropInvariantCategory) {
                cat = 1;
            }
        }

        double meanRate = 0;
        double sumWeight = 0;
	    for (int i = 0; i < rateParameter.getDimension(); i++) {
	        categoryRates[i + cat] = rateParameter.getValue(i);
	        categoryProportions[i + cat] = weightParameter.getValue(i);
	        sumWeight += categoryProportions[i + cat];
	        meanRate += categoryRates[i + cat] * categoryProportions[i + cat];
	    }
	
	    meanRate = (propVariable * meanRate) / sumWeight;
	
	    for (int i = 0; i < categoryRates.length; i++) {
	        categoryRates[i + cat] /= meanRate;
	        categoryProportions[i] *= propVariable/sumWeight;
	    }

        ratesKnown = true;
    }


    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        super.store();
    } // no additional state needs storing

    @Override
    public void restore() {
        super.restore();
        ratesKnown = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        // do explicit check whether any of the non-substitution model parameters changed
        if (categoryCount > 1) {
            if (rateParameter != null && rateParameter.somethingIsDirty() ||
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

    protected boolean ratesKnown;

    protected int categoryCount;

    protected double[] categoryRates;

    protected double[] categoryProportions;

}
