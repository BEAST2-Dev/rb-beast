package rbbeast.evolution.sitemodel;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;

// Soubrier J, Steel M, Lee MS, Der Sarkissian C, Guindon S, Ho SY, Cooper A. The influence of rate heterogeneity among sites on the time dependence of molecular rates. Molecular biology and evolution. 2012 Nov 1;29(11):3345-58.
// https://doi.org/10.1093/molbev/mss140
// Le SQ, Dang CC, Gascuel O. Modeling protein evolution with several amino acid replacement matrices depending on site rates. Molecular biology and evolution. 2012 Oct 1;29(10):2921-36.
// https://academic.oup.com/mbe/article/29/10/2921/1027701
@Description("Bayesian variant of free rate model ")
public class FreeRateModel extends VariableCategorySiteModel {

    public Input<RealParameter> muParameterInput = new Input<RealParameter>("mutationRate", "mutation rate (defaults to 1.0)");
    public Input<RealParameter> weightInput =
            new Input<RealParameter>("weights", "weights of the various categories, should sum to 1", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);
    public Input<RealParameter> invarParameterInput =
            new Input<RealParameter>("proportionInvariant", "proportion of sites that is invariant: should be between 0 (default) and 1");
    public Input<Boolean> incrementalInput = new Input<>("incremental", "whether to add rates incermentally so the rate categories remain in order", false);
    public Input<Integer> categoryCountInput =
            new Input<>("categoryCount", "number of categories used. If set, dimension of weight and rate parameters will be updated and values of rates and weights set to 1/count, "
            		+ "otherwise the parameter dimensions are used as category count and values are left unchanged", -1);
    
    RealParameter muParameter;
    RealParameter rateParameter;
    RealParameter weightParameter;
    RealParameter invarParameter;
    boolean incremental;
    
    @Override
    public void initAndValidate() {
    	incremental = incrementalInput.get();
        muParameter = muParameterInput.get();
        if (muParameter == null) {
            muParameter = new RealParameter("1.0");
        }
        rateParameter = rateParameterInput.get();
        weightParameter = weightInput.get();

    	int dim = categoryCountInput.get() > 0 ? categoryCountInput.get() : 
    			Math.max(rateParameter.getDimension(), weightParameter.getDimension());
        if (rateParameter.getDimension() != dim) {
        	Log.warning.println("Warning: setting dimesnio of " + rateParameter.getID() + " to " + dim);
        	rateParameter.setDimension(dim);
        }
        if (weightParameter.getDimension() != dim) {
        	Log.warning.println("Warning: setting dimesnio of " + weightParameter.getID() + " to " + dim);
        	weightParameter.setDimension(dim);
        }

        if (categoryCountInput.get() > 0) {
	    	for (int i = 0; i < dim; i++) {
	    		rateParameter.setValue(i, 1.0/dim);
	    		weightParameter.setValue(i, 1.0/dim);
	    	}
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
        int categoryCount = rateParameter.getDimension();

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
        return rateParameter.getDimension();
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
        int categoryCount = rateParameter.getDimension();

        if (/*invarParameter != null && */invarParameter.getValue() > 0) {
            if (invarParameter.getValue() >= 1.0) {
            	throw new RuntimeException("Wrong value for parameter " + invarParameter.getID() +
            			". Proportion invariant should be in bewteen 0 and 1 (exclusive)");
            }
            if (hasPropInvariantCategory) {
                categoryCount += 1;
            }
        }
    	if (categoryRates.length != categoryCount) {
    		categoryRates = new double[categoryCount];
            categoryProportions = new double[categoryCount];
    	}
    	
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
        
        if (rateParameter.getDimension() != weightParameter.getDimension()) {
        	int h  =3;
        	h++;
        }

        double meanRate = 0;
        double sumWeight = 0;
	    for (int i = 0; i < rateParameter.getDimension(); i++) {
	        categoryRates[i + cat] = rateParameter.getValue(i);
	        if (incremental && i > 0) {
	        	categoryRates[i + cat] += categoryRates[i + cat - 1];
	        }
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
        if (rateParameter != null && rateParameter.somethingIsDirty() ||
        		weightParameter.somethingIsDirty() ||
                muParameter.somethingIsDirty() ||
                invarParameter.somethingIsDirty()) {
            ratesKnown = false;
        } else {
            if (muParameter.somethingIsDirty() || !hasPropInvariantCategory && invarParameter.somethingIsDirty()) {
                ratesKnown = false;
            }
        }

        // we only get here if something is dirty in its inputs, so always return true
        return true;
    }

    protected boolean ratesKnown;

    protected double[] categoryRates;

    protected double[] categoryProportions;

}
