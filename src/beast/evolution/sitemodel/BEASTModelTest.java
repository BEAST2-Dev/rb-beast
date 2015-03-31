package beast.evolution.sitemodel;

import java.util.Arrays;

import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.evolution.tree.Node;

@Description("Site model that jumps between with and without gamma sites, as well as with and without invariant sites")
public class BEASTModelTest extends SiteModel {

	public Input<BooleanParameter> hasGammaRatesInput = new Input<BooleanParameter>("hasGammaRates", "flag indicating whether gamma rate heterogeneity should be used", Validate.REQUIRED);
	public Input<BooleanParameter> hasInvariantSitesInput = new Input<BooleanParameter>("hasInvariantSites", "flag indicating whether invariant sites should be used", Validate.REQUIRED);

	BooleanParameter hasInvariantSites;
	BooleanParameter hasGammaRates;
	
	@Override
	public void initAndValidate() throws Exception {
		//BooleanParameter dummy = new BooleanParameter("1");
		hasInvariantSites = hasInvariantSitesInput.get();
		//hasInvariantSites.assignFromWithoutID(dummy);
		hasGammaRates = hasGammaRatesInput.get();
		//hasGammaRates.assignFromWithoutID(dummy);
		super.initAndValidate();
		
		// ensure categoryCount = gammaCategories + 1 by checking shape and invar parameters are present
		if (shapeParameterInput.get() == null) {
			throw new Exception("shape parameter must be specified");
		}
		if (shapeParameterInput.get().isEstimatedInput.get() == false) {
			throw new Exception("shape parameter must be estimated");
		}
		if (invarParameterInput.get() == null) {
			throw new Exception("proportionInvariant parameter must be specified");
		}
		if (invarParameterInput.get().isEstimatedInput.get() == false) {
			throw new Exception("proportionInvariant parameter must be estimated");
		}
	}

	@Override
    protected void calculateCategoryRates(final Node node) {
		Arrays.fill(categoryRates, 0.0);
		Arrays.fill(categoryProportions, 0.0);

		double propVariable = 1.0;
        int cat = 0;

        if (/*invarParameter != null && */hasInvariantSites.getValue()) {
            if (hasPropInvariantCategory) {
                categoryRates[0] = 0.0;
                categoryProportions[0] = invarParameter.getValue();
            }
            propVariable = 1.0 - invarParameter.getValue();
            if (hasPropInvariantCategory) {
                cat = 1;
            }
        } else {
            if (hasPropInvariantCategory) {
                categoryProportions[0] = 0.0;
                cat = 1;
            }
        }

        if (hasGammaRates.getValue()) {

            final double a = shapeParameter.getValue();
            double mean = 0.0;
            final int gammaCatCount = categoryCount - cat;

            final GammaDistribution g = new GammaDistributionImpl(a, 1.0 / a);
            for (int i = 0; i < gammaCatCount; i++) {
                try {
                    // RRB: alternative implementation that seems equally good in
                    // the first 5 significant digits, but uses a standard distribution object
                	if (useBeast1StyleGamma) {
                        categoryRates[i + cat] = GammaDistributionQuantile((2.0 * i + 1.0) / (2.0 * gammaCatCount), a, 1.0 / a);
                	} else {
                		categoryRates[i + cat] = g.inverseCumulativeProbability((2.0 * i + 1.0) / (2.0 * gammaCatCount));
                	}

                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Something went wrong with the gamma distribution calculation");
                    System.exit(-1);
                }
                mean += categoryRates[i + cat];

                categoryProportions[i + cat] = propVariable / gammaCatCount;
            }

            mean = (propVariable * mean) / gammaCatCount;

            for (int i = 0; i < gammaCatCount; i++) {
                categoryRates[i + cat] /= mean;
            }
        } else {
            categoryRates[cat] = 1.0 / propVariable;
            categoryProportions[cat] = propVariable;
        }

        // debugging code
//        double sum = 0;
//        for (double r : categoryProportions) {
//        	sum +=r;
//        }
//        if (!hasPropInvariantCategory) {
//        	sum += getProportionInvariant();
//        }
//        if (Math.abs(sum - 1.0) > 1e-10) {
//        	calculateCategoryRates(node);
//        	throw new RuntimeException("Incorrect proportions " + sum + " " + Arrays.toString(categoryProportions));
//        }
//        double meanrate = 0;
//        for (int i = 0; i < categoryProportions.length; i++) {
//        	meanrate += categoryProportions[i] * categoryRates[i];
//        }
//        if (Math.abs(meanrate - 1.0) > 1e-10) {
//        	throw new RuntimeException("Incorrect mean rate");
//        }

        ratesKnown = true;
    }

	
	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty = false;
		if (hasInvariantSites.somethingIsDirty() || hasGammaRates.somethingIsDirty()) {
			isDirty = true;
            ratesKnown = false;
		}
		if (super.requiresRecalculation()) {
			isDirty = true;
		}
		return isDirty;
	}
	
	@Override
    public double getProportionInvariant() {
        if (hasInvariantSites.getValue()) {
        	return invarParameter.getValue();
        } else {
        	return 0.0;
        }
    }
	
}