package rbbeast.evolution.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import rbbeast.evolution.sitemodel.FreeRateModel;

@Description("Proposal to expands or contracts rates and weights for free rate model")
public class FreeRateDimensionMover extends Operator {
    public Input<RealParameter> weightInput =
            new Input<RealParameter>("weights", "weights of the various categories, should sum to 1", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);
    public Input<FreeRateModel> modelInput = new Input<>("model", "free rate model where the weights and rates are used", Validate.REQUIRED);
	

    private int maxCategoryCount;
    private RealParameter rates;
    private RealParameter weights;
    
	@Override
	public void initAndValidate() {
		maxCategoryCount = modelInput.get().maxCategoryCountInput.get();
		rates = rateParameterInput.get();
		weights = weightInput.get();
	}

	@Override
	public double proposal() {
		double HR = 0;

		if (Randomizer.nextBoolean()) {
			// expand # categories
			HR = expand();
		} else {
			// reduce # categories
			HR = contract();
		}
		return HR;
	}

	private double contract() {
		int dim = rates.getDimension();
		if (dim == 1) {
			return Double.NEGATIVE_INFINITY;
		}
		rates.setDimension(dim - 1);
		weights.setDimension(dim - 1);
		
		normalise(weights);
		normalise(rates);
		return 0;
	}

	private void normalise(RealParameter weights) {
		int dim = weights.getDimension();
		double sum = 0;
		for (int i = 0; i < dim; i++) {
			sum += weights.getValue(i);
		}
		for (int i = 0; i < dim; i++) {
			weights.setValue(i, weights.getValue(i)/sum);
		}
	}

	private double expand() {
		int dim = rates.getDimension();
		if (dim == maxCategoryCount) {
			return Double.NEGATIVE_INFINITY;
		}
		rates.setDimension(dim + 1);
		weights.setDimension(dim + 1);

		weights.setValue(dim, 1.0/(dim));
		normalise(weights);
		rates.setValue(dim, 1.0/(dim));
		normalise(rates);
		return 0;
	}

}
