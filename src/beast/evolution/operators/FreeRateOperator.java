package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

@Description("Generates proposals for free rate model")
public class FreeRateOperator extends Operator {
    public Input<RealParameter> weightInput =
            new Input<RealParameter>("weights", "weights of the various categories, should sum to 1", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);
    public final Input<Double> deltaInput = new Input<Double>("delta", "Magnitude of change for two randomly picked values.", 1.0);
	

	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		// enforce sum weight*rate = constant
		if (Randomizer.nextBoolean()) {
			// Weighted Delta Exchange on weights
			return doNormalisedDeltaExchange(weightInput.get(), rateParameterInput.get());
		} else {
			// Weighted Delta Exchange on weights
			return doDeltaExchange(rateParameterInput.get(), weightInput.get());
		}
	}

	private double doNormalisedDeltaExchange(RealParameter realparameter, RealParameter rates) {
        final int dim = realparameter.getDimension();
        if (dim <= 1) {
        	// it is impossible to select two distinct entries in this case, so there is nothing to propose 
        	return 0.0;
        }

        final int dim1 = Randomizer.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = Randomizer.nextInt(dim);
        }

        // operate on real parameter
        double scalar1 = realparameter.getValue(dim1);
        double scalar2 = realparameter.getValue(dim2);

        // exchange a random delta
        final double d = Randomizer.nextDouble() * deltaInput.get();
        scalar1 -= d;
        scalar2 += d;
        if (scalar1 < realparameter.getLower() || scalar1 > realparameter.getUpper() ||
                scalar2 < realparameter.getLower() || scalar2 > realparameter.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        } else {
            realparameter.setValue(dim1, scalar1);
            realparameter.setValue(dim2, scalar2);

            // normalise
        	double sum = 0;
        	for (int i = 0; i < dim; i++) {
        		sum += realparameter.getValue(i) * rates.getValue(i);
        	}
        	for (int i = 0; i < dim; i++) {
        		double newRate = rates.getValue(i) /sum;
        		if (newRate < rates.getLower() || newRate > rates.getUpper()) {
        			return Double.NEGATIVE_INFINITY;
        		}
      			rates.setValue(i, newRate);
         	}        	
        	
        	
        }
		return 0;
	}

	private double doDeltaExchange(RealParameter realparameter, RealParameter parameterWeights) {

        final int dim = realparameter.getDimension();
        if (dim <= 1) {
        	// it is impossible to select two distinct entries in this case, so there is nothing to propose 
        	return 0.0;
        }

        final int dim1 = Randomizer.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = Randomizer.nextInt(dim);
        }

        // operate on real parameter
        double scalar1 = realparameter.getValue(dim1);
        double scalar2 = realparameter.getValue(dim2);

        // exchange a random delta
        final double d = Randomizer.nextDouble() * deltaInput.get();
        scalar1 -= d;
        scalar2 += d * (double) parameterWeights.getValue(dim1) / (double) parameterWeights.getValue(dim2);

        if (scalar1 < realparameter.getLower() || scalar1 > realparameter.getUpper() ||
                scalar2 < realparameter.getLower() || scalar2 > realparameter.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        } else {
            realparameter.setValue(dim1, scalar1);
            realparameter.setValue(dim2, scalar2);
        }
		return 0;
	}

}
