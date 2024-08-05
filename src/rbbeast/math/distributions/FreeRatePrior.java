package rbbeast.math.distributions;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;

@Description("Prior that ensures rates repel each other in order to ensure identifiability")
public class FreeRatePrior extends Distribution {
    public Input<RealParameter> rateParameterInput =
            new Input<>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);

    public Input<Double> thresholdInput =
            new Input<>("threshold", "fraction below which each of pairs incurs a penalty", 1.1);
    
    public Input<Double> penaltyInput =
            new Input<>("penalty", "size of penalty when fraction below which threshold", 1.0);
    
    private double threshold, penalty;
    
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	threshold = thresholdInput.get();
    	penalty = penaltyInput.get();
    }
    
    @Override
    public double calculateLogP() {
    	double [] rates = rateParameterInput.get().getDoubleValues();
    	logP = 0;
    	for (int i = 0; i < rates.length; i++) {
    		for (int j = i+1; j < rates.length; j++) {
    			double r = Math.max(rates[i], rates[j])/Math.min(rates[i], rates[j]);
    			if (r < threshold) {
    				logP -= penalty;
    			}
    		}
    	}
    	
    	return logP;
    }
    
	
	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

}
