package beast.evolution.operators;


import org.apache.commons.math.MathException;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.distribution.Exponential;
import beast.base.inference.distribution.ParametricDistribution;
import beast.math.distributions.RBPrior;
import beast.base.util.Randomizer;



@Description("Operator for RB-SubstitutionModel to jump between states in the hierarchy")
public class RBOperator extends Operator {
	public Input<RealParameter> rateInput = new Input<RealParameter>("rates","rate parameter containing rates for hierarchical subst model", Validate.REQUIRED);
	public Input<Parameter> countInput = new Input<Parameter>("count","count parameter indicating the nr of rates to use", Validate.REQUIRED);
	
	private RealParameter rates;
	private Parameter counts;
//	private Double scaleFactor;
	//final static int [] countmap = {-1, -1, -1, 0, 2, 0}; 
	//final static int [] countmap = {-1, 0, 1, 2, 3, 4}; 
	//final static int [] countmap = {-1, 0, 0, 0, 0, 0}; 
	final static int [] countmap = {-1, -1, -1, -1, -1, -1}; 
	ParametricDistribution distr;
	
	@Override
	public void initAndValidate() {
		rates = rateInput.get();
		distr = new Exponential();
		for (Object plugin : rates.getOutputs()) {
			if (plugin instanceof RBPrior) {
				RBPrior prior = (RBPrior) plugin;
				distr = prior.distInput.get();
			}
		}
		counts = countInput.get();
	}

	@Override
	public double proposal() {
		int count = (int) counts.getArrayValue();
		double oldvalue = (countmap[count] >= 0 ? rates.getArrayValue(countmap[count]) : 1.0);
		double scale = 1;// = (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

		
		if (Randomizer.nextBoolean()) {
			// increase nr of rates
			if (count == rates.getDimension()	) {
				// cannot increase any further
				return Double.NEGATIVE_INFINITY;
			}
			double p = Randomizer.nextDouble();
			try {
				scale = distr.inverseCumulativeProbability(p);
			} catch (MathException e) {
				e.printStackTrace();
			}
			double newValue = oldvalue * scale;
			if (newValue < rates.getLower() || newValue > rates.getUpper()) {
				return Double.NEGATIVE_INFINITY;
			}
			rates.setValue(count, oldvalue * scale);
			counts.setValue(count + 1);
			//double logHR = -Math.log(1/scaleFactor - scaleFactor) + Math.log(oldvalue);
			double logHR = -distr.logDensity(scale) + Math.log(oldvalue);
			return logHR;
		} else {
			// decrease nr of rates
			if (count == 0) {
				// cannot decrease any further
				return Double.NEGATIVE_INFINITY;
			}
			counts.setValue(count - 1);
			scale = rates.getValue(count-1); 
			//double logHR = Math.log(1/scaleFactor - scaleFactor) - Math.log(oldvalue);
			double logHR = distr.logDensity(scale) - Math.log(oldvalue);
			return logHR;
		}
	}

}
