package beast.evolution.operators;


import org.apache.commons.math.MathException;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.Exponential;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.RBPrior;
import beast.util.Randomizer;



@Description("Operator for RB-SubstitutionModel to jump between states in the hierarchy")
public class RBOperator extends Operator {
	public Input<RealParameter> rateInput = new Input<RealParameter>("rates","rate parameter containing rates for hierarchical subst model", Validate.REQUIRED);
	public Input<IntegerParameter> countInput = new Input<IntegerParameter>("count","count parameter indicating the nr of rates to use", Validate.REQUIRED);
	
	private RealParameter rates;
	private IntegerParameter counts;
//	private Double scaleFactor;
	//final static int [] countmap = {-1, -1, -1, 0, 2, 0}; 
	//final static int [] countmap = {-1, 0, 1, 2, 3, 4}; 
	//final static int [] countmap = {-1, 0, 0, 0, 0, 0}; 
	final static int [] countmap = {-1, -1, -1, -1, -1, -1}; 
	ParametricDistribution distr;
	
	@Override
	public void initAndValidate() throws Exception {
		rates = rateInput.get();
		distr = new Exponential();
		for (BEASTObject plugin : rates.outputs) {
			if (plugin instanceof RBPrior) {
				RBPrior prior = (RBPrior) plugin;
				distr = prior.distInput.get();
			}
		}
		counts = countInput.get();
	}

	@Override
	public double proposal() {
		int count = counts.getValue();
		double oldvalue = (countmap[count] >= 0 ? rates.getArrayValue(countmap[count]) : 1.0);
		double scale = 1;// = (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

		
		if (Randomizer.nextBoolean()) {
			// increase nr of rates
			if (count == 5) {
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
