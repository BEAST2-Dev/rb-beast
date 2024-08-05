package rbbeast.math.distributions;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.distribution.Prior;

@Description("Prior for Reversible-jump Based (RB) substitution model, applies prior only to the rates that are in use")
public class RBPrior extends Prior {
	public Input<Parameter> countInput = new Input<Parameter>("count","count parameter indicating the nr of rates to use", Validate.REQUIRED);

	private Parameter counts;
	
	@Override
	public void initAndValidate() {
		counts = countInput.get();
		super.initAndValidate();
	}
	
	@Override
	public double calculateLogP() {
		Function x = m_x.get();
		int dim = (int) counts.getArrayValue();
		double fOffset = dist.offsetInput.get();
		logP = 0;
		for (int i = 0; i < dim; i++) {
			double fX = x.getArrayValue(i) - fOffset;
			//fLogP += Math.log(density(fX));
			logP += dist.logDensity(fX);
		}
		return logP;
	}


}
