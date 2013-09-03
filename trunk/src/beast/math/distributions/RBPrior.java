package beast.math.distributions;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.math.distributions.Prior;

@Description("Prior for Reversible-jump Based (RB) substitution model, applies prior only to the rates that are in use")
public class RBPrior extends Prior {
	public Input<IntegerParameter> countInput = new Input<IntegerParameter>("count","count parameter indicating the nr of rates to use", Validate.REQUIRED);

	private IntegerParameter counts;
	
	@Override
	public void initAndValidate() throws Exception {
		counts = countInput.get();
		super.initAndValidate();
	}
	
	@Override
	public double calculateLogP() throws Exception {
		Function x = m_x.get();
		int dim = counts.getValue();
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
