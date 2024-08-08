package rbbeast.evolution.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

@Description("Returns the dimension of a parameter")
public class GetDimension extends CalculationNode implements Function {
	final public Input<RealParameter> parameterInput = new Input<>("parameter", "parameter for which to get the dimension", Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		return parameterInput.get().getDimension();
	}

}
