package rbbeast.evolution.sitemodel;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel.Base;

@Description("One component of a Mixture Site Model")
public class MixtureSubstModel extends CalculationNode {
	public Input<RealParameter> weightInput = new Input<RealParameter>("weight","weight of this proportion to site model", Validate.REQUIRED);
	public Input<RealParameter> rateInput = new Input<RealParameter>("rate", "rate for this proportion to site model", Validate.REQUIRED);
	public Input<SubstitutionModel.Base> substModelInput = new Input<SubstitutionModel.Base>("substModel","substitution model for this proportionof to site model", Validate.REQUIRED);
	
	
	RealParameter weight;
	RealParameter rate;
	SubstitutionModel.Base substitutionModel;
	
	@Override
	public void initAndValidate() {
		weight = weightInput.get();
		rate = rateInput.get();
		substitutionModel = substModelInput.get();
	}

}
