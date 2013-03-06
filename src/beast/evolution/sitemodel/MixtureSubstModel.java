package beast.evolution.sitemodel;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.substitutionmodel.SubstitutionModel.Base;

@Description("One component of a Mixture Site Model")
public class MixtureSubstModel extends CalculationNode {
	public Input<RealParameter> weightInput = new Input<RealParameter>("weight","weight of this proportion to site model", Validate.REQUIRED);
	public Input<RealParameter> rateInput = new Input<RealParameter>("rate", "rate for this proportion to site model", Validate.REQUIRED);
	public Input<SubstitutionModel.Base> substModelInput = new Input<SubstitutionModel.Base>("substModel","substitution model for this proportionof to site model", Validate.REQUIRED);
	
	
	RealParameter weight;
	RealParameter rate;
	SubstitutionModel.Base substitutionModel;
	
	@Override
	public void initAndValidate() throws Exception {
		weight = weightInput.get();
		rate = rateInput.get();
		substitutionModel = substModelInput.get();
	}

}
