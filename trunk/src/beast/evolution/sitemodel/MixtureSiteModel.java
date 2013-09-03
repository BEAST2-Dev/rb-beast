package beast.evolution.sitemodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;




@Description("General Mixture Model from Stephane Guindon")
public class MixtureSiteModel extends DefaultMixtureSiteModel {
	public Input<List<MixtureSubstModel>> componentInput = new Input<List<MixtureSubstModel>>("component","one component of the mixture site model", new ArrayList<MixtureSubstModel>(), Validate.REQUIRED);
	
	List<MixtureSubstModel> components;
	double [] weights;
	boolean updateWeights;
	double [] rates;
	boolean updateRates;
	
	boolean [] isDirtySiteModel;
	
	@Override
	public void initAndValidate() throws Exception {
		components = componentInput.get();
		weights = new double[components.size()];
		updateWeights = true;
		rates = new double[components.size()];
		updateRates = true;
		isDirtySiteModel = new boolean[components.size()];
	}
	
	@Override
	public int getCategoryCount() {
		return components.size();
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		if (updateWeights) {
			updateWeights();
		}
		return weights;
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		return getCategoryProportions(node)[category];
	}
	
	
	@Override
	public double getProportionForCategory(int category, Node node) {
		if (updateWeights) {
			updateWeights();
		}
		return weights[category];
	}

	private void updateWeights() {
		int k = 0;
		for (MixtureSubstModel model : components) {
			weights[k++] = model.weight.getValue();
		}
		double sum = 0;
		for (double f : weights) {
			sum += f;
		}
		for (k = 0; k < weights.length; k++) {
			weights[k] /= sum;
		}
		updateWeights = false;
	}
			
	@Override
	public double[] getCategoryRates(Node node) {
		if (updateRates || updateWeights) {
			updateRates();
		}
		return rates;
	}

	private void updateRates() {
		if (updateWeights) {
			updateWeights();
		}
		int k = 0;
		for (MixtureSubstModel model : components) {
			rates[k++] = model.rate.getValue();
		}
		double sum = 0;
		k = 0;
		for (double f : weights) {
			sum += f * rates[k++];
		}
		for (k = 0; k < rates.length; k++) {
			rates[k] /= sum;
		}
		updateRates = true;
	}

	@Override
	public double[] getFrequencies(int iClass) {
		return components.get(iClass).substitutionModel.getFrequencies();
	}
	
	@Override
	public void getTransitionProbabilities(Node node, double height, double height2, double jointBranchRate,
			double[] m_fProbabilities, int iClass) {
		components.get(iClass).substitutionModel.getTransitionProbabilities(node, height, height2, jointBranchRate, m_fProbabilities);
	}
		
	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty = false;
		Arrays.fill(isDirtySiteModel, false);
		int k = 0;
		for (MixtureSubstModel model : components) {
			if (model.weight.somethingIsDirty()) {
				updateWeights = true;
			}
			if (model.rate.somethingIsDirty()) {
				updateRates = true;
			}
			if (model.substitutionModel.isDirtyCalculation()) {
				isDirty = true;
				isDirtySiteModel[k] = true;
			}
			k++;
		}
		return updateWeights || updateRates || isDirty;
	}
	
	@Override
	public boolean hasDirtySubstModel(int iClass) {
		return components.get(iClass).substitutionModel.isDirtyCalculation();
	}
	
//	@Override
//	public List<Integer> getDirtyCategories() {
//		List<Integer> classes = new ArrayList<Integer>();
//		for (int k = 0; k < components.size(); k++) {
//			if (components.get(k).substitutionModel.isDirtyCalculation()) {
//				classes.add(k);
//			}
//		}
//		return classes;
//	}
}
