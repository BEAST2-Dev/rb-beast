package beast.evolution.sitemodel;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;

@Description("Gamma site model -- for testing Mixture treelikelihood")
public class DefaultMixtureSiteModel extends SiteModel {

	public void getTransitionProbabilities(Node node, double height, double height2, double jointBranchRate,
			double[] m_fProbabilities, int iClass) {
		substModelInput.get().getTransitionProbabilities(node, height, height2, jointBranchRate, m_fProbabilities);
	}

	public double [] getFrequencies(int iClass) {
		return substModelInput.get().getFrequencies();
	}
	
	/** return whether the substitution model for category iClass is dirty **/ 
	public boolean hasDirtySubstModel(int iClass) {
		return ((CalculationNode) substModelInput.get()).isDirtyCalculation();
	}

//	/** return list of classes/categories that are dirty **/
//	public List<Integer> getDirtyCategories() {
//		return MixtureTreeLikelihood.ALL;
//	}
	
}
