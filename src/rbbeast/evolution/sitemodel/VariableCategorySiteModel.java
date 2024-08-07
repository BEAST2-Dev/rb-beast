package rbbeast.evolution.sitemodel;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.sitemodel.SiteModelInterface.Base;
import beast.base.evolution.tree.Node;

public class VariableCategorySiteModel extends Base {
	final public Input<Integer> maxCategoryCountInput = new Input<>("maxCategoryCount", "maximum number of categories allowed", Validate.REQUIRED); 

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean integrateAcrossCategories() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public int getCategoryCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getCategoryOfSite(int site, Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getCategoryRates(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getProportionForCategory(int category, Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

}
