package rbbeast.evolution.sitemodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;


@Description("Site model representing a single category of a multi-catogry site model -- useful for VariableCategoryTreeLikelihood")
public class SingleCategorySiteModel extends SiteModelInterface.Base {
	public Input<SiteModelInterface> variableCategorySiteModelInput = new Input<>("variableCategorySiteModel", "the sitemodel that this model represents a single site for", Validate.REQUIRED);
	public Input<Integer> catecoryInput = new Input<>("category" , "category in the variableCategorySiteModel", 0);

	private SiteModel.Base siteModel;
	private int category;
	
	@Override
	public void initAndValidate() {
		this.category = catecoryInput.get();
		this.siteModel = (SiteModel.Base) variableCategorySiteModelInput.get();
		substModelInput.setValue(siteModel.substModelInput.get(), this);
	}

	@Override
	public void setDataType(DataType dataType) {
	}


	@Override
	public boolean integrateAcrossCategories() {
		return true;
	}

	@Override
	public int getCategoryCount() {
		return 1;
	}

	@Override
	public int getCategoryOfSite(int site, Node node) {
		return 0;
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		if (siteModel.getCategoryCount() > this.category) {
			return siteModel.getRateForCategory(this.category, node);
		} else {
			return 1.0;
		}
	}

	@Override
	public double[] getCategoryRates(Node node) {
		if (siteModel.getCategoryCount() > this.category) {
			return new double[] {siteModel.getRateForCategory(this.category, node)};
		} else {
			return new double[] {1.0};
		}
	}

	@Override
	public double getProportionForCategory(int category, Node node) {
		return 1;
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		return new double[] {1};
	}

}
