package rbbeast.evolution.sitemodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;


@Description("Site model representing a single category of a multi-catogry site model -- useful for VariableCategoryTreeLikelihood")
public class SingleCategorySiteModel extends SiteModelInterface.Base {

	private SiteModel siteModel;
	private int category;
	
	public SingleCategorySiteModel(int category, SiteModelInterface siteModelInterface) {
		this.category = category;
		this.siteModel = (SiteModel) siteModelInterface;
	}

	@Override
	public void initAndValidate() {
	}

	@Override
	public void setDataType(DataType dataType) {
	}


	@Override
	public boolean integrateAcrossCategories() {
		return false;
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
		return siteModel.getRateForCategory(this.category, node);
	}

	@Override
	public double[] getCategoryRates(Node node) {
		return new double[] {siteModel.getRateForCategory(this.category, node)};
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
