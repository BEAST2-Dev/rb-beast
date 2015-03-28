package test.beast.evolution.sitemodel;

import org.junit.Test;

import beast.evolution.sitemodel.BEASTModelTest;
import beast.evolution.substitutionmodel.JukesCantor;
import junit.framework.TestCase;

public class BEASTModelTestTest extends TestCase {

	@Test
	public void testCategories() throws Exception {
		BEASTModelTest sitemodel = new BEASTModelTest();
		sitemodel.initByName("gammaCategoryCount", 4,
				"shape", "1.0",
				"proportionInvariant", "0.1",
				"hasInvariantSites", "false",
				"hasGammaRates", "false",
				"substModel", new JukesCantor());
		assertEquals(5, sitemodel.getCategoryCount());
		testRates(sitemodel, new double[]{0.0, 1.0, 0.0, 0.0, 0.0});
	}

	@Test
	public void testCategories2() throws Exception {
		BEASTModelTest sitemodel = new BEASTModelTest();
		sitemodel.initByName("gammaCategoryCount", 4,
				"shape", "1.0",
				"proportionInvariant", "0.1",
				"hasInvariantSites", "true",
				"hasGammaRates", "false",
				"substModel", new JukesCantor());
		assertEquals(5, sitemodel.getCategoryCount());
		testRates(sitemodel, new double[]{0.1, 0.9, 0.0, 0.0, 0.0});
	}

	@Test
	public void testCategories3() throws Exception {
		BEASTModelTest sitemodel = new BEASTModelTest();
		sitemodel.initByName("gammaCategoryCount", 4,
				"shape", "1.0",
				"proportionInvariant", "0.1",
				"hasInvariantSites", "false",
				"hasGammaRates", "true",
				"substModel", new JukesCantor());
		assertEquals(5, sitemodel.getCategoryCount());
		testRates(sitemodel, new double[]{0.0, 0.25, 0.25, 0.25, 0.25});
	}

	@Test
	public void testCategories4() throws Exception {
		BEASTModelTest sitemodel = new BEASTModelTest();
		sitemodel.initByName("gammaCategoryCount", 4,
				"shape", "1.0",
				"proportionInvariant", "0.6",
				"hasInvariantSites", "true",
				"hasGammaRates", "true",
				"substModel", new JukesCantor());
		assertEquals(5, sitemodel.getCategoryCount());
		testRates(sitemodel, new double[]{0.6, 0.1, 0.1, 0.1, 0.1});
	}
	
	@Test
	public void testCategories4b() throws Exception {
		BEASTModelTest sitemodel = new BEASTModelTest();
		sitemodel.initByName("gammaCategoryCount", 4,
				"shape", "1.0",
				"proportionInvariant", "0.6",
				"hasInvariantSites", "true",
				"hasGammaRates", "true",
				"substModel", new JukesCantor());
		sitemodel.setPropInvariantIsCategory(false);
		assertEquals(4, sitemodel.getCategoryCount());
		testRates(sitemodel, new double[]{0.1, 0.1, 0.1, 0.1});
		assertEquals(sitemodel.getProportionInvariant(), 0.6, 1e-13);
	}
	
	void testRates(BEASTModelTest sitemodel, double [] expectedprops) {
		double [] props = sitemodel.getCategoryProportions(null);
		double [] rates = sitemodel.getCategoryRates(null);
		for (int i = 0; i < expectedprops.length; i++) {
			assertEquals(expectedprops[i], props[i], 1e-13);
		}
		
		if (sitemodel.hasPropInvariantCategory) {
			assertEquals(0.0, rates[0]);
		}
		
		for (int i = (sitemodel.hasPropInvariantCategory? 1 : 0); i < expectedprops.length; i++) {
			assertEquals(expectedprops[i] < 1e-13, rates[i] < 1e-13);
		}
	}
	

}
