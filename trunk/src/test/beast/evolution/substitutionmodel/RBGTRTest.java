package test.beast.evolution.substitutionmodel;


import org.junit.Test;

import beast.core.Description;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.RB;

import test.beast.evolution.substmodel.GTRTest;

@Description("Test GTR matrix exponentiation")
public class RBGTRTest extends GTRTest {

	Instance[] all = {test2, test1, test0};

	@Test
    public void testGTR() throws Exception {
        for (Instance test : all) {

            RealParameter f = new RealParameter(test.getPi());

            Frequencies freqs = new Frequencies();
            freqs.initByName("frequencies", f, "estimate", false);

            RB rbmodel = new RB();
            Double [] rates = test.getRates();
            Double [] rates2 = new Double[5];
            for (int i = 0; i < 4; i++) {
            	rates2[i] = rates[i] / rates[4];
            }
            rates2[4] = rates[5]/rates[4];
            rbmodel.initByName("count", new IntegerParameter("1"), 
            		"rates", new RealParameter(rates2), 
            		"frequencies", freqs);

            double distance = test.getDistance();

            double[] mat = new double[4 * 4];
            rbmodel.getTransitionProbabilities(null, distance, 0, 1, mat);
            final double[] result = test.getExpectedResult();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 1e-10);
                System.out.println(k + " : " + (mat[k] - result[k]));
            }
        }
    }
}