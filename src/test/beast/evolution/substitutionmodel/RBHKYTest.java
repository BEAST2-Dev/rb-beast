package test.beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.RB;
import test.beast.evolution.substmodel.HKYTest;
import static org.junit.jupiter.api.Assertions.assertEquals;

@Description("Test HKY matrix exponentiation")
public class RBHKYTest extends HKYTest {

	Instance[] all = {test2, test1, test0};

    public void testHKY() throws Exception {
        for (Instance test : all) {

            RealParameter f = new RealParameter(test.getPi());

            Frequencies freqs = new Frequencies();
            freqs.initByName("frequencies", f, "estimate", false);

            RB rbmodel = new RB();
            Double [] rates = new Double[5];
            rates[0] = 1.0/test.getKappa();
            rbmodel.initByName("count", new IntegerParameter("1"), 
            		"rates", new RealParameter(rates), 
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