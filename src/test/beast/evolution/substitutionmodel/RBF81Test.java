package test.beast.evolution.substitutionmodel;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.Frequencies;
import rbbeast.evolution.substitutionmodel.RB;
import junit.framework.TestCase;

/**
 * Test F81 matrix exponentiation
 *
 * @author Walter Xie
 */
public class RBF81Test extends TestCase {

    interface Instance {
        Double[] getPi();

        double getDistance();

        double[] getExpectedResult();
    }

    /*
     * Results obtained by running the following scilab code,

      piQ = diag([.2, .3, .25, .25]) ; d = 0.05 ;

      % Q matrix with zeroed diagonal
      XQ = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0];

      xx = XQ * piQ ;

      % fill diagonal and normalize by total substitution rate
      q0 = (xx + diag(-sum(xx,2))) / sum(piQ * sum(xx,2)) ;
      format('v', 15);
      expm(q0 * d)

     */
    Instance test0 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.25, 0.25, 0.25, 0.25};
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.906379989282, 0.031206670239, 0.031206670239, 0.031206670239,
                    0.031206670239, 0.906379989282, 0.031206670239, 0.031206670239,
                    0.031206670239, 0.031206670239, 0.906379989282, 0.031206670239,
                    0.031206670239, 0.031206670239, 0.031206670239, 0.906379989282
            };
        }
    };

    Instance test1 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.50, 0.20, 0.2, 0.1};
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.929702430444, 0.028119027822, 0.028119027822, 0.014059513911,
                    0.070297569556, 0.887523888711, 0.028119027822, 0.014059513911,
                    0.070297569556, 0.028119027822, 0.887523888711, 0.014059513911,
                    0.070297569556, 0.028119027822, 0.028119027822, 0.873464374800
            };
        }
    };

    Instance test2 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.20, 0.30, 0.25, 0.25};
        }

        public double getDistance() {
            return 0.05;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.948070805840, 0.019473447810, 0.016227873175, 0.016227873175,
                    0.012982298540, 0.954561955110, 0.016227873175, 0.016227873175,
                    0.012982298540, 0.019473447810, 0.951316380475, 0.016227873175,
                    0.012982298540, 0.019473447810, 0.016227873175, 0.951316380475
            };
        }
    };

    Instance test3 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.20, 0.30, 0.25, 0.25};
        }

        public double getDistance() {
            return 0.2;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.811647020254, 0.070632367405, 0.058860306171, 0.058860306171,
                    0.047088244936, 0.835191142722, 0.058860306171, 0.058860306171,
                    0.047088244936, 0.070632367405, 0.823419081488, 0.058860306171,
                    0.047088244936, 0.070632367405, 0.058860306171, 0.823419081488
            };
        }
    };

    Instance[] all = {test0, test1, test2, test3};

    public void testGF81() throws Exception {
        for (Instance test : all) {
            RealParameter f = new RealParameter(test.getPi());
            Frequencies freqs = new Frequencies();
            freqs.initByName("frequencies", f, "estimate", false);

            RB rbmodel = new RB();
            Double [] rates = new Double[5];
            rbmodel.initByName("count", new IntegerParameter("0"), 
            		"rates", new RealParameter(rates), 
            		"frequencies", freqs);

            double distance = test.getDistance();

            double[] mat = new double[4 * 4];
            rbmodel.getTransitionProbabilities(null, distance, 0, 1, mat);
            final double[] result = test.getExpectedResult();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 1e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }
        }
    }
}
