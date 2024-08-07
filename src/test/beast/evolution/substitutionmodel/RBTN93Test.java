package test.beast.evolution.substitutionmodel;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.Frequencies;
import rbbeast.evolution.substitutionmodel.RB;
import junit.framework.TestCase;

/**
 * Test TN93 matrix exponentiation
 *
 * @author Joseph Heled
 *         Date: 7/11/2007
 */
public class RBTN93Test extends TestCase {

    interface Instance {
        Double[] getPi();
        double   getKappa1();
        double   getKappa2();
        double   getDistance();
        double[] getExpectedResult();
    }

    /*
     * Results obtained by running the following scilab code,
     *
     * ACGT
     * 
     * k1 = 2 ; k2 = 3; piQ = diag([.25, .25, .25, .25]) ; d = 0.1 ;
     * % Q matrix with zeroed diagonal
     * XQ = [0  1  k1    1; 1 0 1 k2 ; k1 1 0 1 ; 1 k2 1 0]
     *
     * xx = XQ * piQ ;
     *
     * % fill diagonal and normalize by total substitution rate
     * q0 = (xx + diag(-sum(xx,2))) / sum(piQ * sum(xx,2)) ;
     * expm(q0 * d)
     */
    Instance test0 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.25, 0.25, 0.25,0.25};
        }

        public double getKappa1() {
            return 2;
        }

        public double getKappa2() {
            return 3;
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.916323466703981  ,  0.021263192817492 ,   0.041150147661034 ,   0.021263192817492  ,
                    0.021263192817492  ,  0.897301022862890 ,   0.021263192817493 ,   0.060172591502126  ,
                    0.041150147661034  ,  0.021263192817492 ,   0.916323466703982 ,   0.021263192817492  ,
                    0.021263192817492  ,  0.060172591502126 ,   0.021263192817492 ,   0.897301022862889  ,

            };
        }
    };

    Instance test1 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.1, 0.2, 0.3, 0.4};
        }

        public double getKappa1() {
            return 2;
        }

        public double getKappa2() {
            return 3;
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.895550254199242, 0.017687039418335, 0.051388627545752, 0.035374078836670,
                    0.008843519709168, 0.865344657365451, 0.026530559127503, 0.099281263797879,
                    0.017129542515251, 0.017687039418335, 0.929809339229744, 0.035374078836670,
                    0.008843519709168, 0.049640631898940, 0.026530559127503, 0.914985289264390,
            };
        }
    };

     Instance test2 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.1, 0.2, 0.3, 0.4};
        }

        public double getKappa1() {
            return 1;
        }

        public double getKappa2() {
            return 3;
        }

        public double getDistance() {
            return 0.1;
        }

         public double[] getExpectedResult() {
             return new double[]{

                     0.915952014632886, 0.018677330081581, 0.028015995122371, 0.037354660163162,
                     0.009338665040790, 0.858207194100208, 0.028015995122371, 0.104438145736630,
                     0.009338665040790, 0.018677330081581, 0.934629344714466, 0.037354660163162,
                     0.009338665040790, 0.052219072868315, 0.028015995122371, 0.910426266968523,

             };
         }
    };

    Instance[] all = {test2, test1, test0};

    public void testTN93() throws Exception {
        for( Instance test : all ) {
            Double[] pi = test.getPi();

            RealParameter f = new RealParameter(pi);
            Frequencies freqs = new Frequencies();
            freqs.initByName("frequencies", f, "estimate", false);
            Double [] rates = new Double[5];
            rates[0] = 1.0/test.getKappa2();
            rates[1] = test.getKappa1()/test.getKappa2();
            RB rbmodel = new RB();
            rbmodel.initByName("count", new IntegerParameter("2"), 
            		"rates", new RealParameter(rates), 
            		"frequencies", freqs);

            double distance = test.getDistance();

            double[] mat = new double[4*4];
            rbmodel.getTransitionProbabilities(null, distance, 0, 1, mat);
            final double[] result = test.getExpectedResult();

            for(int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 1e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }
        }
    }
}