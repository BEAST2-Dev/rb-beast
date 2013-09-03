package beast.evolution.likelihood;

import java.util.Arrays;
import java.util.List;

import beast.evolution.likelihood.BeerLikelihoodCore;


public class MixtureBeerLikelihoodCore extends BeerLikelihoodCore implements MixtureLikelihoodCore {

	public MixtureBeerLikelihoodCore(int nStateCount) {
		super(nStateCount);
	}
	
	@Override
	public void integratePartialsMixture(int iNodeIndex, double[] fProportions, double[] fOutPartials,
			double[][] frequencies, double[] fPatternLogLikelihoods) {

		double[] fInPartials = partials[currentPartialsIndex[iNodeIndex]][iNodeIndex];
        int u = 0;
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {
        	double [] freqs = frequencies[0];
            for (int i = 0; i < nrOfStates; i++) {
                fOutPartials[u] = fInPartials[v] * fProportions[0] * freqs[i];
                u++;
                v++;
            }
        }

        for (int l = 1; l < nrOfMatrices; l++) {
        	double [] freqs = frequencies[l];
            u = 0;
            for (int k = 0; k < nrOfPatterns; k++) {
                for (int i = 0; i < nrOfStates; i++) {
                    fOutPartials[u] += fInPartials[v] * fProportions[l] * freqs[i];
                    u++;
                    v++;
                }
            }
        }
	
        v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {
            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {
                sum += fOutPartials[v];
                v++;
            }
            fPatternLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
	} // integratePartialsMixture
	
	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	@Override
	public void calculateStatesStatesPruning(int[] iStates1, double[] fMatrices1,
												int[] iStates2, double[] fMatrices2,
												double[] fPartials3, List<Integer> classes)
	{

		for (int l : classes) {
			int v = nrOfStates * l * nrOfPatterns;

			for (int k = 0; k < nrOfPatterns; k++) {

				int state1 = iStates1[k];
				int state2 = iStates2[k];

				int w = l * matrixSize;

                if (state1 < nrOfStates && state2 < nrOfStates) {

					for (int i = 0; i < nrOfStates; i++) {

						fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];

						v++;
						w += nrOfStates;
					}

				} else if (state1 < nrOfStates) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < nrOfStates; i++) {

						fPartials3[v] = fMatrices1[w + state1];

						v++;
						w += nrOfStates;
					}
				} else if (state2 < nrOfStates) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < nrOfStates; i++) {

						fPartials3[v] = fMatrices2[w + state2];

						v++;
						w += nrOfStates;
					}
				} else {
					// both children have a gap or unknown state so set partials to 1

					for (int j = 0; j < nrOfStates; j++) {
						fPartials3[v] = 1.0;
						v++;
					}
				}
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	@Override
	public  void calculateStatesPartialsPruning(	int[] iStates1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3, List<Integer> classes)
	{

		double sum, tmp;


		for (int l : classes) {
			int v = nrOfStates * l * nrOfPatterns;
			int u = v;
			for (int k = 0; k < nrOfPatterns; k++) {

				int state1 = iStates1[k];

                int w = l * matrixSize;

				if (state1 < nrOfStates) {


					for (int i = 0; i < nrOfStates; i++) {

						tmp = fMatrices1[w + state1];

						sum = 0.0;
						for (int j = 0; j < nrOfStates; j++) {
							sum += fMatrices2[w] * fPartials2[v + j];
							w++;
						}

						fPartials3[u] = tmp * sum;
						u++;
					}

					v += nrOfStates;
				} else {
					// Child 1 has a gap or unknown state so don't use it

					for (int i = 0; i < nrOfStates; i++) {

						sum = 0.0;
						for (int j = 0; j < nrOfStates; j++) {
							sum += fMatrices2[w] * fPartials2[v + j];
							w++;
						}

						fPartials3[u] = sum;
						u++;
					}

					v += nrOfStates;
				}
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	@Override
	public void calculatePartialsPartialsPruning(double[] fPartials1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3, List<Integer> classes)
	{
		double sum1, sum2;

		for (int l : classes) {
			int v = nrOfStates * l * nrOfPatterns;
			int u = v;

			for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

				for (int i = 0; i < nrOfStates; i++) {

					sum1 = sum2 = 0.0;

					for (int j = 0; j < nrOfStates; j++) {
						sum1 += fMatrices1[w] * fPartials1[v + j];
						sum2 += fMatrices2[w] * fPartials2[v + j];

						w++;
					}

					fPartials3[u] = sum1 * sum2;
					u++;
				}
				v += nrOfStates;
			}
		}
	}

	@Override
    public void calculatePartials(int iNodeIndex1, int iNodeIndex2, int iNodeIndex3, List<Integer> classes) {
        currentPartialsIndex[iNodeIndex3] = 1 - storedPartialsIndex[iNodeIndex3];
        
        System.arraycopy(partials[1-currentPartialsIndex[iNodeIndex3]][iNodeIndex3], 0,
        		partials[currentPartialsIndex[iNodeIndex3]][iNodeIndex3], 0, partials[1-currentPartialsIndex[iNodeIndex3]][iNodeIndex3].length);
        
        if (states[iNodeIndex1] != null) {
            if (states[iNodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[iNodeIndex1], matrices[currentMatrixIndex[iNodeIndex1]][iNodeIndex1],
                        states[iNodeIndex2], matrices[currentMatrixIndex[iNodeIndex2]][iNodeIndex2],
                        partials[currentPartialsIndex[iNodeIndex3]][iNodeIndex3], classes);
            } else {
                calculateStatesPartialsPruning(states[iNodeIndex1], matrices[currentMatrixIndex[iNodeIndex1]][iNodeIndex1],
                        partials[currentPartialsIndex[iNodeIndex2]][iNodeIndex2], matrices[currentMatrixIndex[iNodeIndex2]][iNodeIndex2],
                        partials[currentPartialsIndex[iNodeIndex3]][iNodeIndex3], classes);
            }
        } else {
            if (states[iNodeIndex2] != null) {
                calculateStatesPartialsPruning(states[iNodeIndex2], matrices[currentMatrixIndex[iNodeIndex2]][iNodeIndex2],
                        partials[currentPartialsIndex[iNodeIndex1]][iNodeIndex1], matrices[currentMatrixIndex[iNodeIndex1]][iNodeIndex1],
                        partials[currentPartialsIndex[iNodeIndex3]][iNodeIndex3], classes);
            } else {
                calculatePartialsPartialsPruning(partials[currentPartialsIndex[iNodeIndex1]][iNodeIndex1], matrices[currentMatrixIndex[iNodeIndex1]][iNodeIndex1],
                        partials[currentPartialsIndex[iNodeIndex2]][iNodeIndex2], matrices[currentMatrixIndex[iNodeIndex2]][iNodeIndex2],
                        partials[currentPartialsIndex[iNodeIndex3]][iNodeIndex3], classes);
            }
        }

        if (useScaling) {
            scalePartials(iNodeIndex3);
        }
    }

}
