package beast.evolution.likelihood;

import java.util.Arrays;
import java.util.List;

public class MixtureBeerLikelihoodCore extends BeerLikelihoodCore implements MixtureLikelihoodCore {

	public MixtureBeerLikelihoodCore(int nStateCount) {
		super(nStateCount);
	}
	
	@Override
	public void integratePartialsMixture(int iNodeIndex, double[] fProportions, double[] fOutPartials,
			double[][] frequencies, double[] fPatternLogLikelihoods) {

		double[] fInPartials = m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex];
        int u = 0;
        int v = 0;
        for (int k = 0; k < m_nPatterns; k++) {
        	double [] freqs = frequencies[0];
            for (int i = 0; i < m_nStates; i++) {
                fOutPartials[u] = fInPartials[v] * fProportions[0] * freqs[i];
                u++;
                v++;
            }
        }

        for (int l = 1; l < m_nMatrices; l++) {
        	double [] freqs = frequencies[l];
            u = 0;
            for (int k = 0; k < m_nPatterns; k++) {
                for (int i = 0; i < m_nStates; i++) {
                    fOutPartials[u] += fInPartials[v] * fProportions[l] * freqs[i];
                    u++;
                    v++;
                }
            }
        }
	
        v = 0;
        for (int k = 0; k < m_nPatterns; k++) {
            double sum = 0.0;
            for (int i = 0; i < m_nStates; i++) {
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
			int v = m_nStates * l * m_nPatterns;

			for (int k = 0; k < m_nPatterns; k++) {

				int state1 = iStates1[k];
				int state2 = iStates2[k];

				int w = l * m_nMatrixSize;

                if (state1 < m_nStates && state2 < m_nStates) {

					for (int i = 0; i < m_nStates; i++) {

						fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];

						v++;
						w += m_nStates;
					}

				} else if (state1 < m_nStates) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < m_nStates; i++) {

						fPartials3[v] = fMatrices1[w + state1];

						v++;
						w += m_nStates;
					}
				} else if (state2 < m_nStates) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < m_nStates; i++) {

						fPartials3[v] = fMatrices2[w + state2];

						v++;
						w += m_nStates;
					}
				} else {
					// both children have a gap or unknown state so set partials to 1

					for (int j = 0; j < m_nStates; j++) {
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
			int v = m_nStates * l * m_nPatterns;
			int u = v;
			for (int k = 0; k < m_nPatterns; k++) {

				int state1 = iStates1[k];

                int w = l * m_nMatrixSize;

				if (state1 < m_nStates) {


					for (int i = 0; i < m_nStates; i++) {

						tmp = fMatrices1[w + state1];

						sum = 0.0;
						for (int j = 0; j < m_nStates; j++) {
							sum += fMatrices2[w] * fPartials2[v + j];
							w++;
						}

						fPartials3[u] = tmp * sum;
						u++;
					}

					v += m_nStates;
				} else {
					// Child 1 has a gap or unknown state so don't use it

					for (int i = 0; i < m_nStates; i++) {

						sum = 0.0;
						for (int j = 0; j < m_nStates; j++) {
							sum += fMatrices2[w] * fPartials2[v + j];
							w++;
						}

						fPartials3[u] = sum;
						u++;
					}

					v += m_nStates;
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
			int v = m_nStates * l * m_nPatterns;
			int u = v;

			for (int k = 0; k < m_nPatterns; k++) {

                int w = l * m_nMatrixSize;

				for (int i = 0; i < m_nStates; i++) {

					sum1 = sum2 = 0.0;

					for (int j = 0; j < m_nStates; j++) {
						sum1 += fMatrices1[w] * fPartials1[v + j];
						sum2 += fMatrices2[w] * fPartials2[v + j];

						w++;
					}

					fPartials3[u] = sum1 * sum2;
					u++;
				}
				v += m_nStates;
			}
		}
	}

	@Override
    public void calculatePartials(int iNodeIndex1, int iNodeIndex2, int iNodeIndex3, List<Integer> classes) {
        m_iCurrentPartials[iNodeIndex3] = 1 - m_iStoredPartials[iNodeIndex3];
        
        System.arraycopy(m_fPartials[1-m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], 0,
        		m_fPartials[m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], 0, m_fPartials[1-m_iCurrentPartials[iNodeIndex3]][iNodeIndex3].length);
        
        if (m_iStates[iNodeIndex1] != null) {
            if (m_iStates[iNodeIndex2] != null) {
                calculateStatesStatesPruning(
                        m_iStates[iNodeIndex1], m_fMatrices[m_iCurrentMatrices[iNodeIndex1]][iNodeIndex1],
                        m_iStates[iNodeIndex2], m_fMatrices[m_iCurrentMatrices[iNodeIndex2]][iNodeIndex2],
                        m_fPartials[m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], classes);
            } else {
                calculateStatesPartialsPruning(m_iStates[iNodeIndex1], m_fMatrices[m_iCurrentMatrices[iNodeIndex1]][iNodeIndex1],
                        m_fPartials[m_iCurrentPartials[iNodeIndex2]][iNodeIndex2], m_fMatrices[m_iCurrentMatrices[iNodeIndex2]][iNodeIndex2],
                        m_fPartials[m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], classes);
            }
        } else {
            if (m_iStates[iNodeIndex2] != null) {
                calculateStatesPartialsPruning(m_iStates[iNodeIndex2], m_fMatrices[m_iCurrentMatrices[iNodeIndex2]][iNodeIndex2],
                        m_fPartials[m_iCurrentPartials[iNodeIndex1]][iNodeIndex1], m_fMatrices[m_iCurrentMatrices[iNodeIndex1]][iNodeIndex1],
                        m_fPartials[m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], classes);
            } else {
                calculatePartialsPartialsPruning(m_fPartials[m_iCurrentPartials[iNodeIndex1]][iNodeIndex1], m_fMatrices[m_iCurrentMatrices[iNodeIndex1]][iNodeIndex1],
                        m_fPartials[m_iCurrentPartials[iNodeIndex2]][iNodeIndex2], m_fMatrices[m_iCurrentMatrices[iNodeIndex2]][iNodeIndex2],
                        m_fPartials[m_iCurrentPartials[iNodeIndex3]][iNodeIndex3], classes);
            }
        }

        if (m_bUseScaling) {
            scalePartials(iNodeIndex3);
        }
    }

}
