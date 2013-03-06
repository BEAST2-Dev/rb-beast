package beast.evolution.likelihood;

import java.util.List;

public class MixtureBeerLikelihoodCore4 extends BeerLikelihoodCore4 implements MixtureLikelihoodCore {

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
			int v = 4 * l * m_nPatterns;

			for (int k = 0; k < m_nPatterns; k++) {

				int state1 = iStates1[k];
				int state2 = iStates2[k];

				int w = l * m_nMatrixSize;

				if (state1 < 4 && state2 < 4) {

					fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];
					v++;	w += 4;

				} else if (state1 < 4) {
					// child 2 has a gap or unknown state so don't use it

					fPartials3[v] = fMatrices1[w + state1];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1];
					v++;	w += 4;
					fPartials3[v] = fMatrices1[w + state1];
					v++;	w += 4;

				} else if (state2 < 4) {
					// child 2 has a gap or unknown state so don't use it
					fPartials3[v] = fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices2[w + state2];
					v++;	w += 4;
					fPartials3[v] = fMatrices2[w + state2];
					v++;	w += 4;

				} else {
					// both children have a gap or unknown state so set partials to 1
					fPartials3[v] = 1.0;
					v++;
					fPartials3[v] = 1.0;
					v++;
					fPartials3[v] = 1.0;
					v++;
					fPartials3[v] = 1.0;
					v++;
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

		double sum;//, tmp;


		for (int l : classes) {
			int v = 4 * l * m_nPatterns;
			int u = v;
			for (int k = 0; k < m_nPatterns; k++) {

				int state1 = iStates1[k];

                int w = l * m_nMatrixSize;

				if (state1 < 4) {


					sum =	fMatrices2[w] * fPartials2[v];
					sum +=	fMatrices2[w + 1] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 2] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 3] * fPartials2[v + 3];
					fPartials3[u] = fMatrices1[w + state1] * sum;	u++;

					sum =	fMatrices2[w + 4] * fPartials2[v];
					sum +=	fMatrices2[w + 5] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 6] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 7] * fPartials2[v + 3];
					fPartials3[u] = fMatrices1[w + 4 + state1] * sum;	u++;

					sum =	fMatrices2[w + 8] * fPartials2[v];
					sum +=	fMatrices2[w + 9] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 10] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 11] * fPartials2[v + 3];
					fPartials3[u] = fMatrices1[w + 8 + state1] * sum;	u++;

					sum =	fMatrices2[w + 12] * fPartials2[v];
					sum +=	fMatrices2[w + 13] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 14] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 15] * fPartials2[v + 3];
					fPartials3[u] = fMatrices1[w + 12 + state1] * sum;	u++;

					v += 4;

				} else {
					// Child 1 has a gap or unknown state so don't use it


					sum =	fMatrices2[w] * fPartials2[v];
					sum +=	fMatrices2[w + 1] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 2] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 3] * fPartials2[v + 3];
					fPartials3[u] = sum;	u++;

					sum =	fMatrices2[w + 4] * fPartials2[v];
					sum +=	fMatrices2[w + 5] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 6] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 7] * fPartials2[v + 3];
					fPartials3[u] = sum;	u++;

					sum =	fMatrices2[w + 8] * fPartials2[v];
					sum +=	fMatrices2[w + 9] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 10] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 11] * fPartials2[v + 3];
					fPartials3[u] = sum;	u++;

					sum =	fMatrices2[w + 12] * fPartials2[v];
					sum +=	fMatrices2[w + 13] * fPartials2[v + 1];
					sum +=	fMatrices2[w + 14] * fPartials2[v + 2];
					sum +=	fMatrices2[w + 15] * fPartials2[v + 3];
					fPartials3[u] = sum;	u++;

					v += 4;
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
			int v = 4 * l * m_nPatterns;
			int u = v;

			for (int k = 0; k < m_nPatterns; k++) {

                int w = l * m_nMatrixSize;

				sum1 = fMatrices1[w] * fPartials1[v];
				sum2 = fMatrices2[w] * fPartials2[v];
				sum1 += fMatrices1[w + 1] * fPartials1[v + 1];
				sum2 += fMatrices2[w + 1] * fPartials2[v + 1];
				sum1 += fMatrices1[w + 2] * fPartials1[v + 2];
				sum2 += fMatrices2[w + 2] * fPartials2[v + 2];
				sum1 += fMatrices1[w + 3] * fPartials1[v + 3];
				sum2 += fMatrices2[w + 3] * fPartials2[v + 3];
				fPartials3[u] = sum1 * sum2; u++;

				sum1 = fMatrices1[w + 4] * fPartials1[v];
				sum2 = fMatrices2[w + 4] * fPartials2[v];
				sum1 += fMatrices1[w + 5] * fPartials1[v + 1];
				sum2 += fMatrices2[w + 5] * fPartials2[v + 1];
				sum1 += fMatrices1[w + 6] * fPartials1[v + 2];
				sum2 += fMatrices2[w + 6] * fPartials2[v + 2];
				sum1 += fMatrices1[w + 7] * fPartials1[v + 3];
				sum2 += fMatrices2[w + 7] * fPartials2[v + 3];
				fPartials3[u] = sum1 * sum2; u++;

				sum1 = fMatrices1[w + 8] * fPartials1[v];
				sum2 = fMatrices2[w + 8] * fPartials2[v];
				sum1 += fMatrices1[w + 9] * fPartials1[v + 1];
				sum2 += fMatrices2[w + 9] * fPartials2[v + 1];
				sum1 += fMatrices1[w + 10] * fPartials1[v + 2];
				sum2 += fMatrices2[w + 10] * fPartials2[v + 2];
				sum1 += fMatrices1[w + 11] * fPartials1[v + 3];
				sum2 += fMatrices2[w + 11] * fPartials2[v + 3];
				fPartials3[u] = sum1 * sum2; u++;

				sum1 = fMatrices1[w + 12] * fPartials1[v];
				sum2 = fMatrices2[w + 12] * fPartials2[v];
				sum1 += fMatrices1[w + 13] * fPartials1[v + 1];
				sum2 += fMatrices2[w + 13] * fPartials2[v + 1];
				sum1 += fMatrices1[w + 14] * fPartials1[v + 2];
				sum2 += fMatrices2[w + 14] * fPartials2[v + 2];
				sum1 += fMatrices1[w + 15] * fPartials1[v + 3];
				sum2 += fMatrices2[w + 15] * fPartials2[v + 3];
				fPartials3[u] = sum1 * sum2; u++;

				v += 4;
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
