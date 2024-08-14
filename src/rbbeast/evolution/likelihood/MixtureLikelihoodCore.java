package rbbeast.evolution.likelihood;

import java.util.List;

public interface MixtureLikelihoodCore {
	public abstract void integratePartialsMixture(int nr, double[] proportions, double[] m_fRootPartials,
			double[][] frequencies, double[] m_fPatternLogLikelihoods);

	public void calculateStatesPartialsPruning(	int[] iStates1, double[] fMatrices1,
			double[] fPartials2, double[] fMatrices2,
			double[] fPartials3, List<Integer> classes);

	public void calculateStatesStatesPruning(int[] iStates1, double[] fMatrices1,
			int[] iStates2, double[] fMatrices2,
			double[] fPartials3, List<Integer> classes);

	public void calculatePartialsPartialsPruning(double[] fPartials1, double[] fMatrices1,
			double[] fPartials2, double[] fMatrices2,
			double[] fPartials3, List<Integer> classes);

    public void calculatePartials(int iNodeIndex1, int iNodeIndex2, int iNodeIndex3, List<Integer> classes);

	public abstract void integratePartialsMixture(int nr, double[] proportions, double[] m_fRootPartials,
			double[] frequencies, double[] m_fPatternLogLikelihoods, int classCount);

	public void calculateStatesPartialsPruning(	int[] iStates1, double[] fMatrices1,
			double[] fPartials2, double[] fMatrices2,
			double[] fPartials3, int classCount);

	public void calculateStatesStatesPruning(int[] iStates1, double[] fMatrices1,
			int[] iStates2, double[] fMatrices2,
			double[] fPartials3, int classCount);

	public void calculatePartialsPartialsPruning(double[] fPartials1, double[] fMatrices1,
			double[] fPartials2, double[] fMatrices2,
			double[] fPartials3, int classCount);

    public void calculatePartials(int iNodeIndex1, int iNodeIndex2, int iNodeIndex3, int classCount);
}
