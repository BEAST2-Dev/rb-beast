package rbbeast.evolution.util;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beastfx.app.seqgen.SimulatedAlignment;


@Description("Generates sequences with rate heterogeneity drawing rates from a *continuous* gamma distribution")
public class ContinuousGammaSimulatedAlignment extends SimulatedAlignment {
	final public Input<RealParameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of gamma distribution.", Validate.REQUIRED);
	
	private int stateCount;	
	private int [][] seq;
	
	@Override
	public void initAndValidate() {
		stateCount = m_data.get().getDataType().getStateCount();
		super.initAndValidate();
	}
	
	@Override
    public void simulate() {
        Node root = m_tree.getRoot();
        seq = new int[m_tree.getLeafNodeCount()][m_sequenceLength];
        
        // Sample a double from a Gamma distribution with a mean of
        // alpha/lambda and a variance of alpha/lambda^2.
        // alpha = shape, lambda = rate
        double alpha = shapeParameterInput.get().getArrayValue();
        double lambda = alpha;

        double[] rate = new double[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            rate[i] = Randomizer.nextGamma(alpha, lambda);
        }

        double[] frequencies = m_siteModel.getSubstitutionModel().getFrequencies();

        for (int site = 0; site < m_sequenceLength; site++) {
            int seq = Randomizer.randomChoicePDF(frequencies);
        	traverse(root, seq, rate, site);
        }

        for (int i = 0; i < m_tree.getLeafNodeCount(); i++) {
        	sequenceInput.setValue(intArray2ASequence(seq[i], m_tree.getNode(i)), this);
        }

    } // simulate

	private void traverse(Node node, int parentState, double[] rate, int site) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());

            getTransitionProbabilities(m_tree, child, rate[site], m_probabilities[0]);
            double[] cProb = new double[stateCount];
            System.arraycopy(m_probabilities[0], parentState * stateCount, cProb, 0, stateCount);
            int state = Randomizer.randomChoicePDF(cProb);

            if (child.isLeaf()) {
            	this.seq[child.getNr()][site] = state;
            } else {
                traverse(child, state, rate, site);
            }
        }
    } // traverse
	
	
    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    private Sequence intArray2ASequence(int[] seq, Node node) {
        DataType dataType = m_data.get().getDataType();
        String seqString = dataType.encodingToString(seq);
        
        // Find taxon with same name if tree is labelled
        int taxonNum = node.getNr();
        if (node.getID() != null && !node.getID().isEmpty()) {
	        for (int i = 0; i < m_data.get().getTaxaNames().size(); i ++) {
	        	if (m_data.get().getTaxaNames().get(i).equals(node.getID())) {
	        		taxonNum=i;
	        		break;
	        	}
	        }
        }
        
        String taxon = m_data.get().getTaxaNames().get(taxonNum);
        
        
        return new Sequence(taxon, seqString);
    } // intArray2Sequence

	
    void getTransitionProbabilities(Tree tree, Node node, double rate, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(node));
        branchRate *= rate;
        m_siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);

    } // getTransitionProbabilities

	
}
