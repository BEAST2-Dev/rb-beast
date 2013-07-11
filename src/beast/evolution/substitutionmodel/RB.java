package beast.evolution.substitutionmodel;


import java.util.Arrays;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Valuable;
import beast.core.parameter.IntegerParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;

@Description("Substitution model for nucleotides that changes where the count input " +
		"determines the number of parameters used in a hierarchy of models")
@Citation("Bouckaert, Remco, M—nica V. Alvarado-Mora, and J. R. Pinho. Evolutionary rates and HBV: issues of rate estimation with Bayesian molecular methods. Antiviral therapy 18.3 Pt B (2012): 497-503.")
public class RB extends GeneralSubstitutionModel {
	public Input<IntegerParameter> m_countInput = new Input<IntegerParameter>("count", "model number used 0 = JC, 5 and higher if GTR (default 0)", Validate.REQUIRED);
	
	Valuable m_rate;
	IntegerParameter m_count;
	
	@Override
	public void initAndValidate() throws Exception {
		m_count = m_countInput.get();
		m_rate = m_rates.get();
		if (m_rate.getDimension() != 5) {
			throw new Exception("rate input must have dimension 5");
		}
		
    	m_frequencies = frequenciesInput.get();
        updateMatrix = true;
        m_nStates = m_frequencies.getFreqs().length;
    	if (m_nStates != 4) {
    		throw new Exception("Frequencies has wrong size. Expected 4, but got " + m_nStates);
    	}

        eigenSystem = createEigenSystem();
        m_rateMatrix = new double[m_nStates][m_nStates];
        relativeRates = new double[m_nStates * (m_nStates - 1)];
        storedRelativeRates = new double[m_nStates * (m_nStates - 1)];
	}

	@Override
    protected void setupRelativeRates() {
		switch (m_count.getValue()) {
		case 0: // JC96
			Arrays.fill(relativeRates, 1.0);
			break;
		case 1: // HKY
	    	relativeRates[0] = m_rate.getArrayValue(0); // A->C
	    	relativeRates[1] = 1; // A->G
	    	relativeRates[2] = m_rate.getArrayValue(0); // A->T

	    	relativeRates[4] = m_rate.getArrayValue(0); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = m_rate.getArrayValue(0); // G->T
			break;
		case 2: // K3P
	    	relativeRates[0] = m_rate.getArrayValue(0); // A->C
	    	relativeRates[1] = m_rate.getArrayValue(1); // A->G
	    	relativeRates[2] = m_rate.getArrayValue(0); // A->T

	    	relativeRates[4] = m_rate.getArrayValue(0); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = m_rate.getArrayValue(0); // G->T
			break;
		case 3: // TIM
	    	relativeRates[0] = m_rate.getArrayValue(0); // A->C
	    	relativeRates[1] = m_rate.getArrayValue(1); // A->G
	    	relativeRates[2] = m_rate.getArrayValue(2); // A->T

	    	relativeRates[4] = m_rate.getArrayValue(2); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = m_rate.getArrayValue(0); // G->T
			break;
		case 4: // new model 
	    	relativeRates[0] = m_rate.getArrayValue(0); // A->C
	    	relativeRates[1] = m_rate.getArrayValue(1); // A->G
	    	relativeRates[2] = m_rate.getArrayValue(2); // A->T

	    	relativeRates[4] = m_rate.getArrayValue(3); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = m_rate.getArrayValue(0); // G->T
			break;
		default: // GTR
	    	relativeRates[0] = m_rate.getArrayValue(0); // A->C
	    	relativeRates[1] = m_rate.getArrayValue(1); // A->G
	    	relativeRates[2] = m_rate.getArrayValue(2); // A->T

	    	relativeRates[4] = m_rate.getArrayValue(3); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = m_rate.getArrayValue(4); // G->T
			break;
		}

    	relativeRates[3] = relativeRates[0]; // C->A

    	relativeRates[6] = relativeRates[1]; // G->A
    	relativeRates[7] = relativeRates[4]; // G->C

    	relativeRates[9] = relativeRates[2]; // T->A
    	relativeRates[10] = 1; //T->C
    	relativeRates[11] = relativeRates[8]; //T->G
	}
	
	@Override
	public boolean canHandleDataType(DataType dataType) throws Exception {
		if (dataType instanceof Nucleotide) {
			return true;
		}
		throw new Exception("Can only handle nucleotide data");
	}
	
}
