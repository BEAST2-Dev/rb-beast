package beast.evolution.substitutionmodel;

import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

@Description("Substitution model for nucleotides that changes where the count input " +
		"determines the number of parameters used in a hierarchy of models")
@Citation("Alexei J. Drummond and Remco R. Bouckaert. Bayesian evolutionary analysis with BEAST 2. CUP. 2014")
public class VS extends GeneralSubstitutionModel {
	public Input<IntegerParameter> countInput = new Input<IntegerParameter>("count", "model number used 0 = JC, 5 and higher if GTR", Validate.REQUIRED);
	
	Function rate;
	IntegerParameter count;
	
	@Override
	public void initAndValidate() {
		count = countInput.get();
		rate = ratesInput.get();
		if (rate.getDimension() != 5) {
			throw new IllegalArgumentException("rate input must have dimension 5");
		}
		
    	frequencies = frequenciesInput.get();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
    	if (nrOfStates != 4) {
    		throw new IllegalArgumentException("Frequencies has wrong size. Expected 4, but got " + nrOfStates);
    	}

        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e);
		}
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)];
	}

	@Override
    public void setupRelativeRates() {
		switch (count.getValue()) {
		case 0: // JC96
			Arrays.fill(relativeRates, 1.0);
			break;
		case 1: // HKY
	    	relativeRates[0] = rate.getArrayValue(0); // A->C
	    	relativeRates[1] = 1; // A->G
	    	relativeRates[2] = rate.getArrayValue(0); // A->T

	    	relativeRates[4] = rate.getArrayValue(0); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = rate.getArrayValue(0); // G->T
			break;
		case 2: // K3P
	    	relativeRates[0] = rate.getArrayValue(0); // A->C
	    	relativeRates[1] = rate.getArrayValue(1); // A->G
	    	relativeRates[2] = rate.getArrayValue(0); // A->T

	    	relativeRates[4] = rate.getArrayValue(0); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = rate.getArrayValue(0); // G->T
			break;
		case 3: // TIM
	    	relativeRates[0] = rate.getArrayValue(0); // A->C
	    	relativeRates[1] = rate.getArrayValue(1); // A->G
	    	relativeRates[2] = rate.getArrayValue(2); // A->T

	    	relativeRates[4] = rate.getArrayValue(2); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = rate.getArrayValue(0); // G->T
			break;
		case 4: // new model 
	    	relativeRates[0] = rate.getArrayValue(0); // A->C
	    	relativeRates[1] = rate.getArrayValue(1); // A->G
	    	relativeRates[2] = rate.getArrayValue(2); // A->T

	    	relativeRates[4] = rate.getArrayValue(3); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = rate.getArrayValue(0); // G->T
			break;
		default: // GTR
	    	relativeRates[0] = rate.getArrayValue(0); // A->C
	    	relativeRates[1] = rate.getArrayValue(1); // A->G
	    	relativeRates[2] = rate.getArrayValue(2); // A->T

	    	relativeRates[4] = rate.getArrayValue(3); // C->G
	    	relativeRates[5] = 1; // C->T

	    	relativeRates[8] = rate.getArrayValue(4); // G->T
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
	public boolean canHandleDataType(DataType dataType) {
		if (dataType instanceof Nucleotide) {
			return true;
		}
		//throw new Exception("Can only handle nucleotide data");
		return false;
	}
	
}
