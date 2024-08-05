package rbbeast.evolution.util;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.HeapSort;

@Description("logger that logs rates in increasing order + logs associated weights")
public class FreeRateLogger extends BEASTObject implements Loggable {
    public Input<RealParameter> weightInput =
            new Input<RealParameter>("weights", "weights of the various categories, should sum to 1", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rates", "rates for each of the categories, will be weighted and normalised to 1", Validate.REQUIRED);

    private RealParameter rates;
    private RealParameter weights;
    
	@Override
	public void initAndValidate() {
		weights = weightInput.get();
    	rates = rateParameterInput.get();
	}

	@Override
	public void init(PrintStream out) {
		int n = rates.getDimension();
		String partition = rates.getID();
		if (partition.lastIndexOf('.') > 0) {
			partition = partition.substring(partition.lastIndexOf('.'));
		}
		for (int i = 0; i < n ; i++) {
			out.print("sortedRates" + partition + (i+1) + "\t");
		}
		for (int i = 0; i < n ; i++) {
			out.print("sortedWeight" + partition + (i+1) + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		int n = rates.getDimension();
		double [] r = rates.getDoubleValues();
		int [] index = new int[n];
		for (int i = 0; i < n; i++) {
			index[i]= i;
		}
		HeapSort.sort(r, index);
		double sum = 0;
		for (double d : r) {
			sum += d;
		}
		for (int i = 0; i < n ; i++) {
			out.print(rates.getArrayValue(index[i])/sum + "\t");
		}
		for (int i = 0; i < n ; i++) {
			out.print(weights.getArrayValue(index[i]) + "\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
	}

}
