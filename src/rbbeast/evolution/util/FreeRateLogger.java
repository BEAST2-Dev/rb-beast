package rbbeast.evolution.util;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.HeapSort;
import rbbeast.evolution.sitemodel.FreeRateModel;

@Description("logger that logs rates in increasing order + logs associated weights")
public class FreeRateLogger extends BEASTObject implements Loggable {
    public Input<FreeRateModel> modelInput =
            new Input<>("model", "free rate model so parameters can be accessed", Validate.REQUIRED);
	
    private RealParameter rates;
    private RealParameter weights;
    private boolean incremental;
    
	@Override
	public void initAndValidate() {
		weights = modelInput.get().weightInput.get();
    	rates = modelInput.get().rateParameterInput.get();
    	incremental = modelInput.get().incrementalInput.get();
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
		if (incremental) {
			for (int i = 1; i < r.length; i++) {
				r[i] += r[i-1];
			}
		}
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
			out.print(r[index[i]]/sum + "\t");
		}
		for (int i = 0; i < n ; i++) {
			out.print(weights.getArrayValue(index[i]) + "\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
	}

}
