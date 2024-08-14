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
    final public Input<FreeRateModel> modelInput =
            new Input<>("model", "free rate model so parameters can be accessed", Validate.REQUIRED);
    final public Input<Boolean> categoryCountOnlyInput =
            new Input<>("categoryCountOnly", "only log the category count of the model, not the actual rates", false);
    final public Input<Integer> maxCategoryCountInput = new Input<>("maxCategoryCount", "maximum number of categories allowed");
	
    private RealParameter rates;
    private RealParameter weights;
    private boolean incremental;
    private int maxCategoryCount;
    
	@Override
	public void initAndValidate() {
		weights = modelInput.get().weightInput.get();
    	rates = modelInput.get().rateParameterInput.get();
    	incremental = modelInput.get().incrementalInput.get();
    	if (!categoryCountOnlyInput.get() && maxCategoryCountInput.get() == null) {
    		throw new IllegalArgumentException("maxCategoryCount must be specified if categoryCountOnly=true");
    	}
    	if (!categoryCountOnlyInput.get()) {
    		maxCategoryCount = maxCategoryCountInput.get();
    	}
	}

	@Override
	public void init(PrintStream out) {
		String partition = rates.getID();
		if (partition.lastIndexOf('.') > 0) {
			partition = partition.substring(partition.lastIndexOf('.'));
		}
		if (categoryCountOnlyInput.get()) {
			out.print("categoryCount" + partition + "\t");
		} else {
			out.print("categoryCount" + partition + "\t");
			int n = maxCategoryCount;
			if (weights.isEstimated()) {
				for (int i = 0; i < n ; i++) {
					out.print("sortedWeight" + partition + (i+1) + "\t");
				}
			}
			for (int i = 0; i < n ; i++) {
				out.print("sortedRates" + partition + (i+1) + "\t");
			}
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		int n = rates.getDimension();
		if (categoryCountOnlyInput.get()) {
			out.print(n + "\t");
		} else {
			out.print(n + "\t");
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
			n = Math.min(n, maxCategoryCount);
			if (weights.isEstimated()) {
				for (int i = 0; i < n ; i++) {
					out.print(weights.getArrayValue(index[i]) + "\t");
				}
				for (int i = n; i < maxCategoryCount; i++) {
					out.print(0.0 + "\t");
				}
			}
			for (int i = 0; i < n ; i++) {
				out.print(r[index[i]]/sum + "\t");
			}
			for (int i = n; i < maxCategoryCount; i++) {
				out.print(0.0 + "\t");
			}
		}
	}

	@Override
	public void close(PrintStream out) {
	}

}
