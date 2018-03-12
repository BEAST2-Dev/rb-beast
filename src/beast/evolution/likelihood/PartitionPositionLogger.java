package beast.evolution.likelihood;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;



@Description("Logs the position of partion boundaries based on partition lengths")
public class PartitionPositionLogger extends BEASTObject implements Loggable {
    public Input<IntegerParameter> partitionLengthsInput = new Input<IntegerParameter>("partitionLengths", "consecutive length of partitions", Validate.REQUIRED);
    
    IntegerParameter partitionLengths;
    
    @Override
    public void initAndValidate() {
    	partitionLengths = partitionLengthsInput.get();
    }

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i< partitionLengths.getDimension() - 1; i++) {
			out.append(partitionLengths.getID() + ".pos." + (i+1) + "\t");
		}
	}

	@Override
	public void log(long nSample, PrintStream out) {
		int sum = 0;
		for (int i = 0; i< partitionLengths.getDimension() - 1; i++) {
			sum += partitionLengths.getValue(i);
			out.append(sum + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

}
