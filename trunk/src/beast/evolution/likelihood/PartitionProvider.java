package beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;

@Description("Provides dynamic partitioning of an alignment, splitting it in several parts with its own substitution models")
public class PartitionProvider extends CalculationNode {
    public Input<Alignment> alignment= new Input<Alignment>("alignment", "sequence alignment data to be partitioned", Validate.REQUIRED);
    public Input<List<SiteModel.Base>> m_pSiteModel = new Input<List<SiteModel.Base>>("siteModel", "site model for leafs in the beast.tree, one for each partition", new ArrayList<SiteModel.Base>(), Validate.REQUIRED);

    public Input<IntegerParameter> partitionLengthsInput = new Input<IntegerParameter>("partitionLengths", "consecutive length of partitions", Validate.REQUIRED); 
    
    int [][] weights;
    int siteCount;
    int patternCount;
    int [] patternIndex;
    List<Integer>[] patternIndicators;
    IntegerParameter partitionLengths;
    
    boolean needsUpdate = true;
    
    @Override
    public void initAndValidate() throws Exception {
    	partitionLengths = partitionLengthsInput.get();
    	if (partitionLengths.getDimension() != m_pSiteModel.get().size()) {
    		throw new Exception("nr of site models must match nr of partitions defined trough partitionLengths");
    	}
    			
    	siteCount = alignment.get().getSiteCount();
    	patternIndex = new int[siteCount];
    	patternCount = alignment.get().getPatternCount();
    	for (int i = 0; i < siteCount; i++) {
        	patternIndex[i] = alignment.get().getPatternIndex(i);
    	}
    	
    	patternIndicators = new ArrayList[partitionLengths.getDimension()];
    	for (int i = 0; i < patternIndicators.length; i++) {
    		patternIndicators[i] = new ArrayList<Integer>();
    	}

    	// sanity check: each site must be accounted for
    	int sum = 0;
    	for (int i = 0; i < patternIndicators.length; i++) {
    		sum += partitionLengths.getValue(i);
    	}
    	if (sum != siteCount) {
    		throw new Exception("nr of sites does not match partitions lengths (sum=" + sum + " != siteCount=" + siteCount+")");
    	}
    	
//    	weights = new int[patternIndicators.length][patternCount];
//    	int [] patternWeights = alignment.get().getWeights();
//    	for (int i = 0; i < patternCount; i++) {
//    		patternIndicators[i%2].add(i);
//    		weights[i%2][i] = patternWeights[i];
//    	}

    	update();
    }
    
    void update() {
    	for (int i = 0; i < patternIndicators.length; i++) {
    		patternIndicators[i].clear();
    	}
    	weights = new int[patternIndicators.length][patternCount];
    	int partition = 0;
    	int upper = partitionLengths.getValue(partition);
    	for (int i = 0; i < siteCount; i++) {
        	while (i >= upper) {
        		partition++;
        		upper += partitionLengths.getValue(partition);
        	}
    		int pattern = patternIndex[i];
    		if (weights[partition][pattern] == 0) {
    			patternIndicators[partition].add(pattern);
    		}
    		weights[partition][pattern]++;
    	}
    	needsUpdate = false;
    }
    
    public int getPartitionCount() {
    	if (needsUpdate) {
    		update();
    	}
    	return patternIndicators.length;
    }

    public SiteModel.Base getSiteModel(int partition) {
    	if (needsUpdate) {
    		update();
    	}
    	return m_pSiteModel.get().get(partition);
    }
    
    public List<Integer> getPatternIndicators(int partition) {
    	if (needsUpdate) {
    		update();
    	}
    	return patternIndicators[partition];
    }
    
    public int [] getPatternWeights(int partition) {
    	if (needsUpdate) {
    		update();
    	}
    	return weights[partition];
    }

    public boolean aPartitionChanged() {
    	if (needsUpdate) {
    		update();
    	}
    	return true;
    }

    public enum PartitionChange {None, NonZeroWeights_Only, WeightsFromZeroToOne, WeightsFromOneToZero, Added, Deleted} 
	public PartitionChange partitionChanged(int partition) {
    	if (needsUpdate) {
    		update();
    	}
		return PartitionChange.Added;
	}
	
	
	@Override
	protected void store() {
		needsUpdate = true;
		super.store();
	}
	
	@Override
	protected void restore() {
		needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (partitionLengthsInput.isDirty()) {
			needsUpdate = true;
		}
		return true;
	}
}
