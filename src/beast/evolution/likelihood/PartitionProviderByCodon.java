package beast.evolution.likelihood;

import java.util.ArrayList;

import beast.base.core.Description;



@Description("Splits alignment into partititions by lengths, and subdivides these partitions by codon position")
public class PartitionProviderByCodon extends PartitionProvider {

	
	public void initAndValidate() {
    	partitionLengths = partitionLengthsInput.get();
    	if (partitionLengths.getDimension() * 3 != m_pSiteModel.get().size()) {
    		throw new IllegalArgumentException("nr of site models must be 3 times nr of partitions defined trough partitionLengths");
    	}
    			
    	siteCount = alignment.get().getSiteCount();
    	patternIndex = new int[siteCount];
    	patternCount = alignment.get().getPatternCount();
    	for (int i = 0; i < siteCount; i++) {
        	patternIndex[i] = alignment.get().getPatternIndex(i);
    	}
    	
    	patternIndicators = new ArrayList[m_pSiteModel.get().size()];
    	for (int i = 0; i < patternIndicators.length; i++) {
    		patternIndicators[i] = new ArrayList<Integer>();
    	}

    	// sanity check: each site must be accounted for
    	int sum = 0;
    	for (int i = 0; i < partitionLengths.getDimension(); i++) {
    		sum += partitionLengths.getValue(i);
    	}
    	if (sum * 3 != siteCount + (siteCount % 3 == 0? 0 : 3-siteCount % 3)) {
    		throw new IllegalArgumentException("nr of sites does not match partitions lengths (sum=" + sum + "*3 != siteCount=" + siteCount+")");
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
    	int upper = partitionLengths.getValue(partition) * 3;
    	for (int i = 0; i < siteCount; i++) {
        	while (i >= upper) {
        		partition+=3;
        		upper += partitionLengths.getValue(partition/3) * 3;
        	}
    		int pattern = patternIndex[i];
    		if (weights[partition + i % 3][pattern] == 0) {
    			patternIndicators[partition + i % 3].add(pattern);
    		}
    		weights[partition + i % 3][pattern]++;
    	}
    	needsUpdate = false;
    }
}
