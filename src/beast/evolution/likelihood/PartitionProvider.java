package beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.app.beauti.BeautiDoc;
import beast.app.beauti.PartitionContext;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.sitemodel.SiteModel;
import beast.util.Randomizer;



@Description("Provides dynamic partitioning of an alignment, splitting it in several parts with its own substitution models")
public class PartitionProvider extends CalculationNode implements StateNodeInitialiser {
    public Input<Alignment> alignment= new Input<Alignment>("alignment", "sequence alignment data to be partitioned", Validate.REQUIRED);
    public Input<List<SiteModel.Base>> m_pSiteModel = new Input<List<SiteModel.Base>>("siteModel", "site model for leafs in the beast.tree, one for each partition", new ArrayList<SiteModel.Base>());

    public Input<IntegerParameter> partitionLengthsInput = new Input<IntegerParameter>("partitionLengths", "consecutive length of partitions", Validate.REQUIRED); 
    public Input<Integer> numPartitions = new Input<Integer>("numPartitions", "number of partitions to be used, if positive and there is only one sitemodel, the " +
    		"sitemodel will be duplicated (default 0, i.e., ignored)", 0);
    

    
    int [][] weights;
    int siteCount;
    int patternCount;
    int [] patternIndex;
    List<Integer>[] patternIndicators;
    IntegerParameter partitionLengths;
    
    boolean needsUpdate = true;
	boolean needsInitilising = false;
    
    @Override
    public void initAndValidate() throws Exception {    	
    	if (numPartitions.get() > 0) {
    	    /** we have to wait to initialise till everything is parsed, 
    	     * so that the outputs of all plug-ins are properly set up */
    		if (m_pSiteModel.get().size() == 0 && numPartitions.get() > 1) {
    			needsInitilising = true;
    			return;
    		}
    	}
    	
    	
    	partitionLengths = partitionLengthsInput.get();
    	if (partitionLengths.getDimension() != m_pSiteModel.get().size()) {
    		partitionLengths.setDimension(m_pSiteModel.get().size());
    		//throw new Exception("nr of site models must match nr of partitions defined trough partitionLengths");
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
    		Integer [] partitionLengths2 = new Integer[patternIndicators.length];
    		int delta = siteCount / patternIndicators.length;
    		Arrays.fill(partitionLengths2, delta);
    		partitionLengths2[patternIndicators.length - 1] = siteCount - delta * (patternIndicators.length - 1);
    		IntegerParameter dummy = new IntegerParameter(partitionLengths2);
    		dummy.setLower(partitionLengths.getLower());
    		dummy.setUpper(siteCount);
    		partitionLengths.assignFromWithoutID(dummy);
    		//throw new Exception("nr of sites does not match partitions lengths (sum=" + sum + " != siteCount=" + siteCount+")");
    	}
    	
//    	weights = new int[patternIndicators.length][patternCount];
//    	int [] patternWeights = alignment.get().getWeights();
//    	for (int i = 0; i < patternCount; i++) {
//    		patternIndicators[i%2].add(i);
//    		weights[i%2][i] = patternWeights[i];
//    	}

    	//update();
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

	public boolean delayInitialisation() {
		return needsInitilising;
	}

    /** we have to wait to initialise till everything is parsed, 
     * so that the outputs of all plug-ins are properly set up.
     * The StateNodeInitialiser interface can be abused for this purpose
     */
	@Override
	public void initStateNodes() throws Exception {
    	System.err.println(Randomizer.getSeed());
		PartitionedTreeLikelihood likelihood = null;
		for (BEASTInterface plugin : getOutputs()) {
			if (plugin instanceof PartitionedTreeLikelihood) {
				likelihood = (PartitionedTreeLikelihood) plugin;
			}
		}
		if (likelihood == null) {
			throw new Exception("PartitionProvider must have PartitionedTreeLikelihood as output");
		}
		SiteModel.Base siteModel = (SiteModel.Base) likelihood.siteModelInput.get();
		
		if (siteModel instanceof SiteModel) {
			RealParameter mu = ((SiteModel) siteModel).muParameterInput.get();
			if (!mu.isEstimatedInput.get()) {
				System.err.println("Warning: sitemodel mutation rates are NOT estimated. This is probably not what you want");
			}
		}
		
		String sPartition = BeautiDoc.parsePartition(getID());
		Set<BEASTInterface> ancestors = new HashSet<BEASTInterface>();
		BeautiDoc.collectAncestors(this, ancestors, new HashSet<BEASTInterface>());
		MCMC mcmc = null;
		for (BEASTInterface plugin : ancestors) {
			if (plugin instanceof MCMC) {
				mcmc = (MCMC) plugin;
				break;
			}
		}
		m_pSiteModel.get().add(siteModel);
		List<BEASTInterface> tabuList = new ArrayList<BEASTInterface>();
		tabuList.add(this);
		for (int i = 1; i < numPartitions.get(); i++) {
			BeautiDoc doc = new BeautiDoc();
			PartitionContext context = new PartitionContext(sPartition+i);
			BEASTInterface plugin = BeautiDoc.deepCopyPlugin(siteModel, this, mcmc, context, doc, tabuList);
			m_pSiteModel.get().add((SiteModel.Base) plugin);
		}
		if (siteModel instanceof SiteModel) {
			RealParameter mu = ((SiteModel) siteModel).muParameterInput.get();
			mu.isEstimatedInput.setValue(false, mu);
			mu.initByName("lower", ""+mu.getValue(), "upper" , "" + mu.getValue());
		}
		State state = mcmc.startStateInput.get();
		state.initialise();
        state.setPosterior(mcmc.posteriorInput.get());
		needsInitilising = false;
		initAndValidate();

		Set<BEASTInterface> plugins = new HashSet<BEASTInterface>();
		for (BEASTInterface plugin : mcmc.listActivePlugins()) {
			reinitialise(plugin, plugins);
		}
//		
//		for (Plugin plugin : ancestors) {
//			if (plugin instanceof PartitionedTreeLikelihood) {
//				((PartitionedTreeLikelihood) plugin).initAndValidate();
//				break;
//			}
//		}
	}

	private void reinitialise(BEASTInterface plugin, Set<BEASTInterface> plugins) throws Exception {
		for (BEASTInterface plugin2 : plugin.listActivePlugins()) {
			if (!plugins.contains(plugin2)) {
				plugins.add(plugin2);
				reinitialise(plugin2, plugins);
				// somehow, the process now produces a disconnected empty TaxonSet, which needs to be initialised
				// TODO: find out how this TaxonSet came about.
				if (plugin2 instanceof TaxonSet) {
					TaxonSet taxa = (TaxonSet) plugin2;
					if (taxa.alignmentInput.get() == null && taxa.taxonsetInput.get().size() == 0 && taxa.getOutputs().size()==0) {
					} else {
						plugin2.initAndValidate();
					}
				} else if (plugin2 instanceof MCMC) {
					MCMC mcmc = (MCMC) plugin2;
					if (mcmc.sampleFromPriorInput.get()) {
						// likelihood was removed the first time around
						// need to add it back into posterior
						CompoundDistribution posterior = (CompoundDistribution) mcmc.posteriorInput.get();
						CompoundDistribution likelihood = new CompoundDistribution();
						likelihood.setID("likelihood");
						posterior.pDistributions.get().add(likelihood);
					}
					plugin2.initAndValidate();
				} else {
					plugin2.initAndValidate();
				}
			}
		}		
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		SiteModel siteModel = (SiteModel) m_pSiteModel.get().get(0);
		stateNodes.add(siteModel.muParameterInput.get());
	}
}
