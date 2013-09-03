package test.beast.evolution.likelihood;

import java.io.File;

import org.junit.Test;

import beast.core.MCMC;
import beast.util.Randomizer;
import beast.util.XMLParser;



import junit.framework.TestCase;

public class PartitionedLikelihoodTest extends TestCase {
	
	@Test public void testPartitionedLikelihood() {
		Randomizer.setSeed(123);
		XMLParser parser = new XMLParser();
		try {
			new File("test.123.trees").delete();
			new File("test.123.log").delete();
			
			
			MCMC mcmc = (MCMC) parser.parseFile(new File("examples/testHKY.xml"));
			double logP = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
			System.err.println(logP);
			assertEquals(-4053.296578683104, logP, 1e-10);
			
			mcmc.chainLengthInput.setValue(1000, mcmc);
			mcmc.run();
			logP = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
			System.err.println(logP);
			assertEquals(-1817.0926689847693, logP, 1e-10);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
