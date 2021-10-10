package recombination;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.TreeLikelihood;
import junit.framework.TestCase;
import recombination.likelihood.NetworkLikelihood;
import recombination.network.RecombinationNetwork;

public abstract class NetworklikelihoodTestCase extends TestCase {
	
	public NetworklikelihoodTestCase() {
		super();
	}

	   protected NetworkLikelihood newNetworkLikelihood() {
	    	System.setProperty("java.only","true");
	        return new NetworkLikelihood();
	    }
	    
	    protected TreeLikelihood newTreeLikelihood() {
	    	System.setProperty("java.only","true");
	        return new TreeLikelihood();
	    }
	    
	    static public Alignment getAlignmentShort() throws Exception {
	        Sequence t0 = new Sequence("t0", "AGAN");
	        Sequence t1 = new Sequence("t1", "AGAT");
	        Sequence t2 = new Sequence("t2", "AGAA");

	        Alignment data = new Alignment();
	        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
	                "dataType", "nucleotide"
	        );
	        return data;
	    }
	    
	    static public Alignment getAlignmentPosition(int i) throws Exception {
	    	String t0_st = "AGAN";
	    	String t1_st = "AGAT";
	    	String t2_st = "AGAA";
	        Sequence t0 = new Sequence("t0", t0_st.substring(i, i+1));
	        Sequence t1 = new Sequence("t1", t1_st.substring(i, i+1));
	        Sequence t2 = new Sequence("t2", t2_st.substring(i, i+1));

	        Alignment data = new Alignment();
	        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
	                "dataType", "nucleotide"
	        );
	        return data;
	    }
	    
	    
	    static public RecombinationNetwork getNetworkShort() throws Exception {
	        RecombinationNetwork network = new RecombinationNetwork(
	        		"((#H0[&split={0-0},loci={0-0},length=4]:0.4860702162314561,((#H2[&split={1-3},loci={1-3},length=3]:0.036380107508342974,(t1[&loci={0-3},length=4]:0.29041085418573037,t0[&loci={0-3},length=4]:0.29041085418573037)[&loci={0-3},length=4]:0.1528435079144485)[&loci={0-3},length=4]:0.47005038739308824)#H1[&split={0-0},loci={0-0},length=1]:1.5817575814045295)[&loci={0-0},length=1]:0.35688825595463936,(((t2[&loci={0-3},length=4]:0.4068742545918359)#H2[&split={0-0},loci={0-0},length=1]:0.604218532359179,#H1[&split={1-3},loci={1-3},length=3]:0.09778803745774778)[&loci={0-3},length=4]:0.9978993277153256)#H0[&split={1-3},loci={1-3},length=3]:0.8429584721860954)[&loci={0-3},length=4]:0.0;");
	        return network;
	    }
}
