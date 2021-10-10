package recombination.util;

import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import recombination.NetworklikelihoodTestCase;
import recombination.likelihood.NetworkLikelihood;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import test.beast.BEASTTestCase;


public class LikelihoodGainTest extends NetworklikelihoodTestCase {
	
	public LikelihoodGainTest() {
		super();
	}

	@Test
	public void testLikelihoodGain() throws Exception {
    	Alignment data = getAlignmentShort();
		RecombinationNetwork network = getNetworkShort();
	    JukesCantor JC = new JukesCantor();
	    JC.initAndValidate();

	    SiteModel siteModel = new SiteModel();
	    siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

	    NetworkLikelihood likelihood = newNetworkLikelihood();
	    likelihood.initByName("data", data, "recombinationNetwork", network, "siteModel", siteModel);
	    double logP = likelihood.calculateLogP();   

		LikelihoodGain lhGain = new LikelihoodGain();
		lhGain.initByName("recombinationNetwork", network, "likelihood", likelihood);
		
		RecombinationNetworkEdge edgeToRemove = network.getRootEdge().childNode.getChildEdges().get(0).childNode.getChildEdges().get(1);

		//Network with the edge removed:
		String newickString = "(((#H1[&loci={0-0},pr={0-0},length=1]:0.6042185323591791,((t2[&loci={0-3},length=4]:0.4068742545918358)#H1[&loci={1-3},pr={1-3},length=3]:0.036380107508342974,(t1[&loci={0-3},length=4]:0.2904108541857302,t0[&loci={0-3},length=4]:0.2904108541857302)[&loci={0-3},length=4]:0.15284350791444856)[&loci={0-3},length=4]:0.5678384248508361)[&loci={0-3},length=4]:0.9978993277153256)#H0[&loci={1-3},pr={1-3},length=3]:0.8429584721860954,#H0[&loci={0-0},pr={0-0},length=1]:0.8429584721860954)[&loci={0-3},length=4]:0.0;";
		RecombinationNetwork networkEdgeRemoved = new RecombinationNetwork(newickString);
		NetworkLikelihood likelihoodEdgeRemoved = newNetworkLikelihood();
	    likelihoodEdgeRemoved.initByName("data", data, "recombinationNetwork", networkEdgeRemoved, "siteModel", siteModel);
	    double logPEdgeRemoved = likelihoodEdgeRemoved.calculateLogP();
	   	
	    assertEquals(logP-logPEdgeRemoved, lhGain.getLikelihoodGain(edgeToRemove), BEASTTestCase.PRECISION);

	}
}
