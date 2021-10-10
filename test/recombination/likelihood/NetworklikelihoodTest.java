package recombination.likelihood;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import junit.framework.TestCase;
import recombination.NetworklikelihoodTestCase;
import recombination.alignment.RecombinationAlignment;
import recombination.network.RecombinationNetwork;
import test.beast.BEASTTestCase;

public class NetworklikelihoodTest extends NetworklikelihoodTestCase {

    public NetworklikelihoodTest() {
        super();
    }

    @Test
    public void testJC69Likelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	Alignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
                       
        
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("data", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        
        
        System.out.println("f;dl;fdlk;d");
        
        
        // compute the tree likelihoods for each positions individually
        double treeLog = 0.0;
        for (int i = 0; i < 4; i++) {
            Alignment data_pos = getAlignmentPosition(i);
            Node root = network.getLocusChildren(network.getRootEdge().childNode, i);
            TreeParser t = new TreeParser();

            t.initByName("taxa", data_pos,
                    "newick", root.toNewick(false),
                    "IsLabelledNewick", true);            
           
            TreeLikelihood treelikelihood = newTreeLikelihood();
            treelikelihood.initByName("data", data_pos, "tree", t, "siteModel", siteModel);
            treeLog += treelikelihood.calculateLogP();
        }        
        
        assertEquals(logP, treeLog, BEASTTestCase.PRECISION);
        assertEquals(logP, likelihood.calculateLogP(), BEASTTestCase.PRECISION);
    }
    
    @Test
    public void testHKYLikelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	Alignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
                       
       
        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data);

        
        HKY hky = new HKY();
        hky.initByName("kappa", "29.739445", "frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("data", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        
        
        // compute the tree likelihoods for each positions individually
        double treeLog = 0.0;
        for (int i = 0; i < 4; i++) {
            Alignment data_pos = getAlignmentPosition(i);
            Node root = network.getLocusChildren(network.getRootEdge().childNode, i);
            TreeParser t = new TreeParser();

            t.initByName("taxa", data_pos,
                    "newick", root.toNewick(false),
                    "IsLabelledNewick", true);            
           
            TreeLikelihood treelikelihood = newTreeLikelihood();
            treelikelihood.initByName("data", data_pos, "tree", t, "siteModel", siteModel);
            treeLog += treelikelihood.calculateLogP();
        }    
        
        
        assertEquals(logP, treeLog, BEASTTestCase.PRECISION);
    }
    
    @Test
    public void testHKYGammaLikelihood() throws Exception {
    	Randomizer.setSeed(1);
    	
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	Alignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
        
        System.out.println(network);

        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data);
        
        HKY hky = new HKY();
        hky.initByName("kappa", "2.739445", "frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4, "shape", "10.0", "substModel", hky);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("data", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        
        
        // compute the tree likelihoods for each positions individually
        double treeLog = 0.0;
        for (int i = 0; i < 4; i++) {
            Alignment data_pos = getAlignmentPosition(i);
            Node root = network.getLocusChildren(network.getRootEdge().childNode, i);
            TreeParser t = new TreeParser();

            t.initByName("taxa", data_pos,
                    "newick", root.toNewick(false),
                    "IsLabelledNewick", true);            
           
            TreeLikelihood treelikelihood = newTreeLikelihood();
            treelikelihood.initByName("data", data_pos, "tree", t, "siteModel", siteModel);
            treeLog += treelikelihood.calculateLogP();
        }    
        
        System.out.println(logP);
        assertEquals(logP, treeLog, BEASTTestCase.PRECISION);
    }
    
}
