package recombination.likelihood;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Tree;
import junit.framework.TestCase;
import recombination.alignment.RecombinationAlignment;
import recombination.network.RecombinationNetwork;
import test.beast.BEASTTestCase;

public class NetworklikelihoodTest extends TestCase {

    public NetworklikelihoodTest() {
        super();
    }

    protected NetworkLikelihood newNetworkLikelihood() {
    	System.setProperty("java.only","true");
        return new NetworkLikelihood();
    }

    @Test
    public void testJC69Likelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	RecombinationAlignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
        
        
        

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("data", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);

//        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
    }
    
    
    static public RecombinationAlignment getAlignment() throws Exception {
        Sequence human = new Sequence("human", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA");
        Sequence chimp = new Sequence("chimp", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA");
        Sequence bonobo = new Sequence("bonobo", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTTAAATCCCCTTATTTCTACTAGGACTATGAGAGTCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCTCTCAGTAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAGC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTTGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCCCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCAACACCTCTTTACAGTGA");
        Sequence gorilla = new Sequence("gorilla", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA");
        Sequence orangutan = new Sequence("orangutan", "AGAAATATGTCTGACAAAAGAGTTACTTTGATAGAGTAAAAAATAGAGGTCTAAATCCCCTTATTTCTACTAGGACTATGGGAATTGAACCCACCCCTGAGAATCCAAAATTCTCCGTGCCACCCATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTA--CACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTCA-CAGCACTTAATTTCTGTAAGGACTGCAAAACCCCACTTTGCATCAACTGAGCGCAAATCAGCCACTTTAATTAAGCTAAGCCCTCCTAGACCGATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAT-TGGCTTCAGTCCAAAGCCCCGGCAGGCCTTAAAGCTGCTCCTTCGAATTTGCAATTCAACATGACAA-TCACCTCAGGGCTTGGTAAAAAGAGGTCTGACCCCTGTTCTTAGATTTACAGCCTAATGCCTTAACTCGGCCATTTTACCGCAAAAAAGGAAGGAATCGAACCTCCTAAAGCTGGTTTCAAGCCAACCCCATAACCCCCATGACTTTTTCAAAAGGTACTAGAAAAACCATTTCGTAACTTTGTCAAAGTTAAATTACAGGTC-AGACCCTGTGTATCTTA-CATTGCAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACCAGCCTCTCTTTGCAATGA");
        Sequence siamang = new Sequence("siamang", "AGAAATACGTCTGACGAAAGAGTTACTTTGATAGAGTAAATAACAGGGGTTTAAATCCCCTTATTTCTACTAGAACCATAGGAGTCGAACCCATCCTTGAGAATCCAAAACTCTCCGTGCCACCCGTCGCACCCTGTTCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCATACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTAACAAAACTTAATTTCTGCAAGGGCTGCAAAACCCTACTTTGCATCAACCGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCGATGGGACTTAAACCCATAAAAATTTAGTTAACAGCTAAACACCCTAAACAACCTGGCTTCAATCTAAAGCCCCGGCAGA-GTTGAAGCTGCTTCTTTGAACTTGCAATTCAACGTGAAAAATCACTTCGGAGCTTGGCAAAAAGAGGTTTCACCTCTGTCCTTAGATTTACAGTCTAATGCTTTA-CTCAGCCACTTTACCACAAAAAAGGAAGGAATCGAACCCTCTAAAACCGGTTTCAAGCCAGCCCCATAACCTTTATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATCACAGGTCCAAACCCCGTATATCTTATCACTGTAGAGCTAGACCAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACTACCGCCTCTTTACAGTGA");

        RecombinationAlignment data = new RecombinationAlignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public RecombinationAlignment getAlignmentShort() throws Exception {
        Sequence t0 = new Sequence("t0", "AGAA");
        Sequence t1 = new Sequence("t1", "AGAT");
        Sequence t2 = new Sequence("t2", "AGAA");

        RecombinationAlignment data = new RecombinationAlignment();
        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public RecombinationNetwork getNetworkShort() throws Exception {
        RecombinationNetwork network = new RecombinationNetwork(
        		"((#H0[&split={0-929},loci={0-929},length=930]:0.4860702162314561,((#H2[&split={7633-9999},loci={7633-9999},length=2367]:0.036380107508342974,(t1[&loci={0-9999},length=10000]:0.29041085418573037,t0[&loci={0-9999},length=10000]:0.29041085418573037)[&loci={0-9999},length=10000]:0.1528435079144485)[&loci={0-9999},length=10000]:0.47005038739308824)#H1[&split={0-5172},loci={0-5172},length=5173]:1.5817575814045295)[&loci={0-5172},length=5173]:0.35688825595463936,(((t2[&loci={0-9999},length=10000]:0.4068742545918359)#H2[&split={0-7632},loci={0-7632},length=7633]:0.604218532359179,#H1[&split={5173-9999},loci={5173-9999},length=4827]:0.09778803745774778)[&loci={0-9999},length=10000]:0.9978993277153256)#H0[&split={930-9999},loci={930-9999},length=9070]:0.8429584721860954)[&loci={0-9999},length=10000]:0.0;");
        return network;
    }



}
