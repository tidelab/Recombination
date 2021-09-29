package recombination.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.util.NexusParser;
import recombination.likelihood.NetworkLikelihood;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.operators.AddRemoveRecombination;

public class RankRecombinantEdgesFromNetwork {
	
	private RecombinationNetwork network;
	private NetworkLikelihood likelihood;
	private AddRemoveRecombination removeOperator;
	private double logP_full;
	
	public RankRecombinantEdgesFromNetwork(RecombinationNetwork network, NetworkLikelihood likelihood) {
		this.network = network;
		this.likelihood = likelihood;
        this.removeOperator = new AddRemoveRecombination();
    	this.removeOperator.initByName("alpha", 1.0, "network", network, "weight", 1.0);
    	
    	logP_full = likelihood.calculateLogP();
	}
	
	private static class Options {
        File networkFile, alignmentFile, outputFile;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: network_newick alignment_nexus output_file");
        System.exit(exitCode);
    }

    /**
     * Process command line arguments.
     *
     * @param args list of arguments given in the command line
     * @return an Option variable
     */
    public static Options processArguments(String[] args) {

        Options options = new Options();

        if (args.length !=3)
            printUsageAndExit(0);
        int i = 0;
        options.networkFile = new File(args[i++]);
        options.alignmentFile = new File(args[i++]);
        options.outputFile = new File(args[i++]);


        return options;
    }
    
    public static RecombinationNetwork networkFromNewick(File networkFile) throws IOException {
    	BufferedReader reader = new BufferedReader(new FileReader(networkFile));
    	return(new RecombinationNetwork(reader.readLine()));
    }
    
   
    public double getLikelihoodGain(RecombinationNetworkEdge recombinantEdge) {
    	removeOperator.removeRecombinationEdge(recombinantEdge);
    	double logP = likelihood.calculateLogP();
    	network.restore();
    	return(logP_full-logP);
    }
    
    public double getMinimumLikelihoodGain(RecombinationNetworkNode recombinantNode) {
    	final int nodeID = recombinantNode.ID;
    	RecombinationNetworkNode node = network.getNodes().stream().filter(e -> e.ID == nodeID).findFirst().get();
    	RecombinationNetworkEdge edge = node.getParentEdges().get(0);
    	List<RecombinationNetworkEdge> removableEdges = removeOperator.getRemovableEdges();
    	double likelihoodGainLeft = Double.POSITIVE_INFINITY;
    	if (removableEdges.contains(edge)) {
    		likelihoodGainLeft = getLikelihoodGain(edge);
    	}
        removableEdges = removeOperator.getRemovableEdges();
        node = network.getNodes().stream().filter(e -> e.ID == nodeID).findFirst().get();
        edge = node.getParentEdges().get(1);
    	double likelihoodGainRight = Double.POSITIVE_INFINITY;
    	if (removableEdges.contains(edge)) {
    		likelihoodGainRight = getLikelihoodGain(edge);
    	}
    	return(Math.min(likelihoodGainLeft, likelihoodGainRight));
    }
    
    
    
    
    public static void main(String[] args) throws IOException {
    	
    	Options options = processArguments(args);

    	RecombinationNetwork network = networkFromNewick(options.networkFile);
    	Alignment alignment;
        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(options.alignmentFile);
        alignment = nexusParser.m_alignment;
        
        HKY HKY = new HKY();
        RealParameter freqValues = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25});
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", freqValues, "estimate", false);
        HKY.initByName("kappa", "4.0", "frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1", "gammaCategoryCount", 1, "substModel", HKY);

        NetworkLikelihood likelihood = new NetworkLikelihood();
        likelihood.initByName("data", alignment, "recombinationNetwork", network, "siteModel", siteModel);
        
    	try (PrintStream ps = new PrintStream(options.outputFile)) {
        	RankRecombinantEdgesFromNetwork f = new RankRecombinantEdgesFromNetwork(network, likelihood);
        	
        	System.out.println( "-----------------------------------------");

        	
        	List<RecombinationNetworkNode> recombinantNodes = network.getNodes().stream()
                    .filter(e -> e.isRecombination())
                    .collect(Collectors.toList());
        	System.out.println("found " + recombinantNodes.size() + " recombinant nodes");
        	System.out.println(f.removeOperator.getRemovableEdges().size()+" edges can be removed");
        	List<Double> lhGains = new ArrayList<Double>();
        	for (RecombinationNetworkNode n : recombinantNodes) {
        		lhGains.add(f.getMinimumLikelihoodGain(n));
        	}
        	recombinantNodes = network.getNodes().stream()
                    .filter(e -> e.isRecombination())
                    .collect(Collectors.toList());
        	int index = 0;
        	for (RecombinationNetworkNode n : recombinantNodes) {

        		if(n.getMetaData()!=null) {
        			n.setMetaData(n.getMetaData()+",minimumLikelihoodGain="+lhGains.get(index++));
        		} else {
        			n.setMetaData(",minimumLikelihoodGain="+lhGains.get(index++));
        		}
      		}
        	System.out.println("Annotated Tree: "+network.toString());
        	ps.append(network.toString());
        	ps.close();
        }
        
    }
}
