package recombination.util;

import java.util.List;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import recombination.likelihood.NetworkLikelihood;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.operators.AddRemoveRecombination;

public class LikelihoodGain extends BEASTObject {

	
	private RecombinationNetwork network;
	private RecombinationNetwork networkCopy;
	private NetworkLikelihood likelihood;
	private NetworkLikelihood likelihoodCopy;
	private AddRemoveRecombination removeOperator;
	
	private double logP_full;
	
	final public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	final public Input<NetworkLikelihood> likelihoodInput = new Input<>("likelihood", "likelihood function", Validate.REQUIRED);


	@Override
	public void initAndValidate() {
		network = networkInput.get();
		likelihood = likelihoodInput.get();
	   	logP_full = likelihood.calculateLogP();
	   	createCopy();
	}
	
	
	/* Creates a copy of the current network (and likelihood calculation). 
	 * Edge removals are done on the copy only and leave the original network intact.
	 */
	private void createCopy() {
		networkCopy = network.copy();
		likelihoodCopy = new NetworkLikelihood();
		likelihoodCopy.initByName("useAmbiguities", likelihood.m_useAmbiguities.get(), "useTipLikelihoods", likelihood.m_useTipLikelihoods.get(), "implementation", likelihood.implementationInput.get(), "recombinationNetwork", networkCopy, "scaling", likelihood.scaling.get(), "data", likelihood.dataInput.get(), "tree", likelihood.treeInput.get(), "siteModel", likelihood.siteModelInput.get(), "branchRateModel", likelihood.branchRateModelInput.get());
		removeOperator = new AddRemoveRecombination();
	   	removeOperator.initByName("alpha", 1.0, "network", networkCopy, "weight", 1.0);
	}
	
    public double getLikelihoodGain(RecombinationNetworkEdge recombinantEdge) {
    	createCopy();
    	logP_full = likelihood.calculateLogP();
    	RecombinationNetworkEdge copyEdge = networkCopy.getEdges().stream().filter(e -> e.ID.equals(recombinantEdge.ID)).findFirst().get();
    	List<RecombinationNetworkEdge> removableEdges = removeOperator.getRemovableEdges();
    	if (removableEdges.contains(copyEdge)) {
    		removeOperator.removeRecombinationEdge(copyEdge);
    		double logP = likelihoodCopy.calculateLogP();
    		networkCopy.restore();
    		return(logP_full-logP);
    	} else {
    		return(Double.POSITIVE_INFINITY);
    	}
    }
    
    public double getMinimumLikelihoodGain(RecombinationNetworkNode recombinantNode) {
    	List<RecombinationNetworkEdge> parentEdges = recombinantNode.getParentEdges();
    	if(parentEdges.size()!=2) {
    		Log.warning.println("The recombinant node " + recombinantNode.toString()+" has not 2 parent edges");
    	}
    	double likelihoodGainLeft = getLikelihoodGain(parentEdges.get(0));
    	double likelihoodGainRight = getLikelihoodGain(parentEdges.get(0));
    	double minLHGain = Math.min(likelihoodGainLeft, likelihoodGainRight);
    	return (minLHGain!=Double.POSITIVE_INFINITY ? minLHGain : 0.0);
    }
    
    public String getAnnotatedExtendedNewick() {
    	// TODO: implement method. Maybe move to another class?
    	return("");
    }
    

}
