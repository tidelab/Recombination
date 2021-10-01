package recombination.util;

import beast.core.BEASTObject;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import beast.evolution.tree.Tree;
import java.io.PrintStream;

public class MAPnetLogger extends RecombinationNetwork {

    public Input<RecombinationNetwork> networkInput = new Input<>(
            "recombinationNetwork",
            "recombination network state to maximize  posterior",
            Validate.REQUIRED);


    public Input<Distribution> posteriorInput = new Input<>(
            "posterior",
            "Posterior used to identify MAP network",
            Validate.REQUIRED);

    RecombinationNetwork currentMAPNetwork;
    double maxPosterior;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        String currentMAPNetworkNewick = networkInput.get().toString();
        currentMAPNetwork = new RecombinationNetwork(currentMAPNetworkNewick);
        currentMAPNetwork.initAndValidate();
        maxPosterior = Double.NEGATIVE_INFINITY;
    }

    @Override
    public void init(PrintStream out) {
        currentMAPNetwork.init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {
        if (posteriorInput.get().getCurrentLogP()>maxPosterior) {
            maxPosterior = posteriorInput.get().getCurrentLogP();
            currentMAPNetwork.assignFrom(networkInput.get());
        }
        currentMAPNetwork.log(nSample, out);
    }

    @Override
    public void close(PrintStream out) {
        currentMAPNetwork.close(out);
    }

}
