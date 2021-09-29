package recombination.util;


import beast.core.Input;
import beast.core.Runnable;
import beast.evolution.sitemodel.SiteModel;
import recombination.network.RecombinationNetwork;

public class NetworkAlignmentSimulator extends Runnable {

	
    public Input<RecombinationNetwork> recombinationNetworkInput = new Input<>(
            "recombinationNetwork",
            "Recombination network down which to simulate evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "site model for leafs in the beast.tree",
            Input.Validate.REQUIRED);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");

    public Input<Boolean> useNexusInput = new Input<>(
            "useNexus",
            "Use Nexus instead of FASTA format to write alignment file.",
            false);

    @Override
    public void initAndValidate() { }

    @Override
    public void run() throws Exception {
        SimulatedNetworkAlignment alignment = new SimulatedNetworkAlignment();
        alignment.initByName(
                "recombinationNetwork", recombinationNetworkInput.get(),
                "siteModel", siteModelInput.get(),
                "outputFileName", outputFileNameInput.get());
    }

}
