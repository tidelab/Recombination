package recombination.operator;

import beast.util.Randomizer;
import coalre.CoalReTestClass;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkNode;
import recombination.operators.DivertLociOperator;

import org.junit.Assert;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class DivertLociTest extends CoalReTestClass {

    String networkString = "((((#H1[&split={7357-9999},loci={7357-9999},length=2643.0]:0.029028801150026373,(((t7[&loci={0-9999},length=10000.0]:0.3408331510142315,(t8[&loci={0-9999},length=10000.0]:0.14470426889934196,((t4[&loci={0-9999},length=10000.0]:0.008921363131720263,t3[&loci={0-9999},length=10000.0]:0.008921363131720263)[&loci={0-9999},length=10000.0]:0.07368958619017191,(t5[&loci={0-9999},length=10000.0]:0.05651985106014486,((t0[&loci={0-9999},length=10000.0]:0.01371370813629106,t2[&loci={0-9999},length=10000.0]:0.01371370813629106)[&loci={0-9999},length=10000.0]:0.01770276986295235,t1[&loci={0-9999},length=10000.0]:0.03141647799924341)[&loci={0-9999},length=10000.0]:0.02510337306090145)[&loci={0-9999},length=10000.0]:0.02609109826174731)[&loci={0-9999},length=10000.0]:0.06209331957744979)[&loci={0-9999},length=10000.0]:0.19612888211488955)[&loci={0-9999},length=10000.0]:0.2532640786850209)#H2[&split={499-9999},loci={499-9999},length=9501.0]:0.2994016137375858)#H1[&split={0-7356},loci={499-7356},length=6858.0]:0.029028801150026373)[&loci={499-9999},length=9501.0]:0.34436738009841983)#H0[&split={3613-9999},loci={3613-9999},length=6387.0]:0.16198962047105492,#H0[&split={0-3612},loci={499-3612},length=3114.0]:0.16198962047105492)[&loci={499-9999},length=9501.0]:0.2105449461319744,(#H3[&split={0-397},loci={0-397},length=398.0]:0.3596040657302466,(#H2[&split={0-498},loci={0-498},length=499.0]:0.36308767696719646,((t6[&loci={0-9999},length=10000.0]:0.08009733219216959,t9[&loci={0-9999},length=10000.0]:0.08009733219216959)[&loci={0-9999},length=10000.0]:0.703681773516321)#H3[&split={398-9999},loci={398-9999},length=9602.0]:0.1734058009579582)[&loci={0-9999},length=10000.0]:0.1861982647722884)[&loci={0-9999},length=10000.0]:0.4960464198495764)[&loci={0-9999},length=10000.0]:0.0;";
    String networkString2 = "((((((#H2[&split={9240-9999},loci={9240-9999},length=760.0]:0.3579693014979668,(#H3[&split={6573-9999},loci={7302-7756},length=455.0]:0.4198547464044533,#H4[&split={0-7717},loci={6648-7717},length=1070.0]:0.25563256204649165)[&loci={6648-7756},length=1109.0]:1.1619129921459237)[&loci={6648-7756,9240-9999},length=1869.0]:0.34286205275424475)#H1[&split={7372-9999},loci={7372-7756,9240-9999},length=1145.0]:0.23453108442206716,((((((#H7[&split={3116-9999},loci={3116-6068},length=2953.0]:0.04387002757521419,((#H8[&split={0-990},loci={0-902},length=903.0]:0.20086820188684906,((((#H10[&split={6130-9999},loci={6130-7425},length=1296.0]:0.021261671248241076,(((#H12[&split={0-4308},loci={1942-4308},length=2367.0]:0.017015062714914875)#H11[&split={0-3255},loci={1942-3255},length=1314.0]:0.09607673081432438,(((((((t7[&loci={0-9999},length=10000.0]:0.19678534848648216,(t3[&loci={0-9999},length=10000.0]:0.05322159367975174,t6[&loci={0-9999},length=10000.0]:0.05322159367975174)[&loci={0-9999},length=10000.0]:0.14356375480673042)[&loci={0-9999},length=10000.0]:0.05151361444577246)#H16[&split={2275-9999},loci={2275-9999},length=7725.0]:0.07715170782013003)#H15[&split={5532-9999},loci={5532-9999},length=4468.0]:0.018553890999346345,(((t4[&loci={0-9999},length=10000.0]:0.17649068406003138)#H17[&split={0-7301},loci={0-7301},length=7302.0]:0.04039091402995995,(#H18[&split={0-1941},loci={0-1941},length=1942.0]:0.0013853304961968516,(t9[&loci={0-9999},length=10000.0]:0.1857901138507514)#H19[&split={0-7244},loci={0-7244},length=7245.0]:0.015775281812972886)[&loci={0-7244},length=7245.0]:0.01531620242626705)[&loci={0-7301},length=7302.0]:0.10731442948560113,#H20[&split={8739-9999},loci={8739-9999},length=1261.0]:0.16478852729506854)[&loci={0-7301,8739-9999},length=8563.0]:0.019808534176138537)[&loci={0-9999},length=10000.0]:0.06666315890844926,#H21[&split={8920-9999},loci={8920-9999},length=1080.0]:0.017738595522330758)[&loci={0-9999},length=10000.0]:0.01401677843270166)#H14[&split={0-7425},loci={0-7425},length=7426.0]:0.028992203555866236)#H13[&split={2098-9999},loci={2098-7425},length=5328.0]:0.07487331027241306)[&loci={1942-7425},length=5484.0]:0.32310917034218034)#H10[&split={0-6129},loci={1942-6129},length=4188.0]:0.021261671248241076)[&loci={1942-7425},length=5484.0]:0.033310006810186255,#H22[&split={0-6338},loci={903-2097,3256-4308},length=2248.0]:0.04657187850738409)[&loci={903-7425},length=6523.0]:0.13453814563512445,((((((#H19[&split={7245-9999},loci={7245-9999},length=2755.0]:0.0146711225873781,(#H23[&split={6330-9999},loci={6330-9999},length=3670.0]:0.09734477025221858)#H20[&split={0-8738},loci={6330-8738},length=2409.0]:0.041053736157605575)[&loci={6330-9999},length=3670.0]:0.2719204176664495,#H14[&split={7426-9999},loci={7426-9999},length=2574.0]:0.04769715501169708)[&loci={6330-9999},length=3670.0]:0.19505659579144785,(((((((t5[&loci={0-9999},length=10000.0]:0.04149497715939229,t1[&loci={0-9999},length=10000.0]:0.04149497715939229)[&loci={0-9999},length=10000.0]:0.02347000940186421,((t2[&loci={0-9999},length=10000.0]:0.0059770132872838345)#H26[&split={0-7028},loci={0-7028},length=7029.0]:0.005486880263416766)#H25[&split={1773-9999},loci={1773-7028},length=5256.0]:0.0535010930105559)[&loci={0-9999},length=10000.0]:0.11378162648069744,(t8[&loci={0-9999},length=10000.0]:0.06206273002830534)#H23[&split={0-6329},loci={0-6329},length=6330.0]:0.1166838830136486)[&loci={0-9999},length=10000.0]:0.1884032938158804,#H16[&split={0-2274},loci={0-2274},length=2275.0]:0.11885094392557971)[&loci={0-9999},length=10000.0]:0.02577921828001517)#H21[&split={0-8919},loci={0-8919},length=8920.0]:0.2317721230555393,((((t0[&loci={0-9999},length=10000.0]:0.20018006516752743)#H18[&split={1942-9999},loci={1942-9999},length=8058.0]:0.15685647892100132)#H28[&split={0-9487},loci={1942-9487},length=7546.0]:0.058421675303393206)#H12[&split={4309-9999},loci={4309-9487},length=5179.0]:0.14856286090476079)#H27[&split={5095-9999},loci={5095-9487},length=4393.0]:0.06068016789670605)[&loci={0-9487},length=9488.0]:0.024804534978505877)#H24[&split={0-6647},loci={0-6647},length=6648.0]:0.01793246672413218)[&loci={0-9999},length=10000.0]:0.037937854132429116,#H29[&split={9676-9999},loci={9676-9999},length=324.0]:0.012880858126575134)[&loci={0-9999},length=10000.0]:0.0021093878245596054,(((#H25[&split={0-1772},loci={0-1772},length=1773.0]:0.16650906911031738,#H26[&split={7029-9999},loci={7029-9999},length=2971.0]:0.17199594937373414)[&loci={0-1772,7029-9999},length=4744.0]:0.375668936340345)#H31[&split={0-8695},loci={0-1772,7029-8695},length=3440.0]:0.1500478922241023)#H30[&split={0-8483},loci={0-1772,7029-8483},length=3228.0]:0.0037957006275503025)[&loci={0-9999},length=10000.0]:0.01793577252195855,#H30[&split={8484-9999},loci={8484-8695},length=212.0]:0.02173147314950885)[&loci={0-9999},length=10000.0]:0.3153477425819192)[&loci={0-9999},length=10000.0]:0.13327931242701974)#H9[&split={0-6068},loci={0-6068},length=6069.0]:0.5322482481945647)[&loci={0-6068},length=6069.0]:0.480143525875959)#H7[&split={0-3115},loci={0-3115},length=3116.0]:0.04387002757521419)[&loci={0-6068},length=6069.0]:0.379855344642841,(((((((((#H33[&split={0-902},loci={0-902},length=903.0]:0.30734918020033053,#H34[&split={6993-9999},loci={7302-7756},length=455.0]:0.19488281544390607)[&loci={0-902,7302-7756},length=1358.0]:0.06246296606424706)#H3[&split={0-6572},loci={0-902},length=903.0]:0.16521924978447866,((((#H37[&split={4696-9999},loci={9488-9999},length=512.0]:0.03355267416335239,((#H27[&split={0-5094},loci={4309-5094},length=786.0]:0.12142467648671407,#H31[&split={8696-9999},loci={8696-9999},length=1304.0]:0.13180385778203385)[&loci={4309-5094,8696-9999},length=2090.0]:0.03215337284480024,((((#H15[&split={0-5531},loci={2275-5531},length=3257.0]:0.050861127364038916,#H17[&split={7302-9999},loci={7302-9999},length=2698.0]:0.1998211140563922)[&loci={2275-5531,7302-9999},length=5955.0]:0.1844491790963526)#H38[&split={0-7314},loci={2275-5531,7302-7314},length=3270.0]:0.1107829578241889,#H38[&split={7315-9999},loci={7315-9999},length=2685.0]:0.1107829578241889)[&loci={2275-5531,7302-9999},length=5955.0]:0.02095131086491575)#H29[&split={0-9675},loci={2275-5531,7302-9675},length=5631.0]:0.025103883726316223)[&loci={2275-5531,7302-9999},length=5955.0]:0.039060252127992356)[&loci={2275-5531,7302-9999},length=5955.0]:0.16350514165220453)#H36[&split={0-7756},loci={2275-5531,7302-7756},length=3712.0]:0.014311106406990404)#H34[&split={0-6992},loci={2275-5531},length=3257.0]:0.007412004015289075)#H35[&split={3294-9999},loci={3294-5531},length=2238.0]:0.4151530272773427)[&loci={0-902,3294-5531},length=3141.0]:0.046807884912948294,((#H36[&split={7757-9999},loci={7757-9999},length=2243.0]:0.1231585867741507,((((((#H28[&split={9488-9999},loci={9488-9999},length=512.0]:0.11833012494094669,#H11[&split={3256-9999},loci={3256-4308},length=1053.0]:0.04289338692263861)[&loci={3256-4308,9488-9999},length=1565.0]:0.046333283917160184,#H13[&split={0-2097},loci={0-2097},length=2098.0]:0.06802325029788747)[&loci={0-2097,3256-4308,9488-9999},length=3663.0]:0.2014067546462014)#H37[&split={0-4695},loci={0-2097,3256-4308},length=3151.0]:0.0651053177459967,#H24[&split={6648-9999},loci={6648-9487},length=2840.0]:0.13870624216693905)[&loci={0-2097,3256-4308,6648-9487},length=5991.0]:0.03379723972012616)#H33[&split={903-9999},loci={903-2097,3256-4308,6648-9487},length=5088.0]:0.037649717755424916)#H22[&split={6339-9999},loci={6648-9487},length=2840.0]:0.18366412736815985)[&loci={6648-9999},length=3352.0]:0.3127204854989545)#H4[&split={7718-9999},loci={7718-9999},length=2282.0]:0.0478049503394653)[&loci={0-902,3294-5531,7718-9999},length=5423.0]:0.10157981967066432)#H8[&split={991-9999},loci={3294-5531,7718-9999},length=4520.0]:0.11373471539674518,#H9[&split={6069-9999},loci={6069-9999},length=3931.0]:0.44511476170446085)[&loci={3294-5531,6069-9999},length=6169.0]:0.012728616236560075)#H32[&split={6233-9999},loci={6233-9999},length=3767.0]:0.5280998221021458,(#H32[&split={0-6232},loci={3294-5531,6069-6232},length=2402.0]:0.14509037171222605,#H35[&split={0-3293},loci={2275-3293},length=1019.0]:0.8350944352064866)[&loci={2275-5531,6069-6232},length=3421.0]:0.38300945038991974)[&loci={2275-5531,6069-9999},length=7188.0]:0.2556283289488679)#H2[&split={0-9239},loci={2275-5531,6069-9239},length=6428.0]:0.1945456172965443)[&loci={0-9239},length=9240.0]:0.3456122897294751)#H6[&split={0-8604},loci={0-8604},length=8605.0]:0.06276462209584555)#H5[&split={0-7661},loci={0-7661},length=7662.0]:0.03841238216566989,#H6[&split={8605-9999},loci={8605-9239},length=635.0]:0.10117700426151544)[&loci={0-7661,8605-9239},length=8297.0]:0.1841877682963844,#H1[&split={0-7371},loci={6648-7371},length=724.0]:0.12469132533170768)[&loci={0-7661,8605-9239},length=8297.0]:0.10983975909035948)[&loci={0-7756,8605-9999},length=9152.0]:0.015448252837242826)#H0[&split={0-6768},loci={0-6768},length=6769.0]:0.1191130035471124,#H0[&split={6769-9999},loci={6769-7756,8605-9999},length=2383.0]:0.1191130035471124)[&loci={0-7756,8605-9999},length=9152.0]:0.13688107726724486,#H5[&split={7662-9999},loci={7662-8604},length=943.0]:0.6038822432040138)[&loci={0-9999},length=10000.0]:0.0;";
    

    @Test
    public void testAddRemoveSegment() {

        RecombinationNetwork network = new RecombinationNetwork(networkString2);
        RecombinationNetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(2);

        DivertLociOperator operator = new DivertLociOperator();
        
        List<Integer> listToAdd = new ArrayList<Integer>();
        BreakPoints lociToAdd = new BreakPoints();
        
        listToAdd.add(500);
        listToAdd.add(800);
        lociToAdd.init(listToAdd);
        
        
        String net1 = network.toString();
        
        double logPremove = operator.removeLociFromAncestors(
              leafNode.getParentEdges().get(0), lociToAdd);
        

        double logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
                
        Assert.assertEquals(network.toString(), net1);
    }

//    @Test
//    public void testAddRemoveReassortmentEdge() {
//        Network network = getContempNetwork(2, 8, 0.0);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("alpha", 1.0,
//                "network", network,
//                "weight", 1.0);
//
//        NetworkNode origRoot = network.getRootEdge().childNode;
//
//        NetworkEdge sourceEdge = network.getRootEdge().childNode.getChildEdges().get(0);
//        double sourceTime = sourceEdge.getLength()/2.0;
//        NetworkEdge destEdge = network.getRootEdge();
//        double destTime = destEdge.childNode.getHeight() + 1.0;
//
//        double logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
//
//        NetworkEdge edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0).parentNode == origRoot
//                ? sourceEdge.parentNode.getParentEdges().get(1)
//                : sourceEdge.parentNode.getParentEdges().get(0);
//
//        double logP2 = operator.removeReassortmentEdge(edgeToRemove);
//
//        Assert.assertEquals(logP1, -logP2, 1e-10);
//
//        sourceEdge = network.getRootEdge().childNode.getChildEdges().get(1);
//        sourceTime = sourceEdge.getLength()/4.0;
//        destEdge = sourceEdge;
//        destTime = sourceEdge.getLength()*3.0/4.0;
//
//        logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
//
//        edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0);
//
//        logP2 = operator.removeReassortmentEdge(edgeToRemove);
//
//        Assert.assertEquals(logP1, -logP2, 1e-10);
//    }
//
//    @Test
//    public void testRemoveReassortment() {
//        // TODO Flesh out this test
//
//        Network network = new Network(networkString);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("network", network,
//                "alpha", 1.0,
//                "weight", 1.0);
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        operator.removeReassortment();
//
//        System.out.println(network.getExtendedNewickVerbose());
//    }
//
//    @Test
//    public void testAddReassortment() {
//        // TODO Flesh out this test
//
//        Network network = new Network(networkString);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("network", network,
//                "alpha", 1.0,
//                "weight", 1.0);
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        double logHR = operator.addReassortment();
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        System.out.println(logHR);
//    }
}