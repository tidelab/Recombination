package recombination.network;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

import org.junit.Assert;

public class SerializationDeserializationTest {
	
	String networkString = "((#H0[&loci={0-0},pr={0-0},length=1]:0.4860702162314561,#H1[&loci={0-0},pr={0-0},length=1]:1.5817575814045295)[&loci={0-0},length=1]:0.35688825595463936,((#H2[&loci={0-0},pr={0-0},length=1]:0.6042185323591791,(((t2[&loci={0-3},length=4]:0.4068742545918358)#H2[&loci={1-3},pr={1-3},length=3]:0.036380107508342974,(t1[&loci={0-3},length=4]:0.2904108541857302,t0[&loci={0-3},length=4]:0.2904108541857302)[&loci={0-3},length=4]:0.15284350791444856)[&loci={0-3},length=4]:0.47005038739308835)#H1[&loci={1-3},pr={1-3},length=3]:0.09778803745774778)[&loci={0-3},length=4]:0.9978993277153256)#H0[&loci={1-3},pr={1-3},length=3]:0.8429584721860954)[&loci={0-3},length=4]:0.0;";


	@Test
	public void extendedNewickTest() {
		
		RecombinationNetwork network = new RecombinationNetwork(networkString);
		assertEquals(networkString, network.getExtendedNewick());
		
		network = new RecombinationNetwork();
		network.initByName("extendedNewick", networkString);
		assertEquals(networkString, network.getExtendedNewick());
	}
}
