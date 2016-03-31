package common;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class Stats {
	
	
	public static void main(String[] args) {
		double p = poisson(1.0,4);
		System.out.println(p);
	}
	
	
	public static double poisson(double mean, int value){
		PoissonDistribution pdist = new PoissonDistribution(mean);
		double probability = pdist.cumulativeProbability(value);
		return(probability);
	}
	

	
	

}
