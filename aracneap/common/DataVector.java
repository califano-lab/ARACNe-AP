package common;

import java.util.ArrayList;

public class DataVector {
	// Variable Declaration
	public boolean [] NAs;
	public double [] values;
	public String geneName;
	
	// Constructor
	public DataVector(boolean [] NAs, double[] values, String geneName){
		this.NAs=NAs;
		this.values=values;
		this.geneName=geneName;
	}
	
	
	// Methods

	/*
	 * This Method provides two vectors (for not NAs samples) and three counts:
	 * Nr. of samples with NAs in both vectors
	 * Nr. of samples with NAs in vector 1
	 * Nr. of samples with NAs in vector 2
	 * 
	*/
	public ArrayList<double[]> getQuadrants(DataVector x){
		ArrayList<double[]> output = new ArrayList<double[]>();

		// Quadrant counts
		double[] bothNAs = new double[1];
		double[] NAs1 = new double[1];
		double[] NAs2 = new double[1];
		bothNAs[0] = NAs1[0] = NAs2[0] = 0;
		
		// First loop does the counting
		int lengthOfNotNAsInBoth = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				if (x.NAs[i]==false){
					lengthOfNotNAsInBoth++;
				} else {
					NAs2[0]++;
				}
			}  else {
				if (x.NAs[i]==false){
					NAs1[0]++;
				} else {
					bothNAs[0]++;
				}
			}
		}
	
		// Second loop generates the not NAs vectors
		double [] notNA1 = new double[lengthOfNotNAsInBoth];
		double [] notNA2 = new double[lengthOfNotNAsInBoth];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				notNA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		output.add(notNA1);
		output.add(notNA2);
		output.add(bothNAs);
		output.add(NAs1);
		output.add(NAs2);

		return output;
	}
	
	
	public double[] getNotNAvalues() {
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				length++;
			}
		}
		
		double [] notNAvalues = new double[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				notNAvalues[j] = values[i];
				j++;		
			}
		}
		return notNAvalues;
	}
	
	public ArrayList<double[]> getNotNAsInBoth(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				length++;
			}
		}
	
		double [] notNA1 = new double[length];
		double [] notNA2 = new double[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				notNA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<double[]> output = new ArrayList<double[]>();
		output.add(notNA1);
		output.add(notNA2);
		return output;
	}

	
	public ArrayList<double[]> getNotNAsInFirst(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==true){
				length++;
			}
		}
	
		double [] notNA1 = new double[length];
		double [] NA2 = new double[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==true){
				notNA1[j] = values[i];
				NA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<double[]> output = new ArrayList<double[]>();
		output.add(notNA1);
		output.add(NA2);
		return output;
	}

	public ArrayList<double[]> getNotNAsInSecond(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==true & x.NAs[i]==false){
				length++;
			}
		}
	
		double [] NA1 = new double[length];
		double [] notNA2 = new double[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==true & x.NAs[i]==false){
				NA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<double[]> output = new ArrayList<double[]>();
		output.add(NA1);
		output.add(notNA2);
		return output;
	}
	
}
