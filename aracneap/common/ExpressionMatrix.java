package common;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class ExpressionMatrix {
	// Variable declaration
	double[][] data;
	boolean[][] naBoolean;
	ArrayList<String> genes;
	ArrayList<String> samples;
	private int[] bootsamples;

	// Constructor (inherited objects)
	public ExpressionMatrix(double[][] data, ArrayList<String> genes, ArrayList<String> samples){
		this.data=data;
		this.genes=genes;
		this.samples=samples;
	}
	
	// Constructor (inherited objects) with the addition of boolean information about NA status of values
	public ExpressionMatrix(double[][] data, boolean[][] naBoolean, ArrayList<String> genes, ArrayList<String> samples){
		this.data=data;
		this.naBoolean=naBoolean;
		this.genes=genes;
		this.samples=samples;
	}

	// Constructor (parses a file)
	public ExpressionMatrix(File file) throws NumberFormatException, IOException, Exception{
		//System.out.println("Loading file "+file);
		// Open the file that is the first
		// command line parameter
		//FileInputStream fstream0 = new FileInputStream(file);
		// Get the object of DataInputStream
		//DataInputStream in0 = new DataInputStream(fstream0);
		BufferedReader br0 = new BufferedReader(new FileReader(file));
		String strLine0;
		// Count number of lines (features, variables) and columns (arrays, samples)
		int rows = -1; // We won't count the first line
		int columns = 0;
		while ((strLine0 = br0.readLine()) != null) {
			if (strLine0.length()>0) { // This statement will skip empty lines.
				rows++;
				if (rows == 0) {
					String[] splitter = strLine0.split("\t");
					for (int j = 1; j<splitter.length; j++) {
						if (splitter[j].trim().length() > 0) {
							columns++;
						}
					}
				}
			}
		}
		br0.close();

		// Now reloop on the file to load it
	
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		// Read File Line By Line
		data = new double[rows][columns];
		genes = new ArrayList<String>();
		samples = new ArrayList<String>();

		int i = 0;
		while ((strLine = br.readLine()) != null) {
			String[] splitter = strLine.split("\t");
			if (i == 0) { // Fill the column names
				for (int j = 1; j<splitter.length; j++) {
					samples.add(splitter[j]);
				}
			} else {
				String gene = null;
				for (int j = 0; j<splitter.length; j++) {
					if (j == 0) { // If this is the first column
						gene=splitter[j];
						genes.add(gene);
					} else {
						if(!splitter[j].equals("NA")){
							data[i-1][j-1] = Double.parseDouble(splitter[j]);
						} else {
							data[i-1][j-1] = 0;
						}
					}
				}
			}
			i++;
		}
		br.close();
	}

	// Constructor (parses a file) with missing values
	public ExpressionMatrix(File file, boolean _na) throws NumberFormatException, IOException, Exception{
		//System.out.println("Loading file "+file);
		// Open the file that is the first
		// command line parameter
		
		
		BufferedReader br0 = new BufferedReader(new FileReader(file));
		String strLine0;
		// Count number of lines (features, variables) and columns (arrays, samples)
		int rows = -1; // We won't count the first line
		int columns = 0;
		while ((strLine0 = br0.readLine()) != null) {
			if (strLine0.length()>0) { // This statement will skip empty lines.
				rows++;
				if (rows == 0) {
					String[] splitter = strLine0.split("\t");
					for (int j = 1; j<splitter.length; j++) {
						if (splitter[j].trim().length() > 0) {
							columns++;
						}
					}
				}
			}
		}
		br0.close();

		// Now reloop on the file to load it
	
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		// Read File Line By Line
		data = new double[rows][columns];
		naBoolean = new boolean[rows][columns];
		genes = new ArrayList<String>();
		samples = new ArrayList<String>();

		int i = 0;
		while ((strLine = br.readLine()) != null) {
			String[] splitter = strLine.split("\t");
			if (i == 0) { // Fill the column names
				for (int j = 1; j<splitter.length; j++) {
					samples.add(splitter[j]);
				}
			} else {
				String gene = null;
				for (int j = 0; j<splitter.length; j++) {
					if (j == 0) { // If this is the first column
						gene=splitter[j];
						genes.add(gene);
					} else {
						if(!splitter[j].toUpperCase().equals("NA") && !splitter[j].toUpperCase().equals("NAN") && !splitter[j].toUpperCase().equals("")){
							data[i-1][j-1] = Double.parseDouble(splitter[j]);
							naBoolean[i-1][j-1] = false;
						} else {
							data[i-1][j-1] = 0;
							naBoolean[i-1][j-1] = true;
						}
					}
				}
			}
			i++;
		}
		br.close();
	}


	
	//// Methods
	
	// Rank with white noise
	public HashMap<String,short[]> rank(Random random){
		HashMap<String,short[]> rankData = new HashMap<String, short[]>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			// Add white noise
			for(int ii = 0; ii<inputVector.length; ii++){
				 inputVector[ii] = inputVector[ii] + random.nextDouble() / 100000;
			}
			
			short[] values = Methods.rankVector(inputVector);
			rankData.put(gene, values);
			i++;
		}
		return(rankData);
	}
	
	// Rank without white nois
	public HashMap<String,short[]> rank(){
		HashMap<String,short[]> rankData = new HashMap<String, short[]>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			short[] values = Methods.rankVector(inputVector);
			rankData.put(gene, values);
			i++;
		}
		return(rankData);
	}


	public ExpressionMatrix bootstrap(Random random){

		// Which samples to get
		int[] bootsamples = new int[samples.size()];
		for(int i = 0; i<bootsamples.length; i++){
			bootsamples[i] = random.nextInt(samples.size()-1);
		}
		this.bootsamples=bootsamples;
		// Generate a bootstrapped data
		double[][] newdata = new double[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				newdata[i][j] = data [i][bootsamples[j]] + random.nextDouble() / 100000;// WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
			}
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, genes, samples);
		return(em);
	}
	
	
	/**
	 * Method to bootstrap NA containing data
	 */
	public ExpressionMatrix bootstrapNA(Random random){

		// Which samples to get
		int[] bootsamples = new int[samples.size()];
		for(int i = 0; i<bootsamples.length; i++){
			bootsamples[i] = random.nextInt(samples.size()-1);
		}
		this.bootsamples=bootsamples;
		// Generate a bootstrapped data
		double[][] newdata = new double[genes.size()][samples.size()];
		boolean[][] newdataNA = new boolean[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				newdata[i][j] = data[i][bootsamples[j]] + random.nextDouble() / 100000;// WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
				newdataNA[i][j] = naBoolean[i][bootsamples[j]];// WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
			}
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, newdataNA, genes, samples);
		return(em);
	}

	public int[] getBootsamples() {
		return bootsamples;
	}


	public ExpressionMatrix bootstrap(int[] bootsamples){
		Random noiseGenerator = new Random(); // WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
		// Generate a bootstrapped data
		double[][] newdata = new double[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				newdata[i][j] = data [i][bootsamples[j]] + noiseGenerator.nextDouble() / 100000;
			}		
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, genes, samples);
		return(em);
	}

	/**
	 * get bootstraps with precomputed bootsamples
	 */
	public ExpressionMatrix bootstrapNA(int[] bootsamples){
		Random noiseGenerator = new Random(); // WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
		// Generate a bootstrapped data
		double[][] newdata = new double[genes.size()][samples.size()];
		boolean[][] newdataNA = new boolean[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				newdata[i][j] = data [i][bootsamples[j]] + noiseGenerator.nextDouble() / 100000;
				newdataNA[i][j] = naBoolean[i][bootsamples[j]];
			}		
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, newdataNA, genes, samples);
		return(em);
	}

	public double[][] getData() {
		return data;
	}

	public ArrayList<String> getSamples() {
		return samples;
	}


	public ArrayList<String> getGenes() {
		return genes;
	}


	// This method generates a new expressionMatrix which is a subset of the current object
	public ExpressionMatrix selectSamples(ArrayList<String> common) {
		double[][] newdata = new double [genes.size()][];
		int igene = 0;
		for (String gene : genes){
			double vector[] = new double[common.size()];
			int i = 0;
			for (String sample : common){
				int j = samples.indexOf(sample);
				vector[i] = data[igene][j];
				i++;
			}
			newdata[igene]=vector;
			igene++;
		}
		ExpressionMatrix newEm = new ExpressionMatrix(newdata, genes, common);
		return(newEm);
	}

	// This method generates a new expressionMatrix which is a subset of the current object supporting NA values
	public ExpressionMatrix selectSamplesNA(ArrayList<String> common) {
		double[][] newdata = new double [genes.size()][];
		boolean[][] newdataNA = new boolean [genes.size()][];
		int igene = 0;
		for (String gene : genes){
			double vector[] = new double[common.size()];
			int i = 0;
			for (String sample : common){
				int j = samples.indexOf(sample);
				vector[i] = data[igene][j];
				i++;
			}
			newdata[igene]=vector;
			igene++;
		}
		ExpressionMatrix newEm = new ExpressionMatrix(newdata, newdataNA, genes, common);
		return(newEm);
	}

}