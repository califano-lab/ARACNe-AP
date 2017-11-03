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

public class DataMatrix {

	public HashMap<String, DataVector> data;
	public ArrayList<String> genes;
	public ArrayList<String> samples;

	
	// Constructor (inherited objects)
	public DataMatrix(HashMap<String, DataVector> data, ArrayList<String> genes, ArrayList<String> samples){
		this.data=data;
		this.genes=genes;
		this.samples=samples;
	}
	
	// Constructor (parses a file) with missing values
	public DataMatrix(File file) throws NumberFormatException, IOException, Exception{
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
		
		data = new HashMap<String, DataVector>();

		// Now reloop on the file to load it
	
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		// Read File Line By Line
		
		genes = new ArrayList<String>();
		samples = new ArrayList<String>();

		int i = 0;
		while ((strLine = br.readLine()) != null) {
			
			String gene = "";
			double[] tempData = new double[columns];
			boolean[] tempNAinfo = new boolean[columns];

			String[] splitter = strLine.split("\t");
			if (i == 0) { // Fill the column names
				for (int j = 1; j<splitter.length; j++) {
					samples.add(splitter[j]);
				}
			} else {
				
				for (int j = 0; j<splitter.length; j++) {
					if (j == 0) { // If this is the first column
						gene=splitter[j];
						genes.add(gene);
					} else {
						if(!splitter[j].toUpperCase().equals("NA") && !splitter[j].toUpperCase().equals("NAN") && !splitter[j].toUpperCase().equals("")){
							tempData[j-1] = Double.parseDouble(splitter[j]);
							tempNAinfo[j-1] = false;
						} else {
							tempData[j-1] = 0;
							tempNAinfo[j-1] = true;
						}
					}
				}
				
				DataVector tempDataVector = new DataVector(tempNAinfo, tempData, gene);
				data.put(gene,tempDataVector);
			}

			i++;
		}
		br.close();
	}
	
	/**
	 * Method to bootstrap NA containing data
	 */
	public void bootstrap(Random random){

		// Which samples to get
		int[] bootsamples = new int[samples.size()];
		for(int i = 0; i<bootsamples.length; i++){
			bootsamples[i] = random.nextInt(samples.size()-1);
		}

		// Generate a bootstrapped data
		for(String gene : genes){
			double[] newdata = new double[samples.size()];
			boolean[] newdataNA = new boolean[samples.size()];
			for(int j = 0; j<samples.size(); j++){
				newdata[j] = data.get(gene).values[bootsamples[j]] + random.nextDouble() / 100000;// WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
				newdataNA[j] = data.get(gene).NAs[bootsamples[j]];// WHITE NOISE AGAIN! FOR BREAKING THE TIES IN THE BOOTSTRAP!
			}
			data.get(gene).values = newdata;
			data.get(gene).NAs = newdataNA;
		}
	}
	

	public HashMap<String,DataVector> getData() {
		return data;
	}

	public ArrayList<String> getSamples() {
		return samples;
	}


	public ArrayList<String> getGenes() {
		return genes;
	}

	
}
