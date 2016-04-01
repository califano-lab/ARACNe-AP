package aracne;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class BootstrapConsolidator {

	private double[] poissonLibrary = new double[100];
	
	private HashMap<String, Integer> edgesOccurrences = new HashMap<String, Integer>();
	private HashMap<String, Double> mi = new HashMap<String, Double>();
	
	private HashSet<String> tfs = new HashSet<String>();
	private HashSet<String> targets = new HashSet<String>();
	
	private int maxCount = 0;

	private boolean nobonferroni;

	
	public BootstrapConsolidator(boolean nobonferroni) {
		this.nobonferroni=nobonferroni;
	}

	public void mergeFiles(File folder) throws IOException{
		FilenameFilter filter = new FilenameFilter() {
		    @Override
		    public boolean accept(File dir, String name) {
		        return name.matches("bootstrapNetwork_.*");
		    }
		};
		File[] listOfFiles = folder.listFiles(filter);
	
		// Calculate the occurrence of each edge (count object)
		// And sum the global mi for each edge
		System.out.println("Integrating "+listOfFiles.length+" bootstraps...");
		for (File file : listOfFiles) {
		    if (file.isFile()) {
		        System.out.println(file.getName());
		        try {
		        	BufferedReader br = new BufferedReader(new FileReader(file));
		        	String line = "";
		        	boolean firstline = true;
		        	while((line = br.readLine()) != null){
		        		if(firstline){
		        			firstline=false;
		        		} else {
			        		String[] sp = line.split("\t");
			        		String key = sp[0]+"#"+sp[1]; // edge key
			        		tfs.add(sp[0]);
			        		targets.add(sp[1]);
			        		if(edgesOccurrences.containsKey(key)){
			        			edgesOccurrences.put(key, edgesOccurrences.get(key)+1);
			        			mi.put(key, mi.get(key)+Double.parseDouble(sp[2]));
			        		}
			        		else{
			        			edgesOccurrences.put(key, 1);
			        			mi.put(key, Double.parseDouble(sp[2]));
			        		}
		        		}
		        	}
		        	br.close();
		        }
		        catch(Exception e){
		        	e.printStackTrace();
		        }
		    }
		}
		
		// Initialize the poisson distribution based on the average number of edge appearance
		long allCount = 0;
		//FileWriter fw = new FileWriter(new File("file.txt"));
		for(String key : edgesOccurrences.keySet()){
			// Calculate the average mi of each edge
			int edgeOccurrence=edgesOccurrences.get(key);
			mi.put(key, mi.get(key)/edgeOccurrence);
			allCount += edgeOccurrence; // Total number of edges in all the bootstrap files
			if(edgeOccurrence > maxCount){
				maxCount = edgeOccurrence; // Edge with the most count
			}
			//fw.write(edgeOccurrence+"\n");
		}
		//fw.close();
		// Lambda
		//long allObservedEdges = edgesOccurrences.size();
		//double meanEdgeNumber = allCount/allObservedEdges; // Wrong way: observed edges in the denominator
		double meanEdgeNumber = (allCount*1.0)/(targets.size()*tfs.size()); // Right way: all possible edges
		generatePoissonPvalues(meanEdgeNumber);
	}
	
	public void writeSignificant(String outFile, double poissonPvalue){
		int totalCount = tfs.size()*targets.size();
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outFile)));

			// Header
			bw.write("Regulator\tTarget\tMI\tpvalue\n");

			String[] keys = edgesOccurrences.keySet().toArray(new String[0]);
			for(String key : keys){
				int occurrence = edgesOccurrences.get(key);
				// Bonferroni correction
				// In the C++ version, the correction was done using the observed number of edges
				// Here, the correction is done based on ALL the possible combinations
				double thisPvalue;
				if(nobonferroni) {
					thisPvalue = poissonLibrary[occurrence];
				} else {
					thisPvalue = poissonLibrary[occurrence]*totalCount;	
				}
				
				if(thisPvalue>1){
					thisPvalue=1;
				}
				if(thisPvalue < poissonPvalue){
					String[] sp = key.split("#");
					bw.write(sp[0]+"\t"+sp[1]+"\t"+mi.get(key)+"\t"+thisPvalue+"\n");
				}
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void generatePoissonPvalues(double mean){
		poissonLibrary = new double[maxCount+1];
		PoissonDistribution pdist = new PoissonDistribution(mean);
		for(int i=0; i<maxCount+1; i++){
			poissonLibrary[i] = 1.0-pdist.cumulativeProbability(i);
		}
	}
	
}








