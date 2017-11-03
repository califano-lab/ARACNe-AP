package common;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeMap;

/**
 * Class defining static methods to load specific data formats
 */
public class DataParser {
	static ArrayList<String> genes;
	private int[] bootRandomSamples = null;
	private int sampleNumber = 0;

	public String[] getGenes() {
		String[] genearray = new String[genes.size()];
		genearray = genes.toArray(genearray);
		return genearray;
	}
	/**
	 * Method to read a file with a single number
	 * @throws IOException 
	 */
	public static Double readValue(File inFile) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		String string = br.readLine();
		Double value = Double.valueOf(string);
		br.close();
		return(value);
	}
	/**
	 * Method to write a single value into a file
	 */
	public static void writeValue(Double value, File outFile) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		bw.write(String.valueOf(value));
		bw.close();
	}

	/**
	 * Method to read a gene list
	 */
	public static HashSet<String> readGeneSet(File genelistFile){
		try{
			BufferedReader br = new BufferedReader(new FileReader(genelistFile));
			String l = "";
			HashSet<String> tfset = new HashSet<String>();
			while((l = br.readLine()) != null){
				tfset.add(l);
			}
			br.close();
			return(tfset);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * Method to write a two columns array into a file
	 * If several identical MIs are present, the max CMI will be selected
	 * @throws IOException 
	 */
	public static void writeCMIthresholds(double[][] cmiThresholds, File outFile) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		TreeMap<Double,Double> treshes = new TreeMap<Double,Double>();
		for (int i = 0; i < cmiThresholds.length; i++) {
			double mi = cmiThresholds[i][0];
			double cmi = cmiThresholds[i][1];
			if(treshes.containsKey(mi)){
				if(cmi>treshes.get(mi)){
					treshes.put(mi, cmi);
				}
			} else {
				treshes.put(mi, cmi);
			}
		}
		
		for(Double mi : treshes.keySet()){
			String str = mi + "\t" + treshes.get(mi) + "\n";
			bw.write(str);
		}
		bw.close();
	}

	/**
	 * Method to read a two columns array from a file
	 * @throws IOException 
	 */
	public static double[][] readCMIthresholds(File inFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		String line = "";
		int nrlines = 0;
		while((line = br.readLine()) != null){
			nrlines++;
		}
		br.close();

		double[][]output = new double[nrlines][2];
		BufferedReader br2 = new BufferedReader(new FileReader(inFile));
		int i = 0;
		while((line = br2.readLine()) != null){
			String[] splittedLine = line.split("\t");
			output[i][0]=Double.parseDouble(splittedLine[0]);
			output[i][1]=Double.parseDouble(splittedLine[1]);
			i++;
		}
		br2.close();
		return output;
	}
	
	
	/**
	 * Method to read an integer list
	 */
	public static ArrayList<Integer> readIntegerArray(File inFile) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(inFile));
			String l = "";
			ArrayList<Integer> triplets = new ArrayList<Integer>();
			while((l = br.readLine()) != null){
				triplets.add(Integer.parseInt(l));
			}
			br.close();
			return(triplets);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}

	
	/**
	 * Method to write an integer list
	 * @throws IOException 
	 */
	public static void writeIntegerList(ArrayList<Integer> input,  File outFile) throws IOException{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			for(Integer i : input){
				String str = String.valueOf(i)+"\n";
				bw.write(str);
			}
			bw.close();
	}

	/**
	 * Method to write a string list
	 * @throws IOException 
	 */
	public static void writeStringArray(String[] input, File outFile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		for(int i = 0; i<input.length; i++){
			String str = input[i]+"\n";
			bw.write(str);
		}
		bw.close();
	}
	
	/**
	 * Method to write an integer list
	 * @throws IOException 
	 */
	public static void writeArraylist(ArrayList<String> input, File outFile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		for(int i = 0; i<input.size(); i++){
			String str = input.get(i)+"\n";
			bw.write(str);
		}
		bw.close();
	}
	
	
	/**
	 * Method to get the last line of a file
	 */
	public static String tail( File file ) {
	    RandomAccessFile fileHandler = null;
	    try {
	        fileHandler = new RandomAccessFile( file, "r" );
	        long fileLength = fileHandler.length() - 1;
	        StringBuilder sb = new StringBuilder();

	        for(long filePointer = fileLength; filePointer != -1; filePointer--){
	            fileHandler.seek( filePointer );
	            int readByte = fileHandler.readByte();

	            if( readByte == 0xA ) {
	                if( filePointer == fileLength ) {
	                    continue;
	                }
	                break;

	            } else if( readByte == 0xD ) {
	                if( filePointer == fileLength - 1 ) {
	                    continue;
	                }
	                break;
	            }

	            sb.append( ( char ) readByte );
	        }

	        String lastLine = sb.reverse().toString();
	        return lastLine;
	    } catch( java.io.FileNotFoundException e ) {
	        e.printStackTrace();
	        return null;
	    } catch( java.io.IOException e ) {
	        e.printStackTrace();
	        return null;
	    } finally {
	        if (fileHandler != null )
	            try {
	                fileHandler.close();
	            } catch (IOException e) {
	                /* ignore */
	            }
	    }
	}

	
}