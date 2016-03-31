package paracne;

import jargs.gnu.CmdLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;

import aracne.BootstrapConsolidator;
import aracne.DPI;
import aracne.MI;

import common.DataParser;
import common.ExpressionMatrix;


public class Paracne {
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static float simCut = 0; // In the original paper this was a parameter. E.g. if two TFs have a very high MI, DPI is not calculated
	final static int geneNumberInMIthresholding = 1000000;
	static Random random = new Random();
	private static boolean singlemode = false;

	// Main Method
	public static void main(String[] args) throws Exception {
		Locale.setDefault(Locale.US);

		//// Parse arguments
		CmdLineParser parser = new CmdLineParser();
		CmdLineParser.Option optionExpressionFile1 = parser.addStringOption('e', "expfile_upstream");
		CmdLineParser.Option optionTranscriptionFactorsFile = parser.addStringOption('t', "tfs");
		CmdLineParser.Option optionOutputFolder = parser.addStringOption('o', "output");
		CmdLineParser.Option optionConsolidate = parser.addBooleanOption('c',"consolidate");
		CmdLineParser.Option optionCalculateThreshold = parser.addBooleanOption('t',"calculateThreshold");
		CmdLineParser.Option optionPvalue = parser.addDoubleOption('p',"pvalue");
		CmdLineParser.Option optionSeed = parser.addIntegerOption('s',"seed");
		CmdLineParser.Option optionThreaded = parser.addIntegerOption('m',"threads");
		CmdLineParser.Option optionNodpi = parser.addBooleanOption('n',"nodpi");

		try {
			parser.parse(args);
		}
		catch ( CmdLineParser.OptionException e ) {
			System.err.println(e.getMessage());
			System.exit(2);
		}

		File outputFolder = new File((String)parser.getOptionValue(optionOutputFolder));


		boolean isConsolidate = false;
		if ((Boolean)parser.getOptionValue(optionConsolidate)!=null) {
			isConsolidate = true;
		}

		boolean isThreshold = false;
		if ((Boolean)parser.getOptionValue(optionCalculateThreshold)!=null) {
			isThreshold = true;
		}

		boolean noDPI = false;
		if ((Boolean)parser.getOptionValue(optionNodpi)!=null) {
			noDPI = true;
		}


		// Here the program forks
		// You can calculate the MI threshold
		// You can run a single bootstrap
		// You can consolidate bootstraps
		if(noDPI){
			String exp1File = (String)parser.getOptionValue(optionExpressionFile1);
			File expressionFile1 = new File(exp1File);
			File transcriptionFactorsFile = new File((String)parser.getOptionValue(optionTranscriptionFactorsFile));

			Integer threadCount = (Integer)parser.getOptionValue(optionThreaded);

			singlemode = true;

			if(threadCount == null){
				threadCount = 1;
			}
			
			Double miPvalue = new Double((Double) parser.getOptionValue(optionPvalue));
			runNoDPI(expressionFile1,transcriptionFactorsFile,outputFolder,miPvalue,threadCount);
			System.exit(0);
		}

		if(isThreshold){
			File expressionFile = new File((String)parser.getOptionValue(optionExpressionFile1));
			outputFolder.mkdir();
			Double miPvalue = new Double((Double) parser.getOptionValue(optionPvalue));
			Integer seed = (Integer)parser.getOptionValue(optionSeed);
			if(seed!=null){
				random = new Random(seed);
			} else{
				random = new Random();
			}
			runThreshold(expressionFile,outputFolder,miPvalue,seed);
		}
		else if(!isConsolidate){
			String exp1File = (String)parser.getOptionValue(optionExpressionFile1);
			File expressionFile1 = new File(exp1File);
			
			File transcriptionFactorsFile = new File((String)parser.getOptionValue(optionTranscriptionFactorsFile));
			Integer seed = (Integer)parser.getOptionValue(optionSeed);
			Integer threadCount = (Integer)parser.getOptionValue(optionThreaded);

			singlemode = true;

			if(threadCount == null){
				threadCount = 1;
			}

			if(seed!=null){
				random = new Random(seed);
			} else{
				random = new Random();
			}
			String processId = new BigInteger(130, random).toString(32);
			Double miPvalue = new Double((Double) parser.getOptionValue(optionPvalue));
			runBootsrap(expressionFile1,transcriptionFactorsFile,outputFolder,processId,miPvalue,threadCount,singlemode);
		} else {
			runConsolidate(outputFolder);
		}
	}


	// Calculate Threshold mode
	private static void runThreshold(File _expressionFile, File _outputFolder, double _miPvalue, int _seed) throws NumberFormatException, Exception{
		// Read expression matrix and transcription factor lists
		DataMatrix em = new DataMatrix(_expressionFile);
		HashMap<String, DataVector> data = em.data;
		
		// Check if the sample size
		int sampleNumber = em.samples.size();
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}
		
		//// Calculate threshold for the required p-value
		// Don't if a threshold file already exists
		File miThresholdFile = new File(_outputFolder+"/miThreshold_p"+formatter.format(_miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;

		if(miThresholdFile.exists()){
			System.out.println("MI threshold file was already there, but I am recalculating it.");
		}
		MI miCPU = new MI(em);
		
		System.out.println(sampleNumber);
		miThreshold = miCPU.calibrateMIThresholdNA(data,geneNumberInMIthresholding,_miPvalue,_seed);
		DataParser.writeValue(miThreshold, miThresholdFile);
	}

	// Single bootstrap mode
	private static void runBootsrap(
			File expressionFile1,
			File transcriptionFactorsFile,
			File outputFolder, 
			String processId,
			Double miPvalue,
			Integer threadCount, boolean singlemode
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();
		// Read expression matrices and bootstrap them
		DataMatrix em1 = new DataMatrix(expressionFile1);
		System.out.println("Bootstrapping input matrix 1 with "+em1.getGenes().size()+" genes and "+em1.getSamples().size()+" samples");
		
		//compute a random selection of samples to build bootstrapped data
		em1.bootstrap(random);

		HashMap<String, DataVector> data = em1.data;
		
		// Check if the sample size is less than the short limit
		int sampleNumber = em1.getSamples().size();
		
		if(sampleNumber > 32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		// TF list
		HashSet<String> tfSet = DataParser.readGeneSet(transcriptionFactorsFile);
		String[] tfList = tfSet.toArray(new String[0]);
		Arrays.sort(tfList);

		if(tfList.length==0){
			System.err.println("The transcription factor file is badly formatted or empty");
			System.exit(1);
		}

		// Check if the threshold file exists
		File miThresholdFile = new File(outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;
		if(!miThresholdFile.exists()){
			System.err.println("MI threshold file is not present.");
			System.err.println("Please run ARACNE in --calculateThreshold mode first");
			System.exit(1);
		}
		System.out.println("MI threshold file is present");
		miThreshold = DataParser.readValue(miThresholdFile);
		//miThreshold = 0;
		// Calculate a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		
		System.out.println("Calculate network from: "+expressionFile1);
		
		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		MI miCPU = new MI(data,tfList,miThreshold,threadCount);
		finalNetwork = miCPU.getFinalNetwork();
		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then calculate DPI
		long time2 = System.currentTimeMillis();

		HashMap<String, HashSet<String>> removedEdges = dpi(finalNetwork, threadCount);
		System.out.println("DPI time elapsed: "+(System.currentTimeMillis() - time2)/1000+" sec");

		// And write out the single bootstrap network
		File outputFile = new File(outputFolder.getAbsolutePath()+"/bootstrapNetwork_"+processId+".txt");
		writeFinal(outputFile,removedEdges,finalNetwork);

		long finalTime = System.currentTimeMillis();
		System.out.println("Total time elapsed: "+(finalTime - initialTime)/1000+" sec");
	}

	// This method consolidates the several bootstraps
	private static void runConsolidate(File outputFolder) throws IOException {
		BootstrapConsolidator c = new BootstrapConsolidator(false);
		c.mergeFiles(outputFolder);

		// P-value for the Poisson distribution. Aka how many times an edge has to appear in the bootstraps to be kept.
		// Hard-coded to 0.3 in the original ARACNe
		double pvalue = 0.05; 

		c.writeSignificant(outputFolder+"/finalNetwork_4col.tsv", pvalue);

		System.out.println("\n        :");
		System.out.println("       :");
		System.out.println("        :       Cool! All done!");
		System.out.println("    /\\('')/\\");
		System.out.println("    \\      /");
	}


	// Method to read the expression file
	// Method to load the TF list
	// DPI method
	private static HashMap<String, HashSet<String>> dpi(HashMap<String, HashMap<String, Double>> finalNet, int _threadNumber){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet, _threadNumber);
	}
	
	
	private static void runNoDPI(File expressionFile1,
			File transcriptionFactorsFile, File outputFolder,
			Double miPvalue, int threadCount) 
					throws NumberFormatException, IOException, Exception {
		long initialTime = System.currentTimeMillis();
		// Read expression matrix and do not bootstrap it
		DataMatrix em1 = new DataMatrix(expressionFile1);
		HashMap<String, DataVector> data = em1.data;

		// Check if the sample size is less than the short limit
		int sampleNumber = em1.getSamples().size();
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		// TF list
		HashSet<String> tfSet = DataParser.readGeneSet(transcriptionFactorsFile);
		String[] tfList = tfSet.toArray(new String[0]);
		Arrays.sort(tfList);

		if(tfList.length==0){
			System.err.println("The transcription factor file is badly formatted or empty");
			System.exit(1);
		}

		// Check if the threshold file exists
//		File miThresholdFile = new File(outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
//		double miThreshold;
//		if(!miThresholdFile.exists()){
//			System.err.println("MI threshold file is not present.");
//			System.err.println("Please run ARACNE in --calculateThreshold mode first");
//			System.exit(1);
//		}
//		System.out.println("MI threshold file is present");
//		miThreshold = DataParser.readValue(miThresholdFile);

		// Calculate a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		System.out.println("Calculate network from: "+expressionFile1);

		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		
		MI miCPU = new MI(data,tfList,(double)-1,threadCount);
		
		finalNetwork = miCPU.getFinalNetwork();
		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then DO NOT calculate DPI
		HashMap<String, HashSet<String>> removedEdges = new HashMap<String,HashSet<String>>();

		// And write out the single bootstrap network
		File outputFile = new File(outputFolder.getAbsolutePath()+"/aracneOutput_3col.txt");
		writeFinal(outputFile,removedEdges,finalNetwork);

		long finalTime = System.currentTimeMillis();
		System.out.println("Total time elapsed: "+(finalTime - initialTime)/1000+" sec");

	}


	public static void writeFinal(File finalDPIfile, HashMap<String, HashSet<String>> removedEdges, HashMap<String, HashMap<String, Double>> finalNet){
		try{
			int left = 0;
			int removed = 0;

			BufferedWriter bw = new BufferedWriter(new FileWriter(finalDPIfile));
			for(String k : finalNet.keySet()){
				HashSet<String> tr = null;
				if(removedEdges.containsKey(k)){
					tr = removedEdges.get(k);
				}
				else{
					tr = new HashSet<String>();
				}

				for(String kk : finalNet.get(k).keySet()){
					if(!tr.contains(kk)){
						bw.write(k+"\t"+kk+"\t"+finalNet.get(k).get(kk)+"\n");
						left++;
					} else {
						removed++;
					}
				}
			}
			bw.close();
			System.out.println("Edges removed by DPI:\t"+removed);
			System.out.println("Final Network size:\t"+left);
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}







