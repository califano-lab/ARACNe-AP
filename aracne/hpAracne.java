package aracne;

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

import common.DataMatrix;
import common.DataParser;
import common.DataVector;
import common.ExpressionMatrix;


public class hpAracne{
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static float simCut = 0; // In the original paper this was a parameter. E.g. if two TFs have a very high MI, DPI is not calculated
	final static int geneNumberInMIthresholding = 1000;
	static Random random = new Random();

	// Main Method
	public static void main(String[] args) throws Exception {
		Locale.setDefault(Locale.US);

		//// Parse arguments
		CmdLineParser parser = new CmdLineParser();
		CmdLineParser.Option optionExpressionFile = parser.addStringOption('e', "expfile");
		CmdLineParser.Option optionTranscriptionFactorsFile = parser.addStringOption('t', "tfs");
		CmdLineParser.Option optionOutputFolder = parser.addStringOption('o', "output");
		CmdLineParser.Option optionConsolidate = parser.addBooleanOption('c',"consolidate");
		CmdLineParser.Option optionCalculateThreshold = parser.addBooleanOption('j',"calculateThreshold");
		CmdLineParser.Option optionPvalue = parser.addDoubleOption('p',"pvalue");
		CmdLineParser.Option optionSeed = parser.addIntegerOption('s',"seed");
		CmdLineParser.Option optionThreaded = parser.addIntegerOption('m',"threads");
		CmdLineParser.Option optionNodpi = parser.addBooleanOption('n',"nodpi");
		CmdLineParser.Option optionNobootstrap = parser.addBooleanOption('b',"nobootstrap");
		CmdLineParser.Option optionNobonferroni = parser.addBooleanOption('r',"nobonferroni");
		CmdLineParser.Option optionConsolidatepvalue = parser.addDoubleOption('v',"consolidatepvalue");

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

		boolean nobootstrap = false;
		if ((Boolean)parser.getOptionValue(optionNobootstrap)!=null) {
			nobootstrap = true;
		}
		
		boolean nobonferroni = false;
		if ((Boolean)parser.getOptionValue(optionNobonferroni)!=null) {
			nobonferroni = true;
		}
		
		Double consolidatePvalue;
		try {
			consolidatePvalue = new Double((Double)parser.getOptionValue(optionConsolidatepvalue));
		} catch (Exception e) {
			consolidatePvalue = 0.05;
		}
		

		// Here the program forks
		// You can calculate the MI threshold
		// You can run a single bootstrap (or non bootstrap)
		// You can consolidate bootstraps

		if(isThreshold){
			File expressionFile = new File((String)parser.getOptionValue(optionExpressionFile));
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
			String expFile = (String)parser.getOptionValue(optionExpressionFile);
			File expressionFile = new File(expFile);
			File transcriptionFactorsFile = new File((String)parser.getOptionValue(optionTranscriptionFactorsFile));
			Integer seed = (Integer)parser.getOptionValue(optionSeed);
			Integer threadCount = (Integer)parser.getOptionValue(optionThreaded);

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
			runAracne(
					expressionFile,
					transcriptionFactorsFile,
					outputFolder,
					processId,
					miPvalue,
					threadCount,
					noDPI,
					nobootstrap
					);
		} else {
			runConsolidate(outputFolder,nobonferroni,consolidatePvalue);
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
		
		miThreshold = miCPU.calibrateMIThresholdNA(data,geneNumberInMIthresholding,_miPvalue,_seed);
		DataParser.writeValue(miThreshold, miThresholdFile);
	}

	// Single run mode
	private static void runAracne(
			File expressionFile,
			File transcriptionFactorsFile,
			File outputFolder, 
			String processId,
			Double miPvalue,
			Integer threadCount,
			boolean noDPI, // Do not use DPI
			boolean nobootstrap // Do not use bootstrap
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();
		// Read expression matrices 
		DataMatrix dm = new DataMatrix(expressionFile);
		
		// Don't bootstrap them
		if(!nobootstrap){
			System.out.println("Bootstrapping input matrix 1 with "+dm.getGenes().size()+" genes and "+dm.getSamples().size()+" samples");
			dm.bootstrap(random);
		}
		HashMap<String, DataVector> data = dm.data;

		// Check if the sample size is less than the short limit
		int sampleNumber = dm.getSamples().size();
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
		File miThresholdFile = new File(outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;
		if(!miThresholdFile.exists()){
			System.err.println("MI threshold file is not present.");
			System.err.println("Please run ARACNE in --calculateThreshold mode first");
			System.exit(1);
		}
		System.out.println("MI threshold file is present");
		miThreshold = DataParser.readValue(miThresholdFile);

		// Calculate a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		System.out.println("Calculate network from: "+expressionFile);
		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		MI miCPU = new MI(data,tfList,miThreshold,threadCount);
		finalNetwork = miCPU.getFinalNetwork();
		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then calculate DPI
		long time2 = System.currentTimeMillis();
		HashMap<String, HashSet<String>> removedEdges;
		if(noDPI){
			removedEdges = new HashMap<String, HashSet<String>>();
		}else {
			removedEdges = dpi(finalNetwork, threadCount);
			System.out.println("DPI time elapsed: "+(System.currentTimeMillis() - time2)/1000+" sec");
		}

		// And write out the single bootstrap network
		File outputFile;
		if(nobootstrap){
			outputFile = new File(outputFolder.getAbsolutePath()+"/nobootstrap_network.txt");
		} else {
			outputFile = new File(outputFolder.getAbsolutePath()+"/bootstrapNetwork_"+processId+".txt");
		}
		writeFinal(outputFile,removedEdges,finalNetwork);

		long finalTime = System.currentTimeMillis();
		System.out.println("Total time elapsed: "+(finalTime - initialTime)/1000+" sec");
	}

	// This method consolidates the several bootstraps
	private static void runConsolidate(File outputFolder, boolean nobonferroni, Double consolidatePvalue) throws IOException {
		BootstrapConsolidator c = new BootstrapConsolidator(nobonferroni);
		c.mergeFiles(outputFolder);
		String outputFile = outputFolder+"/network.txt";

		
		// consolidatePvalue is the P-value for the Poisson distribution. Aka how many times an edge has to appear in the bootstraps to be kept.
		// Hard-coded to 0.3 in the original ARACNe
		c.writeSignificant(outputFile, consolidatePvalue);

		System.out.println("\n        :");
		System.out.println("       :");
		System.out.println("        :       Cool! All done!");
		System.out.println("    /\\('')/\\");
		System.out.println("    \\      /");
	}


	// Method to read the expression file
	// Method to load the TF list
	// DPI method
	private static HashMap<String, HashSet<String>> dpi(
			HashMap<String,
			HashMap<String, Double>> finalNet,
			int threadNumber
		){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet, threadNumber);
	}

	public static void writeFinal(File finalDPIfile, HashMap<String, HashSet<String>> removedEdges, HashMap<String, HashMap<String, Double>> finalNet){
		try{
			int left = 0;
			int removed = 0;

			BufferedWriter bw = new BufferedWriter(new FileWriter(finalDPIfile));
			
			// Header
			bw.write("Regulator\tTarget\tMI\n");

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






