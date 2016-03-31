package aracne;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import common.DataMatrix;
import common.DataVector;
import common.Methods;

public class MI {
	// Variables Declaration
	private String[] genes;
	private int sampleNumber;

	private HashMap<String, HashMap<String, Double>> finalNetwork;

	// Main method to test the calculation
	public static void main(String[] args) throws IOException{
		printOrderedMIdistribution();
		//testOnePair(args);
	}

	private static void printOrderedMIdistribution() {
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(new File("expdata/coad_fix.exp")));
			String line = br.readLine(); // Discard header

			// Fill an object which we know the size of
			int nrsamples = 466;
			int nrgenes = 20318;
			short data[][] = new short[nrgenes][nrsamples];

			for(int i=0;i<nrgenes;i++){
				line = br.readLine();
				String[] split = line.split("\t");
				double[] vector = new double[split.length-2];
				for(int j=2; j<split.length; j++){
					vector[j-2]=Double.parseDouble(split[j]);
				}
				data[i] = Methods.rankVector(vector); 
			}
			br.close();
			System.out.println("Data loaded");

			long time1 = System.currentTimeMillis();

			/// Now loop over triplets (not all possible triplets, but incrementally increasing
			double [] mivalues = new double[nrgenes-1];
			for(int i = 0;i<nrgenes-1;i++){
				short[]vectorX = data[i];
				short[]vectorY = data[i+1];
				//	mivalues[i] = computeMI(vectorX, vectorY);
			}
			PrintWriter writer = new PrintWriter("output/distros/mivalues_java.txt", "UTF-8");
			for(int i = 0;i<mivalues.length;i++){
				writer.println(mivalues[i]);
			}			
			writer.close();

			long time2 = System.currentTimeMillis();

			long time = time2-time1;
			System.out.println("Time elapsed: "+(time));


		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unused")
	private static void testOnePair(String[] args) throws IOException {
		FileReader fr = new FileReader(new File("expdata/coad_fix.exp"));
		BufferedReader br = new BufferedReader(fr);
		String line = br.readLine();

		line = br.readLine();
		String[] split = line.split("\t");

		double [] d1 = new double[split.length-2];
		double [] d2 = new double[split.length-2];

		for(int i=2; i<split.length; i++){
			d1[i-2] = Double.parseDouble(split[i])+i/1000000;
		}

		line = br.readLine();
		split = line.split("\t");
		for(int i=2; i<split.length; i++){
			d2[i-2] = Double.parseDouble(split[i])+i/1000000;
		}

		br.close();

		// Rank Transform
		short[] s1 = Methods.rankVector(d1);
		short[] s2 = Methods.rankVector(d2);

		//double mi = computeMI(s1,s2); 

		//	System.out.println(mi);
	}


	// Constructor
	public MI(
			HashMap<String, short[]> rankData1,
			HashMap<String, short[]> rankData2,
			String[] tfList,
			Double miThreshold,
			Integer threadCount
			) {
		this.genes = rankData2.keySet().toArray(new String[0]);
		Arrays.sort(genes);
		this.setSampleNumber(rankData1.get(genes[0]).length);

		List<String> templist = Arrays.asList(tfList);
		HashSet<String> temptf = new HashSet<String>(templist);
		temptf.retainAll(new HashSet<String>(Arrays.asList(genes)));
		tfList = temptf.toArray(new String[0]);
		Arrays.sort(tfList);

		// Loop to check which edges are kept. It will generate the finalNetwork HashMap
		finalNetwork = new HashMap<String, HashMap<String, Double>>();
		for(int i=0; i<tfList.length; i++){
			HashMap<String, Double> tt = new HashMap<String, Double>();
			finalNetwork.put(tfList[i], tt);
		}

		// multi threading here, run threadCount many parallel MI calculations
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		for(int i=0; i<tfList.length; i++){
			MIThread mt = new MIThread(tfList, genes, rankData1, rankData2, miThreshold, i);
			executor.execute(mt);
		}

		executor.shutdown();
		while (!executor.isTerminated()) {
			// do not continue before all threads finish
		}

		System.out.println("TFs processed: "+tfList.length);
		//System.out.println("Edges passing Mi threshold: "+pass);
		//System.out.println("Edges discarded: "+notpass);
	}

	class MIThread extends Thread {

		private String[] tfList;
		private String[] genes;
		private HashMap<String, short[]> rankData1;
		private HashMap<String, short[]> rankData2;
		private double miThreshold;
		private int tfNumber;

		MIThread(String[] tfs, String[] genes, HashMap<String, short[]> rd1, HashMap<String, short[]> rd2, double mi, int tfNumber) {
			this.tfList = tfs;
			this.genes = genes;
			this.rankData1 = rd1;
			this.rankData2 = rd2;
			this.miThreshold = mi;
			this.tfNumber = tfNumber;
		}

		public void run() {
			for(int j=0; j<genes.length; j++){
				if(!genes[j].equals(tfList[tfNumber])){
					short[] vectorX = rankData1.get(tfList[tfNumber]);
					short[] vectorY = rankData2.get(genes[j]);

					double mi = computeMI(vectorX,vectorY);		// mutithread

					// adaptive.miThresh supplants micut
					if(mi >= miThreshold){
						setMI(tfList[tfNumber], genes[j], mi);
					}
				}
			}
		}

		private synchronized void setMI(String _tf, String _gene, double _mi){
			finalNetwork.get(_tf).put(_gene, _mi);
		}
	}

	/**
	 * get bootstraps with pre-computed bootsamples for Paracne.
	 *
	 * @param data the data
	 * @param tfList the tf list
	 * @param miThreshold the mi threshold
	 * @param threadCount the thread count
	 */
	public MI(
			HashMap<String, DataVector> data,
			String[] tfList,
			Double miThreshold,
			Integer threadCount
			) {
		this.genes = data.keySet().toArray(new String[0]);
		Arrays.sort(genes);
		this.setSampleNumber(data.get(genes[0]).values.length);

		List<String> templist = Arrays.asList(tfList);
		HashSet<String> temptf = new HashSet<String>(templist);
		temptf.retainAll(new HashSet<String>(Arrays.asList(genes)));
		tfList = temptf.toArray(new String[0]);
		Arrays.sort(tfList);

		// Loop to check which edges are kept. It will generate the finalNetwork HashMap
		finalNetwork = new HashMap<String, HashMap<String, Double>>();
		for(int i=0; i<tfList.length; i++){
			HashMap<String, Double> tt = new HashMap<String, Double>();
			finalNetwork.put(tfList[i], tt);
		}

		// multi threading here, run threadCount many parallel MI calculations
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		for(int i=0; i<tfList.length; i++){
			MIThreadNA mt = new MIThreadNA(tfList, genes, data, miThreshold, i);
			executor.execute(mt);
		}

		executor.shutdown();
		while (!executor.isTerminated()) {
			// do not continue before all threads finish
		}

		System.out.println("TFs processed: "+tfList.length);
		//System.out.println("Edges passing Mi threshold: "+pass);
		//System.out.println("Edges discarded: "+notpass);
	}

	/**
	 * 
	 * MI thread to compute MI for data with NA values
	 * 
	 */
	class MIThreadNA extends Thread {

		private String[] tfList;
		private String[] genes;
		private HashMap<String, DataVector> data;
		private double miThreshold;
		private int tfNumber;

		/**
		 * Instantiates a new MI thread na.
		 *
		 * @param _tfs the _tfs
		 * @param _genes the _genes
		 * @param _rd1 the _rd1
		 * @param _mi the _mi
		 * @param _tfNumber the _tf number
		 */
		MIThreadNA(String[] _tfs, String[] _genes, HashMap<String, DataVector> _rd1, double _mi, int _tfNumber) {
			tfList = _tfs;
			genes = _genes;
			data = _rd1;
			miThreshold = _mi;
			tfNumber = _tfNumber;
		}

		public void run() {
			for(int j=0; j<genes.length; j++){
				if(!genes[j].equals(tfList[tfNumber])){
					DataVector vectorX = data.get(tfList[tfNumber]);
					DataVector vectorY = data.get(genes[j]);

					double mi = computeMI(vectorX,vectorY);		// mutithread

					// adaptive.miThresh supplants micut
					if(mi >= miThreshold){
						setMI(tfList[tfNumber], genes[j], mi);
					}
				}
			}
		}

		private synchronized void setMI(String _tf, String _gene, double _mi){
			finalNetwork.get(_tf).put(_gene, _mi);
		}
	}

	/**
	 * computeMI for DataVectors when data contains NA values Paracne
	 * 
	 * @param vectorX
	 * @param vectorY
	 * @return
	 */
	public static double computeMI(DataVector vectorX, DataVector vectorY){

		ArrayList<double[]> splitQuad = vectorX.getQuadrants(vectorY);
		short[] valuesX = Methods.rankVector(splitQuad.get(0));
		short[] valuesY = Methods.rankVector(splitQuad.get(1));

		//get the counts for na samples, perform first split by using the counts and the recursive mi call for sample pairs that have both values 
		int bothNACount = (int)splitQuad.get(2)[0];
		int na1 = (int)splitQuad.get(3)[0];
		int na2 = (int)splitQuad.get(4)[0];

		boolean firstLoop = true;
		//double mi = computeMI(valuesX,valuesY);

		//mi=(mi*valuesX.length)-Math.log(valuesX.lenghth);

		double mi = 0;

		if(valuesX.length < 8){
			mi = getInformation(valuesX.length, (valuesX.length+na1), (valuesX.length+na2));
		}
		else{
			mi = computeMI(valuesX,valuesY,(short)0,(short)(valuesX.length),(short)0,(short)(valuesY.length),firstLoop);
		}

		int numberOfSamples = vectorX.values.length;
		mi += getInformation(na1, na1+bothNACount, na1+valuesY.length);
		mi += getInformation(na2, na2+bothNACount, na2+valuesY.length);
		mi += getInformation(bothNACount, bothNACount+na1, bothNACount+na2);

		mi = mi/numberOfSamples + Math.log(numberOfSamples); // TODO why this (proof on whiteboard picture by Yishai, late January 2015)

		return mi;
	}

	// This constructor is to initialize the object for computeMiThreshold
	/**
	 * Instantiates a new mi.
	 *
	 * @param rankData the rank data
	 */
	public MI(
			HashMap<String, short[]> rankData
			) {
		this.genes = rankData.keySet().toArray(new String[0]);
		this.setSampleNumber(rankData.get(genes[0]).length);
	}

	// This constructor is to initialize the object for computeMiThreshold
	/**
	 * Instantiates a new mi.
	 *
	 * @param data the data
	 */
	public MI(DataMatrix data) {
		this.genes = data.genes.toArray(new String[0]);
		this.setSampleNumber(data.samples.size());
	}


	// This constructor is to calculate the full MI (as in --nodpi mode)
	/**
	 * Instantiates a new mi.
	 *
	 * @param rankData1 the rank data1
	 * @param tfList the tf list
	 * @param miThreshold the mi threshold
	 */
	public MI(HashMap<String, short[]> rankData1, String[] tfList, double miThreshold) {
		// Loop to check which edges are kept. It will generate the finalNetwork HashMap
		// All edges between TFs in the provided list will be kept
		finalNetwork = new HashMap<String, HashMap<String, Double>>();
		for(int i=0; i<tfList.length; i++){
			HashMap<String, Double> tt = new HashMap<String, Double>();
			finalNetwork.put(tfList[i], tt);
		}

		for(int i=0; i<tfList.length; i++){
			for(int j=i+1; j<tfList.length; j++){
				short[] vectorX = rankData1.get(tfList[i]);
				short[] vectorY = rankData1.get(tfList[j]);
				double mi = computeMI(vectorX,vectorY);

				if(mi >= miThreshold){
					finalNetwork.get(tfList[i]).put(tfList[j], mi);
				}

			}
		}
		System.out.println("TFs processed: "+tfList.length);
	}

	//// Methods
	public HashMap<String, HashMap<String, Double>> getFinalNetwork() {
		return finalNetwork;
	}

	// Elegant method to find the MI threshold (based on permutation)
	public double calibrateMIThreshold(int numberOfSamples, int numberOfGenes, double miPvalue, int _seed){
		System.out.println("Finding threshold for "+numberOfSamples+" samples");
		short[][] x = new short[numberOfGenes][numberOfSamples];
		Random r = new Random(_seed);
		for(int i=0; i<x.length; i++){
			for(int j=0; j<x[0].length; j++){
				x[i][j] = (short)j;
			}
		}
		for(int i=0; i<x.length; i++){
			for(int j=0; j<x[0].length*10; j++){
				int r1 = r.nextInt(x[0].length);
				int r2 = r.nextInt(x[0].length);
				short temp = x[i][r1];
				x[i][r1] = x[i][r2];
				x[i][r2] = temp;
			}
		}
		ArrayList<Double> mit = new ArrayList<Double>();
		for(int i=0; i<x.length; i++){
			for(int j=(short)(i+1); j<x.length; j++){
				boolean firstLoop = true;
				mit.add(computeMI(x[i],x[j],(short)0,(short)x[i].length,(short)0,(short)x[j].length, firstLoop)/x[i].length + Math.log(x[j].length));
			}
		}
		double[] mis = new double[mit.size()];
		for(int i=0; i<mit.size(); i++){
			mis[i] = mit.get(i);
		}
		double[] rr = fitNull(mis,20);
		double miThreshold = rr[1]*Math.log(miPvalue)+rr[0]; // This looks a lot like the config_kernel.txt parameters
		System.out.println("Parameters for fitted threshold function: "+Arrays.toString(fitNull(mis,100)));
		System.out.println("MI threshold: "+miThreshold);
		return(miThreshold);
	}

	// Fit a null distribution method
	private double[] fitNull(double[] _x, int _tailLength){
		Arrays.sort(_x);
		double[] y = new double[_x.length];
		for(int i=0; i<_x.length; i++){
			y[i] = (_x.length-i)*1.0/_x.length;
		}
		double[] tailx = new double[_tailLength];
		double[] tailMI = new double[_tailLength];
		System.arraycopy(y, y.length-_tailLength, tailx, 0, _tailLength);
		System.arraycopy(_x, _x.length-_tailLength, tailMI, 0, _tailLength);
		for(int i=0; i<tailx.length; i++){
			tailx[i] = Math.log(tailx[i]);
		}
		double[] temp = tailx;
		tailx = tailMI;
		tailMI = temp;
		double sumx = 0.0;
		double sumy = 0.0;
		for(int i=0; i<tailx.length; i++){
			sumx  += tailMI[i];
			sumy  += tailx[i];
		}
		double xbar = sumx / tailx.length;
		double ybar = sumy / tailx.length;
		double xxbar = 0.0;
		double xybar = 0.0;
		for (int i = 0; i < tailx.length; i++) {
			xxbar += (tailMI[i] - xbar) * (tailMI[i] - xbar);
			xybar += (tailMI[i] - xbar) * (y[i] - ybar);
		}
		double beta1 = xybar / xxbar;
		double beta0 = ybar - beta1 * xbar;
		double[] re = new double[2];
		re[0] = beta0;
		re[1] = beta1;
		return re;
	}

	// This method takes into consideration NAs
	public double calibrateMIThresholdNA(HashMap<String, DataVector> _data, int randomPaircount, double miPvalue, int _seed){
		System.out.println("Finding threshold for NA "+randomPaircount+" gene pairs");

		HashMap<String, DataVector> tempData = (HashMap<String, DataVector>)_data.clone();
		Random r = new Random(_seed);

		String[] kset = _data.keySet().toArray(new String[0]);

		ArrayList<Double> mit = new ArrayList<Double>();

		for(int i=0; i<randomPaircount; i++){

			String r1 = kset[r.nextInt(kset.length)];
			String r2 = kset[r.nextInt(kset.length)];

			for(int j=0; j<tempData.get(r1).values.length; j++){

				int randPos = r.nextInt(tempData.get(r1).values.length);
				double temp = tempData.get(r1).values[j];
				boolean temp2 =  tempData.get(r1).NAs[j];

				tempData.get(r1).values[j] =  tempData.get(r1).values[randPos];
				tempData.get(r1).NAs[j] =  tempData.get(r1).NAs[randPos];

				tempData.get(r1).values[randPos] = temp;
				tempData.get(r1).NAs[randPos] = temp2;


				randPos = r.nextInt(tempData.get(r2).values.length);
				temp = tempData.get(r2).values[j];
				temp2 =  tempData.get(r2).NAs[j];

				tempData.get(r2).values[j] =  tempData.get(r2).values[randPos];
				tempData.get(r2).NAs[j] =  tempData.get(r2).NAs[randPos];

				tempData.get(r2).values[randPos] = temp;
				tempData.get(r2).NAs[randPos] = temp2;
			}

			mit.add(computeMI(tempData.get(r1), tempData.get(r2)));
		}

		double[] mis = new double[mit.size()];
		for(int i=0; i<mit.size(); i++){
			mis[i] = mit.get(i);
		}
		double[] rr = fitNull(mis,(int)Math.round(randomPaircount*0.05));
		double miThreshold = rr[1]*Math.log(miPvalue)+rr[0]; // This looks a lot like the config_kernel.txt parameters
		System.out.println("Parameters for fitted threshold function: "+Arrays.toString(fitNull(mis,100)));
		System.out.println("MI threshold: "+miThreshold);
		return(miThreshold);
	}


	//// This method calculates MI through the adaptive partitioning technique
	// This is ideally the one to be turned into a GPU marvel
	public static double computeMI(short[] vectorX, short[] vectorY){
		boolean firstLoop = true;
		return computeMI(vectorX,vectorY,(short)0,(short)vectorX.length,(short)0,(short)vectorY.length,firstLoop)
				/vectorX.length 
				+ Math.log(vectorX.length); // TODO why this (proof on whiteboard picture by Yishai, late January 2015)
	}

	private static double computeMI(
			short[] vectorX,
			short[] vectorY,
			short limX1,
			short limX2,
			short limY1,
			short limY2,
			boolean firstLoop
			){
		double mi = 0;
		short pointsInQuadrant1 = 0;
		short pointsInQuadrant2 = 0;
		short pointsInQuadrant3 = 0;
		short pointsInQuadrant4 = 0;

		// Pivots on half quadrants
		short splitX = (short)((limX2+limX1)/2);
		short splitY = (short)((limY2+limY1)/2);

		short[] libraryX1 = new short[vectorX.length];
		short[] libraryX2 = new short[vectorX.length];
		short[] libraryX3 = new short[vectorX.length];
		short[] libraryX4 = new short[vectorX.length];
		short[] libraryY1 = new short[vectorX.length];
		short[] libraryY2 = new short[vectorX.length];
		short[] libraryY3 = new short[vectorX.length];
		short[] libraryY4 = new short[vectorX.length];


		for(int i=0; i<vectorX.length; i++){
			if(vectorX[i] <= splitX){
				if(vectorY[i] <= splitY){
					libraryX1[pointsInQuadrant1] = vectorX[i];
					libraryY1[pointsInQuadrant1] = vectorY[i];
					pointsInQuadrant1++;
				}
				else{
					libraryX2[pointsInQuadrant2] = vectorX[i];
					libraryY2[pointsInQuadrant2] = vectorY[i];
					pointsInQuadrant2++;
				}
			}
			else{
				if(vectorY[i] > splitY){
					libraryX3[pointsInQuadrant3] = vectorX[i];
					libraryY3[pointsInQuadrant3] = vectorY[i];
					pointsInQuadrant3++;
				}
				else{
					libraryX4[pointsInQuadrant4] = vectorX[i];
					libraryY4[pointsInQuadrant4] = vectorY[i];
					pointsInQuadrant4++;
				}
			}
		}

		// Chi-square uniformity test
		double expected = vectorX.length/4.0; 
		double tst = (Math.pow(pointsInQuadrant1-expected,2) + Math.pow(pointsInQuadrant2-expected,2) + Math.pow(pointsInQuadrant3-expected,2) + Math.pow(pointsInQuadrant4-expected,2))/expected;

		// 7.815 is the hard-coded threshold, below which you pass the uniformity Chi-Square test with confidence >95% 

		if(tst > 7.815 || firstLoop){
			firstLoop = false;

			// TODO System.arraycopy may be faster for hundreds of samples

			// Operations to put points in different quadrants
			// TODO the splitX-X1 can be precomputed for future optimization
			if(pointsInQuadrant1 > 2){
				short[] x1 = Arrays.copyOfRange(libraryX1, 0, pointsInQuadrant1);
				short[] y1 = Arrays.copyOfRange(libraryY1, 0, pointsInQuadrant1);

				mi += computeMI(x1, y1, limX1, splitX, limY1, splitY, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant1, splitX-limX1, splitY-limY1);
			}

			if(pointsInQuadrant2 > 2){
				short[] x2 = Arrays.copyOfRange(libraryX2, 0, pointsInQuadrant2);
				short[] y2 = Arrays.copyOfRange(libraryY2, 0, pointsInQuadrant2);

				mi += computeMI(x2, y2, limX1, splitX, splitY, limY2, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant2, splitX-limX1, limY2-splitY);
			}

			if(pointsInQuadrant3 > 2){
				short[] x3 = Arrays.copyOfRange(libraryX3, 0, pointsInQuadrant3);
				short[] y3 = Arrays.copyOfRange(libraryY3, 0, pointsInQuadrant3);

				try{
					mi += computeMI(x3, y3, splitX, limX2, splitY, limY2, firstLoop);
				}
				catch(Exception e){
					System.out.println(e);
				}
			}
			else{
				mi += getInformation(pointsInQuadrant3, limX2-splitX, limY2-splitY);
			}

			if(pointsInQuadrant4 > 2){
				short[] x4 = Arrays.copyOfRange(libraryX4, 0, pointsInQuadrant4);
				short[] y4 = Arrays.copyOfRange(libraryY4, 0, pointsInQuadrant4);

				mi += computeMI(x4, y4, splitX, limX2, limY1, splitY, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant4, limX2-splitX, splitY-limY1);
			}
		}
		else{
			mi += getInformation(vectorX.length, limX2-limX1, limY2-limY1);
		}

		return mi;
	}

	private static double getInformation(double pXY, double pX, double pY){
		if(pXY == 0){
			return  0;
		}
		else{
			// This little formula contains the entire MI calculation. Wow
			return pXY*Math.log( pXY/(pX*pY));
		}
	}

	public int getSampleNumber() {
		return sampleNumber;
	}

	public void setSampleNumber(int sampleNumber) {
		this.sampleNumber = sampleNumber;
	}

	public short median(short[] input){
		short median;
		Arrays.sort(input);
		median = (short) input[input.length/2];
		return(median);
	}
}

