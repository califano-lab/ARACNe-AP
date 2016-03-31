package aracne;

import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class DPI {

	public DPI(){
		
	}
	
	public HashMap<String, HashSet<String>> dpi(
			HashMap<String, HashMap<String, Double>> finalNet,
			int threadNumber
			){
			
			// Keys are hubs
			String[] transcriptionFactors = finalNet.keySet().toArray(new String[0]);

			// finalNet contains TFs vs. targets MIs

			// transcriptionFactorsMI: first level keys are TFs, second level keys are TFs, the doubles are MI values
			HashMap<String, HashMap<String, Double>> tftfNetwork = new HashMap<String, HashMap<String, Double>>();
			HashMap<String, HashSet<String>> removedEdges = new HashMap<String,HashSet<String>>();

			for(String tf : transcriptionFactors){
				HashSet<String> temp = new HashSet<String>();
				removedEdges.put(tf, temp);
			}
			
			// This loop populates the TF-TF Mi Hashmap (tftfNetwork)
			for (int i = 0; i < transcriptionFactors.length; i++) {
				for (int j = 0; j < transcriptionFactors.length; j++) {
					if (finalNet.get(transcriptionFactors[i]).containsKey(transcriptionFactors[j])) {
						//if (finalNet.get(keys[i]).get(keys[j]) >= simCut) { // So, since simCut is not used, this check is unnecessary
						
						if (tftfNetwork.containsKey(transcriptionFactors[i])) {
							tftfNetwork.get(transcriptionFactors[i]).put(transcriptionFactors[j], finalNet.get(transcriptionFactors[i]).get(transcriptionFactors[j]));
						} else {
							HashMap<String, Double> dtemp = new HashMap<String, Double>();
							dtemp.put(transcriptionFactors[j], finalNet.get(transcriptionFactors[i]).get(transcriptionFactors[j]));
							tftfNetwork.put(transcriptionFactors[i], dtemp);
						}
						//}
					}
				}
			}

			// This block checks common targets between TFs and operates the DPI check on all TF-target-TF clique triplets
			try{
				ExecutorService executor = Executors.newFixedThreadPool(threadNumber);
				
				for (int i = 0; i < transcriptionFactors.length; i++) {
					// targetsOfI is all genes having a significant edge with the TF i
					// If the first TF is in our tftfMI...
					if (tftfNetwork.containsKey(transcriptionFactors[i])) {
						DPIThread mt = new DPIThread(i, transcriptionFactors, finalNet, tftfNetwork, removedEdges);
						executor.execute(mt);
					}
				}
				
				executor.shutdown();
		        while (!executor.isTerminated()) {
		        	// do not continue before all threads finish
		        }
				
				return(removedEdges);
			}
			catch(Exception e){
				e.printStackTrace();
			}
			return(null);
	}
	

	/**
	 * The Class DPIThread.
	 */
	public class DPIThread extends Thread{
		
		private int i = 0;
		private String[] transcriptionFactors;
		private HashMap<String, HashMap<String, Double>> finalNet;
		HashMap<String, HashMap<String, Double>> tftfNetwork;
		HashMap<String, HashSet<String>> removedEdges;
		
		public DPIThread(int _i, String[] _tfs, HashMap<String, HashMap<String, Double>> _finalNet, HashMap<String, HashMap<String, Double>> _tftfNetwork, HashMap<String, HashSet<String>> _removedEdges){
			i = _i;
			transcriptionFactors = _tfs;
			finalNet = _finalNet;
			tftfNetwork = _tftfNetwork;
			removedEdges = _removedEdges;
		}
		
		public void run(){
			
			HashSet<String> targetsOfI = new HashSet<String>(finalNet.get(transcriptionFactors[i]).keySet());
			HashMap<String, Double> fin1 = finalNet.get(transcriptionFactors[i]);
			HashMap<String, Double> tft1 = tftfNetwork.get(transcriptionFactors[i]);
			
			HashSet<String> rem1 = removedEdges.get(transcriptionFactors[i]);
			
			for (int j = i + 1; j < transcriptionFactors.length; j++) {
				// And if the second TF has an edge with the first TF...
				if (tft1.containsKey(transcriptionFactors[j])) {
					
					HashMap<String, Double> fin2 = finalNet.get(transcriptionFactors[j]);
					HashSet<String> rem2 = removedEdges.get(transcriptionFactors[j]);
					
					// targetsOfJ is all genes having a significant edge with the TF j
					HashSet<String> targetsOfJ = new HashSet<String>(fin2.keySet());
					
					// Intersection: targetsOfJ now contains only targets in common between the two TFs
					targetsOfJ.retainAll(targetsOfI);		// this eats most of the time needed in the DPI calculation
					
					// Obtain the TF-TF MI
					double tftfMI = tft1.get(transcriptionFactors[j]);

					// Loop over the common targets
					for (String target : targetsOfJ){
	
						double v1 = fin1.get(target);
						double v2 = fin2.get(target);
						
						if (v1 < tftfMI && v1 < v2) {
							synchronized(rem1){
								rem1.add(target);
							}
						} else if (v2 < tftfMI && v2 < v1) {
							synchronized(rem2){
								rem2.add(target);
							}
						}
					}
				}
			}
		}
	}
}
