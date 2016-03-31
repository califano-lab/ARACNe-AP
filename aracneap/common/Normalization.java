package common;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;



// This class apparently rank-transforms an expression matrix
public class Normalization {
	public static HashMap<String, short[]> rankTransform(HashMap<String, HashMap<String, Double>> inputData, double noise){
		Random noiseGenerator = new Random();
		String[] profiles = inputData.keySet().toArray(new String[0]);
		Arrays.sort(profiles);
		String[] genes = inputData.get(profiles[0]).keySet().toArray(new String[0]);
		Arrays.sort(genes);
		
		HashMap<String, short[]> rankTransform = new HashMap<String, short[]>();
		
		for(String g : genes){
			double[] values = new double[profiles.length];
			short[] rankValue = new short[profiles.length];
			for(int i=0; i<profiles.length; i++){
				// White noise is added here
				values[i] = inputData.get(profiles[i]).get(g)+noiseGenerator.nextDouble()*noise;
			}
			for(int i=0; i<values.length; i++){
				short counter = 1;
				for(int j=0; j<values.length; j++){
					if(values[i] < values[j]){
						counter++;
					}
				}
				rankValue[i] = counter;
			}
			rankTransform.put(g, rankValue);
		}
		return rankTransform;
	}
}
