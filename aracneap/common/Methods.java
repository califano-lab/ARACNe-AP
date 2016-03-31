package common;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.ArrayUtils;

public class Methods {
	public static short median(short[] input){
		short median;
		Arrays.sort(input);
		median = (short) input[input.length/2];
		return(median);
	}

	public static short mean(short[] input){
		int intmean = 0;
		for (short i = 0; i<input.length; i++){
			intmean += input[i]; 
		}
		short mean = (short)(Math.round(intmean / input.length));
		return(mean);
	}

	public static double mean(int[] input){
		double mean = 0.0;
		for (int i = 0; i<input.length; i++){
			mean += input[i]; 
		}
		mean = mean / input.length;
		return(mean);
	}

	public static double mean(ArrayList<Integer> input){
		double mean = 0.0;
		for (int i = 0; i<input.size(); i++){
			mean += input.get(i); 
		}
		mean = mean / input.size();
		return(mean);
	}


	public static double[][] transposeMatrix(double [][] m){
		double[][] temp = new double[m[0].length][m.length];
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[0].length; j++)
				temp[j][i] = m[i][j];
		return temp;
	}

	public static short[] rankVector(double[] inputVector){
		short[] rankVector = new short[inputVector.length];
		for(int i=0; i<inputVector.length; i++){
			int counter = 1;
			for(int j=0; j<inputVector.length; j++){
				if(inputVector[i] < inputVector[j]){
					counter++;
				}
			}
			rankVector[i] = (short)counter;
		}
		return rankVector;
	}

	public static short[] union(short[] shortA, short[] shortB){
		Short[] inputA = ArrayUtils.toObject(shortA);
		Short[] inputB = ArrayUtils.toObject(shortB);
		Set<Short> a = new TreeSet<Short>(Arrays.asList(inputA));
		Set<Short> b = new TreeSet<Short>(Arrays.asList(inputB));

		Set<Short> c = new TreeSet<Short>(a);
		c.addAll(b);

		short[] output = ArrayUtils.toPrimitive(c.toArray(new Short[c.size()]));
		return (output);
	}

	public static short[] intersect(short[] shortA, short[] shortB){
		Short[] inputA = ArrayUtils.toObject(shortA);
		Short[] inputB = ArrayUtils.toObject(shortB);
		Set<Short> a = new TreeSet<Short>(Arrays.asList(inputA));
		Set<Short> b = new TreeSet<Short>(Arrays.asList(inputB));

		Set<Short> c = new TreeSet<Short>(a);
		c.retainAll(b);

		short[] output = ArrayUtils.toPrimitive(c.toArray(new Short[c.size()]));
		return (output);
	}

	public static short[] difference(short[] shortA, short[] shortB){
		Short[] inputA = ArrayUtils.toObject(shortA);
		Short[] inputB = ArrayUtils.toObject(shortB);
		Set<Short> a = new TreeSet<Short>(Arrays.asList(inputA));
		Set<Short> b = new TreeSet<Short>(Arrays.asList(inputB));

		Set<Short> c = new TreeSet<Short>(a);
		c.removeAll(b);

		short[] output = ArrayUtils.toPrimitive(c.toArray(new Short[c.size()]));
		return (output);
	}


	public static short[] reverse(short[] shortA){
		Short[] inputA = ArrayUtils.toObject(shortA);
		ArrayList<Short> list = (ArrayList<Short>) Arrays.asList(inputA);
		Collections.reverse(list);

		short[] output = ArrayUtils.toPrimitive(list.toArray(new Short[list.size()]));
		return (output);
	}


	public static short[] generateRandomVector(int length, Random random) {
		short [] randomArray = new short[length];
		for(short i = 0; i< length; i++){
			randomArray[i]=i;
		}
		shuffleArray(randomArray, random);
		return(randomArray);

	}


	public static void shuffleArray(short[] array, Random random) {
		int index;
		short temp;
		for (int i = array.length - 1; i > 0; i--) {
			index = random.nextInt(i + 1);
			temp = array[index];
			array[index] = array[i];
			array[i] = temp;
		}
	}

	public static ArrayList<String> intersect(ArrayList<String> list1, ArrayList<String> list2){
		ArrayList<String> list = new ArrayList<String>();
		for (String obj : list1) {
			if(list2.contains(obj)) {
				list.add(obj);
			}
		}
		return list;
	}


}


















