package common;

public class Sorter {

	private double[] numbers;
	private String[] gene;
	
	public double[] getVaues(){
		return numbers;
	}
	
	public String[] getGenes(){
		return gene;
	}
	
	public void sort(double[] exp, String[] _gene, boolean _desc) {
		// Check for empty or null array
		if (exp ==null || exp.length==0){
			return;
		}
		
		this.numbers = exp;
		int number = exp.length;
		gene = _gene;
		
		quicksort(0, number - 1);
		
		if(_desc){
			
			for(int i=0; i<exp.length/2; i++){
				double td = 0;
				String tg = "";
				
				td = exp[i];
				exp[i] = exp[exp.length-1-i];
				exp[exp.length-1-i] = td;
				
				tg = gene[i];
				gene[i] = gene[gene.length-1-i];
				gene[gene.length-1-i] = tg;
			}
		}
	}

	private void quicksort(int low, int high) {
		int i = low, j = high;

		double pivot = numbers[low + (high-low)/2];

		while (i <= j) {

			while (numbers[i] < pivot) {
				i++;
			}

			while (numbers[j] > pivot) {
				j--;
			}

			if (i <= j) {
				exchange(i, j);
				i++;
				j--;
			}
		}

		if (low < j)
			quicksort(low, j);
		if (i < high)
			quicksort(i, high);
	}

	private void exchange(int i, int j) {
		double temp = numbers[i];
		numbers[i] = numbers[j];
		numbers[j] = temp;
		
		String tg = gene[i];
		gene[i] = gene[j];
		gene[j] = tg;
	}
	
}
