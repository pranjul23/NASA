package org.nasa;

public class Util {

	public static double[] getRandomDistributionForSubsetValues(int start, int end, int totalLength){
		double randomProbabilityArray[] = new double[totalLength];
		double randomProbabilityArrayTemp[] = getRandomDistribution(end-start+1);
		int j=0;
		for(int i= start;i<=end;i++){
			randomProbabilityArray[i] = randomProbabilityArrayTemp[j];
			j++;
		}
		return randomProbabilityArray;
	}
	
	public static int getCorrectIndexForProbabilityDistribution(double d[]){
		double randomNumber = Math.random();
		double sum = 0;
		int i=0;
		
		for(i=0;i<d.length;i++){
			if(randomNumber<sum){
				break;
			}
			sum = sum + d[i];
		}
		
		return i-1;
	}
	
	public static int getRandomNumberFromZeroToLength(int length){
		return  (int)Math.floor(Math.random()*length);
	}
	
	public static int getRandomNumberWithinARange(int min, int max){
		return min +   (int)Math.floor(Math.random()*(max-min));
	}
	
	
	public static double[] getRandomDistribution(int length){
		double randomProbabilityArray[] = new double[length];
		double sum = 0;
		
		for(int i=0;i<length;i++){
			randomProbabilityArray[i] = Math.random();
			sum = sum + randomProbabilityArray[i];
		}
		
		for(int j=0;j<length;j++){
			randomProbabilityArray[j]  = randomProbabilityArray[j] /sum;
		}
		return randomProbabilityArray;
	}
	
	public static void main(String args[]){	
		double d[] = {0.3, 0.3, 0.3};
		System.out.println(  getCorrectIndexForProbabilityDistribution(d));
	}
}
