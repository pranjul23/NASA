package org.nasa;

public class InitialHiddenStateProbabilities {

	double hiddenStateProbabilities[] = null;
	
	public InitialHiddenStateProbabilities(int smallestState, int highestState, int maximumState){
		hiddenStateProbabilities = new double[maximumState];
		hiddenStateProbabilities = Util.getRandomDistributionForSubsetValues(smallestState, highestState, maximumState);		
		
	}
	
	
}
