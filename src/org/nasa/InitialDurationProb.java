package org.nasa;

/*
 * Class to store the init Duration Probabilities
 */

public class InitialDurationProb {

	//Assume that Duration starts from Time Stamp 0 until maxDuration -1
	double[] durationProbabilities = null;
	
	public InitialDurationProb(int minDuration , int maxDuration , int totalDuration) {
		durationProbabilities = Util.getRandomDistributionForSubsetValues(minDuration, maxDuration, totalDuration);		
	}
	
}
