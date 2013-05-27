package org.nasa;
import java.util.*;

public class ProbDurationGivenState {

	double probDurationGivenState[][] ;
	
	public ProbDurationGivenState(int duration, int numStates) {
		probDurationGivenState = new double[duration][numStates];
		
	}
	
	public void initializeProbDurationGivenState(NormalClassDataGeneration normalClassDataGeneration){
		
		List<HiddenState> hiddenStatesList = normalClassDataGeneration.hiddenStateList;
		int stateCounter = 0;
		for(HiddenState hiddenState : hiddenStatesList){
		  Duration minDuration = hiddenState.minDuration;
		  Duration maxDuration = hiddenState.maxDuration;
		  double probs[] = Util.getRandomDistributionForSubsetValues(minDuration.timeStamp, maxDuration.timeStamp, normalClassDataGeneration.totalDuration);
		  
		  for(int j=0;j<probs.length;j++){
			  probDurationGivenState[j][stateCounter] = probs[j];
		  }
		  stateCounter++;
		}
		
	}
	
}
