package org.nasa;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;

public class ProbofObservationGivenState {

	double probObsGivenState[][] = null;
	
	ProbofObservationGivenState(int numberOfObservations, int numberOfStates){
		probObsGivenState = new double[numberOfObservations][numberOfStates];
	}
	
	public int getIndexOfTheObservation(List<String> allObservations, String observation){
		int counter = 0;
		for(String obs : allObservations){
		  if(obs.equals(observation)){
			  return counter;
		  }
		  counter++;
		}
		
		return 0;
		
	}
	
	public void initializeProbabilityTable(NormalClassDataGeneration normalClassDataGeneration){	  
		int stateCounter = 0;
		int probabilityCounter = 0;	
	
		for(HiddenState hiddenState : normalClassDataGeneration.hiddenStateList ){
			List<String> observationList = hiddenState.observationsList ;
			double probablities[] = Util.getRandomDistribution(observationList.size());
			probabilityCounter = 0;	
			for(String observation : observationList){
				int index =  this.getIndexOfTheObservation(normalClassDataGeneration.possibleObservationsList, observation);
				probObsGivenState[index][stateCounter] = probablities[probabilityCounter];
				probabilityCounter++;
			  }			
			stateCounter++;
		}
		
		
	}
	
}
