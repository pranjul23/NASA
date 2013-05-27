package org.nasa;
import java.util.*;

public class NormalClassDataGeneration {

	//Stores All Possible Observations In All The Sequences.
	Set<String> possibleObservationsSet = new HashSet<String>();
	List<String> possibleObservationsList = new ArrayList<String>();
	
	//The initial State of the System
	int initialState ;
	
	//The initial Duration of the System
	int initialDuration;
	
	//The number flight samples the model can generate.
	int maxSequenceLength ;
	
	//The smallest timestamp for this observation will continue.
	int minDuration;
	
	//The largest timestamp  for which this observation will continue.
	int maxDuration;
	
	// The maximum timestamp for which any observation can continue.
	int totalDuration;
	
	// The smallest Possible Initial State.
	int lowestPossibleInitialState;
	
	// The highest Possible Initial State.
	int highestPossibleInitialState;
	
	// The maximum number of Initial State.
	int totalPossibleInitialState;
	
	//The number of observations one flight can generate.
	int maxObservationLength;
	
	//Possible Set of hidden state a system can be in.
	List<HiddenState> hiddenStateList = new ArrayList<HiddenState>();
	
	//Generating the initial set of probabilities for the duration
	InitialDurationProb initDuratiProb = null;
	
	//Generating the initial set of probabilities for the Hidden State the system would be in.
	InitialHiddenStateProbabilities initHidProb	= null;
	ProbofObservationGivenState probofObservationGivenState = null;
	ProbOfNewStateGivenOldState probOfNewStateGivenOldState = null;
	ProbDurationGivenState probDurationGivenState = null;
	
	public void initlaizeProbabilites(){
	
		//For Now let the lowestPossibleInitialState be zero
	 lowestPossibleInitialState = 0;
	 highestPossibleInitialState = hiddenStateList.size()-1;
	 totalPossibleInitialState = hiddenStateList.size();
	
		// Generating A Random Initial State Distribution.
	 initHidProb = new InitialHiddenStateProbabilities(lowestPossibleInitialState, highestPossibleInitialState , totalPossibleInitialState);
	    //Generating A Random Initial State From it.
	 initialState = Util.getCorrectIndexForProbabilityDistribution(initHidProb.hiddenStateProbabilities);
	 
	    //Getting the Durations From this State.
	 int minDuration = hiddenStateList.get(initialState).minDuration.timeStamp;
	 int maxDuration = hiddenStateList.get(initialState).maxDuration.timeStamp;
	 
       //Total Duration is already set by the InputDataLoader.
	 
	    //Getting the Duration Probabilites...Igor needs it for latter verifying from EM model.
	 initDuratiProb = new InitialDurationProb(minDuration, maxDuration,totalDuration);
	 initialDuration = Util.getCorrectIndexForProbabilityDistribution(initDuratiProb.durationProbabilities);
	 
	 //At this point Initial Duration and Initial State have been already Set.
	 //Now Initialzing the probability of observation given state.
	 probofObservationGivenState = new ProbofObservationGivenState(possibleObservationsSet.size(), hiddenStateList.size());
	 probofObservationGivenState.initializeProbabilityTable(this);
	 
	 probOfNewStateGivenOldState = new ProbOfNewStateGivenOldState(this.hiddenStateList.size());
	 probOfNewStateGivenOldState.initializeProbabilityTables();
	 
	 probDurationGivenState = new ProbDurationGivenState(this.totalDuration, this.hiddenStateList.size());
	 probDurationGivenState.initializeProbDurationGivenState(this);
	 	
	}
	
	
	
	
}
