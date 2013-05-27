package org.nasa;
import java.util.*;
/*
 * One Hidden state will have many a limit of duration i.e the maximum amount of time
 * the plane can stay in that state and a set of observations that is possible set of 
 * observations one can  see from that state.
 */
public class HiddenState {

	Duration minDuration ;
	Duration maxDuration ;
	List<String> observationsList = new ArrayList<String>();
	
	HiddenState(List<String> observationsList , Duration minDuration, Duration maxDuration){
		this.minDuration = minDuration;
		this.maxDuration = maxDuration;
		this.observationsList = observationsList;
	}	
	
}
