package org.nasa;

public class ConstructDurationGivenOldStateAndOldDuration {

	public static double[][][] generateNewDurationGivenOldDurationAndOldState(NormalClassDataGeneration normal)
	{
		double durationBasedOnStateAndDuration[][][] = new double[normal.totalDuration+1][normal.hiddenStateList.size()][normal.totalDuration+1];
		
		for(int newDuration=0;newDuration<normal.totalDuration+1;newDuration++){
			for(int oldState=0;oldState<normal.hiddenStateList.size();oldState++){
				for(int oldDuration=0;oldDuration<normal.totalDuration+1;oldDuration++){
					
					if(oldDuration!=0){
						durationBasedOnStateAndDuration[newDuration][oldState][newDuration] = 0;
					}else{
						durationBasedOnStateAndDuration[newDuration][oldState][newDuration] =
								normal.probDurationGivenState.probDurationGivenState[newDuration][oldState];
					}
				}
			}
		}
		
		return durationBasedOnStateAndDuration;
		
	}
	
	public static double[][][] generateNewStateGivenOldDurationAndOldState(NormalClassDataGeneration normal)
	{
		double stateBasedOnStateAndDuration[][][] = new double[normal.hiddenStateList.size()][normal.hiddenStateList.size()][normal.totalDuration+1];
		
		for(int newState=0;newState<normal.hiddenStateList.size();newState++){
			for(int oldState=0;oldState<normal.hiddenStateList.size();oldState++){
				for(int oldDuration=0;oldDuration<normal.totalDuration+1;oldDuration++){
					
					if(oldDuration!=0){
						stateBasedOnStateAndDuration[newState][oldState][newState] = 0;
					}else{
						stateBasedOnStateAndDuration[newState][oldState][newState] =
								normal.probOfNewStateGivenOldState.probNewStateGivenOldState[newState][oldState];
					}
				}
			}
		}
		
		return stateBasedOnStateAndDuration;
		
	}
}
