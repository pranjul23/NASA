package org.nasa;

public class ProbOfNewStateGivenOldState {

	double probNewStateGivenOldState[][] = null;
	int numStates = 0;
	
	public ProbOfNewStateGivenOldState(int numStates) {
		probNewStateGivenOldState = new double[numStates][numStates];
		this.numStates = numStates;
		
	}
	
	public void initializeProbabilityTables(){
		for(int i=0;i<numStates;i++){
			double d[] = Util.getRandomDistribution(numStates);
			for(int j=0;j<numStates;j++){
				probNewStateGivenOldState[j][i] = d[j];
			}
		}
	}
}
