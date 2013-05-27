package org.nasa;
import java.io.*;
import java.util.*;

public class InputDataLoader {

	public static void main(String args[]) throws Exception{

		
		NormalClassDataGeneration normalData = getOneDataDistribution();
		NormalClassDataGeneration anamolousData = getOneDataDistribution();
		
		 debugMethod(normalData);
		 int numberOfObservationSteps = 200;
		 int trainingData = 500;
		 int normalTestData = 400;
		 int anamolousTestData = 100;
		 SequenceGenerator.GenerateDataFromAnamolies(anamolousTestData, normalTestData, numberOfObservationSteps, anamolousData);		
		 SequenceGenerator.GenerateTrainDataWithNormalDataInTestFile(trainingData, normalTestData, numberOfObservationSteps, normalData);
		 CreateFactorGraphForHSMM.createFactorGraph(normalData);

	}
	
	public static NormalClassDataGeneration getOneDataDistribution() throws Exception{
		FileInputStream fstream = new FileInputStream("./InputData");
		// Get the object of DataInputStream
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;


		strLine = br.readLine();
		strLine = br.readLine();
		//This will contain the maximum number of durations.
		NormalClassDataGeneration normalData = new NormalClassDataGeneration();
		normalData.totalDuration = Integer.parseInt(strLine);

		//Reading one more Comment Line
		strLine = br.readLine();
        normalData.hiddenStateList = new ArrayList<HiddenState>();
        
		while ((strLine = br.readLine()) != null) {
			//Data Reading Starts
			String data[] = strLine.split(" ");
            Duration minDuration = new Duration(Integer.parseInt(data[0]));
			Duration maxDuration = new Duration(Integer.parseInt(data[1]));
			List<String> observationsList = new ArrayList<String>();

			for(int i=2;i<data.length;i++){
				
				observationsList.add(data[i]);
				
				if(!normalData.possibleObservationsSet.contains(data[i])){
					normalData.possibleObservationsSet.add(data[i]);
					normalData.possibleObservationsList.add(data[i]);
				}
				
			}


			HiddenState hiddenState = new HiddenState(observationsList, minDuration, maxDuration);
			normalData.hiddenStateList.add(hiddenState);
           //Data Loading Comes To An End.
		}
		
		
		// Initializing the probability Distributions
		 normalData.initlaizeProbabilites();
		 
		 return normalData;
	}
	
	public static void debugMethod(NormalClassDataGeneration normal){
		System.out.println("Total Duration:" + normal.totalDuration);
		System.out.println("Initial State:" + normal.initialState);
		
		//State Transition Probabilities
		System.out.println("State Transition Probabilites");
		for(int i=0;i<normal.hiddenStateList.size();i++){
			for(int j=0;j<normal.hiddenStateList.size();j++){
			System.out.println("Old State:: " + i + " ::New State:: " + j + " ::Prob:: " + normal.probOfNewStateGivenOldState.probNewStateGivenOldState[j][i]);	
			}
			System.out.println();
		}
		
		//Observation Given State Probabilities
		System.out.println("Observation Given State Probabilities");
		for(int i=0;i<normal.hiddenStateList.size();i++){
			for(int j=0;j<normal.possibleObservationsList.size() ;j++){
				System.out.println("State:: " + i + " ::Observation:: " + normal.possibleObservationsList.get(j) + " ::Prob:: " + normal.probofObservationGivenState.probObsGivenState[j][i]);	
			}
		   System.out.println();
		}
		
		System.out.println("Duration Given State Probabilities");
		for(int i=0;i<normal.hiddenStateList.size();i++){
			for(int j=0;j<normal.totalDuration ;j++){
				System.out.println("State:: " + i + " ::Duration:: " + j + " ::Prob:: " + normal.probDurationGivenState.probDurationGivenState[j][i]);	
			}
		   System.out.println();
		}
		
	  System.out.println("Initial State Probabilities");
	  for(int i=0;i<normal.hiddenStateList.size();i++){
		 System.out.println("State:: " + i + " :Has Initial Probability: " + normal.initHidProb.hiddenStateProbabilities[i]);  
	  }
	  
	  System.out.println("Initial Duration Probabilities");
	  for(int i=0;i<normal.totalDuration;i++){
		 System.out.println("Duration:: " + i + " :Has Initial Probability: " + normal.initDuratiProb.durationProbabilities[i]);  
	  }
	  
	}
}
