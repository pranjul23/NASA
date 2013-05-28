package org.nasa;

import java.io.*;

public class SequenceGenerator {
	
	public static void GenerateDataFromAnamoliesForOtherModels(int anamolousTestData, int normalTestData,  int length,
			NormalClassDataGeneration normal) throws Exception {
		
// Writing in the Training Data
		FileWriter fstream = new FileWriter("HMMtesting.txt");
		BufferedWriter out = new BufferedWriter(fstream);

	
		
		// Duration 1 means Actually Duration 0
		int j = 0;
		while (j < anamolousTestData) {
			String result = "";
			int i = 0;
			int oldState = normal.initialState;
			int oldDuration = normal.initialDuration;

			while (i < length) {
				int nextState = getNewState(oldState, oldDuration, normal);
				int nextDuration = getNewDuration(oldState, oldDuration, normal);
				String observationSymbol = getObservationSymbol(nextState,
						normal);

				oldState = nextState;
				oldDuration = nextDuration;
				result = result + test.getInteger(observationSymbol) + " ";
				i++;

			}
			System.out.println(result);
			result = result.substring(0, result.length() - 1);
			out.write(result);
			out.write("\n");
			j++;

		}
		out.close();
	}

	public static void GenerateDataFromAnamolies(int anamolousTestData, int normalTestData,  int length,
			NormalClassDataGeneration normal) throws Exception {
		
// Writing in the Training Data
		FileWriter fstream = new FileWriter("HSMMtesting.txt");
		BufferedWriter out = new BufferedWriter(fstream);

		out.write(""+ (anamolousTestData+normalTestData));
		out.write("\n");
		
		String demoLine = "";
		int j =4;
		for(int i=0;i<length;i++){
			demoLine = demoLine +   j + " ";
			j = j + 3;
		}
		
		// Duration 1 means Actually Duration 0
		j = 0;
		while (j < anamolousTestData) {
			String result = "";
			int i = 0;
			int oldState = normal.initialState;
			int oldDuration = normal.initialDuration;

			while (i < length) {
				int nextState = getNewState(oldState, oldDuration, normal);
				int nextDuration = getNewDuration(oldState, oldDuration, normal);
				String observationSymbol = getObservationSymbol(nextState,
						normal);

				oldState = nextState;
				oldDuration = nextDuration;
				result = result + test.getInteger(observationSymbol) + " ";
				i++;

			}
			System.out.println(result);
			result = result.substring(0, result.length() - 1);
			out.write(""+length);
			out.write("\n");
			out.write(demoLine);
			out.write("\n");
			out.write(result);
			out.write("\n");
			j++;

		}
		out.close();
	}
	
	public static void GenerateTrainDataWithNormalDataInTestFile(int numTrainingSequences, int normalDataInTest,  int length,
			NormalClassDataGeneration normal) throws Exception {
		
// Writing in the Training Data
		FileWriter fstream = new FileWriter("HSMMtraining.txt");
		BufferedWriter out = new BufferedWriter(fstream);
		
		String demoLine = "";
		int j =4;
		for(int i=0;i<length;i++){
			demoLine = demoLine +   j + " ";
			j = j + 3;
		}
		
		out.write(""+numTrainingSequences);
		out.write("\n");

		// Duration 1 means Actually Duration 0
		 j = 0;
		while (j < numTrainingSequences) {
			String result = "";
			int i = 0;
			int oldState = normal.initialState;
			int oldDuration = normal.initialDuration;

			while (i < length) {
				int nextState = getNewState(oldState, oldDuration, normal);
				int nextDuration = getNewDuration(oldState, oldDuration, normal);
				String observationSymbol = getObservationSymbol(nextState,
						normal);

				oldState = nextState;
				oldDuration = nextDuration;
				result = result + test.getInteger(observationSymbol) + " ";
				i++;

			}
			System.out.println(result);
			result = result.substring(0, result.length() - 1);
			out.write(""+length);
			out.write("\n");
			out.write(demoLine);
			out.write("\n");
			out.write(result);
			out.write("\n");
			j++;

		}
		out.close();
		
 // Writing Normal Data in the Test File
		// Writing in the Training Data
				 fstream = new FileWriter("Test",true);
				 out = new BufferedWriter(fstream);

				// Duration 1 means Actually Duration 0
			       j = 0;
				while (j < normalDataInTest) {
					String result = "";
					int i = 0;
					int oldState = normal.initialState;
					int oldDuration = normal.initialDuration;

					while (i < length) {
						int nextState = getNewState(oldState, oldDuration, normal);
						int nextDuration = getNewDuration(oldState, oldDuration, normal);
						String observationSymbol = getObservationSymbol(nextState,
								normal);

						oldState = nextState;
						oldDuration = nextDuration;
						result = result + test.getInteger(observationSymbol) + " ";
						i++;

					}
					System.out.println(result);
					result = result.substring(0, result.length() - 1);
					out.write(""+length);
					out.write("\n");
					out.write(demoLine);
					out.write("\n");
					out.write(result);
					out.write("\n");
					j++;

				}
				out.close();

	}
	
	public static void GenerateTrainDataWithNormalDataInTestFileForOtherModels(int numTrainingSequences, int normalDataInTest,  int length,
			NormalClassDataGeneration normal) throws Exception {
		
// Writing in the Training Data
		FileWriter fstream = new FileWriter("HMMtraining.txt");
		BufferedWriter out = new BufferedWriter(fstream);
		
		// Duration 1 means Actually Duration 0
		 int j = 0;
		while (j < numTrainingSequences) {
			String result = "";
			int i = 0;
			int oldState = normal.initialState;
			int oldDuration = normal.initialDuration;

			while (i < length) {
				int nextState = getNewState(oldState, oldDuration, normal);
				int nextDuration = getNewDuration(oldState, oldDuration, normal);
				String observationSymbol = getObservationSymbol(nextState,
						normal);

				oldState = nextState;
				oldDuration = nextDuration;
				result = result + test.getInteger(observationSymbol) + " ";
				i++;

			}
			System.out.println(result);
			result = result.substring(0, result.length() - 1);
			out.write(result);
			out.write("\n");
			j++;

		}
		out.close();
		
 // Writing Normal Data in the Test File
		// Writing in the Training Data
				 fstream = new FileWriter("HMMtesting.txt",true);
				 out = new BufferedWriter(fstream);

				// Duration 1 means Actually Duration 0
			       j = 0;
				while (j < normalDataInTest) {
					String result = "";
					int i = 0;
					int oldState = normal.initialState;
					int oldDuration = normal.initialDuration;

					while (i < length) {
						int nextState = getNewState(oldState, oldDuration, normal);
						int nextDuration = getNewDuration(oldState, oldDuration, normal);
						String observationSymbol = getObservationSymbol(nextState,
								normal);

						oldState = nextState;
						oldDuration = nextDuration;
						result = result + test.getInteger(observationSymbol) + " ";
						i++;

					}
					System.out.println(result);
					result = result.substring(0, result.length() - 1);
					out.write(result);
					out.write("\n");
					j++;

				}
				out.close();

	}

	public static String getObservationSymbol(int nextState,
			NormalClassDataGeneration normal) {

		double d[] = new double[normal.possibleObservationsList.size()];
		for (int i = 0; i < d.length; i++) {
			d[i] = normal.probofObservationGivenState.probObsGivenState[i][nextState];
		}

		int index = Util.getCorrectIndexForProbabilityDistribution(d);
		return normal.possibleObservationsList.get(index);

	}

	public static int getNewDuration(int oldState, int oldDuration,
			NormalClassDataGeneration normal) {

		if (oldDuration > 0) {
			return oldDuration - 1;
		} else {
			double d[] = new double[normal.totalDuration];
			for (int i = 0; i < d.length; i++) {
				d[i] = normal.probDurationGivenState.probDurationGivenState[i][oldState];
			}
			return Util.getCorrectIndexForProbabilityDistribution(d);
		}
	}

	public static int getNewState(int oldState, int oldDuration,
			NormalClassDataGeneration normal) {

		if (oldDuration > 0) {
			return oldState;
		} else {
			double d[] = new double[normal.hiddenStateList.size()];
			for (int i = 0; i < d.length; i++) {
				d[i] = normal.probOfNewStateGivenOldState.probNewStateGivenOldState[i][oldState];
			}
			return Util.getCorrectIndexForProbabilityDistribution(d);
		}
	}
}
