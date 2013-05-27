package org.nasa;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

public class CreateFactorGraphForHSMM {

	public static void main(String args[]) {

	}

	public static void createFactorGraph(NormalClassDataGeneration normal)
			throws Exception {

		FileWriter fstream = new FileWriter("hsmm_factor_graph_init.fg");
		BufferedWriter out = new BufferedWriter(fstream);

		out.write("5");
		out.write("\n");out.write("\n");

		List<String> list = getDO(normal);
		for (String s : list) {
			out.write(s);
			out.write("\n");
		}
		out.write("\n");

		list = getAO(normal);
		for (String s : list) {
			out.write(s);
			out.write("\n");
		}
		out.write("\n");

		list = getD1GivenA1D0(normal);
		for (String s : list) {
			out.write(s);
			out.write("\n");
		}
		out.write("\n");

		list = getA1GivenA0D0(normal);
		for (String s : list) {
			out.write(s);
			out.write("\n");
		}
		out.write("\n");

		list = getO1GivenA1(normal);
		for (String s : list) {
			out.write(s);
			out.write("\n");
		}
		out.write("\n");
		out.close();

	}

	public static double[] createOneDimRandomMatrix(int dim) {
		double d[] = new double[dim];

		double sum = 0;

		for (int i = 0; i < dim; i++) {
			d[i] = Math.random();
			sum = sum + d[i];
		}

		for (int i = 0; i < dim; i++) {
			d[i] = d[i] / sum;

		}
		return d;
	}

	public static double[][][] createThreeDimRandomMatrix(int dimA, int dimB,
			int dimC, int flag) {

		double d[][][] = new double[dimA][dimB][dimC];

		if (flag == 0) {
			// This is duration given old duration and old state
			// dimA is the new duration....dimB is the old state and dimC is the
			// old duration
			for (int i = 0; i < dimA; i++) {
				for (int j = 0; j < dimB; j++) {
					for (int k = 0; k < dimC; k++) {
						// old duration is not zero
						if (k != 0) {
                          // new duration is one minus the old duration
							if((i-1)==k){
								d[i][j][k] = 1;
							}else{
								d[i][j][k] = 0;
							}
							
						}else{
							d[i][j][k] = Math.random();
						}
						
					}
				}
			}
		}

		if (flag == 1) {
			// This is new state given old state and old duration
			// dimA is the new state....dimB is the old state and dimC is the
			// old duration
			for (int i = 0; i < dimA; i++) {
				for (int j = 0; j < dimB; j++) {
					for (int k = 0; k < dimC; k++) {
						// if the old duration is not zero
						if (k != 0) {
							// if the old state and new state are same
							if (i == j)
								{
								d[i][j][k] = 1;
								}
							else
								d[i][j][k] = 0;
						} else {
							d[i][j][k] = Math.random();
						}
					}
				}
			}
		}

		double sum = 0;
		for (int j = 0; j < dimB; j++) {
			for (int k = 0; k < dimC; k++) {
				sum = 0;
				for (int i = 0; i < dimA; i++) {
					sum = sum + d[i][j][k];
				}

				for (int i = 0; i < dimA; i++) {
					if(sum!=0)
					d[i][j][k] = d[i][j][k] / sum;

				}
			}
		}

		return d;
	}

	public static double[][] createTwoDimRandomMatrix(int dimA, int dimB) {
		double d[][] = new double[dimA][dimB];
		double sum = 0;

		// Just Generating the data
		for (int i = 0; i < dimA; i++) {
			for (int j = 0; j < dimB; j++) {
				d[i][j] = Math.random();
			}
		}

		// Making sure it sums upto one
		for (int i = 0; i < dimB; i++) {
			sum = 0;
			for (int j = 0; j < dimA; j++) {
				sum = sum + d[j][i];
			}

			for (int j = 0; j < dimA; j++) {
				d[j][i] = d[j][i] / sum;
			}
		}

		return d;
	}

	public static int getNumberofNonZerosInOneDim(double d[]) {
		int count = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] != 0) {
				count++;
			}
		}

		return count;
	}

	public static int getNumberofNonZerosInTwoDim(double d[][]) {
		int count = 0;

		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[0].length; j++) {
				if (d[i][j] != 0) {
					count++;
				}
			}
		}

		return count;
	}

	public static int getNumberofNonZerosInThreeDim(double d[][][]) {
		int count = 0;

		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[0].length; j++) {
				for (int k = 0; k < d[0][0].length; k++) {
					if (d[i][j][k] != 0) {
						count++;
					}
				}
			}
		}

		return count;
	}

	public static List<String> getNonZeroElementsWithIndexInOneDim(double d[]) {
		List<String> data = new ArrayList<String>();

		for (int i = 0; i < d.length; i++) {
			if (d[i] != 0) {
				data.add(i + "\t" + d[i]);
			}
		}

		return data;
	}

	public static List<String> getNonZeroElementsWithIndexInTwoDim(double d[][]) {
		List<String> data = new ArrayList<String>();
		int counter = 0;

		for (int i = 0; i < d[0].length; i++) {
			for (int j = 0; j < d.length; j++) {
				if (d[j][i] != 0) {
					data.add(counter + "\t" + d[j][i]);
				}
				counter++;
			}
		}
		return data;
	}

	public static List<String> getNonZeroElementsWithIndexInThreeDim(
			double d[][][]) {
		List<String> data = new ArrayList<String>();
		int counter = 0;

		for (int k = 0; k < d[0][0].length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d.length; j++) {
					if (d[j][i][k] != 0) {
						data.add(counter + "\t" + d[j][i][k]);
					}
					counter++;
				}
			}
		}
		return data;
	}

	public static List<String> getDO(NormalClassDataGeneration normal) {
		List<String> outputData = new ArrayList<String>();
		outputData.add("1");
		outputData.add("0");
		outputData.add("" + (normal.totalDuration + 1));

		double d[] = createOneDimRandomMatrix(normal.totalDuration + 1);
		outputData.add("" + getNumberofNonZerosInOneDim(d));

		outputData.addAll(getNonZeroElementsWithIndexInOneDim(d));
		return outputData;
	}

	public static List<String> getAO(NormalClassDataGeneration normal) {
		List<String> outputData = new ArrayList<String>();
		outputData.add("1");
		outputData.add("1");
		outputData.add("" + normal.hiddenStateList.size());

		double d[] = createOneDimRandomMatrix(normal.hiddenStateList.size());
		outputData.add("" + getNumberofNonZerosInOneDim(d));

		outputData.addAll(getNonZeroElementsWithIndexInOneDim(d));
		return outputData;
	}

	public static List<String> getD1GivenA1D0(NormalClassDataGeneration normal) {
		List<String> outputData = new ArrayList<String>();
		outputData.add("3");
		outputData.add("2\t3\t0");
		outputData.add((normal.totalDuration + 1) + "\t"
				+ normal.hiddenStateList.size() + "\t"
				+ (normal.totalDuration + 1));

		double d[][][] = createThreeDimRandomMatrix(normal.totalDuration + 1,
				normal.hiddenStateList.size(), normal.totalDuration + 1,0);
		outputData.add("" + getNumberofNonZerosInThreeDim(d));

		outputData.addAll(getNonZeroElementsWithIndexInThreeDim(d));
		return outputData;
	}

	public static List<String> getA1GivenA0D0(NormalClassDataGeneration normal) {
		List<String> outputData = new ArrayList<String>();
		outputData.add("3");
		outputData.add("3\t1\t0");
		outputData.add((normal.hiddenStateList.size()) + "\t"
				+ normal.hiddenStateList.size() + "\t"
				+ (normal.totalDuration + 1));

		double d[][][] = createThreeDimRandomMatrix(
				normal.hiddenStateList.size(), normal.hiddenStateList.size(),
				normal.totalDuration + 1,1);
		outputData.add("" + getNumberofNonZerosInThreeDim(d));

		outputData.addAll(getNonZeroElementsWithIndexInThreeDim(d));
		return outputData;
	}

	public static List<String> getO1GivenA1(NormalClassDataGeneration normal) {
		List<String> outputData = new ArrayList<String>();
		outputData.add("2");
		outputData.add("4\t3");
		outputData.add((normal.possibleObservationsList.size()) + "\t"
				+ normal.hiddenStateList.size());

		double d[][] = createTwoDimRandomMatrix(
				normal.possibleObservationsList.size(),
				normal.hiddenStateList.size());
		outputData.add("" + getNumberofNonZerosInTwoDim(d));

		outputData.addAll(getNonZeroElementsWithIndexInTwoDim(d));
		return outputData;
	}

}
