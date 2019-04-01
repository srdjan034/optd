import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.json.JSONException;
import org.json.JSONObject;
import org.moeaframework.Executor;
import org.moeaframework.analysis.plot.Plot;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Solution;

public class PredictOptimalParameterValues {
	public static void main(String[] args) {
		try {
			
			int nNode = Integer.parseInt(args[0]);
			int ppn = Integer.parseInt(args[1]);
			int nParticle  = Integer.parseInt(args[2]);
			
			// Read config file
			JSONObject configJson = new JSONObject(new String(Files.readAllBytes(Paths.get("conf.json"))));
			String AlgorithmGA = configJson.getString("AlgorithmGA");
			int nEvaluation = configJson.getInt("NumberOfEvulations");
			int PopulationSize = configJson.getInt("PopulationSize");
			double sbxRate = configJson.getDouble("sbx.rate");
			double sbxDistributionIndex = configJson.getDouble("sbx.distributionIndex");
			double pmRate = configJson.getDouble("pm.rate");
			double pmDistributionIndex = configJson.getDouble("pm.distributionIndex");
			String QueueOverflowPredictionModelPath = configJson.getString("QueueOverflowPredictionModelPath");
			String ExecutionTimePredictionModelPath = configJson.getString("ExecutionTimePredictionModelPath");

			double[] bestSolutionValuesY = new double[nEvaluation / PopulationSize + 1];

			// Run
			List<NondominatedPopulation> results = new Executor()
					.withProblemClass(LPT_Problem.class, nNode * ppn, ExecutionTimePredictionModelPath,
							QueueOverflowPredictionModelPath, PopulationSize, bestSolutionValuesY)
					.withProperty("populationSize", PopulationSize).withAlgorithm(AlgorithmGA)
					.withProperty("sbx.rate", sbxRate).withProperty("sbx.distributionIndex", sbxDistributionIndex)
					.withProperty("pm.rate", pmRate).withProperty("pm.distributionIndex", pmDistributionIndex)
					.withMaxEvaluations(nEvaluation).runSeeds(1);

			// Print best solution and create configuration file
			for (NondominatedPopulation result : results)
				for (Solution solution : result) {
					int nConsumer = (int) (Double.parseDouble(solution.getVariable(1).toString()));
					int nChunk = (int) (Double.parseDouble(solution.getVariable(2).toString()));
					int pathLength = (int) (Double.parseDouble(solution.getVariable(3).toString()));
					try {
						writeConfFile(nNode, ppn, nParticle, nConsumer, nChunk, pathLength);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
		} catch (JSONException e1) {
			e1.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	private static void writeConfFile(int nNode, int ppn, int nParticle, int nConsumer, int nChunk, int pathLength)
			throws IOException {
		int nProc = nNode * ppn;

		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("PTD_optimal_values.json")));

		bw.write("{\n");
		bw.write("\t\"nodeCount\": " + nNode + ",\n");
		bw.write("\t\"processorsPerNode\": " + ppn + ",\n");
		bw.write("\t\"producerCount\": " + (nProc - nConsumer - 1) + ",\n");
		bw.write("\t\"consumerCount\": " + nConsumer + ",\n");
		bw.write("\t\"partialTrajectoryLength\": " + pathLength + ",\n");
		bw.write("\t\"chunkSize\": " + nChunk + ",\n");
		bw.write("\t\"particleCount\": " + nParticle + "\n");
		bw.write("}\n");

		bw.close();

	}
}
