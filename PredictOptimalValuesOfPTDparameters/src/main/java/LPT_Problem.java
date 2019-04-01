import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;

import weka.classifiers.Classifier;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ConverterUtils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.gui.experiment.GeneratorPropertyIteratorPanel;
import weka.core.Attribute;

public class LPT_Problem extends AbstractProblem {

	Classifier executionTimeClassifier = null;
	Classifier queueOverflowClassifier = null;

	Instance executionTimeInstance = null;
	Instance queueOverflowInstance = null;
	
	int nProc;
	int populationIndex = 0;
	int populationSize;
	
	double minPredictionTimeInGeneration = Double.MAX_VALUE;
	double[] bestSolutionValues;

	Attribute nProcAttribute = new Attribute("nProc");
	Attribute nConsumersAttribute = new Attribute("nConsumers");
	Attribute nChunkAttribute = new Attribute("nChunk");
	Attribute TrajectoryLengthAttribute = new Attribute("nTrajectoryLength");
	Attribute timeAttribute = new Attribute("time");
	Attribute queueOverflowAttribute = null;

	DataSource trainigSet;
	Instance instance;

	PredictIdealSimulationTime idealTimePrediction = null;

	private LPT_Problem(int numberOfVariables, int numberOfObjectives, double[] bestSolutionValues) {
		super(numberOfVariables, numberOfObjectives);

		try {

		} catch (Exception e) {
			e.printStackTrace();
		}

		this.bestSolutionValues = bestSolutionValues;
		initialiseInstance();

		try {
			idealTimePrediction = new PredictIdealSimulationTime("mlModel\\IdealTimes.csv");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private void initialiseInstance() {
		Instances data = null;
		ArrayList<Attribute> attributeList1 = new ArrayList<Attribute>();

		attributeList1.add(nProcAttribute);
		attributeList1.add(nConsumersAttribute);
		attributeList1.add(nChunkAttribute);
		attributeList1.add(TrajectoryLengthAttribute);
		attributeList1.add(timeAttribute);

		data = new Instances("ExecutionTimePredictionInstances", attributeList1, 5);
		data.setClassIndex(data.numAttributes() - 1);

		executionTimeInstance = new DenseInstance(data.numAttributes());
		executionTimeInstance.setDataset(data);

		List<String> queueOverflow_values = new ArrayList<String>(2);
		queueOverflow_values.add("No");
		queueOverflow_values.add("Yes");

		queueOverflowAttribute = new Attribute("queueOverflow", queueOverflow_values);

		ArrayList<Attribute> attributeList2 = new ArrayList<Attribute>();
		attributeList2.add(nProcAttribute);
		attributeList2.add(nConsumersAttribute);
		attributeList2.add(nChunkAttribute);
		attributeList2.add(TrajectoryLengthAttribute);
		attributeList2.add(queueOverflowAttribute);

		data = new Instances("QueueOverflowPredictionInstances", attributeList2, 0);
		data.setClassIndex(data.numAttributes() - 1);

		queueOverflowInstance = new DenseInstance(data.numAttributes());
		queueOverflowInstance.setDataset(data);
	}

	public LPT_Problem(int nProc, String ExecutionTimePredictionModelFilePath,
			String QueueOverflowPredictionModelFilePath, int populationSize, double[] bestSolutionValues)
					throws Exception {
		this(4, 1, bestSolutionValues);

		this.nProc = nProc;

		File file = new File(ExecutionTimePredictionModelFilePath);

		if (!file.exists())
			throw new FileNotFoundException(ExecutionTimePredictionModelFilePath + "path is not valid!");

		file = new File(QueueOverflowPredictionModelFilePath);

		if (!file.exists())
			throw new FileNotFoundException(QueueOverflowPredictionModelFilePath + "path is not valid!");

		executionTimeClassifier = (Classifier) weka.core.SerializationHelper.read(ExecutionTimePredictionModelFilePath);
		queueOverflowClassifier = (Classifier) weka.core.SerializationHelper.read(QueueOverflowPredictionModelFilePath);

		this.populationSize = populationSize;

	}

	public void evaluate(Solution solution) {

		double[] values = EncodingUtils.getReal(solution);
		double[] f = new double[numberOfObjectives];
		double predictedTime = Double.MAX_VALUE;

		for (int nProc = 15; nProc < 260; nProc += 5) {

			executionTimeInstance.setValue(nProcAttribute, nProc);
			executionTimeInstance.setValue(nConsumersAttribute, 1);
			executionTimeInstance.setValue(nChunkAttribute, 100000);
			executionTimeInstance.setValue(TrajectoryLengthAttribute, 100000);

			try {
				predictedTime = executionTimeClassifier.classifyInstance(executionTimeInstance);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		// Set solution values to weka instance
		queueOverflowInstance.setValue(nProcAttribute, values[0]);
		queueOverflowInstance.setValue(nConsumersAttribute, values[1]);
		queueOverflowInstance.setValue(nChunkAttribute, values[2]);
		queueOverflowInstance.setValue(TrajectoryLengthAttribute, values[3]);

		double predictedQueueOverflow = 1.0;

		// Queue overflow prediction
		try {
			synchronized (queueOverflowClassifier) {
				predictedQueueOverflow = queueOverflowClassifier.classifyInstance(queueOverflowInstance);
			}
		} catch (Exception e) {
			System.out.println(e.getMessage());
		}

		// Total simulation time prediction
		if (predictedQueueOverflow < 1) {

			executionTimeInstance.setValue(nProcAttribute, values[0]);
			executionTimeInstance.setValue(nConsumersAttribute, values[1]);
			executionTimeInstance.setValue(nChunkAttribute, values[2]);
			executionTimeInstance.setValue(TrajectoryLengthAttribute, values[3]);

			try {
				synchronized (executionTimeClassifier) {
					predictedTime = executionTimeClassifier.classifyInstance(executionTimeInstance);

					if (predictedTime < idealTimePrediction.getPredictedTime(values[0])) {
						predictedTime = Double.MAX_VALUE;
					}
				}

			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		f[0] = predictedTime;
		solution.setObjectives(f);

		populationIndex++;

		if (populationIndex % populationSize == 0) {
			bestSolutionValues[populationIndex / populationSize - 1] = minPredictionTimeInGeneration;
		} else {
			if (minPredictionTimeInGeneration > predictedTime) {
				minPredictionTimeInGeneration = predictedTime;
			}
		}
	}

	public Solution newSolution() {
		Solution solution = new Solution(getNumberOfVariables(), getNumberOfObjectives());

		solution.setVariable(0, new RealVariable(1, nProc));
		solution.setVariable(1, new RealVariable(1, 5));
		solution.setVariable(2, new RealVariable(256, 262000));
		solution.setVariable(3, new RealVariable(1000, 400000));

		return solution;
	}

}
