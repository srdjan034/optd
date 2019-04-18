# OPTD (Optimized Partial Trajectory Decomposition)

The OPTD introduces generic, scalable and adaptive load balancing parallel Lagrangian particle tracking (LPT) approach, especially suited for the problems where particles are characterized by highly variable computation time. The algorithm introduces master/slave model dubbed Partial Trajectory Decomposition (PTD), in which a certain number of processes produce partial trajectories and put them into the shared queue, while the remaining processes simulate actual particle motion using previously generated partial trajectories. Our approach also introduces meta-heuristics for determining the optimal values of partial trajectory length, chunk size and the number of processes acting as producers/consumers, for the given total number of participating processes (Optimized Partial Trajectory Decomposition, OPTD). The optimization process employs a surrogate model to estimate the simulation time. The surrogate is based on historical data and uses a coupled machine learning model, consisting of classification and regression phases.

# MPI implementation

We have implemented OPTD in C, using standard MPI for message passing and benchmarked on a model of <sup>220</sup>Rn progeny behaviour in the diffusion chamber, where particle motion is characterized by an exponential lifetime distribution and Maxwell velocity distribution.

### Run MPI implementation on cluster

It is necessary to start the cluster job by sumbiting the PBS script <i>run_ptd_job.sub</i>. 

```
#!/bin/sh
#PBS -N OPTD
#PBS -q batch
#PBS -l nodes=1:ppn=16

cd $PBS_O_WORKDIR

mpicc -O2 ptd.c -o ptd -lm

time mpirun ./ptd
```

Number of nodes and processors per nodes can be set in file <i>run_ptd_job.sub</i>. Run the script using <i>qsub run_ptd_job.sub</i>


All PTD parameter values can be set in file <i>PTD_optimal_values.json</i>.

```
{
    "nodeCount" : 1,
    "processorsPerNode" : 16,
    "producerCount" : 14,
    "consumerCount" : 1,
    "partialTrajectoryLength" : 10000,
    "chunkSize" : 10000,
    "particleCount" : 1024
}
```


### Predicting optimal PTD parametter values in file <i>PTD_optimal_values.json</i>

We carried out testing on the cluster containing 22 worker nodes, each equipped with dual Intel Xeon E5-2670 @ 2.60GHz (16 cores per node) with Infiniband QDR interconnection. Based od historical testing data which contains simulation time for each historical PTD parameter we have created surrogat model based on machine learning. The surrogate model learns from the acquired historical data to predict the simulation time for various PTD parameter combinations. In  order  to  benchmark  our  novel  OPTD  approach,  we have  employed  the  EA  optimization  with  elitism  enabled, using  MOEA  optimization  framework  by  Hadka  (2016), to  determine  the  optimal  PTD  parameters. The surrogate model was used to evaluate various combinations of PTD parameters. 

In order to predict optimal PTD parametter values we have created JAVA maven project <i>PredictOptimalValuesOfPTDparameters</i>. All EA optimization parameter values can be set in file <i>conf.json</i>. When starting the program, the following input arguments should be sent: <i>(i)</i> number od nodes, <i>(ii)</i> processors per node, <i>(iii)</i> number of particles. The result of program execution is a file <i>PTD_optimal_values.json</i>.

# Licence and contributing

All code is licensed under the MIT license. We enthusiastically welcome contributions and feedback.
