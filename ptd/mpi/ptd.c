#include "ptd.h"

int main(int argc, char *argv[]) 
{
    int rank, size; 

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Conf * conf = readConf("PTD_optimal_values.json");
    
    BUFFER_SIZE = conf->chunkSize / size;
   
    if(conf->consumerCount + conf->producerCount + 1 == size)
    {
        initialize();
        
        if(rank == 0) // Run Shared Queue
        {   
            runQueue(conf, conf->consumerCount, conf->producerCount);
        }
        else if(rank <= conf->consumerCount) // Run Particle Simulators
        {
            simulateParticles(rank, conf);
        }
        else if(conf->consumerCount < rank) // Run Partial Trajectory Producers
        {
            generatePathTrajektories(rank, conf->partialTrajectoryLength, conf);
        }
    }
    else if(rank == 0)
    {
        printf("Error in file PTD_optimal_values.json");
    }
    
    MPI_Finalize();
    
    return 0;
}
