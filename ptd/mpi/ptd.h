#include "model.h"

void runQueue(Conf * conf, int consumerCount, int producerCount);
void simulateParticles( int rank, Conf * conf);
void generatePathTrajektories(int rank,int partialTrajectoryLength, Conf * conf);
void createMpiTypeForPartialTrajectory(MPI_Datatype  * dt_partial_trajectory);
void createMpiTypeForParticle(MPI_Datatype  * dt_particle);

void generatePathTrajektories(int rank, int partialTrajectoryLength, Conf * conf)
{   
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    struct drand48_data buffer_seed;
    srand48_r(rank, &buffer_seed);
    
    MPI_Status status;

    int gauss_kon = 1;
    double gauss_y = 0.0;
    int endOfSimulationFlag = 0;
    int chunkSize = BUFFER_SIZE;
    int rankSource;
    int i;
    
    Partial_trajectory partial_trajectory;
    
    while(endOfSimulationFlag == 0)
    {   
        MPI_Recv(&rankSource, 1, MPI_INT, MPI_ANY_SOURCE, NEW_CHUNK, MPI_COMM_WORLD, &status);
        
        for(i = 0; i < chunkSize && endOfSimulationFlag == 0; i++)
        {
            generatePartialTrajectory(&partial_trajectory, &gauss_kon, &gauss_y, 
                    partialTrajectoryLength, &buffer_seed);

            MPI_Send(&partial_trajectory, 1, dt_partial_trajectory, 0, 
                    PARTIAL_TRAJECTORY, MPI_COMM_WORLD); 

            MPI_Iprobe(0, SIMULATION_END, MPI_COMM_WORLD, &endOfSimulationFlag, &status);
        }
    }
}

void simulateParticles( int rank, Conf * conf)
{ 
    Partial_trajectory * partialTrajectory = (Partial_trajectory *)malloc(sizeof(Partial_trajectory));
    Particle particle;
    MPI_Request request;
    MPI_Status status;
    
    MPI_Datatype dt_particle;         
    createMpiTypeForParticle(&dt_particle);
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory); 
    
    int nAir, nBottom, nWall, nTop;
    nAir = nBottom = nWall = nTop = 0;
    
    int chunkSize = BUFFER_SIZE * conf->producerCount;
    
    while(1)
    {
        //Request for new particle.
        MPI_Send(&rank, 1, MPI_INT, 0, PARTICLE_REQUEST, MPI_COMM_WORLD);
        
        // Receive particle.
        MPI_Recv(&particle, 1, dt_particle, 0, PARTICLE, MPI_COMM_WORLD, &status);
        
        if(particle.status == NONE) 
        {
            break;
        }

        while(1)
        {       
            // Send request for partial trajectory.
            MPI_Send(&rank, 1, MPI_INT, 0, PARTIAL_TRAJECTORY_REQUEST, 
                    MPI_COMM_WORLD);

            MPI_Recv(partialTrajectory, 1, dt_partial_trajectory, 0, 
                    PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &status);
            
            if( checkBoundingBox(&particle, partialTrajectory) != AIR)
            {   
                struct drand48_data seed;
                
                initializeSeed(&seed, partialTrajectory);
                
                if(checkDepositionPlace(&particle, conf->partialTrajectoryLength, &seed) == 1)
                {
                    switch(particle.status)
                    {
                        case TOP:
                            nTop++;
                        break;
                        case BOTTOM:
                            nBottom++;
                        break;
                        case WALL:
                            nWall++;
                        break;
                        case AIR:
                            nAir++;
                        break;
                    }
                    
                    break;
                }
                else
                {
                    particle.life_time += partialTrajectory->lifeTimeStepSum;
                    particle.x += partialTrajectory->xLast;
                    particle.y += partialTrajectory->yLast;
                    particle.z += partialTrajectory->zLast;
                }
            }
            else
            {
                particle.life_time += partialTrajectory->lifeTimeStepSum;
                particle.x += partialTrajectory->xLast;
                particle.y += partialTrajectory->yLast;
                particle.z += partialTrajectory->zLast;
            }
        }
    }
    
    sendDepositionStatistics(nAir, nBottom, nTop, nWall);
}

void runQueue(Conf * conf, int consumerCount, int producerCount)
{  
    clock_t begin = clock();
    
    int running = 1;
    int i, j;
    MPI_Request request;
    int indexProducer = -1;
    int indexConsumer = -1;
    int processedPathCount = 0;
    
    MPI_Datatype dt_particle;
    createMpiTypeForParticle(&dt_particle);
    
    int particleCountTmp = conf->particleCount;
    int consumerCountTmp = conf->consumerCount;
    
    int nextParticleIndex = 0;
    
    int bufferParticles[conf->consumerCount];
    MPI_Request particleRequests[conf->consumerCount];
    int outCountParticles;
    int array_of_indices_particles[conf->consumerCount];
    MPI_Status array_of_statuses_particles[conf->consumerCount]; 
    
    struct drand48_data buffer_seed;
    srand48_r(time(NULL), &buffer_seed);
    

    // Generate particles in chamber
    Particle particles[conf->particleCount]; 
    for(i = 0; i < conf->particleCount; i++)
        initialize_particle(&particles[i], i, tau, r, &buffer_seed);
    
    // Subscribe for particle requests.
    for(i = 0; i < conf->consumerCount; i++)
    {   
        MPI_Irecv(&bufferParticles[i], 1, MPI_INT, i+1, PARTICLE_REQUEST, 
                MPI_COMM_WORLD, &particleRequests[i]);
    } 
    
    /*
        Prijavi se da prihvatas zahteve za parcijalnim putanjama.
     */
    
    int bufferConsumerTrajectoryRequest[consumerCount];
    MPI_Request consumerRequests[consumerCount];
    MPI_Status  consumerStatuses[consumerCount];
    int outCountConsumers = 0;
    int array_of_indices_consumers[consumerCount];
    
    for(i = 0; i < consumerCount; i++)
    {   
        MPI_Irecv(&bufferConsumerTrajectoryRequest[i], 1, MPI_INT, i + 1, 
                PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[i]);
    }

    /* 
        Prijavi se da prihvatas parcijalne putanje od proizvodjaca.
    */
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    Partial_trajectory * partial_trajectory;
    
    Partial_trajectory bufferProducer[producerCount];
    MPI_Request producerRequests[producerCount];
    MPI_Status  producerStatuses[producerCount];
    int outCountProducers = 0;
    int array_of_indices_producers[producerCount];
    
    for(i = consumerCount, j = 0; i < consumerCount + producerCount; i++, j++)
    {   
        MPI_Irecv(&bufferProducer[j], 1, dt_partial_trajectory, i + 1, 
                PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[j]);
    }
    
    int chunkSize = BUFFER_SIZE * conf->producerCount;
    int countTrajectories = chunkSize;
    
    int n = 1;
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Send(&n, 1, MPI_INT, i, NEW_CHUNK, MPI_COMM_WORLD);
    
    /*
        Shared Queue simulation.
     */
    int count;
    MPI_Status status;
    while(running)
    {
        MPI_Testsome(producerCount, producerRequests, &outCountProducers, 
                array_of_indices_producers, producerStatuses); 
        
        // Send all partial trajectories
        while(outCountProducers > 0 && running)
        {
            if(outCountConsumers == 0 && running)
                MPI_Testsome(consumerCount, consumerRequests, &outCountConsumers, 
                        array_of_indices_consumers, consumerStatuses);

            // If there are requests for partial trajectories
            if(outCountConsumers > 0)
            {
                // Handle all requests for partial trajectories.
                while(outCountProducers && outCountConsumers && running)
                {   
                    outCountProducers--;
                    outCountConsumers--;

                    indexProducer = array_of_indices_producers[outCountProducers];
                    indexConsumer = array_of_indices_consumers[outCountConsumers];

                    partial_trajectory = &bufferProducer[indexProducer];

                    processedPathCount++;
                    
                    countTrajectories--;
            
                    if(countTrajectories == 0)
                    {
                        countTrajectories = chunkSize;

                        for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
                        {
                            MPI_Send(&n, 1, MPI_INT, i, NEW_CHUNK, MPI_COMM_WORLD);
                        }
                    }

                    MPI_Send(partial_trajectory, 1, dt_partial_trajectory, 
                            consumerStatuses[outCountConsumers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY, MPI_COMM_WORLD);

                    MPI_Irecv(&bufferConsumerTrajectoryRequest[indexConsumer], 
                            1, MPI_INT, consumerStatuses[outCountConsumers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY_REQUEST, 
                            MPI_COMM_WORLD, &consumerRequests[indexConsumer]);

                    MPI_Irecv(&bufferProducer[indexProducer], 1, dt_partial_trajectory, 
                            producerStatuses[outCountProducers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY, MPI_COMM_WORLD, 
                            &producerRequests[indexProducer]);
                }
            }

            // Check if there are requests for particles.
            MPI_Testsome(conf->consumerCount, particleRequests, &outCountParticles, 
                        array_of_indices_particles, array_of_statuses_particles);

            // Handle all requests for particles.
            for(i = 0; i < outCountParticles; i++)
            { 
                 // Check if there is more particles in chamber.
                if(nextParticleIndex < particleCountTmp)
                {
                    MPI_Send(&(particles[nextParticleIndex]), 1, dt_particle, 
                                array_of_statuses_particles[i].MPI_SOURCE, PARTICLE, 
                                MPI_COMM_WORLD);

                    MPI_Irecv(&bufferParticles[array_of_indices_particles[i]], 1, MPI_INT, 
                                array_of_statuses_particles[i].MPI_SOURCE, PARTICLE_REQUEST, 
                                MPI_COMM_WORLD, &particleRequests[array_of_indices_particles[i]]);

                    nextParticleIndex++;
                }
                else // If all particles are sent to simulators, send status NONE.
                {
                    Particle p;
                    p.status = NONE;
                    MPI_Send(&p, 1, dt_particle, array_of_statuses_particles[i].MPI_SOURCE, 
                               PARTICLE, MPI_COMM_WORLD);

                    consumerCountTmp--;

                    // If all particles are deposited.
                    if(consumerCountTmp == 0)
                        running = 0;
                }
            }
        }
    }

    // Send messages for termination of simulation.
    n = 1;
    
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Send(&n, 1, MPI_INT, i, SIMULATION_END, MPI_COMM_WORLD);
    
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
    {
        MPI_Send(&n, 1, MPI_INT, i, NEW_CHUNK, MPI_COMM_WORLD);
    }
    
    getDepositionStatistics(conf);
    
    printf("\nTotal simulation time = %lf s.\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    
}

void createMpiTypeForParticle(MPI_Datatype  * dt_particle)
{
    int          blocklengths[3];
    MPI_Aint     offsets[3], intExtent, doubleExtent;
    MPI_Datatype types[3];
    
    offsets[0] = 0;
    types[0] = MPI_DOUBLE;
    blocklengths[0] = 5;
    
    MPI_Type_extent(MPI_DOUBLE, &doubleExtent);
    offsets[1] = 5 * doubleExtent;
    types[1] = MPI_INT;
    blocklengths[1] = 1;

    MPI_Type_extent(MPI_INT, &intExtent);
    offsets[2] = intExtent +  5 * doubleExtent;
    types[2] = MPI_SHORT_INT;
    blocklengths[2] = 1;

    MPI_Type_create_struct(3, blocklengths, offsets, types, dt_particle);
    MPI_Type_commit(dt_particle);  
}

void createMpiTypeForPartialTrajectory(MPI_Datatype  * dt_partial_trajectory)
{
    int          blocklengths[3];
    MPI_Aint     offsets[3], extent;
    MPI_Datatype types[3];
    
    types[0] = MPI_DOUBLE;
    blocklengths[0] = 10;
    offsets[0] = 0;
    
    types[1] = MPI_UNSIGNED_LONG_LONG;
    blocklengths[1] = 1;
    MPI_Type_extent(MPI_DOUBLE, &extent); 
    offsets[1] = 10 * extent;
    
    types[2] = MPI_UNSIGNED_SHORT;
    blocklengths[2] = 8;
    MPI_Type_extent(MPI_UNSIGNED_LONG_LONG, &extent);
    offsets[2] = offsets[1] + extent;

    MPI_Type_create_struct(3, blocklengths, offsets, types, dt_partial_trajectory);
    MPI_Type_commit(dt_partial_trajectory);
}