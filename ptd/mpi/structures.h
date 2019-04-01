
enum DEPORTATION_PLACE
{
    INIT = -1,
    AIR = 0,
    WALL = 1,
    TOP = 2,
    BOTTOM = 3,
    DECOMPOSED = 4,
    NONE = 5
};

enum MESSAGE_TYPE
{
    PARTICLE_REQUEST = 10,
    PARTICLE = 11,
    PARTIAL_TRAJECTORY_REQUEST = 12,
    PARTIAL_TRAJECTORY = 13,
    SIMULATION_END = 14,
    NEW_CHUNK = 15
};

typedef struct
{
    double life_time;
    double half_life;
    double z;
    double x;
    double y;
    
    int id;
    enum DEPORTATION_PLACE status;

} Particle;

typedef struct 
{
    // Minimal bounding box of the partial trajectory
    double xMin;
    double xMax;
    double yMin;
    double yMax;
    double zMin;
    double zMax;

    // Translation vector X
    double xLast;
    double yLast;
    double zLast;

    /* 
    	Partial  trajectory  time,  which  represents  
    	the  timeneeded to cross the partial trajectory.
    */
    double lifeTimeStepSum;
    
    // Random number generator seed - drand48_data
    unsigned long long int __a;
    unsigned short int __c;
    unsigned short int __init;
    unsigned short int __old_x1;
    unsigned short int __old_x2;
    unsigned short int __old_x3;
    unsigned short int __x1;
    unsigned short int __x2;
    unsigned short int __x3;
    
} Partial_trajectory;

typedef struct
{
    int nodeCount;
    int processorsPerNode;
    int producerCount;
    int consumerCount;
    int partialTrajectoryLength;
    int chunkSize;
    int particleCount;

} Conf;