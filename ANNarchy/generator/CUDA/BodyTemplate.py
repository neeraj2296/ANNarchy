body_template = '''
#ifdef __CUDA_ARCH__
/***********************************************************************************/
/*                                                                                 */
/*                                                                                 */
/*          DEVICE - code                                                          */
/*                                                                                 */
/*                                                                                 */
/***********************************************************************************/
#include <curand_kernel.h>

/****************************************
 * init random states                   *
 ****************************************/
__global__ void rng_setup_kernel( int N, curandState* states, unsigned long seed )
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if( tid < N )
    {
        curand_init( seed, tid, 0, &states[ tid ] );
    }
}

/****************************************
 * inline functions                     *
 ****************************************/
__device__ __forceinline__ double positive( double x ) { return (x>0) ? x : 0; }
__device__ __forceinline__ double negative( double x ) { return x<0.0? x : 0.0; }
__device__ __forceinline__ double clip(double x, double a, double b) { return x<a? a : (x>b? b :x); }

/****************************************
 * custom functions                     *
 ****************************************/
%(custom_func)s

/****************************************
 * updating neural variables            *
 ****************************************/
%(pop_kernel)s
 
/****************************************
 * weighted sum kernels                 *
 ****************************************/
%(psp_kernel)s

/****************************************
 * update synapses kernel               *
 ****************************************/
%(syn_kernel)s

/****************************************
 * global operations kernel             *
 ****************************************/
%(glob_ops_kernel)s

#else
#include "ANNarchy.h"
#include <math.h>

// cuda specific header
#include <cuda_runtime_api.h>
#include <curand.h>
#include <float.h>

/***********************************************************************************/
/*                                                                                 */
/*                                                                                 */
/*          HOST - code                                                            */
/*                                                                                 */
/*                                                                                 */
/***********************************************************************************/
// kernel config
%(kernel_config)s

// RNG
__global__ void rng_setup_kernel( int N, curandState* states, unsigned long seed );

void init_curand_states( int N, curandState* states, unsigned long seed ) {
    int numThreads = 64;
    int numBlocks = ceil (double(N) / double(numThreads));

    rng_setup_kernel<<< numBlocks, numThreads >>>( N, states, seed);
}

/*
 * Internal data
 *
 */
double dt;
long int t;

// Recorders 
std::vector<Monitor*> recorders;
void addRecorder(Monitor* recorder){
    recorders.push_back(recorder);
}
void removeRecorder(Monitor* recorder){
    for(int i=0; i<recorders.size(); i++){
        if(recorders[i] == recorder){
            recorders.erase(recorders.begin()+i);
            break;
        }
    }
}

// Populations
%(pop_ptr)s

// Projections
%(proj_ptr)s

template<typename T>
std::vector<int> flattenIdx(std::vector<std::vector<T> > in)
{
    std::vector<T> flatIdx = std::vector<T>();
    typename std::vector<std::vector<T> >::iterator it;

    for ( it = in.begin(); it != in.end(); it++)
    {
        flatIdx.push_back(it->size());
    }

    return flatIdx;
}

template<typename T>
std::vector<int> flattenOff(std::vector<std::vector<T> > in)
{
    std::vector<T> flatOff = std::vector<T>();
    typename std::vector<std::vector<T> >::iterator it;

    int t = 0;
    for ( it = in.begin(); it != in.end(); it++)
    {
        flatOff.push_back(t);
        t+= it->size();
    }

    return flatOff;
}

template<typename T> 
std::vector<T> flattenArray(std::vector<std::vector<T> > in) 
{
    std::vector<T> flatVec = std::vector<T>();
    typename std::vector<std::vector<T> >::iterator it;

    for ( it = in.begin(); it != in.end(); it++)
    {
        flatVec.insert(flatVec.end(), it->begin(), it->end());
    }

    return flatVec;
}

template<typename T> 
std::vector<std::vector<T> > deFlattenArray(std::vector<T> in, std::vector<int> idx) 
{
    std::vector<std::vector<T> > deFlatVec = std::vector<std::vector<T> >();
    std::vector<int>::iterator it;

    int t=0;
    for ( it = idx.begin(); it != idx.end(); it++)
    {
        std::vector<T> tmp = std::vector<T>(in.begin()+t, in.begin()+t+*it);
        t += *it;

        deFlatVec.push_back(tmp);
    }

    return deFlatVec;
}

%(stream_setup)s

void progress(int i, int nbSteps) {
    double tInMs = nbSteps * dt;
    if ( tInMs > 1000.0 )
        std::cout << "\\rSimulate " << (int)(tInMs/1000.0) << " s: " << (int)( (double)(i+1)/double(nbSteps) * 100.0 )<< " finished.";
    else
        std::cout << "\\rSimulate " << tInMs << " ms: " << (int)( (double)(i+1)/double(nbSteps) * 100.0 )<< " finished.";
    std::flush(std::cout);
}

/**
 *  Implementation remark (27.02.2015, HD) to: run(int), step() and single_step()
 *
 *  we have two functions in ANNarchy to run simulation code: run(int) and step(). The latter one to
 *  run several steps at once, the other one just a single step. On CUDA I face the problem, that I
 *  propably need to update variables before step() and definitly changed variables after step().
 *  run(int) calls step() normally, if I add transfer codes in step(), run(N) would end up in N
 *  back transfers from GPUs, whereas we only need the last one.
 *  As solution I renamed step() to single_step(), the interface behaves as OMP side and I only
 *  transfer at begin and at end, as planned.
 */
void single_step(); // function prototype

// Simulate the network for the given number of steps
void run(int nbSteps) {
%(host_device_transfer)s

    stream_assign();

    // simulation loop
    for(int i=0; i<nbSteps; i++) {
        single_step();
        //progress(i, nbSteps);
    }
    
    //std::cout << std::endl;
%(device_host_transfer)s
}

void step() {
%(host_device_transfer)s
    single_step();
%(device_host_transfer)s
}

// Initialize the internal data and random numbers generators
void initialize(double _dt, long seed) {

    dt = _dt;
    t = (long int)(0);

%(device_init)s
%(random_dist_init)s
%(delay_init)s
%(spike_init)s
%(globalops_init)s
%(projection_init)s

    // create streams
    stream_setup();
}

%(kernel_def)s

// Step method. Generated by ANNarchy. (analog to step() in OMP)
void single_step()
{

    ////////////////////////////////
    // Presynaptic events
    ////////////////////////////////
%(compute_sums)s
    cudaDeviceSynchronize();

    ////////////////////////////////
    // Reset spikes
    ////////////////////////////////

    ////////////////////////////////
    // Update neural variables
    ////////////////////////////////
%(update_neuron)s    

    ////////////////////////////////
    // Delay outputs
    ////////////////////////////////
%(delay_code)s

    ////////////////////////////////
    // Global operations (min/max/mean)
    ////////////////////////////////
%(update_globalops)s

    cudaDeviceSynchronize();

    ////////////////////////////////
    // Update synaptic variables
    ////////////////////////////////
%(update_synapse)s    

    ////////////////////////////////
    // Postsynaptic events
    ////////////////////////////////

    ////////////////////////////////
    // Recording
    ////////////////////////////////
    for(int i=0; i < recorders.size(); i++){
        recorders[i]->record();
    }
    

    ////////////////////////////////
    // Increase internal time
    ////////////////////////////////
    t++;

}


/*
 * Access to time and dt
 *
*/
long int getTime() {return t;}
void setTime(long int t_) { t=t_;}
double getDt() { return dt;}
void setDt(double dt_) { dt=dt_;}
#endif
'''
