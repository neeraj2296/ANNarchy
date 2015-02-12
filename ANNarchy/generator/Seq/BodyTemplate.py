body_template = '''
#include "ANNarchy.h"

/*
 * Internal data
 *
*/
double dt;
long int t;
std::mt19937  rng;

// Populations
%(pop_ptr)s

// Projections
%(proj_ptr)s

// Global operations
%(glops_def)s

// Simulate the network for the given number of steps
void run(int nbSteps) {

    for(int i=0; i<nbSteps; i++)
    {
        step();
    }

}

int run_until(int steps, std::vector<int> populations, bool or_and)
{

%(run_until)s

}

// Initialize the internal data and random numbers generators
void initialize(double _dt, long int seed) {

    // Internal variables
    dt = _dt;
    t = (long int)(0);

    // Random number generators
    if(seed==-1){
        rng = std::mt19937(time(NULL));
    }
    else{
        rng = std::mt19937(seed);
    }
    
%(random_dist_init)s
%(delay_init)s
%(spike_init)s
%(projection_init)s
%(globalops_init)s
}

// Step method. Generated by ANNarchy.
void step()
{

    double sum;
    int rk_pre, rk_post, i, j, rk_j;
    
    ////////////////////////////////
    // Presynaptic events
    ////////////////////////////////
%(reset_sums)s
%(compute_sums)s

    ////////////////////////////////
    // Update random distributions
    ////////////////////////////////
%(random_dist_update)s

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


    ////////////////////////////////
    // Update synaptic variables
    ////////////////////////////////
%(update_synapse)s    


    ////////////////////////////////
    // Postsynaptic events
    ////////////////////////////////
%(post_event)s

    ////////////////////////////////
    // Structural plasticity
    ////////////////////////////////
%(structural_plasticity)s

    ////////////////////////////////
    // Recording
    ////////////////////////////////
%(record)s

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
'''