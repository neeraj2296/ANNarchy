profile_header=\
"""
#pragma once

#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include <papi.h>
#include <math.h>

#ifdef _DEBUG
    #define debug_cout(x) (std::cout << x << std::endl)
#else
    #define debug_cout(x)
#endif

class Measurement{
    std::vector<double> _raw_data;
    double _mean;
    double _std;
    long_long _start;
    long_long _stop;

public:
    Measurement() {
        debug_cout("Create Measurement object");
        _mean = 0.0;
        _std = 0.0;
        _start = 0.0;
        _stop = 0.0;
        _raw_data = std::vector<double>();
    }

    ~Measurement() {
        debug_cout("Destroy Measurement object");
        _raw_data.clear();
    }

    inline void start_wall_time() {
        _start = PAPI_get_real_usec();
    }

    inline void stop_wall_time() {
        _stop = PAPI_get_real_usec();
        _raw_data.push_back(double(_stop-_start));
    }

    void reset() {
        debug_cout("Reset Measurement object");
        _raw_data.clear();
        _mean = 0.0;
        _std = 0.0;
    }

    void evaluate() {
        if (_raw_data.empty())
            return;

        long_long sum = 0.0;
        int num_elem = _raw_data.size();

        // mean over all times
        for( auto it = _raw_data.begin(); it != _raw_data.end(); it++ ) {
            sum += *it;
        }

        _mean = double(sum) / double(num_elem);

        // variance over all times
        double var = 0;
        for( auto it = _raw_data.begin(); it != _raw_data.end(); it++ ) {
            var += (*it-_mean)*(*it-_mean);
        }

        _std = sqrt(var);
    }

    friend std::ostream& operator << (std::ostream& stream, const Measurement& measure);
};

/**
 *  \brief      Profiling class
 *  \details    Creation and initialization of the class are not thread-safe. 
 *              The class is implemented as singleton ensure uniqueness during
 *              runtime.
 */
class Profiling {
    static std::unique_ptr<Profiling> _instance;    ///< reference to this class, created on first call of Profiling::get_instance()
    std::vector<Measurement*> _datasets;    ///< Instances of measurement objects. The index for correct access is retrieved from Profiling::_identifier
    std::map<std::pair<std::string, std::string>, int > _identifier;    ///< maps a (obj, func) descriptor to index for Profiling::_datasets

    /**
     *  \brief  Constructor
     */
    Profiling() {
        debug_cout("Create Profiling instance.");

        // initialize PAPI
        if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
            exit(1);
    }

public:
    /**
     *  \brief  Destructor
     */
    ~Profiling() {
        debug_cout("Destroy Profiling instance.");

        for(auto it = _datasets.begin(); it != _datasets.end(); it++ )
            delete *it;
    }

    /**
     *  \brief      Retrieve Profiling class instance
     *  \details    First call initialize the class
     *  \return     Reference to Profiling class for further interaction
     */
    static Profiling* get_instance() {
        if ( _instance.get() == nullptr )
            _instance = std::unique_ptr<Profiling>(new Profiling());

        return _instance.get();
    }

    /**
     *  \brief      Add a function to measurement dataset
     *  \details    A measurment is uniquely described by an object (either population or projection name)
     *              and a function name. The provided items will be used as key for the internal data map.
     *  \param[IN]  obj object name as string
     *  \param[IN]  func function name as string
     *  \return     Instance of measurement class. If the function is called multiple times no additional 
     *              object will be created.
     */
    Measurement* register_function(std::string obj, std::string func) {
        auto pair = std::pair<std::string, std::string>(obj, func);

        if ( _identifier.count(pair) == 0) { // not in list
            _identifier.insert(std::pair< std::pair<std::string, std::string>, int >(pair, _datasets.size()));
            _datasets.push_back(new Measurement());

            debug_cout( "(" + pair.first + ", " + pair.second + ") added to dataset." );
        }

        return get_measurement(obj, func);
    }

    /**
     *  \brief      Retrieve registered function from measurement dataset
     *  \details    A measurment is uniquely described by an object (either population or projection name)
     *              and a function name. The provided items will be used as key for the internal data map.
     *  \param[IN]  obj object name as string
     *  \param[IN]  func function name as string
     *  \return     Instance of measurement class. If the measurment were not registered before, no additional 
     *              object will be created. The function returns a 'nullptr' in this case.
     */
    Measurement* get_measurement(std::string obj, std::string func) {

        auto idx = _identifier.find(std::pair<std::string, std::string>(obj, func));
        if ( idx == _identifier.end() )
            return nullptr;

        return _datasets[idx->second];
    }

    /**
     *  \brief      Clear performance data of registered functions.
     *  \details    A measurment is uniquely described by an object (either population or projection name)
     *              and a function name. We iterate over all datasets and clear them seperatly.
     */
    void reset() {
        for( auto it = _datasets.begin(); it != _datasets.end(); it++ )
            (*it)->reset();
    }

    /**
     *  \brief      Evaluate registered functions.
     *  \details    A measurment is uniquely described by an object (either population or projection name)
     *              and a function name. We iterate over all datasets and evaluate them seperatly.
     */
    void evaluate() {
        for( auto it = _datasets.begin(); it != _datasets.end(); it++ )
            (*it)->evaluate();
    }

    friend std::ostream& operator << (std::ostream& stream, const Profiling& profiling);
};

// out stream operators
inline std::ostream& operator << (std::ostream& stream, const Measurement& measure) {
    if ( measure._raw_data.empty() )
        return stream;

    stream << "mean: " << measure._mean << ", "
           << "std: " << measure._std << " ( over " << measure._raw_data.size() << " measurements )";

    return stream;
}

inline std::ostream& operator << (std::ostream& stream, const Profiling& profiling) {
    if ( profiling._datasets.empty() )
        return stream;

    for ( auto it = profiling._identifier.begin(); it != profiling._identifier.end(); it++ ) {
        stream << "  (" << it->first.first << ", " << it->first.second << ") - " << *(profiling._datasets[ it->second ]) << std::endl;
    }

    return stream;
}
"""

profile_template = {
    'include': """//Profiling
#include "Profiling.h"
std::unique_ptr<Profiling> Profiling::_instance(nullptr);
""",
    'init': """
    //initialize profiler, create singleton instance
    auto profiler = Profiling::get_instance();
    profiler->register_function("net", "step");
    """,
    'step_pre': """// before
    auto measure = Profiling::get_instance()->get_measurement("net", "step");
    measure->start_wall_time();
    """,
    'step_post': """// after
    measure->stop_wall_time();
    """,
    'run_pre': """
    Profiling::get_instance()->reset();
    """,
    'run_post': """
    Profiling::get_instance()->evaluate();
    std::cout << "run " << nbSteps << " steps: " << std::endl;
    std::cout << *Profiling::get_instance() << std::endl;
    """,
    #
    # Operations
    'compute_psp': {
        'before' : "measure_psp->start_wall_time();",
        'after' : "measure_psp->stop_wall_time();"
    },
    'update_synapse': {
        'before' : "measure_step->start_wall_time();",
        'after' : "measure_step->stop_wall_time();"
    },
    'update_neuron': {
        'before' : "measure_step->start_wall_time();",
        'after' : "measure_step->stop_wall_time();"
    },
    'spike_prop': {
        'before' : "measure_prop->start_wall_time();",
        'after' : "measure_prop->stop_wall_time();"
    }
}
