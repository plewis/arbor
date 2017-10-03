#pragma once

#include <iostream>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

namespace arbor
{

class Arbor
    {
    public:
                            Arbor();
                            ~Arbor();

        void                processCommandLineOptions(int argc, const char * argv[]);
        void                run();

    private:
        typedef std::pair<unsigned, unsigned>   size_topol_pair_t;
        typedef std::vector<double>             dlb_vect_t;
        typedef std::vector<double>             param_vect_t;
        typedef std::vector<param_vect_t>       sample_vect_t;
        typedef std::pair<unsigned, double>     pair_uint_dbl_t;
        typedef pair_uint_dbl_t &               pair_uint_dbl_ref_t;

        enum ParamType {ptype_edge_length, ptype_gamma_shape, ptype_pinvar, ptype_exchangeability, ptype_frequency};

        // member functions
        void                showUserSpecifiedOptions() const;
        void                calcParamTypesFromColumnHeaders();
        void                readData();

        // data members representing user-configurable options
        std::string                         _data_file_name;            // the name of the sample file to be processed
        unsigned                            _min_sample_size;           // minimum number of samples for tree topology inclusion
        unsigned                            _rnseed;                    // pseudorandom number seed
        double                              _grape_fraction;            // fraction of sample used to define grapes

        // data members that could become user-configurable options in the future
        unsigned                            _skip;                      // number of samples to skip (as burn-in)

        // other data members
        std::vector<std::string>            _column_headers;            // labels in first row of parameter file
        std::map<unsigned, dlb_vect_t>      _log_likelihoods;           // map with key=topology and value=vector of log-likelihood (3rd column of parameter file)
        std::map<unsigned, dlb_vect_t>      _log_priors;                // map with key=topology and value=vector of log-prior (4th column of parameter file)
        std::map<unsigned, sample_vect_t>   _parameters;                // map with key=topology and value=vector of parameter vectors
        unsigned                            _nsubsets;                  // number of data subsets (determined from prefixes in column headers)
        unsigned                            _nparams;                   // length of each vector stored in _parameters
        std::vector<ParamType>              _param_types;               // parameter types in order of appearance in sample vectors
        std::vector<unsigned>               _param_subsets;             // data subset of each parameter in order of appearance in sample vectors
        std::vector<size_topol_pair_t>      _tree_topologies;           // vector storing topology number and frequency pairs

        // static data members
        static std::string  _program_name;              // name of the program
        static std::string  _author;                    // author of program
        static unsigned     _major_version;             // major version number (e.g. "1" in "version 1.2")
        static unsigned     _minor_version;             // minor version number (e.g. "2" in "version 1.2")
        static unsigned     _def_rnseed;                // default value for _rnseed
        static unsigned     _def_minsamplesize;         // default value for _min_sample_size
        static double       _def_grapefraction;         // default value for _grape_fraction

    };

inline Arbor::Arbor()
  : _grape_fraction(_def_grapefraction)
  , _min_sample_size(_def_minsamplesize)
  , _rnseed(_def_rnseed)
  , _skip(1)
    {
    std::cout << boost::str(boost::format("%s %d.%d (written by %s)") % _program_name % _major_version % _minor_version % _author) << std::endl;
    }

inline Arbor::~Arbor()
    {
    }
    
}
