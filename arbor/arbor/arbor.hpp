#pragma once

#include <iostream>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include "grape.hpp"
#include "grapedb.hpp"
#include "lot.hpp"

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
        void                initRandom();
        void                showUserSpecifiedOptions() const;
        void                calcParamTypesFromColumnHeaders();
        void                readData();
        void                showParamTable() const;
        void                calcIndivTopolMargLikes();
        void                showTopoFreq();
        void                showTopoTable();
        double              summationWithFloatingControl(const std::vector<double> & v) const;
        void                margLikeOneTopology(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors);
        void                chooseReferenceSample(sample_vect_t & parameters);
        void                calcMeansAndStdevs(sample_vect_t & parameters);
        void                standardizeSampleAllowingCorrelation(sample_vect_t & parameters);
        void                standardizeSampleAssumingIndependence(sample_vect_t & parameters);
        double              calcDistance(std::vector<double> & a, std::vector<double> & b) const;
        double              calcPutativeRadius(unsigned index, sample_vect_t & parameters, unsigned num_to_keep, std::vector<pair_uint_dbl_t> & neighbors);
        void                createGrapes(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians);
        void                calcDeltaTerms(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians);
        void                calcRatioTerms(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians);
        void                calcOverallMargLike();


        // data members representing user-configurable options
        std::string                             _data_file_name;            // the name of the sample file to be processed
        unsigned                                _min_sample_size;           // minimum number of samples for tree topology inclusion
        unsigned                                _rnseed;                    // pseudorandom number seed
        double                                  _grape_fraction;            // fraction of total sample used as the reference sample
        double                                  _keep_fraction;             // fraction of the reference sample used to define grapes
        double                                  _forced_radius;             // if specified, this will be the radius used for all grapes

        // data members that could become user-configurable options in the future
        unsigned                                _skip;                      // number of samples to skip (as burn-in)
        bool                                    _correlation;               // if true, uses full covariance matrix to standardize parameters

        // other data members
        bool                                     _force_radius;             // true if and only if _forced_radius is specified by user
        Lot::SharedPtr                          _lot;                       // Pseudorandom number generator
        GrapeDatabase                           _db;                        // object that stores results for each tree topology and used for outputting summaries of results
        std::vector<std::string>                _column_headers;            // labels in first row of parameter file
        std::map<unsigned, dlb_vect_t>          _log_likelihoods;           // map with key=topology and value=vector of log-likelihood (3rd column of parameter file)
        std::map<unsigned, dlb_vect_t>          _log_priors;                // map with key=topology and value=vector of log-prior (4th column of parameter file)
        std::map<unsigned, sample_vect_t>       _parameters;                // map with key=topology and value=vector of parameter vectors
        unsigned                                _nsubsets;                  // number of data subsets (determined from prefixes in column headers)
        unsigned                                _nparams;                   // length of each vector stored in _parameters
        std::vector<ParamType>                  _param_types;               // parameter types in order of appearance in sample vectors
        std::vector<unsigned>                   _param_subsets;             // data subset of each parameter in order of appearance in sample vectors
        std::vector<size_topol_pair_t>          _tree_topologies;           // vector storing topology number and frequency pairs
        unsigned                                _N;                         // total number of estimation sample points used over all tree topologies
        unsigned                                _total_sample_size;         // total number of samples used (does not count samples from tree topologies not used)
        double                                  _min_radius;                // minimum radius of all grapes for one particular tree topology
        std::vector<double>                     _log_ratio_terms;           // terms composing numerator of estimator of 1/c
        std::vector<double>                     _delta_terms;               // terms composing denominator of estimator of 1/c

        // - used if assuming independence among parameters
        std::vector<double>                     _means;                     // means of parameters from reference sample
        std::vector<double>                     _stdevs;                    // standard deviations of parameters from reference sample

        // - used if allowing correlations among parameters
        Eigen::VectorXd                         _Means;                     // 1 x p vector of column means
        Eigen::VectorXd                         _primary_eigenvector;       // eigenvector with largest eigenvalue from variance-covariance matrix
        Eigen::MatrixXd                         _Y;                         // n x p matrix of sampled parameter vectors
        Eigen::MatrixXd                         _M;                         // n x p matrix in which each row is vector of column means
        Eigen::MatrixXd                         _S;                         // p x p covariance matrix raised to power 0.5
        Eigen::MatrixXd                         _Sinv;                      // p x p covariance matrix raised to power -0.5
        double                                  _logdetS;                   // log |determinant of _S|

        // data members that are reused for each topology
        unsigned                                _ngrapes;                   // equals _keep_fraction*_reference_indices.size()
        unsigned                                _in_multiple_grapes;        // number of estimation sample points placed in more than one grape
        std::vector<double>                     _placed;                    // vector of absolute differences between log-kernel of placed sample point and the log-kernel of the center of the grape in which it was placed
        std::vector<Grape>                      _grapes;                    // vector of Grape objects
        std::vector<double>                     _log_jacobians;             // log of the jacobian for log-transformation plus standardization
        std::vector<unsigned>                   _reference_indices;         // Indices into _parameters of reference sample points
        std::vector<unsigned>                   _estimation_indices;        // Indices into _parameters of estimation sample points

        // static data members
        static std::string                      _program_name;              // name of the program
        static std::string                      _author;                    // author of program
        static unsigned                         _major_version;             // major version number (e.g. "1" in "version 1.2")
        static unsigned                         _minor_version;             // minor version number (e.g. "2" in "version 1.2")
        static unsigned                         _def_rnseed;                // default value for _rnseed
        static unsigned                         _def_minsamplesize;         // default value for _min_sample_size
        static double                           _def_grapefraction;         // default value for _grape_fraction
        static double                           _def_keepfraction;          // default value for _keep_fraction
        static double                           _def_forcedradius;          // default value for _force_radius

    };

inline Arbor::Arbor()
  : _min_sample_size(_def_minsamplesize)
  , _rnseed(_def_rnseed)
  , _grape_fraction(_def_grapefraction)
  , _keep_fraction(_def_keepfraction)
  , _forced_radius(_def_forcedradius)
  , _skip(1)
  , _correlation(true)
  , _force_radius(false)
    {
    std::cout << boost::str(boost::format("%s %d.%d (written by %s)") % _program_name % _major_version % _minor_version % _author) << std::endl;
    }

inline Arbor::~Arbor()
    {
    }

}
