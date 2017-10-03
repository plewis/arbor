#include <regex>
#include <fstream>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/gamma.hpp>   // lgamma
#include "utilfunc.hpp"
#include "arbor.hpp"
#include "xarbor.hpp"

using namespace arbor;

void Arbor::processCommandLineOptions(int argc, const char * argv[])
    {
    // Get user-supplied configuration settings from command line and/or config file
    boost::program_options::variables_map       vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h",                                                                                                                  "produce help message")
        ("version,v",                                                                                                               "show program version")
        ("datafile,d",  boost::program_options::value<std::string>(&_data_file_name),                                           "name of data file in NEXUS format")
        ("grapefraction,f", boost::program_options::value< double >(&_grape_fraction)->default_value(_def_grapefraction),           "fraction of sample to use in forming grapes")
        ("minsamplesize,m", boost::program_options::value< unsigned >(&_min_sample_size)->default_value(_def_minsamplesize),             "minimum number of sampled points for tree topology to be considered")
        ("seed,s",          boost::program_options::value< unsigned >(&_rnseed)->default_value(_def_rnseed),             "seed for pseudorandom number generation")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    // If user specified --help on command line, output usage summary and quit
    if (vm.count("help") > 0)
        {
        std::cout << desc << "\n";
        std::exit(1);
        }

    // If user specified --version on command line, output version and quit
    if (vm.count("version") > 0)
        {
        std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
        std::exit(1);
        }

    // Ensure that user specified --datafile on command line
    if (vm.count("datafile") == 0)
        {
        std::cout << "Must specify --datafile on command line\n";
        std::cout << desc << std::endl;
        std::exit(1);
        }

    // If user used --seed on command line, perform sanity check to ensure that value specified was greater than 0
    if (vm.count("seed") > 0)
        {
        if (_rnseed <= 0.0)
            {
            std::cout << boost::str(boost::format("Error: specified pseudorandom number seed (%.3f) should be greater than 0.") % _rnseed) << std::endl;
            std::exit(1);
            }
        }
    // If user used --grapefraction on command line, perform sanity check to ensure that value specified was greater than or equal to 0 and less than or equal to 1
    if (vm.count("grapefraction") > 0)
        {
        if (_grape_fraction <= 0.0 || _grape_fraction >= 1.0)
            {
            std::cout << boost::str(boost::format("Error: specified grapefraction (%.3f) should be between 0 and 1.") % _grape_fraction) << std::endl;
            std::exit(1);
            }
        }

    // If user used --minsamplesize on command line, perform sanity check to ensure that value specified was greater than 0
    if (vm.count("minsamplesize") > 0)
        {
        if (_min_sample_size <= 0)
            {
            std::cout << boost::str(boost::format("Error: specified minsamplesize (%d) should be greater than 0.") % _min_sample_size) << std::endl;
            std::exit(1);
            }
        }

    }

void Arbor::showUserSpecifiedOptions() const
    {
    std::cout << "\nUser-configurable settings:" << std::endl;
    std::cout << boost::str(boost::format("  seed:               %d")   % _rnseed) << std::endl;
    std::cout << boost::str(boost::format("  datafile:           %s")   % _data_file_name) << std::endl;
    std::cout << boost::str(boost::format("  minsamplesize:      %d")   % _min_sample_size) << std::endl;
    std::cout << boost::str(boost::format("  grapefraction:      %.3f") % _grape_fraction) << std::endl;
    }

void Arbor::calcParamTypesFromColumnHeaders()
    {
    // Fill in _param_types and _param_subsets vectors using strings in _column_headers vector

    std::string subsetstr;
    std::smatch what;
    std::regex rex_shape("log\\((\\d+)_gamma_shape\\)");
    std::regex rex_pinvar("logit\\((\\d+)_pinvar\\)");
    std::regex rex_exchangeability("log\\((\\d+)_r([ACGT][ACGT])/(\\d+)_rAC\\)");
    std::regex rex_frequency("log\\((\\d+)_freq([ACGT])/(\\d+)_freqA\\)");

    _nsubsets = 0;
    _param_types.clear();
    _param_subsets.clear();
    for (unsigned i = 4; i < _column_headers.size(); ++i)
        {
        std::string & h = _column_headers[i];
        if (std::regex_match(h, what, rex_shape))
            {
            _param_types.push_back(ptype_gamma_shape);
            subsetstr = std::string(what[1].first, what[1].second);
            unsigned subset = boost::lexical_cast<unsigned>(subsetstr);
            _param_subsets.push_back(subset);
            if (subset > _nsubsets)
                _nsubsets = subset;
            }
        else if (std::regex_match(h, what, rex_pinvar))
            {
            _param_types.push_back(ptype_pinvar);
            subsetstr = std::string(what[1].first, what[1].second);
            unsigned subset = boost::lexical_cast<unsigned>(subsetstr);
            _param_subsets.push_back(subset);
            if (subset > _nsubsets)
                _nsubsets = subset;
            }
        else if (std::regex_match(h, what, rex_exchangeability))
            {
            _param_types.push_back(ptype_exchangeability);
            subsetstr = std::string(what[1].first, what[1].second);
            unsigned subset = boost::lexical_cast<unsigned>(subsetstr);
            _param_subsets.push_back(subset);
            if (subset > _nsubsets)
                _nsubsets = subset;
            }
        else if (std::regex_match(h, what, rex_frequency))
            {
            _param_types.push_back(ptype_frequency);
            subsetstr = std::string(what[1].first, what[1].second);
            unsigned subset = boost::lexical_cast<unsigned>(subsetstr);
            _param_subsets.push_back(subset);
            if (subset > _nsubsets)
                _nsubsets = subset;
            }
        else
            {
            _param_types.push_back(ptype_edge_length);
            _param_subsets.push_back(0);
            }
        }

    _nparams = (unsigned)_param_types.size();
    }

void Arbor::readData()
    {
    // Create a temporary map for storing frequency of each topology
    typedef std::map< unsigned, unsigned > topology_freq_t;
    topology_freq_t topology_freq;

    // Clear containers
    _column_headers.clear();
    _log_likelihoods.clear();
    _log_priors.clear();
    _parameters.clear();

    // Open the data file
    std::ifstream inf(_data_file_name.c_str());
    if (!inf)
        throw XArbor(boost::str(boost::format("Could not open file named \"%s\"") % _data_file_name));
        
    // Read contents of file into string
    std::string contents((std::istreambuf_iterator<char>(inf)), std::istreambuf_iterator<char>());

    // Regular expression that splits contents into separate lines
    const char * rex = "[\\r\\n]+";

    // Set up regular expression reader to iterate through rows
    std::regex row_expr(rex);
    std::sregex_token_iterator m(contents.begin(), contents.end(), row_expr, -1);
    std::sregex_token_iterator m2;

    for (unsigned row = 0; m != m2; ++m, ++row)
        {
        std::string rowstr = m->str();
        boost::algorithm::trim(rowstr);
        std::vector< std::string > s;
        boost::split(s, rowstr, boost::is_any_of("\t"), boost::token_compress_on);

        unsigned nelements = (unsigned)s.size();

        if (row == 0) {
            // first row should be column headers
            _column_headers.resize(nelements);
            std::copy(s.begin(), s.end(), _column_headers.begin());

            calcParamTypesFromColumnHeaders();
        }
        else if (row > _skip) {
            // s[0] is the iteration: we ignore this

            // s[1] is the tree topology: this will be used as a key for _log_likelihoods, _log_priors, and _parameters
            unsigned tree_topology = boost::lexical_cast<unsigned>(s[1]);

            // Update the frequency of this tree topology
            topology_freq_t::iterator lower_bound = topology_freq.lower_bound(tree_topology);
            if (lower_bound != topology_freq.end() && !(topology_freq.key_comp()(tree_topology, lower_bound->first)))
                {
                // this tree topology has already been seen, so increment the count
                lower_bound->second += 1;
                }
            else
                {
                // this tree topology has not yet been seen, so create a new entry
                topology_freq.insert(lower_bound, topology_freq_t::value_type(tree_topology, 1));
                }

            // s[2] is the log-likelihood: store this in _log_likelihoods map
            double logL = boost::lexical_cast<double>(s[2]);
            _log_likelihoods[tree_topology].push_back(logL);

            // s[3] is the log-prior: store this in _log_priors map
            double logP = boost::lexical_cast<double>(s[3]);
            _log_priors[tree_topology].push_back(logP);

            // s[4] is the first log-transformed parameter: store all parameters in _parameters map
            std::vector<double> tmp;
            for (unsigned i = 4; i < nelements; ++i)
                {
                double v = boost::lexical_cast<double>(s[i]);
                tmp.push_back(v);
                }
            _parameters[tree_topology].push_back(tmp);

        }   // else row > 0
    } // loop over rows

    // Copy entries of topology_freq map to the vector _tree_topologies, each element of which is a (frequency,topology) pair
    _tree_topologies.clear();
    for (auto p : topology_freq)
        {
        _tree_topologies.push_back(std::make_pair(p.second,p.first));
        }

    // Sort the entries of _tree_topologies so that topologies with highest frequency come first
    std::sort(_tree_topologies.begin(), _tree_topologies.end(), std::greater<size_topol_pair_t>());

    // _nparams should equal the length of the first parameter vector for the most frequently sampled tree topology
    assert(_nparams == _parameters[_tree_topologies[0].second][0].size());

    //unsigned i = 1;
    //std::cout << "\nTopologies found:\n" << std::endl;
    //for (auto p : _tree_topologies)
    //    {
    //    std::cout << boost::str(boost::format("%12d %12d %12d") % i % p.first % p.second) << std::endl;
    //    ++i;
    //    }
    }

void Arbor::initRandom()
    {
    _lot = Lot::SharedPtr(new Lot());
    _lot->setSeed(_rnseed);
    }

void Arbor::showParamTable() const
    {
    std::cout << std::endl;
    std::cout << boost::str(boost::format("  %12s %12s   %s") % "parameter" % "data subset" % "parameter type") << std::endl;
    for (unsigned k = 0; k < _param_types.size(); ++k)
        {
        if (_param_types[k] == 0)
            std::cout << boost::str(boost::format("  %12d %12s   %s") % (k+1) % "NA" % "edge length") << std::endl;
        else if (_param_types[k] == 1)
            std::cout << boost::str(boost::format("  %12d %12d   %s") % (k+1) % _param_subsets[k] % "gamma shape") << std::endl;
        else if (_param_types[k] == 2)
            std::cout << boost::str(boost::format("  %12d %12d   %s") % (k+1) % _param_subsets[k] % "proportion of invariable sites") << std::endl;
        else if (_param_types[k] == 3)
            std::cout << boost::str(boost::format("  %12d %12d   %s") % (k+1) % _param_subsets[k] % "GTR exchangeability") << std::endl;
        else if (_param_types[k] == 4)
            std::cout << boost::str(boost::format("  %12d %12d   %s") % (k+1) % _param_subsets[k] % "state frequency") << std::endl;
        else
            std::cout << boost::str(boost::format("  %12d %12s   %s") % (k+1) % "NA" % "unknown") << std::endl;
        }
    }

void Arbor::chooseGrapeSeeds(sample_vect_t & parameters)
    {
    // Builds two vectors of indices into _parameters vector:
    // _reference_indices holds indices of points used as grape seeds (a random proportion _grape_fraction of all points in _parameters);
    // _estimation_indices holds indices of all remaining points.

    // Clear vectors that will be used to store indices
    _reference_indices.clear();
    _estimation_indices.clear();

    for (unsigned i = 0; i < parameters.size(); ++i)
        {
        if (_lot->uniform() < _grape_fraction)
            _reference_indices.push_back(i);
        else
            _estimation_indices.push_back(i);
        }

    // Add to tally of total number of estimation points used over all topologies
    _N += (unsigned)_estimation_indices.size();
    }

void Arbor::calcMeansAndStdevs(sample_vect_t & parameters)
    {
    // Computes means and standard deviations of each parameter for the points whose indices are in _reference_indices.
    assert(_nparams > 0);

    // Clear and zero vectors that will be used to store means and standard deviations
    _means.resize(_nparams);
    _means.assign(_nparams, 0.0);
    _stdevs.resize(_nparams);
    _stdevs.assign(_nparams, 0.0);

    // Loop through all vectors in reference sample
    unsigned reference_sample_size = (unsigned)_reference_indices.size();
    for (auto i : _reference_indices)
        {
        std::vector<double> & v = parameters[i];

        // Update sum stored in _means and sum-of-squares stored in _stdevs
        for (unsigned j = 0; j < _nparams; ++j)
            {
            _means[j]  += v[j];
            _stdevs[j] += pow(v[j],2);
            }
        }

    // Complete _means and _stdevs calculation
    for (unsigned i = 0; i < _nparams; ++i)
        {
        _means[i] /= reference_sample_size;
        double tmp = (_stdevs[i] - pow(_means[i],2)*reference_sample_size)/(reference_sample_size - 1);
        _stdevs[i] = sqrt(tmp);
        }
    }

void Arbor::standardizeSampleAllowingCorrelation(sample_vect_t & parameters)
    {
    // Computes variance-covariance matrix for the points whose indices are in _reference_indices and standardizes sample.
    assert(_nparams > 0);

    unsigned n = (unsigned)parameters.size();
    unsigned p = _nparams;

    // calculate _means
    _Means.resize(p);
    _Means.setZero();
    for (auto & sample : parameters)
        {
        for (unsigned j = 0; j < p; ++j)
            _Means(j) += sample[j];
        }
    for (unsigned j = 0; j < p; ++j)
        _Means(j) /= n;

    // Store observations in matrix _Y and means in _M
    _Y.resize(n,p);
    _M.resize(n,p);
    unsigned i = 0;
    for (auto & sample : parameters)
        {
        for (unsigned j = 0; j < p; ++j)
            {
            _Y(i,j) = sample[j];
            _M(i,j) = _Means(j);
            }
        ++i;
        }

    // Let M hold the centered observations
    // TODO: should avoid building _M as it repeats the mean vector n times
    Eigen::MatrixXd M = _Y - _M;

    // Compute variance-covariance matrix V
    Eigen::MatrixXd V = M.transpose()*M/(n-1);

    // Compute square root matrix _S = V^{0.5} and inverse square root matrix _Sinv = V^{-0.5}
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(V);
    if (eigensolver.info() != Eigen::Success)
        {
        if (eigensolver.info() == Eigen::NumericalIssue)
            throw XArbor("could not compute eigendecomposition of variance-covariance matrix due to a numerical issue");
        else if (eigensolver.info() == Eigen::NoConvergence)
            throw XArbor("could not compute eigendecomposition of variance-covariance matrix due to no convergence");
        else if (eigensolver.info() == Eigen::InvalidInput)
            throw XArbor("could not compute eigendecomposition of variance-covariance matrix due to invalid input");
        else
            throw XArbor("could not compute eigendecomposition of variance-covariance matrix due to some unknown cause");
        }

    // Compute square root and inverse square root of eigenvalues
    Eigen::VectorXd eval = eigensolver.eigenvalues();
    Eigen::VectorXd eval_sqrt = eval;
    Eigen::VectorXd eval_sqrt_inv = eval;
    unsigned index_of_largest_eval = 0;
    double largest_eval = 0.0;
    for (i = 0; i < p; ++i)
        {
        if (eval[i] > largest_eval) {
            largest_eval = eval[i];
            index_of_largest_eval = i;
        }
        eval_sqrt[i]     = sqrt(eval[i]);
        eval_sqrt_inv[i] = pow(eval[i], -0.5);
        }

    // Obtain eigenvectors
    Eigen::MatrixXd evec = eigensolver.eigenvectors();

    // Save eigenvector with largest eigenvalue
    _primary_eigenvector = evec.col(index_of_largest_eval);

    // Compute _S (square root of variance-covariance matrix)
    _S = evec * eval_sqrt.asDiagonal() * evec.inverse();

    // Compute _Sinv (inverse of square root of variance-covariance matrix)
    _Sinv = _S.inverse();

    // Let matrix X hold the standardized samples
    Eigen::MatrixXd X = M * _Sinv;

    // Let _logdetS hold the log of |determinant of _S|
    // |det(S)| is the Jacobian for the standardization transformation
    _logdetS = log(fabs(_S.determinant()));

    // These vectors will be used to compute the denominator sums for jacobians of dirichlet-like parameters
    std::vector<double> freqdenom(_nsubsets, 1.0);
    std::vector<double> exchgdenom(_nsubsets, 1.0);
    std::vector<double> pinvdenom(_nsubsets, 1.0);

    // Copy values in X back to parameters and compute log jacobian terms
    i = 0;
    _log_jacobians.clear();
    for (auto & sample : parameters)
        {
        freqdenom.assign(_nsubsets, 1.0);
        exchgdenom.assign(_nsubsets, 1.0);
        pinvdenom.assign(_nsubsets, 1.0);

        double logj = 0.0;
        for (unsigned j = 0; j < p; ++j)
            {
            double value = sample[j];

            // log-jacobian
            logj += value;   // term in log-jacobian

            if (_param_types[j] == ptype_pinvar)
                pinvdenom[_param_subsets[j] - 1] += exp(value);
            else if (_param_types[j] == ptype_exchangeability)
                exchgdenom[_param_subsets[j] - 1] += exp(value);
            else if (_param_types[j] == ptype_frequency)
                freqdenom[_param_subsets[j] - 1] += exp(value);

            sample[j] = X(i,j);
            }

        logj += _logdetS;
        for (unsigned k = 0; k < _nsubsets; ++k)
            {
            logj -= 2.0*log(pinvdenom[k]);
            logj -= 4.0*log(freqdenom[k]);
            logj -= 6.0*log(exchgdenom[k]);
            }

        _log_jacobians.push_back(logj);

        ++i;
        }
    }

void Arbor::standardizeSampleAssumingIndependence(sample_vect_t & parameters)
    {
    // Standardizes vectors in parameters by subtracting mean and dividing by standard deviation for each parameter separately.

    // Sanity checks: should have stored sample in _parameters and calculated _means and _stdevs
    assert(_nparams > 0);
    assert(_means.size() == _nparams);
    assert(_stdevs.size() == _nparams);

    double sum_log_sd = std::accumulate(_stdevs.begin(), _stdevs.end(), 0.0, pluslog());
    _log_jacobians.clear();

    // These vectors will be used to compute the denominator sums for jacobians of dirichlet-like parameters
    std::vector<double> freqdenom(_nsubsets, 1.0);
    std::vector<double> exchgdenom(_nsubsets, 1.0);
    std::vector<double> pinvdenom(_nsubsets, 1.0);

    for (auto & v : parameters)
        {
        freqdenom.assign(_nsubsets, 1.0);
        exchgdenom.assign(_nsubsets, 1.0);
        pinvdenom.assign(_nsubsets, 1.0);

        double logj = 0.0;
        for (unsigned i = 0; i < _nparams; ++i)
            {
            double value = v[i];

            // log-jacobian
            logj += value;   // term in log-jacobian

            if (_param_types[i] == ptype_pinvar)
                pinvdenom[_param_subsets[i] - 1] += exp(value);
            else if (_param_types[i] == ptype_exchangeability)
                exchgdenom[_param_subsets[i] - 1] += exp(value);
            else if (_param_types[i] == ptype_frequency)
                freqdenom[_param_subsets[i] - 1] += exp(value);

            // standardize ith parameter vector element
            v[i] = (value - _means[i])/_stdevs[i];
            }

        logj += sum_log_sd;
        for (unsigned k = 0; k < _nsubsets; ++k)
            {
            logj -= 2.0*log(pinvdenom[k]);
            logj -= 4.0*log(freqdenom[k]);
            logj -= 6.0*log(exchgdenom[k]);
            }

        _log_jacobians.push_back(logj);

        }
    }

double Arbor::calcDistance(std::vector<double> & a, std::vector<double> & b) const
    {
    // Returns Euclidean distance betweeen points a and b.
    assert(_nparams == a.size());
    assert(_nparams == b.size());
    double d = 0.0;
    for (unsigned j = 0; j < _nparams; ++j)
        {
        d += pow(a[j] - b[j], 2.0);
        }
    d = sqrt(d);
    return d;
    }

double Arbor::calcPutativeRadius(unsigned index, sample_vect_t & parameters, unsigned num_to_keep, std::vector<pair_uint_dbl_t> & neighbors)
    {
    // Returns smallest radius around _parameters[index] such that num_to_keep reference points are included.
    neighbors.clear();
    for (unsigned i = 0; i < _reference_indices.size(); ++i)
        {
        unsigned k = _reference_indices[i];
        double d = 0.0;
        for (unsigned j = 0; j < _nparams; ++j)
            {
            d += pow(parameters[k][j] - parameters[index][j], 2.0);
            }
        d = sqrt(d);
        neighbors.push_back(std::make_pair(i, d));
        }
    std::sort(neighbors.begin(), neighbors.end(), UintDblPairSmallestFirst);
    neighbors.erase(neighbors.begin() + num_to_keep + 1, neighbors.end());

    return neighbors[num_to_keep].second;
    }

void Arbor::createGrapes(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians)
    {
    // Creates grape objects from points stored in _reference_indices.
    assert(_reference_indices.size() > 0);

    // Sort _reference_indices so that reference points having highest log-kernel are first
    std::vector< pair_uint_dbl_t > tmp;
    for (auto j : _reference_indices)
        {
        double log_kernel = log_likelihoods[j] + log_priors[j] + log_jacobians[j];
        tmp.push_back(std::make_pair(j, log_kernel));
        }
    std::sort(tmp.begin(), tmp.end(), UintDblPairLargestFirst);

    // Recreate _reference_indices using tmp
    _reference_indices.clear();
    for (auto x : tmp)
        {
        _reference_indices.push_back(x.first);
        }

    // Initialize min_max_radius, starting min off at the largest possible radius and starting max at the smallest possible radius
    std::pair<double,double> min_max_radius = std::pair<double,double>(std::numeric_limits<double>().max(),0.0);

    // Each grape is centered at one reference point, and the radius of the grape is set just large enough
    // to enclose the num_to_keep closest reference points
    unsigned num_to_keep = 1;

    _grapes.clear();
    std::vector<pair_uint_dbl_t> neighbors;
    for (unsigned i = 0; i < _reference_indices.size(); ++i)
        {
        unsigned k = _reference_indices[i];
        double r = calcPutativeRadius(k, parameters, num_to_keep, neighbors);
        std::vector<double> & param_vect = parameters[k];
        _grapes.push_back(Grape(k, r, param_vect, neighbors));

        if (r < min_max_radius.first)
            min_max_radius.first = r;
        if (r > min_max_radius.second)
            min_max_radius.second = r;
        }

    // Debug sanity check: make sure all grapes contain at least 2 reference sample points (their center plus 1 other)
    for (auto & g : _grapes)    //temporary!
        {
        unsigned ninside = 0;
        for (unsigned i = 0; i < _reference_indices.size(); ++i)
            {
            unsigned k = _reference_indices[i];
            double d = calcDistance(g._v, parameters[k]);
            if (d <= g._radius + 1.e-8)
                ninside++;
            }
        if (ninside != 2)
            throw XArbor(boost::str(boost::format("at least one grape does not enclose 1 point besides itself (ninside = %d)") % ninside));
        }

    // Trying Yu-Bo's method: all grapes have same radius equal to the smallest radius
    // calculated for any grape
    for (auto & g : _grapes)
        {
        g._radius = min_max_radius.first; //temporary!
        }
    }

void Arbor::calcDeltaTerms(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians)
    {
    // Calculates hyperball and hypercylinder volumes for each grape and appends terms for current topology to to _delta_terms vector.
    double log_const_term = log(M_PI)*_nparams/2 - boost::math::lgamma<double>(0.5*_nparams + 1.0);
    for (unsigned i = 0; i < _grapes.size(); ++i)
        {
        Grape & g = _grapes[i];
        g._log_hyperball_volume     = log_const_term + log(g._radius)*_nparams;
        g._log_kernel               = log_likelihoods[g._index] + log_priors[g._index] + log_jacobians[g._index];
        g._log_hypercylinder_volume = g._log_hyperball_volume + g._log_kernel;
        _delta_terms.push_back(g._log_hypercylinder_volume);
        }
    std::sort(_grapes.begin(), _grapes.end(), std::greater<Grape>());
    }

void Arbor::calcRatioTerms(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors, dlb_vect_t & log_jacobians)
    {
    // Appends terms for current topology to to _log_ratio_terms vector.
    // The tmp_count_num_grapes vector counts the number of grapes to which each estimation sample belongs
    // Only useful if _grape_overlap is true; otherwise, every element will be 1
    std::vector<unsigned> tmp_count_num_grapes(_estimation_indices.size(), 0);

    _placed.clear();

    // Go through _estimation_indices and, for each point, determine if it falls within any grape. If so, append its log kernel minus
    // the representative (grape's) log kernel to the log_ratios vector.
    unsigned k = 0;
    for (auto i : _estimation_indices)
        {
        // Note: grapes are sorted so that those with greatest hypercylinder volume are first, which should result in
        // greater efficiency (unless _grape_overlap is true) because those with larger volumes should capture the
        // most estimation sample points
        for (unsigned j = 0; j < _grapes.size(); ++j)
            {
            Grape & g = _grapes[j];
            double d = g.distPointToGrapeBoundary(parameters[i]);
            if (d < 0.0)
                {
                double log_kernel = log_likelihoods[i] + log_priors[i] + log_jacobians[i];

                //std::cerr << "adding " << log_kernel << " to grape " << j << std::endl; //temporary!

                g._placed_log_kernels.push_back(log_kernel);            // record log_kernel of point placed in grape g so that variance can later be calculated
                _placed.push_back(fabs(g._log_kernel - log_kernel));
                _log_ratio_terms.push_back(g._log_kernel - log_kernel);
                tmp_count_num_grapes[k] += 1;
                }
            }
        k += 1;
        }

    // Count number of estimation sample points belonging to multiple grapes
    _in_multiple_grapes = 0;
    for (k = 0; k < _estimation_indices.size(); ++k)
        {
        if (tmp_count_num_grapes[k] > 1)
            _in_multiple_grapes++;
        }
    }

void Arbor::margLikeOneTopology(sample_vect_t & parameters, dlb_vect_t & log_likelihoods, dlb_vect_t & log_priors)
    {
    // Builds _reference_indices and _estimation_indices for this tree topology
    // and adds number of estimation points to _N)
    chooseGrapeSeeds(parameters);

    // Calculates _means and _stdevs from reference sample only
    calcMeansAndStdevs(parameters);

    // Calculates _log_jacobians (one element for each element of parameters)
    if (_correlation)
        standardizeSampleAllowingCorrelation(parameters);
    else
        standardizeSampleAssumingIndependence(parameters);

    // Builds _grapes from reference sample only; sets _radius and _v for each grape created.
    createGrapes(parameters, log_likelihoods, log_priors, _log_jacobians);

    // Calculates _log_hyperball_volume, _log_kernel, and _log_hypercylinder_volume for each grape, then
    // sorts _grapes from largest to smallest _log_hypercylinder_volume and adds terms to _delta_terms vector
    calcDeltaTerms(parameters, log_likelihoods, log_priors, _log_jacobians);

    // Processes estimation sample, adding terms to _log_ratio_terms vector and _placed vector
    calcRatioTerms(parameters, log_likelihoods, log_priors, _log_jacobians);
    }

void Arbor::calcIndivTopolMargLikes()
    {
    _N = 0;
    std::vector<double> total_placed;
    _total_sample_size = 0;
    _min_radius = std::numeric_limits<double>().max();

    _log_ratio_terms.clear();
    _delta_terms.clear();

    for (auto p : _tree_topologies)
        {
        unsigned topology  = p.second;
        unsigned frequency = p.first;
        GrapeData::SharedPtr topol_record = _db.addTopology(topology, frequency);

        _total_sample_size += frequency;
        if (frequency > _min_sample_size)
            {
            margLikeOneTopology(_parameters[topology], _log_likelihoods[topology], _log_priors[topology]);

            topol_record->_used = true;
            topol_record->_reference_samplesize = (unsigned)_reference_indices.size();
            topol_record->_estimation_samplesize = (unsigned)_estimation_indices.size();
            topol_record->_num_placed = (unsigned)_placed.size();
            topol_record->_radius = _min_radius;
            }
        }
    }

void Arbor::showTopoTable()
    {
    std::cout << "\n  Found " << _tree_topologies.size() << " tree topologies:" << std::endl;
    std::cout << boost::str(boost::format("  %10s %10s %10s %10s %10s %10s %10s\n") % "topology" % "frequency" % "rsample" % "esample" % "placed" % "%placed" % "radius");
    for (auto topol_record : _db._topologies)
        {
        if (topol_record->_used) {
            std::cout << boost::str(boost::format("  %10d %10d %10d %10d %10d %10.3f %10.3f\n") % topol_record->_id % topol_record->_frequency % topol_record->_reference_samplesize % topol_record->_estimation_samplesize % topol_record->_num_placed % topol_record->calcPercentPlaced() % topol_record->_radius);
            }
        else {
            std::cout << boost::str(boost::format("  %10d %10d below cutoff (%d samples required to consider topology)\n") % topol_record->_id % topol_record->_frequency % _min_sample_size);
            }
        }
    }

double Arbor::summationWithFloatingControl(const std::vector<double> & v) const
    {
    // Computes log of the sum (linear scale) of the (log scale) elements in supplied vector v.
    // For example, let v = [log(a), log(b), log(c)] and suppose max({a,b,c}} = c and that you desire log(a + b + c)
    // log(a + b + c) = log(c[a/c + b/c + c/c])
    //                = log(c) + log([a/c + b/c + 1])
    //                = log(c) + log[exp(log(a) - log(c)) + exp(log(b) - log(c)) + exp(log(c) - log(c))]
    assert(v.size() > 0);
    double max_v = *(std::max_element(v.begin(), v.end()));
    double sum_exp_diff = std::accumulate(v.begin(), v.end(), 0.0, plusexp(-max_v));
    return max_v + log(sum_exp_diff);
    }

void Arbor::calcOverallMargLike()
    {
    double log_denominator_eq11 = summationWithFloatingControl(_delta_terms);
    double log_marginal_likelihood = log_denominator_eq11;
    if (_log_ratio_terms.size() > 0)
        {
        double log_numerator_eq11 = summationWithFloatingControl(_log_ratio_terms);
        log_marginal_likelihood = log_denominator_eq11 - (log_numerator_eq11 - log(_N));
        }

    std::cout << boost::str(boost::format("Log marginal likelihood     = %.5f") % log_marginal_likelihood) << std::endl;
    }

void Arbor::run()
    {
    initRandom();
    showUserSpecifiedOptions();
    readData();
    showParamTable();
    calcIndivTopolMargLikes();
    showTopoTable();
    calcOverallMargLike();
    }
