#include <regex>
#include <fstream>
#include <boost/algorithm/string.hpp>
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

    unsigned i = 1;
    std::cout << "\nTopologies found:\n" << std::endl;
    for (auto p : _tree_topologies)
        {
        std::cout << boost::str(boost::format("%12d %12d %12d") % i % p.first % p.second) << std::endl;
        ++i;
        }
    }

void Arbor::run()
    {
    try
        {
        showUserSpecifiedOptions();
        readData();
        }
    catch (XArbor x)
        {
        std::cerr << boost::str(boost::format("%s encountered a problem:\n  %s") % _program_name % x.what()) << std::endl;
        }
    }
