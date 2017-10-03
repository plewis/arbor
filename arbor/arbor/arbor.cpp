#include <stdio.h>
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

void Arbor::run()
    {
    try
        {
        showUserSpecifiedOptions();
        }
    catch (XArbor x)
        {
        std::cerr << boost::str(boost::format("%s encountered a problem:\n  %s") % _program_name % x.what()) << std::endl;
        }
    }
