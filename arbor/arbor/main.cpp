#include <iostream>
#include "arbor.hpp"

using namespace arbor;

// static data member initializations
std::string Arbor::_program_name        = "Arbor";
std::string Arbor::_author              = "Paul O. Lewis";
unsigned    Arbor::_major_version       = 1;
unsigned    Arbor::_minor_version       = 0;
unsigned    Arbor::_def_rnseed          = 1;
unsigned    Arbor::_def_minsamplesize   = 50;
double      Arbor::_def_grapefraction   = 0.1;
double      Arbor::_def_keepfraction    = 0.9;
double      Arbor::_def_forcedradius    = 0.01;

int main(int argc, const char * argv[])
    {
    //std::cerr << "Current working directory: " << getcwd(NULL, 0) << std::endl;

    Arbor arbor;
    try {
        arbor.processCommandLineOptions(argc, argv);
        arbor.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
    }
