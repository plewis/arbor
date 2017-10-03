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

        // member functions
        void showUserSpecifiedOptions() const;

        // data members
        std::string         _data_file_name;            // the name of the sample file to be processed
        unsigned            _min_sample_size;           // minimum number of samples for tree topology inclusion
        unsigned            _rnseed;                    // pseudorandom number seed
        double              _grape_fraction;            // fraction of sample used to define grapes

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
    {
    std::cout << boost::str(boost::format("%s %d.%d (written by %s)") % _program_name % _major_version % _minor_version % _author) << std::endl;
    }

inline Arbor::~Arbor()
    {
    }
    
}
