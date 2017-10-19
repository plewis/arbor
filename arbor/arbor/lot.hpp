#pragma once

#include <ctime>
#include <memory>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace arbor
{

class Lot
    {
    public:
                                Lot();
                                ~Lot();
        
        void                    setSeed(unsigned seed);
        double                  uniform();
        double                  normal();
        double                  logUniform();
        
        typedef std::shared_ptr<Lot> SharedPtr;

    private:
    
        typedef boost::variate_generator<boost::mt19937 &, boost::uniform_real<> > uniform_variate_generator_t;
        typedef boost::variate_generator<boost::mt19937 &, boost::normal_distribution<> > normal_variate_generator_t;

        unsigned                        _seed;
        boost::mt19937                  _generator;
        uniform_variate_generator_t *   _uniform_variate_generator;
        normal_variate_generator_t  *   _normal_variate_generator;
    };

inline Lot::Lot() : _seed(0), _generator(1), _uniform_variate_generator(0)
    {
    _generator.seed(static_cast<unsigned int>(std::time(0)));
    _uniform_variate_generator = new uniform_variate_generator_t(_generator, boost::uniform_real<>(0,1));
    _normal_variate_generator = new normal_variate_generator_t(_generator, boost::normal_distribution<>(0,1));
    }
    
inline Lot::~Lot()
    {
    }
    
inline void Lot::setSeed(unsigned seed)
    {
    _seed = seed;
    _generator.seed(_seed > 0 ? _seed : static_cast<unsigned int>(std::time(0)));
    }
    
inline double Lot::uniform()
    {
    return (*_uniform_variate_generator)();
    }

inline double Lot::normal()
    {
    return (*_normal_variate_generator)();
    }

inline double Lot::logUniform()
    {
    double u = (*_uniform_variate_generator)();
    assert(u > 0.0);
    return std::log(u);
    }

}