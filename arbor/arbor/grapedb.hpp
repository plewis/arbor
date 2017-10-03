#pragma once

#include <map>
#include <memory>

namespace arbor
{

struct GrapeData
    {
                    GrapeData();
        
        unsigned    _id;
        unsigned    _frequency;
        unsigned    _reference_samplesize;
        unsigned    _estimation_samplesize;
        unsigned    _num_placed;
        double      _radius;
        bool        _used;

        double      calcPercentPlaced();

        typedef std::shared_ptr<GrapeData> SharedPtr;
    };

inline GrapeData::GrapeData() :
  _frequency(0),
  _reference_samplesize(0),
  _estimation_samplesize(0),
  _num_placed(0),
  _radius(0.0),
  _used(false)
    {
    }

inline double GrapeData::calcPercentPlaced()
    {
    return (100.0*_num_placed/_estimation_samplesize);
    }

class GrapeDatabase
    {
    public:
                                                    GrapeDatabase();
                                                    ~GrapeDatabase();

        std::vector<GrapeData::SharedPtr>            _topologies;

        GrapeData::SharedPtr                         addTopology(unsigned id, unsigned freq);

        typedef std::shared_ptr<GrapeDatabase> SharedPtr;

    };

inline GrapeDatabase::GrapeDatabase()
    {
    }
    
inline GrapeDatabase::~GrapeDatabase()
    {
    }

inline GrapeData::SharedPtr GrapeDatabase::addTopology(unsigned id, unsigned freq)
    {
    GrapeData::SharedPtr p = GrapeData::SharedPtr(new GrapeData());
    p->_id        = id;
    p->_frequency = freq;
    _topologies.push_back(p);
    return p;
    }


}

