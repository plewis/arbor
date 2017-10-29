#pragma once

namespace arbor
{

class Grape
    {
    public:
                                                Grape(const Grape & other);
                                                Grape(unsigned index, double radius, const std::vector<double> & v, std::vector<std::pair<unsigned, double> > & n);
                                                ~Grape();

        bool                                    operator==(const Grape & other) const;
        bool                                    operator<(const Grape & other) const;
        bool                                    operator>(const Grape & other) const;

        double                                  distPointToGrapeBoundary(std::vector<double> & point) const;

        std::vector<std::pair<unsigned, double> > _neighbors;

        std::vector<double>                     _v;                         // the point serving as the center of this Grape
        double                                  _radius;                    // the radius of this Grape
        unsigned                                _index;                     // index of center point of this Grape
        double                                  _log_kernel;                // height of hypercylinder for which this Grape is base
        double                                  _log_hyperball_volume;      // volume of the hyperball base for this Grape
        double                                  _log_hypercylinder_volume;  // log of the volume of the hypercylinder for this Grape

        std::vector<double>                     _placed_log_kernels;        // log_kernel for every estimation sample point placed in this grape
        std::vector<unsigned>                   _placed_est_indices;        // index of every estimation sample point placed in this grape
    };

inline Grape::Grape(const Grape & other) : _radius(other._radius), _index(other._index), _log_kernel(other._log_kernel), _log_hyperball_volume(other._log_hyperball_volume), _log_hypercylinder_volume(other._log_hypercylinder_volume)
    {
    _v.resize(other._v.size());
    std::copy(other._v.begin(), other._v.end(), _v.begin());
    _neighbors.resize(other._neighbors.size());
    std::copy(other._neighbors.begin(), other._neighbors.end(), _neighbors.begin());
    }

inline Grape::Grape(unsigned index, double radius, const std::vector<double> & v, std::vector<std::pair<unsigned, double> > & n) : _radius(radius), _index(index), _log_kernel(0.0), _log_hyperball_volume(0.0), _log_hypercylinder_volume(0.0)
    {
    _v.resize(v.size());
    std::copy(v.begin(), v.end(), _v.begin());

    _neighbors.resize(n.size());
    std::copy(n.begin(), n.end(), _neighbors.begin());
    }

inline Grape::~Grape()
    {
    _v.clear();
    }

inline bool Grape::operator==(const Grape & other) const
    {
    return (_log_hypercylinder_volume == other._log_hypercylinder_volume);
    }

inline bool Grape::operator<(const Grape & other) const
    {
    return (_log_hypercylinder_volume < other._log_hypercylinder_volume);
    }

inline bool Grape::operator>(const Grape & other) const
    {
    return (_log_hypercylinder_volume > other._log_hypercylinder_volume);
    }

inline double Grape::distPointToGrapeBoundary(std::vector<double> & point) const
    {
    // Calculate euclidean distance between point and the boundary of the grape
    double d = 0.0;
    for (unsigned j = 0; j < _v.size(); ++j)
        {
        d += pow(point[j] - _v[j], 2);
        }
    return sqrt(d) - _radius;
    }

}
