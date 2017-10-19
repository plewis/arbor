#pragma once

namespace arbor
{

class PairDist
    {
    public:
                                            PairDist(unsigned dim);
                                            ~PairDist();

        void                                setPairwiseDistance(unsigned row, unsigned col, float dist);
        float                               getPairwiseDistance(unsigned row, unsigned col) const;

        typedef std::shared_ptr<PairDist>   SharedPtr;

    private:

        unsigned                            getIndex(unsigned row, unsigned col) const;

        unsigned                            _n;
        std::vector<float>                 _d;
    };

inline PairDist::PairDist(unsigned n)
    {
    _n = n;
    unsigned dim = n*(n-1)/2;
    _d.resize(dim, 0.0);
    }

inline PairDist::~PairDist()
    {
    }

inline unsigned PairDist::getIndex(unsigned row, unsigned col) const
    {
    // Access element at row i, column j, for dimension n = 5 when only upper triangle is stored
    //
    //        j=0   j=1   j=2   j=3   j=4
    //      +-----+-----+-----+-----+-----+
    // i=0  |     |  0  |  1  |  2  |  3  |  0 =                       = 0*5 - (0)
    //      +-----+-----+-----+-----+-----+
    // i=1  |     |     |  4  |  5  |  6  |  4 = (5-1)                 = 1*5 - (0 + 1)
    //      +-----+-----+-----+-----+-----+
    // i=2  |     |     |     |  7  |  8  |  7 = (5-1) + (5-2)         = 2*5 - (0 + 1 + 2)
    //      +-----+-----+-----+-----+-----+
    // i=3  |     |     |     |     |  9  |  9 = (5-1) + (5-2) + (5-3) = 3*5 - (0 + 1 + 2 + 3)
    //      +-----+-----+-----+-----+-----+
    // i=4  |     |     |     |     |     |       1st element on row i = n*i - i*(i+1)/2
    //      +-----+-----+-----+-----+-----+
    //
    assert(row < _n);
    assert(col < _n);
    assert(row != col);
    unsigned i = row;
    unsigned j = col;
    if (row > col)
        {
        i = col;
        j = row;
        }
    unsigned k = _n*i - i*(i+1)/2 + j - i - 1;
    return k;
    }

inline float PairDist::getPairwiseDistance(unsigned row, unsigned col) const
    {
    unsigned k = getIndex(row, col);
    assert(k < _d.size());
    return _d[k];
    }

inline void PairDist::setPairwiseDistance(unsigned row, unsigned col, float dist)
    {
    unsigned k = getIndex(row, col);
    assert(k < _d.size());
    _d[k] = dist;
    }

}
