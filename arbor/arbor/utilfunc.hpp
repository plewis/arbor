#pragma once 

namespace arbor
{

struct plusexp
    {
    plusexp() : _a(0.0) {}
    plusexp(double a) : _a(a) {}
    double operator()(double x, double y) const
        {return x + exp(y + _a);}
    double _a;
    };

struct pluslog
    {
    pluslog() : _a(0.0) {}
    pluslog(double a) : _a(a) {}
    double operator()(double x, double y) const
        {return x + log(y + _a);}
    double _a;
    };

template<class T1, class T2>
struct squared_difference : public std::binary_function<T1, T2, T1>
    {
    T1 operator()(T1 first, T2 second) const {
        return std::pow(first - second, 2.0);
        }
    };

bool UintDblPairLargestFirst(const std::pair<unsigned, double> & one, const std::pair<unsigned, double> & two)
    {
    return (one.second > two.second);
    }

bool UintDblPairSmallestFirst(const std::pair<unsigned, double> & one, const std::pair<unsigned, double> & two)
    {
    return (one.second < two.second);
    }

}