#pragma once
#include <stdint.h>
#include <ostream>
#include <boost/function.hpp>

typedef boost::function<double(double)> Function1;
typedef boost::function<double(double, double)> Function2;

struct DomainOfFunction
{
    double xStart, xEnd, yStart, yEnd;
};

struct Problem
{
    DomainOfFunction domain;
    Function2 k, q, F, u;
    Function1 phiXStart, phiXEnd, phiYStart, phiYEnd;
};

inline std::ostream& 
operator << (std::ostream& os, const DomainOfFunction& d)
{
    os << "xStart: " << d.xStart << std::endl;
    os << "xEnd:   " << d.xEnd   << std::endl;
    os << "yStart: " << d.yStart << std::endl;
    os << "yEnd:   " << d.yEnd   << std::endl;

    return os;
}
