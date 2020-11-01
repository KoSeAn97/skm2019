#pragma once
#include <problem/types.h>

class AwesomeFunc
{
public:
    AwesomeFunc(Function2 f, double step)
        : f_(f)
        , coeff_(step * step)
    {}
    double operator () (double x, double y)
    {
        return f_(x, y) / coeff_;
    }

private:
    Function2 f_;
    double coeff_;
};
