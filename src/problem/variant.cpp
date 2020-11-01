#include <cmath>
#include <problem/variant.h>

namespace problem
{
double const
    xStart = -1,   xEnd = 2,
    yStart = -2,   yEnd = 2;

double k(double x, double y)
{
    return 4 + x;
}

double q(double x, double y)
{
    double z = x + y;
    return z * z;
}

double u(double x, double y)
{
    double z = x + y;
    return std::exp(1 - z * z);
}

// 2 e^(1 - (x + y)^2) (-4 - 2 x + 8 x^2 + 2 x^3 - y + 16 x y + 4 x^2 y + 8 y^2 + 2 x y^2)
// +
// 2 e^(1 - x^2 - 2 x y - y^2) (4 + x) (-1 + 2 x^2 + 4 x y + 2 y^2)
// = 
// 2 e^(1 - (x + y)^2) (-8 + x (-3 + 4 x (4 + x)) - y + 8 x (4 + x) y + 4 (4 + x) y^2)
// 
double F(double x, double y)
{
    double x_plus_y = x + y;
    double x_plus_4 = x + 4;
    double laplace = 
        2 * std::exp(1 - x_plus_y * x_plus_y) *
        (-8 + x * (-3 + 4 * x * x_plus_4) - y + 8 * x * x_plus_4 * y + 4 * x_plus_4 * y * y);

    return -laplace + q(x, y) * u(x, y);
}

double phiYEnd(double x) 
{
    return u(x, yEnd);
}

double phiYStart(double x)
{
    return u(x, yStart);
}

double phiXEnd(double y)
{
    return u(xEnd, y);
}

double phiXStart(double y)
{
    return u(xStart, y);
}
}