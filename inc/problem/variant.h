#pragma once

namespace problem
{
extern const double xStart, xEnd, yStart, yEnd;

double k(double x, double y);
double q(double x, double y);
double u(double x, double y);
double F(double x, double y);
double phiXStart(double x);
double phiXEnd(double x);
double phiYStart(double y);
double phiYEnd(double y);
}