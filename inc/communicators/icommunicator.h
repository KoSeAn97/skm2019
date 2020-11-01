#pragma once
#include <utils/matrix.h>
#include <utils/directions.h>

class ICommunicator
{
public:
    virtual ~ICommunicator() {}
    virtual bool IsMaster() const = 0;
    virtual int GetRank() const = 0;
    virtual int CountSubtask(Dimension i) const = 0;
    virtual int GetSubtask(Dimension i) const = 0;

    virtual void UpdateBorders(Matrix& w) const = 0;
    virtual double UpdateTau(double localNumer, double localDenom) const = 0;
    virtual double UpdateNormOfDifference(double localNorm) const = 0;
    virtual double UpdateError(double localError) const = 0;
};
