#pragma once
#include <communicators/icommunicator.h>
#include <utils/directions.h>

class LocalCommunicator
    : public ICommunicator
{
public:
    bool IsMaster() const;
    int GetRank() const;
    int CountSubtask(Dimension i) const;
    int GetSubtask(Dimension i) const;

    void UpdateBorders(Matrix& w) const;
    double UpdateTau(double localNumer, double localDenom) const;
    double UpdateNormOfDifference(double localNorm) const;
    double UpdateError(double localError) const;
};
