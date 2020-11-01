#pragma once
#include <mpi.h>
#include <communicators/icommunicator.h>

class MpiCommunicator
    : public ICommunicator
{
public:
    MpiCommunicator(int xSplit, int ySplit);

    bool IsMaster() const;
    int GetRank() const;
    int CountSubtask(Dimension i) const;
    int GetSubtask(Dimension i) const;

    void UpdateBorders(Matrix& w) const;
    double UpdateTau(double localNumer, double localDenom) const;
    double UpdateNormOfDifference(double localNorm) const;
    double UpdateError(double localError) const;

private:
    void SetUpGridCommunicator(int xSplit, int ySplit);

    MPI_Comm communicator_;
    int rank_;
    int dims_[2];
    int coords_[2];
    int neighbors_[2][2];
};
