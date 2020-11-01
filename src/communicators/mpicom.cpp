#include <stdexcept>
#include <iostream>
#include <utility>
#include <cmath>
#include <communicators/mpicom.h>
#include <utils/directions.h>

std::pair<int, int>
ComputeDimensions(int size);

MpiCommunicator::MpiCommunicator(int xSplit, int ySplit)
{
    SetUpGridCommunicator(xSplit, ySplit);

    MPI_Comm_rank(communicator_, &rank_);
    MPI_Cart_coords(communicator_, rank_, 2, coords_);
    MPI_Cart_shift(communicator_, Vertical,   1, &neighbors_[Vertical][0],   &neighbors_[Vertical][1]  );
    MPI_Cart_shift(communicator_, Horizontal, 1, &neighbors_[Horizontal][0], &neighbors_[Horizontal][1]);

    if (IsMaster())
    {
        std::cerr << "rows: " << dims_[0] << std::endl;
        std::cerr << "cols: " << dims_[1] << std::endl;
    }
}

bool 
MpiCommunicator::IsMaster() const
{
    return rank_ == 0;
}

int 
MpiCommunicator::GetRank() const
{
    return rank_;
}

int 
MpiCommunicator::CountSubtask(Dimension i) const
{
    return dims_[i];
}

int 
MpiCommunicator::GetSubtask(Dimension i) const
{
    return coords_[i];
}

void 
MpiCommunicator::UpdateBorders(Matrix& w) const
{
    Direction directionMap[2][2] =
    {
        { Down, Up    }, // Vertical
        { Left, Right }  // Horizontal
    };

    int sizes[2] =
    { 
        w.GetCols(), // inverted
        w.GetRows()
    };

    std::vector<double> in, out;
    int maxSize = std::max(sizes[0], sizes[1]);
    in.reserve(maxSize);
    out.reserve(maxSize);

    for (int dimIdx = 0; dimIdx < 2; dimIdx++) // Vertical, Horizontal
    {
        for (int j = 0; j < 2; j++) // Torwards, Backwards
        {
            int dirIdx = (coords_[dimIdx] + j) % 2; // odds and evens
            int neighbor = neighbors_[dimIdx][dirIdx];
            Direction d = directionMap[dimIdx][dirIdx];
            
            if (neighbor == MPI_PROC_NULL)
            {
                continue;
            }

            MPI_Status _;
            w.GetBorder(d, out);
            in.resize(sizes[dimIdx]);
            MPI_Sendrecv(
                &out[0], sizes[dimIdx], MPI_DOUBLE, neighbor, 0,
                &in[0],  sizes[dimIdx], MPI_DOUBLE, neighbor, 0,
                communicator_, &_
            );
            w.SetUnwrappedBorder(d, in);
        }
    }
}

double 
MpiCommunicator::UpdateTau(double localNumer, double localDenom) const
{
    double numer, denom;
    MPI_Allreduce(&localNumer, &numer, 1, MPI_DOUBLE, MPI_SUM, communicator_);
    MPI_Allreduce(&localDenom, &denom, 1, MPI_DOUBLE, MPI_SUM, communicator_);

    return numer / denom;
}

double 
MpiCommunicator::UpdateNormOfDifference(double localNorm) const
{
    double norm;
    MPI_Allreduce(&localNorm, &norm, 1, MPI_DOUBLE, MPI_SUM, communicator_);

    return norm;
}

double 
MpiCommunicator::UpdateError(double localError) const
{
    double error;
    MPI_Allreduce(&localError, &error, 1, MPI_DOUBLE, MPI_MAX, communicator_);

    return error;
}

void
MpiCommunicator::SetUpGridCommunicator(int xSplit, int ySplit)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int allDots = (xSplit - 2) * (ySplit - 2);
    std::pair<int, int> d = allDots <= size
        ? std::make_pair(xSplit - 2, ySplit - 2)
        : ComputeDimensions(size);

    int periods[2] = {0, 0};
    dims_[0] = d.first;
    dims_[1] = d.second;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims_, periods, 1, &communicator_);
}

std::pair<int, int>
ComputeDimensions(int size)
{
    int x = std::sqrt(size);
    while (size % x)
    {
        x--;
    }
    return std::make_pair(size / x, x);
}
