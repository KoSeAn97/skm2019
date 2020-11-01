#include <communicators/localcom.h>

bool 
LocalCommunicator::IsMaster() const
{
    return true;
}

int 
LocalCommunicator::GetRank() const
{
    return 0;
}

int 
LocalCommunicator::CountSubtask(Dimension) const
{
    return 1;
}

int 
LocalCommunicator::GetSubtask(Dimension) const
{
    return 0;
}

void 
LocalCommunicator::UpdateBorders(Matrix&) const
{
}

double 
LocalCommunicator::UpdateTau(double localNumer, double localDenom) const
{
    return localNumer / localDenom;
}

double 
LocalCommunicator::UpdateNormOfDifference(double localNorm) const
{
    return localNorm;
}

double 
LocalCommunicator::UpdateError(double localError) const
{
    return localError;
}
