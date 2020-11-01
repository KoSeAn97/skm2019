#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <problem/fem.h>
#include <problem/types.h>
#include <communicators/icommunicator.h>
#include <utils/mpi_holder.h>
#include <utils/matrix.h>

class Solver
{
public:
    Solver(int& argc, char**& argv);

    void SetProblem(Problem p, double eps);
    void Solve();
    void PrintReport();

private:
    void DumpSolution(int row, int col) const;
    void DumpEtalon(int row, int col) const;

    boost::scoped_ptr<MpiHolder> mpi_;
    boost::shared_ptr<ICommunicator> com_;
    boost::shared_ptr<FiniteElementMethod> m_;

    int xSplit_, ySplit_;
    double eps_;
    double normOfDifference_, maxDifference_;
    size_t iterations_;
    double elapsedTime_;
    Matrix solution_;
    bool dump_;
};