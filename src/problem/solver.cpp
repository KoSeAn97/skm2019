#include <iostream>
#include <cmath>
#include <boost/make_shared.hpp>
#include <problem/solver.h>
#include <communicators/mpicom.h>
#include <communicators/localcom.h>
#include <problem/subtasks.h>
#include <fstream>
#include <getopt.h>
#include <cstring>
#include <omp.h>

Solver::Solver(int& argc, char**& argv)
    : xSplit_(0)
    , ySplit_(0)
    , normOfDifference_(0)
    , maxDifference_(0)
    , iterations_(0)
    , elapsedTime_(0)
    , dump_(false)
{
    int opt;
    while (-1 != (opt = getopt(argc, argv, "m:n:d")))
    {
        switch (opt)
        {
        case 'm':
            xSplit_ = std::atoi(optarg);
            break;
        case 'n':
            ySplit_ = std::atoi(optarg);
            break;
        case 'd':
            dump_ = true;
            break;
        default:
            throw std::runtime_error("invalid arg");
            break;
        }
    }
    if (xSplit_ < 3 || ySplit_ < 3)
    {
        throw std::runtime_error("set nonzero splits");
    }

#ifdef CMC_OPENMPI
    mpi_.reset(new MpiHolder(&argc, &argv));
    com_ = boost::make_shared<MpiCommunicator>(xSplit_, ySplit_);
#else
    com_ = boost::make_shared<LocalCommunicator>();
#endif
}

void 
Solver::SetProblem(Problem p, double eps)
{
    SetSubtask(
        p.domain, 
        xSplit_, ySplit_, 
        com_->CountSubtask(Vertical), com_->CountSubtask(Horizontal),
        com_->GetSubtask(Vertical),   com_->GetSubtask(Horizontal)
    );
    m_ = boost::make_shared<FiniteElementMethod>(com_, p, xSplit_, ySplit_);
    eps_ = eps;
}

void 
Solver::Solve()
{
    double startTime = MPI_Wtime();

    Matrix w = m_->FirstApproximation();
    do
    {
        Matrix next = m_->EvaluateNext(w);
        normOfDifference_ = m_->NormOfDifference(w, next);
        w = next;
        iterations_++;
    } while (normOfDifference_ > eps_);

    double endTime = MPI_Wtime();

    elapsedTime_ = endTime - startTime;

    Matrix etalon = m_->Etalon();
    double localError = 0;
    for (int i = 0; i < w.GetRows(); i++)
    {
        for (int j = 0; j < w.GetCols(); j++)
        {
            double diff = std::fabs(w(i, j) - etalon(i, j));
            if (diff > localError)
            {
                localError = diff;
            }
        }
    }
    maxDifference_ = com_->UpdateError(localError);
    solution_ = w.Unwrapped();
}

void 
Solver::PrintReport()
{
    if (dump_)
    {
        int row = com_->GetSubtask(Vertical);
        int col = com_->GetSubtask(Horizontal);
        DumpSolution(row, col);
        DumpEtalon(row, col);
    }

    if (!com_->IsMaster())
    {
        return;
    }

    std::cout << "Max difference: " << maxDifference_ << std::endl;
    std::cout << "Norm:           " << normOfDifference_ << std::endl;
    std::cout << "Iterations:     " << iterations_ << std::endl;
    std::cout << "Elapsed time:   " << elapsedTime_ << std::endl;
}

void
Solver::DumpSolution(int row, int col) const
{
    char filename[32];
    std::sprintf(filename, "solution_%d-%d", row, col);
    std::ofstream f(filename);
    for (int i = 0; i < solution_.GetRows(); i++)
    {
        for (int j = 0; j < solution_.GetCols(); j++)
        {
            f << i << ',' << j << ',' << solution_(i, j) << std::endl;
        }
    }
}

void
Solver::DumpEtalon(int row, int col) const
{
    char filename[32];
    std::sprintf(filename, "etalon_%d-%d", row, col);
    std::ofstream f(filename);

    const Matrix etalon = m_->Etalon().Unwrapped();
    for (int i = 0; i < etalon.GetRows(); i++)
    {
        for (int j = 0; j < etalon.GetCols(); j++)
        {
            f << i << ',' << j << ',' << etalon(i, j) << std::endl;
        }
    }
}