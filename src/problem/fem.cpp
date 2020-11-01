#include <iostream>
#include <fstream>
#include <problem/fem.h>
#include <problem/types.h>
#include <utils/grid_fn.h>
#include <problem/subtasks.h>
#include <fstream>
#include <utils/awesome_fn.h>
#include <cstring>

static double LaplaceConvolution(const Matrix& w, const Matrix& a, const Matrix& b);
Matrix FirstApproximation(const Problem& p, const int xSplit, const int ySplit);

FiniteElementMethod::FiniteElementMethod(boost::shared_ptr<ICommunicator> c, const Problem& p, int xSplit, int ySplit)
    : com_(c)
{
    double xStep = FindStepSize(p.domain.xStart, p.domain.xEnd, xSplit);
    double yStep = FindStepSize(p.domain.yStart, p.domain.yEnd, ySplit);

#ifdef CMC_DEBUG
    {
        char x[20];
        std::sprintf(x, "%d_%d", com_->GetSubtask(Vertical), com_->GetSubtask(Horizontal));
        std::ofstream f(x);
        f << p.domain << std::endl;
        f << xStep << std::endl << yStep << std::endl;
    }
#endif

    double xHalfStep = xStep / 2;
    double yHalfStep = yStep / 2;

    q_ = CreateGrid(p.q, p.domain, xSplit, ySplit);
    F_ = CreateGrid(p.F, p.domain, xSplit, ySplit);
    u_ = CreateGrid(p.u, p.domain, xSplit, ySplit);

    first_ = ::FirstApproximation(p, xSplit, ySplit);

    DomainOfFunction domainA = p.domain;
    domainA.xStart -= xHalfStep;
    domainA.xEnd -= xHalfStep;
    a_ = CreateGrid(AwesomeFunc(p.k, xStep), domainA, xSplit, ySplit);

    DomainOfFunction domainB = p.domain;
    domainB.yStart -= yHalfStep;
    domainB.yEnd -= yHalfStep;
    b_ = CreateGrid(AwesomeFunc(p.k, yStep), domainB, xSplit, ySplit);

#ifdef CMC_DEBUG
    {
        char x[20];
        std::sprintf(x, "%d", com_->GetRank());
        std::ofstream f(x, std::ios::app);

        f << "a:\n" << a_ << std::endl;
        f << "b:\n" << b_ << std::endl;
        f << "q:\n" << q_ << std::endl;
        f << "F:\n" << F_ << std::endl;
    }
#endif
}

Matrix 
FiniteElementMethod::FirstApproximation() const
{
    return first_.DeepCopy();
}

Matrix 
FiniteElementMethod::Etalon() const
{
    return u_.DeepCopy();
}

Matrix 
FiniteElementMethod::EvaluateNext(const Matrix& w) const
{
    Matrix next = w.DeepCopy();
#ifdef CMC_DEBUG
    {
        char x[20];
        std::sprintf(x, "%d", com_->GetRank());
        std::ofstream f(x, std::ios::app);
        f << "copy:\n" << next.Unwrapped() << std::endl;
    }
#endif
    com_->UpdateBorders(next);
#ifdef CMC_DEBUG
    {
        char x[20];
        std::sprintf(x, "%d", com_->GetRank());
        std::ofstream f(x, std::ios::app);
        f << "upd:\n" << next.Unwrapped() << std::endl;
    }
#endif
    Matrix r = EvaluateResidual(next);

    double t = EvaluateTau(r);

    const int nRows = next.GetRows(), nCols = next.GetCols();
#ifdef CMC_OMP
# ifdef BLUEGENE_KLUDGE
#  pragma omp parallel for
# else
#  pragma omp parallel for collapse(2)
# endif
#endif
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            next(i, j) -= t * r(i, j);
        }
    }

    return next;
}

double 
FiniteElementMethod::NormOfDifference(const Matrix& prev, const Matrix& curr) const
{
    const int nRows = prev.GetRows();
    const int nCols = prev.GetCols();
    if ( curr.GetRows() != nRows || curr.GetCols() != nCols)
    {
        throw std::runtime_error("matrix sizes differ");
    }

    double localNorm = 0;

#ifdef CMC_OMP
# ifdef BLUEGENE_KLUDGE
#  pragma omp parallel for reduction(+:localNorm)
# else
#  pragma omp parallel for reduction(+:localNorm) collapse(2)
# endif
#endif
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            double diff = prev(i, j) - curr(i, j);
            localNorm += diff * diff;
        }
    }

    return com_->UpdateNormOfDifference(localNorm);
}

Matrix 
FiniteElementMethod::EvaluateResidual(const Matrix& w) const
{
    Matrix r(w.GetRows(), w.GetCols(), true);

    const int nRows = r.GetRows(), nCols = r.GetCols();
#ifdef CMC_OMP
# ifdef BLUEGENE_KLUDGE
#  pragma omp parallel for
# else
#  pragma omp parallel for collapse(2)
# endif
#endif
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            Matrix ws = w .Window(i, j);
            Matrix as = a_.Window(i, j);
            Matrix bs = b_.Window(i, j);

            r(i, j) = -LaplaceConvolution(ws, as, bs) + q_(i, j) * w(i, j) - F_(i, j);
        }
    }
#ifdef CMC_DEBUG
    {
        char x[20];
        std::sprintf(x, "%d", com_->GetRank());
        std::ofstream f(x, std::ios::app);
        f << "res:\n" << r.Unwrapped() << std::endl;
    }
#endif
    com_->UpdateBorders(r);

#ifdef CMC_DEBUG
    {        
        char x[20];
        std::sprintf(x, "%d", com_->GetRank());
        std::ofstream f(x, std::ios::app);
        f << "updr:\n" << r.Unwrapped() << std::endl;
    }
#endif
    return r;
}

double
FiniteElementMethod::EvaluateTau(const Matrix& r) const
{
    double tau_numer = 0; // scalar_product(Ar, r);
    double tau_denum = 0; // scalar_product(Ar, Ar);

    const int nRows = r.GetRows(), nCols = r.GetCols();
#ifdef CMC_OMP
# ifdef BLUEGENE_KLUDGE
#  pragma omp parallel for reduction(+:tau_numer,tau_denum)
# else
#  pragma omp parallel for reduction(+:tau_numer,tau_denum) collapse(2)
# endif
#endif
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            Matrix rs = r .Window(i, j);
            Matrix as = a_.Window(i, j);
            Matrix bs = b_.Window(i, j);

            double ar = -LaplaceConvolution(rs, as, bs) + q_(i, j) * r(i, j);
            tau_numer += ar * r(i, j);
            tau_denum += ar * ar;
        }
    }

    return com_->UpdateTau(tau_numer, tau_denum);
}

double
LaplaceConvolution(Matrix const& w, Matrix const& a, Matrix const& b)
{
    double k = a(1, 1) + a(2, 1) + b(1, 1) + b(1, 2);
    return
        w(2, 1) * a(2, 1) +
        w(0, 1) * a(1, 1) +
        w(1, 2) * b(1, 2) +
        w(1, 0) * b(1, 1) +
        w(1, 1) * -k;
    return 0;
}

Matrix 
FirstApproximation(const Problem& p, const int xSplit, const int ySplit)
{
    Matrix m(xSplit, ySplit, true);
    double xStep = FindStepSize(p.domain.xStart, p.domain.xEnd, xSplit);
    double yStep = FindStepSize(p.domain.yStart, p.domain.yEnd, ySplit);
    
    {
        std::vector<double> up(ySplit), down(ySplit);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < ySplit; i++) 
        {
            double y = p.domain.yStart + i * yStep;
            up[i] = p.phiXEnd(y);
            down[i] = p.phiXStart(y);
        }
        m.SetUnwrappedBorder(Up, up);
        m.SetUnwrappedBorder(Down, down);
    }
    {
        std::vector<double> left(xSplit), right(xSplit);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < xSplit; i++)
        {
            double x = p.domain.xStart + i * xStep;
            left[i] = p.phiYStart(x);
            right[i] = p.phiYEnd(x);
        }
        m.SetUnwrappedBorder(Left, left);
        m.SetUnwrappedBorder(Right, right);
    }

    return m;
}
