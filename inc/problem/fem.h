#pragma once
#include <boost/shared_ptr.hpp>
#include <communicators/icommunicator.h>
#include <problem/types.h>

class FiniteElementMethod
{
public:
    FiniteElementMethod(boost::shared_ptr<ICommunicator> c, const Problem& p, int xSplit, int ySplit);

    Matrix FirstApproximation() const;
    Matrix Etalon() const;
    Matrix EvaluateNext(const Matrix& w) const;
    double NormOfDifference(const Matrix& prev, const Matrix& curr) const;

private:
    Matrix EvaluateResidual(const Matrix& w) const;
    double EvaluateTau(const Matrix& r) const;

    boost::shared_ptr<ICommunicator> com_;

    Matrix a_, b_, q_, F_, u_, first_;
};

Matrix FirstApproximation(const Problem& p, int xSplit, int ySplit);
