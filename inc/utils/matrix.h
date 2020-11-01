#pragma once
#include <ostream>
#include <vector>
#include <boost/shared_array.hpp>
#include <utils/directions.h>

class Matrix {
public:
    Matrix(int n_rows = 0, int n_cols = 0, bool wrapped = false);
    Matrix(Matrix const& rhs);
    Matrix const& operator = (Matrix const& rhs);

    Matrix DeepCopy() const;

    int GetRows() const;
    int GetCols() const;

    void SetUnwrappedBorder(Direction d, const std::vector<double>& border);
    void GetUnwrappedBorder(Direction d, std::vector<double>& border);
    void SetBorder(Direction d, const std::vector<double>& border);
    void GetBorder(Direction d, std::vector<double>& border);

    double& operator()(int row, int col);
    double operator()(int row, int col) const;

    Matrix const Unwrapped() const;
    Matrix const Window(int row, int col) const;

private:
    void SetSomeBorder(Direction d, const std::vector<double>& border, int shift);
    void GetSomeBorder(Direction d, std::vector<double>& border, int shift);

    boost::shared_array<double> data_;
    int nRows_, nCols_, stride_, pinRow_, pinCol_;
    bool wrapped_;
};

std::ostream& operator<<(std::ostream& out, Matrix const& m);
