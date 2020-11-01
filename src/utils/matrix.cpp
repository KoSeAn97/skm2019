#include <cstring>
#include <sstream>
#include <stdexcept>
#include <utils/matrix.h>

Matrix::Matrix(int nRows, int nCols, bool wrapped)
    : nRows_(nRows)
    , nCols_(nCols)
    , stride_(nCols)
    , pinRow_(0)
    , pinCol_(0)
    , wrapped_(wrapped)
{
    if (wrapped_) 
    {
        stride_ += 2;
        pinRow_ += 1;
        pinCol_ += 1;
    }

    int size = (nCols_ + 2) * (nRows_ + 2);
    if (size) 
    {
        data_.reset(new double[size]);
        std::memset(data_.get(), 0, size * sizeof(double));
    }
}

Matrix::Matrix(Matrix const& rhs)
    : data_(rhs.data_)
    , nRows_(rhs.nRows_)
    , nCols_(rhs.nCols_)
    , stride_(rhs.stride_)
    , pinRow_(rhs.pinRow_)
    , pinCol_(rhs.pinCol_)
    , wrapped_(rhs.wrapped_)
{
}

Matrix const&
Matrix::operator=(Matrix const& rhs)
{
    if (this == &rhs) {
        return *this;
    }

    data_ = rhs.data_;
    nRows_ = rhs.nRows_;
    nCols_ = rhs.nCols_;
    stride_ = rhs.stride_;
    pinRow_ = rhs.pinRow_;
    pinCol_ = rhs.pinCol_;
    wrapped_ = rhs.wrapped_;

    return *this;
}

Matrix
Matrix::DeepCopy() const
{
    int size;
    if (wrapped_) 
    {
        size = (nRows_ + 2) * (nCols_ + 2);
    } 
    else 
    {
        size = nRows_ * nCols_;
    }

    Matrix rv(nRows_, nCols_, wrapped_);
    std::memcpy(rv.data_.get(), data_.get(), size * sizeof(double));

    return rv;
}

int Matrix::GetRows() const
{
    return nRows_;
}

int Matrix::GetCols() const
{
    return nCols_;
}

void Matrix::SetSomeBorder(Direction d, std::vector<double> const& border, int shift)
{
    if (!wrapped_ && !shift) 
    {
        throw std::invalid_argument("out of range");
    }

    const int borderSize = border.size();

    switch (d) 
    {
    case Up:
        if (borderSize != nCols_) 
        {
            std::stringstream ss;
            ss << "invalid border size: " << borderSize << ' ' << nCols_;
            throw std::invalid_argument(ss.str());
        }
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < borderSize; i++) 
        {
            (*this)(nRows_ - shift, i) = border[i];
        }
        break;

    case Down:
        if (borderSize != nCols_) 
        {
            std::stringstream ss;
            ss << "invalid border size: " << borderSize << ' ' << nCols_;
            throw std::invalid_argument(ss.str());
        }
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < borderSize; i++) 
        {
            (*this)(-1 + shift, i) = border[i];
        }
        break;

    case Left:
        if (borderSize != nRows_) 
        {
            std::stringstream ss;
            ss << "invalid border size: " << borderSize << ' ' << nCols_;
            throw std::invalid_argument(ss.str());
        }
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < borderSize; i++) 
        {
            (*this)(i, -1 + shift) = border[i];
        }
        break;

    case Right:
        if (borderSize != nRows_) 
        {
            std::stringstream ss;
            ss << "invalid border size: " << borderSize << ' ' << nCols_;
            throw std::invalid_argument(ss.str());
        }
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < borderSize; i++) 
        {
            (*this)(i, nCols_ - shift) = border[i];
        }
    }
}

void Matrix::GetSomeBorder(Direction d, std::vector<double>& border, int shift)
{
    if (!wrapped_ && !shift) 
    {
        throw std::invalid_argument("out of range");
    }

    switch (d) 
    {
    case Up:
        border.resize(nCols_);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < nCols_; i++) 
        {
            border[i] = (*this)(nRows_ - shift, i);
        }
        break;

    case Down:
        border.resize(nCols_);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < nCols_; i++) 
        {
            border[i] = (*this)(-1 + shift, i);
        }
        break;

    case Left:
        border.resize(nRows_);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < nRows_; i++) 
        {
            border[i] = (*this)(i, -1 + shift);
        }
        break;

    case Right:
        border.resize(nRows_);
#ifdef CMC_OMP
# pragma omp parallel for
#endif
        for (int i = 0; i < nRows_; i++) 
        {
            border[i] = (*this)(i, nCols_ - shift);
        }
    }
}

void Matrix::SetUnwrappedBorder(Direction d, std::vector<double> const& border)
{
    SetSomeBorder(d, border, 0);
}

void Matrix::GetUnwrappedBorder(Direction d, std::vector<double>& border)
{
    GetSomeBorder(d, border, 0);
}

void Matrix::SetBorder(Direction d, std::vector<double> const& border)
{
    SetSomeBorder(d, border, 1);
}

void Matrix::GetBorder(Direction d, std::vector<double>& border)
{
    GetSomeBorder(d, border, 1);
}

double&
Matrix::operator()(int row, int col)
{
    row += pinRow_;
    col += pinCol_;
    return data_.get()[row * stride_ + col];
}

double
Matrix::operator()(int row, int col) const
{
    row += pinRow_;
    col += pinCol_;
    return data_.get()[row * stride_ + col];
}

Matrix const
Matrix::Unwrapped() const
{
    if (!wrapped_) 
    {
        throw std::invalid_argument("matrix isn't wrapped");
    }

    Matrix rv(*this);
    rv.pinRow_ = 0;
    rv.pinCol_ = 0;
    rv.nRows_ += 2;
    rv.nCols_ += 2;
    rv.wrapped_ = false;

    return rv;
}

Matrix const
Matrix::Window(int row, int col) const
{
    if (!wrapped_) 
    {
        throw std::invalid_argument("matrix isn't wrapped");
    }
    Matrix rv(*this);
    rv.pinRow_ = row;
    rv.pinCol_ = col;
    rv.nRows_ = 3;
    rv.nCols_ = 3;
    rv.wrapped_ = false;

    return rv;
}

std::ostream&
operator<<(std::ostream& out, Matrix const& m)
{
    int nRows = m.GetRows();
    int nCols = m.GetCols();

    if (nRows * nCols == 0) 
    {
        out << "empty matrix" << std::endl;
        return out;
    }
    out << std::fixed;
    for (int i = 0; i < nRows; ++i) 
    {
        for (int j = 0; j < nCols; ++j) 
        {
            out << m(i, j) << " ";
        }
        out << std::endl;
    }
    return out;
}
