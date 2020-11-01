#include <utils/grid_fn.h>
#include <problem/subtasks.h>

Matrix 
CreateGrid(const Function2& f, const DomainOfFunction& d, const int xSplit, const int ySplit)
{
    double xStep = FindStepSize(d.xStart, d.xEnd, xSplit);
    double yStep = FindStepSize(d.yStart, d.yEnd, ySplit);

    Matrix m(xSplit, ySplit, true);
#ifdef CMC_OMP
# ifdef BLUEGENE_KLUDGE
#  pragma omp parallel for
# else
#  pragma omp parallel for collapse(2)
# endif
#endif
    for (int i = -1; i <= xSplit; i++)
    {
        for (int j = -1; j <= ySplit; j++)
        {
            m(i, j) = f(d.xStart + i * xStep, d.yStart + j * yStep);
        }
    }

    return m;
}