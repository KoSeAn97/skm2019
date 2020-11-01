#include <iostream>
#include <stdexcept>
#include <problem/subtasks.h>
#include <problem/types.h>

double 
FindStepSize(double start, double end, int split)
{
    if (split < 4)
    {
        throw std::invalid_argument("split must be >= 4");
    }
    return (end - start) / (split - 1);
}

void 
RemoveBorder(DomainOfFunction& d, int& xSplit, int& ySplit)
{
    double xStep = FindStepSize(d.xStart, d.xEnd, xSplit);
    double yStep = FindStepSize(d.yStart, d.yEnd, ySplit);

    d.xStart += xStep;
    d.yStart += yStep;
    d.xEnd -= xStep;
    d.yEnd -= yStep;
    xSplit -= 2;
    ySplit -= 2;
}

void
SetSubtaskLine(double& start, double& end, int& split, int i, int n)
{
    if (n == 1)
    {
        return;
    }

    double step = FindStepSize(start, end, split);
    int splitPerProc = split / n;
    int tail = split - n * splitPerProc;
    split = splitPerProc;
    if (i == n - 1)
    {
        split += tail;
    }

    start += step * splitPerProc * i;
    end = start + step * (split - 1);
}

void
SetSubtask(DomainOfFunction& d, int& xSplit, int& ySplit, int nRows, int nCols, int row, int col)
{
    RemoveBorder(d, xSplit, ySplit);
    SetSubtaskLine(d.xStart, d.xEnd, xSplit, row, nRows);
    SetSubtaskLine(d.yStart, d.yEnd, ySplit, col, nCols);
}
