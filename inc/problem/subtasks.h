#pragma once
#include <problem/types.h>

double FindStepSize(double start, double end, int split);
void RemoveBorder(DomainOfFunction& d, int& xSplit, int& ySplit);
void SetSubtask(DomainOfFunction& d, int& xSplit, int& ySplit, int nRows, int nCols, int row, int col);
