#pragma once
#include <problem/types.h>
#include <utils/matrix.h>

Matrix CreateGrid(const Function2& f, const DomainOfFunction& d, const int xSplit, const int ySplit);
