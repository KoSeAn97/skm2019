#include <iostream>
#include <cstdlib>
#include <problem/types.h>
#include <problem/variant.h>
#include <problem/solver.h>
#include <problem/subtasks.h>
#include <problem/fem.h>
#include <utils/grid_fn.h>
#include <utils/awesome_fn.h>

static double
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

int 
_main(int argc, char** argv)
{
    const DomainOfFunction domain = 
    {
        .xStart = problem::xStart,
        .xEnd   = problem::xEnd,
        .yStart = problem::yStart,
        .yEnd   = problem::yEnd
    };

    const Problem p =
    {
        .domain = domain,
        .k = problem::k,
        .q = problem::q,
        .F = problem::F,
        .u = problem::u,
        .phiXStart = problem::phiXStart,
        .phiXEnd   = problem::phiXEnd,
        .phiYStart = problem::phiYStart,
        .phiYEnd   = problem::phiYEnd
    };

    int rv = EXIT_FAILURE;
    try
    {
        DomainOfFunction d = p.domain;
        int xSplit = 2000;
        int ySplit = 2000;
        
        SetSubtask(d, xSplit, ySplit, 1, 1, 0, 0);

        std::cout << d << std::endl << xSplit << ' ' << ySplit << std::endl;        

        Matrix w = CreateGrid(p.u, d, xSplit, ySplit);

        Matrix q_ = CreateGrid(p.q, d, xSplit, ySplit);
        Matrix F_ = CreateGrid(p.F, d, xSplit, ySplit);

        double xStep = FindStepSize(d.xStart, d.xEnd, xSplit);
        double yStep = FindStepSize(d.yStart, d.yEnd, ySplit);

        double xHalfStep = xStep / 2;
        double yHalfStep = yStep / 2;

        DomainOfFunction domainA = d;
        domainA.xStart -= xHalfStep;
        domainA.xEnd -= xHalfStep;
        Matrix a_ = CreateGrid(AwesomeFunc(p.k, xStep), domainA, xSplit, ySplit);

        DomainOfFunction domainB = d;
        domainB.yStart -= yHalfStep;
        domainB.yEnd -= yHalfStep;
        Matrix b_ = CreateGrid(AwesomeFunc(p.k, yStep), domainB, xSplit, ySplit);

        double maxR = 0;
        for (size_t i = 0; i < w.GetRows(); i++)
        {
            for (size_t j = 0; j < w.GetCols(); j++)
            {
                Matrix ws = w .Window(i, j);
                Matrix as = a_.Window(i, j);
                Matrix bs = b_.Window(i, j);

                double r = -LaplaceConvolution(ws, as, bs) + q_(i, j) * w(i, j) - F_(i, j);
                if (r > maxR)
                {
                    maxR = r;
                }
            }
        }

        std::cout << "Max R: " << maxR << std::endl;

        rv = EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    return rv;
}