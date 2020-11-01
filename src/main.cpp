#include <iostream>
#include <cstdlib>
#include <problem/types.h>
#include <problem/variant.h>
#include <problem/solver.h>
#include <problem/subtasks.h>

const double eps = 0.000001;

int 
main(int argc, char** argv)
{
    const DomainOfFunction domain = 
    {
        problem::xStart,
        problem::xEnd,
        problem::yStart,
        problem::yEnd
    };

    const Problem p =
    {
        domain,
        problem::k, problem::q, problem::F, problem::u,
        problem::phiXStart, problem::phiXEnd, problem::phiYStart, problem::phiYEnd
    };

    int rv = EXIT_FAILURE;
    try
    {
        Solver solver(argc, argv);
        
        solver.SetProblem(p, eps);
        solver.Solve();
        solver.PrintReport();

        rv = EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    return rv;
}