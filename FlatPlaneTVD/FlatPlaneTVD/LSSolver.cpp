#include "LSSolver.h"

std::vector<double> LSSolver::SolveSparseCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
    const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut)
{
    Solver solve(std::tie(nCol, ptr, col, val));

    std::vector<double> x(nCol, 0.0);

    std::tie(nItOut, errorOut) = solve(rhs, x);

    return x;
}

LSSolver::LSSolver()
{}
LSSolver::~LSSolver()
{}