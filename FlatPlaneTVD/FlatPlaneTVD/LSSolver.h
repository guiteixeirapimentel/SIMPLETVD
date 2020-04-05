#pragma once
#pragma once

#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

typedef amgcl::backend::builtin<double> Backend;
typedef amgcl::backend::builtin<double> BackendD;

typedef amgcl::make_solver<
    // Use AMG as preconditioner:
    amgcl::amg<
    Backend,
    amgcl::coarsening::smoothed_aggregation,
    amgcl::relaxation::spai0
    >,
    // And BiCGStab as iterative solver:
    amgcl::solver::bicgstab<Backend>
> Solver;


class LSSolver
{
public:
    LSSolver();
    ~LSSolver();

    std::vector<double> SolveSparseCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
        const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut);
	
private:

};