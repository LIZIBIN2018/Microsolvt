#pragma once
#include "TransferOperator.h"
#include <functional>
#include "Multigrid.h"

enum class SolverType
{
    VCycle,
    FullMultigridVCycle
};

template <int dim>
class CycleSolver
{
public: // ctor & dtor
    CycleSolver(Json::Value &root)
    {
        try
        {
            auto rstOptStr = root["restriction_operator"].asString();
            auto itpOptStr = root["interpolation_operator"].asString();
            if (rstOptStr == "injection")
                rst_opt = new InjectionOperator<dim>();
            else if (rstOptStr == "fullWeighting")
                rst_opt = new FullWeightOperator<dim>();
            else
                throw 1;
            if (itpOptStr == "linear")
                itp_opt = new LinearInterpolationOperator<dim>();
            else if (itpOptStr == "quadratic")
                itp_opt = new QuadraticInterpolationOperator<dim>();
            else
                throw 1;
        }
        catch (...)
        {
            std::cerr << "Invalid transfer opetator" << '\n';
        }

        try
        {
            auto cycle_type = root["cycle"].asString();
            if (cycle_type == "V-cycle")
                solver_type = SolverType::VCycle;
            else if (cycle_type == "FullMultigridVCycle")
                solver_type = SolverType::FullMultigridVCycle;
            else
                throw 1;
        }
        catch (...)
        {
            std::cerr << "Invalid cycle type" << std::endl;
            exit(1);
        }

    }
    virtual ~CycleSolver() { }

    virtual void solve(const Multigrid<dim> &grid,std::function<double(double)> f)
    {

    }
    virtual void solve(const Multigrid<dim> &grid,std::function<double(double, double)> f)
    {
        
    }

private:                                 // data
    RestrictionOperator<dim>   *rst_opt;   // restriction
    InterpolationOperator<dim> *itp_opt; // interpolation
    SolverType                  solver_type; 

    
};
