#pragma once
#include "TransferOperator.h"
#include <functional>
template<int dim>
class Multigrid;



template<int dim>
class CycleSolver
{
public:     //ctor & dtor
    CycleSolver(Multigrid<dim> *grid_ptr_) : grid_ptr(grid_ptr_){ }
    virtual ~CycleSolver() { grid_ptr = nullptr;}
private:    // data
    RestrictionOperator<dim> *rst_opt;   // restriction
    InterpolationOperator<dim> *itp_opt; // interpolation
protected:  // methods
    Multigrid<dim> *grid_ptr;
    virtual void solve(std::function<double(double)>) = 0;
    virtual void solve(std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>) = 0;
};



template<int dim>
class VCycle : public CycleSolver<dim>
{
public:
    VCycle(Multigrid<dim> *grid_ptr_): CycleSolver<dim>(grid_ptr_) { }
    virtual void solve(std::function<double(double)>) override
    {

    }
    virtual void solve(std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>) override
    {

    }
};
};



template<int dim>
class FullMultigridVCycle : public CycleSolver<dim>
{
public:
    FullMultigridVCycle(Multigrid<dim> *grid_ptr_): CycleSolver<dim>(grid_ptr_) { }
    virtual void solve(std::function<double(double)>) override
    {

    }
    virtual void solve(std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>) override
    {

    }
};