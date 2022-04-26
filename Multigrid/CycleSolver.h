#pragma once
#include "TransferOperator.h"
#include "Multigrid.h"
#include <functional>

template<int dim>
class Multigrid;



template<int dim>
class CycleSolver
{
public:     //ctor & dtor
    CycleSolver(Multigrid<dim> *grid_ptr_) : grid_ptr(grid_ptr_)
    { 
        if(grid_ptr->grid_rst_opt_type == RstOptType::injection)
            itp_opt = new InjectionOperator(this);
        else if(grid_ptr->grid_rst_opt_type == RstOptType::fullWeighting)
            itp_opt = new FullWeightOperator(this);
        if(grid_ptr->grid_itp_opt_type == ItpOptType::linear)
            itp_opt = new LinearInterpolationOperator(this);
        else if(grid_ptr->grid_itp_opt_type == ItpOptType::quadratic)
            itp_opt = new QuadraticInterpolationOperator(this);
        
    }
    virtual ~CycleSolver() { grid_ptr = nullptr;}
    virtual void solve(std::function<double(double)>) = 0;
    virtual void solve(std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>) = 0;
private:    // data
    RestrictionOperator<dim>   *rst_opt;   // restriction
    InterpolationOperator<dim> *itp_opt;   // interpolation
protected:  // methods
    Multigrid<dim> *grid_ptr;

};



template<int dim>
class VCycle : public CycleSolver<dim>
{
public:
    VCycle(Multigrid<dim> *grid_ptr_): CycleSolver<dim>(grid_ptr_) { }
    virtual void solve(std::function<double(double)>) override        //1d
    {
    
    }
    virtual void solve(std::function<double(double,double)>,  
                       std::function<double(double,double)>,
                       std::function<double(double,double)>,
                       std::function<double(double,double)>) override  //2d
    {

    }

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