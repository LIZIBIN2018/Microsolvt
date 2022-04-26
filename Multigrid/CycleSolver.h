#pragma once


template<int dim>
class Multigrid;

template<int dim>
class CycleSolver
{
public:
    CycleSolver()
    {
        
    }
private:
    RestrictionOperator *rst_opt;   // restriction
    InterpolationOperator *itp_opt; // interpolation
protected: 
    Multigrid<dim> grid_ptr;
    void solve() = 0;
};

template<int dim>
class VCycle : public CycleSolver
{
public:
    void solve() override
    {

    }
};

template<int dim>
class FullMultigridCycle : public CycleSolver
{
public:
    void solve() override
    {

    }
};