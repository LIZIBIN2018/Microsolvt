#pragma once


template<int dim>
class Multigrid;

template<int dim>
class Cycle
{
protected: 
    Multigrid<dim> grid_ptr;
    void solve() = 0;
};

template<int dim>
class VCycle : public Cycle
{
public:
    void solve() override
    {

    }
};

template<int dim>
class FullMultigridCycle : public Cycle
{
public:
    void solve() override
    {

    }
};