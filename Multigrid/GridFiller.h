#pragma once 

template<int dim>
class Multigrid;

template<int dim>
class GridFiller
{
protected:
    Multigrid<dim> *grid_ptr;

public:
    void fillGrid(); // TODO 用填充器来填充网格
};