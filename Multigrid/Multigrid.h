#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "../Grid.h"
#include "../Array2D.h" //TODO Cmake里解决这个问题
#include "../Array1D.h" //TODO Cmake里解决这个问题

template<int dim>
class Multigrid
{
private:
    Array2D<GridNode>        *grid_node_arr_ptr_2d;      // 2D grid
    Array1D<GridNode>        *grid_node_arr_ptr_1d;      // 1D grid
    GridBdryType              bdryType;      
    double                    grid_length;    
    double                    inv_grid_length;
    unsigned                  grid_size;                 // pow(2,n)
    int                       grid_node_num;             


public:
    Multigrid() { }
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;

    void grid_initialization();
    void grid_solve();
    void grid_output();
};

template<int dim>
Multigrid<dim>::Multigrid(const Json::Value &root)
{
    // 写入属性

    // 初始化
}