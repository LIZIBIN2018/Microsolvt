#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "../Grid.h"
#include "../Array2D.h" //TODO Cmake里解决这个问题
#include "../Array1D.h" //TODO Cmake里解决这个问题
#include "TransferOperator.h"
#include <eigen3/Eigen/Eigen>


template <int dim>
class Multigrid
{
private:
    Eigen::MatrixXd      data;
    GridBdryType         bdryType;
    double               grid_length;     // h
    double               inv_grid_length; // floor(?) = n
    unsigned             grid_size;     // n-1
    int                  grid_node_num;
    int                  max_itr_num; // maximum times for iteration
    double               epsilon;  // tolerance
    

public:
    Multigrid() {}
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;

    void grid_initialization();
    void grid_output(std::string path);

private: //TOOLS
};


// implementation 
template <int dim>
Multigrid<dim>::Multigrid(const Json::Value &root)
{
    // 写入属性
    try
    {
        grid_length = root["grid_h"].asDouble();
        if (grid_length <= 0)
            throw 1;
        inv_grid_length = 1.0 / grid_length;
        grid_size = floor(inv_grid_length) - 1; // n - 1
        grid_length = 1.0 / (grid_size - 1);    // h 重新计算，保证等分
        grid_node_num = grid_size; //TODO
        if (dim == 2)
            grid_node_num *= grid_size;
    }
    catch (...)
    {
        std::cerr << "Invalid parameter grid_h." << std::endl;
        exit(1);
    }


    // 生成网格
    data.resize(grid_size,grid_size); //TOCHECK

}

template <int dim>
void Multigrid<dim>::grid_output(std::string path)
{
    //TODO
}
