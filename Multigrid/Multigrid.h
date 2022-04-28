#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "TransferOperator.h"
#include <eigen3/Eigen/Eigen>


template <int dim>
class Multigrid
{
public:
    Eigen::MatrixXd      data;
    double               grid_length;       // h
    unsigned             grid_size;         // n-1
    unsigned long        grid_node_num;

public:
    Multigrid():grid_length(0),grid_size(0) { }
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;

    void grid_output(std::string path);
};


// implementation 
template <int dim>
Multigrid<dim>::Multigrid(const Json::Value &root) //DONE
{
    // 写入属性
    grid_size = root["grid_n"].asInt() - 1; // n - 1
    if(grid_size <= 0 && (grid_size & (grid_size + 1)) != 0)
        throw "Invalid grid_n";
    grid_length = 1.0 / (grid_size + 1);

    // 生成网格
    grid_node_num = pow(grid_size,dim);
    data = Eigen::MatrixXd::Zero(grid_node_num, 1);
}

template <int dim>
void Multigrid<dim>::grid_output(std::string path)  //DONE
{  
    if(grid_size == 0)
        return;

    path +=std::string("result") + std::to_string(dim) + std::string("D_") + 
           std::to_string(grid_size+1) + std::string(".txt");
    std::ofstream f(path);
    f.precision(17);

    int index = 0;
    int current_col = 0;
    if(dim == 2)
    {
        for(int i = 0; i < grid_size + 2; i++)
            f << 0 << ' ';
        f << std::endl;
    }
    while(index < grid_node_num)
    {
        if(current_col == 0) 
            f << 0 << ' ';
        f << data(index++) << ' ';
        if(++current_col == grid_size)
        {
            f << 0 << std::endl;
            current_col = 0;
        } 
    }
    if(dim == 2)
    {
        for(int i = 0; i < grid_size + 2; i++)
            f << 0 << ' ';
        f << std::endl;
    }
    f.close();
    std::cout << "Job completes with results written to " << path << std::endl;
}
