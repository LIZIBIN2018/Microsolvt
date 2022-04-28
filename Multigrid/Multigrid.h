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
    GridBdryType         bdryType;
    double               grid_length;     // h
    unsigned             grid_size;     // n-1
    

public:
    Multigrid() {}
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;

    void grid_initialization();
    void grid_output(std::string path);
};


// implementation 
template <int dim>
Multigrid<dim>::Multigrid(const Json::Value &root)
{
    // 写入属性
    grid_size = root["grid_n"].asInt() - 1; // n - 1
    if(grid_size <= 0 && (grid_size & (grid_size + 1)) != 0)
        throw "Invalid grid_n";
    grid_length = 1.0 / (grid_size + 1);

    // 生成网格
    if(dim == 1)
        data.resize(grid_size, 1); //TOCHECK
    if(dim == 2)
        data.resize(grid_size * grid_size, 1); //TOCHECK
}

template <int dim>
void Multigrid<dim>::grid_output(std::string path)
{
    std::ofstream f(path);
    f.precision(17);
    if (dim == 1)
    {
        f << 0 << std::endl;
        for (int i = 0; i < grid_size; i++)
        {
            f << data(i) << std::endl;
        }
        f << 0 << std::endl;
    }
    if (dim == 2)
    {
        for (int i = -1; i < grid_size + 1; ++i)
        {
            for (int j = -1; j < grid_size + 1; ++j)
            {
                if (i == -1 || i == grid_size || j == -1 || j == grid_size)
                {
                    f << 0 << ' ';
                }
                else
                {
                    f << data(j + grid_size * i) << ' ';
                }
            }
            f << std::endl;
        }
    }

    f.close();
    std::cout << "Job completes with results written to " << path << std::endl;
}
