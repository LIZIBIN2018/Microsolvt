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
public:
    Eigen::MatrixXd      data;
    GridBdryType         bdryType;
    double               grid_length;     // h
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
};


// implementation 
template <int dim>
Multigrid<dim>::Multigrid(const Json::Value &root)
{
    // 写入属性
    grid_size = root["grid_n"].asInt() - 1; // n - 1
    if(grid_size <= 0 && (grid_size & (grid_size + 1)) != 0)
        throw std::exception("Invalid grid_n");
    grid_length = 1.0 / (grid_size + 1);
    
    grid_node_num = grid_size;
    if (dim == 2)
        grid_node_num *= grid_size;

    // 生成网格
    if(dim == 1)
        data.resize(grid_size,1); //TOCHECK
    if(dim == 2)
        data.resize(grid_size,grid_size); //TOCHECK
}

template <int dim>
void Multigrid<dim>::grid_output(std::string path)
{
    std::ofstream f(path);
    f.precision(17);
    if(dim == 1){
        f << 0 << ' ';
        for(int i = 0;i < grid_size; i++){
            f << data(i) << ' ';
        }
        f<< 0 << ' ';
    }
    if(dim == 2){
        for(int i = 0;i<grid_size + 2;++i){
            for(int j = 0;j<grid_size + 2;++j){
                if(i == 0 || i == grid_size+1 || j == 0 || j == grid_size+1){
                    f << 0 << ' ';
                }
                else{
                    f << data(j + grid_size*i) << ' ';
                }
            }
            f << endl;
        }
    }

    f.close();
    std::cout << "Job completes with results written to " << path << std::endl;
}
