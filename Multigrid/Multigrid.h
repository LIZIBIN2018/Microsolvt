#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "../Grid.h"
#include "../Array2D.h" //TODO Cmake里解决这个问题
#include "../Array1D.h" //TODO Cmake里解决这个问题
#include "TransferOperator.h"
#include "CycleSolver.h"  // 代表cycle的类型

enum class RstOptType // Restriction Operator
{
    injection,
    fullWeighting
};

enum class ItpOptType // Interpolation Operator
{
    linear,
    quadratic
};

template <int dim>
class Multigrid
{
private:
    Array2D<GridNode>   *grid_node_arr_ptr_2d; // 2D grid
    Array1D<GridNode>   *grid_node_arr_ptr_1d; // 1D grid
    GridBdryType         bdryType;
    double               grid_length;     // h
    double               inv_grid_length; // floor(?) = n
    unsigned             grid_size;     // n+1
    int                  grid_node_num;
    RstOptType           grid_rst_opt_type;
    ItpOptType           grid_itp_opt_type;
    int                  max_itr_num; // maximum times for iteration
    double               epsilon;  // tolerance

    CycleSolver<dim>    *grid_solver; 
   
public:
    Multigrid() {}
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;
    ~Multigrid() 
    {
        delete grid_node_arr_ptr_1d;  grid_node_arr_ptr_1d = nullptr;
        delete grid_node_arr_ptr_2d;  grid_node_arr_ptr_2d = nullptr;
        delete grid_solver;           grid_solver = nullptr;
    } 

    void grid_initialization();
    // TODO这丑陋的写法
    void grid_solve(std::function<double(double)> f); 
    void grid_solve(std::function<double(double,double)>f,
                      std::function<double(double,double)>df_dx,
                      std::function<double(double,double)>df_dy,
                      std::function<double(double,double)>neg_laplacian_f);

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
        grid_size = floor(inv_grid_length) + 1; // n + 1
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

    // 设置网格的边界条件 //TODELETE 没必要
    try
    {
        auto bdry_type_str = root["bdry_type"].asString();
        if (bdry_type_str == std::string("Dirichlet"))
            this->bdryType = GridBdryType::Dirichlet;
        else if (bdry_type_str == std::string("Neumann"))
            this->bdryType = GridBdryType::Neumann;
        else if (bdry_type_str == std::string("DirichletNeumann"))
            this->bdryType = GridBdryType::DirichletNeumann;
        else if (bdry_type_str == std::string("NeumannDirichlet"))
            this->bdryType = GridBdryType::NeumannDirichlet;
        else
            throw 1;
    }
    catch (...)
    {
        std::cerr << "Invalid parameter bdry_type." << std::endl;
        exit(1);
    }

    // 生成网格
    if (dim == 2)
        grid_node_arr_ptr_2d = new Array2D<GridNode>(grid_size, grid_size);
    if (dim == 1)
        grid_node_arr_ptr_1d = new Array1D<GridNode>(grid_size);

    // 设置求解器
    try
    {
        auto cycle_type = root["cycle"].asString();
        if(cycle_type == "V-cycle")
            grid_solver = new VCycle<dim>(this);
        else if(cycle_type == "FullMultigridVCycle")
            grid_solver = new FullMultigridVCycle<dim>(this);
        else
            throw 1;
    }
    catch(...)
    {
        std::cerr << "Invalid cycle type" << std::endl;
        exit(1);
    }
}


template <int dim> //
void Multigrid<dim>::grid_solve(
    std::function<double(double)> f)
{
    
}

template<int dim>
void Multigrid<dim>::grid_solve(std::function<double(double,double)>f,
                      std::function<double(double,double)>df_dx,
                      std::function<double(double,double)>df_dy,
                      std::function<double(double,double)>neg_laplacian_f)
{

}

template <int dim>
void Multigrid<dim>::grid_output(std::string path)
{
    //TODO
}
