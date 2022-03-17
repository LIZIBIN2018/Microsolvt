#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "Array2D.h"

// 用于指示点在内部或者外部或者是在边界上
enum class GridNodeType
{
    interior,
    exterior,
    squareBoundary,
    circleBoundary
};

// 用于指示网格的边界条件类型
enum class GridBdryType
{
    Dirichlet,
    DirichletNeumann, // Dirichlet on square border, Neumann on circle border
    NeumannDirichlet, // Neumann on square border, Dirichlet on circle border
    Neumann,
};

// 网格节点
struct GridNode
{
    double val;
    GridNodeType type = GridNodeType::interior;   
};

class Grid
{
    Array2D<GridNode>        *grid_node_arr_ptr;      // 2D grid
    GridBdryType              bdryType;      
    double                    grid_length;    
    double                    inv_grid_length;
    unsigned                  grid_size;
    int                       grid_node_num;

    bool                      circ_exist = true;      // 默认有圆
    double                    circ_radius;
    std::pair<double, double> circ_center;
    
public:  //主要函数
    // ctor & dtor
    Grid() = default;
    Grid(const Json::Value &root);
    Grid(const Grid &) = delete;        // can not copy the grid
    //Grid(Grid && grid) { *this = std::move(grid); } TODO
    ~Grid();
    
    // date manipulation 
    constexpr GridNode &at(const size_t x,const size_t y) 
    {
        return grid_node_arr_ptr->at(x,y);
    }

    // 区域是否符合要求
    bool is_regular();
    
    // 初始化：进行内外点标记，以及尽量直接在边界点写入已知的数值。
    void grid_initialization();

    // 求解
    void grid_solve(std::function<double(double,double)>f,
                    std::function<double(double,double)>df_dx,
                    std::function<double(double,double)>df_dy,
                    std::function<double(double,double)>neg_laplacian_f);
    
    // 输出
    void grid_output(std::string path);
public: // 辅助函数
    void safe_delete(); 
    GridNode *next_node(int &x, int &y); //寻找指定位置的下一个有效节点
    GridNode *next_node(int &idx);
public:
    // TODELETE debug tools 
    void show()
    {
        std::cout<< "grid length = " << grid_length << '\n';
        std::cout<< "grid size = "   << grid_size   << '\n';
        std::cout<< "circ_radius = " << circ_radius << '\n';
        //std::cout<< "circ_center = " << grid_length << '\n';
    }
    void generate_bit_graph(int flag = 0);

};
