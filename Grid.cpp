#include "Grid.h"
#include <cmath>
#include <iostream>
#include <fstream>

Grid::Grid(const Json::Value &root)
: grid_node_arr_ptr(nullptr)
{
    // 计算网格的大小，步长
    grid_length = root["grid_h"].asDouble();
    grid_size   = floor(1.0 / grid_length);
    grid_length = 1.0 / grid_size; //重新计算，保证等分

    // 读取圆属性
    circ_radius = root["circ_r"].asDouble();
    circ_center = std::move(
        std::pair<double, double>{
            root["circ_c"][0].asDouble(),
            root["circ_c"][1].asDouble() });

    // 是否符合要求？（圆至少包含了四个网格点？区域是否连通？）(是不是有圆, 不包含四点我们就认为没有)
    if(! is_regular())  return;

    // 设置网格的边界条件
    if(root["bdry_type"].asString() == std::string("Dirichlet"))
        this->bdryType = GridBdryType::Dirichlet;
    else if(root["bdry_type"].asString() == std::string("Neumann"))
        this->bdryType = GridBdryType::Neumann;
    else
        this->bdryType = GridBdryType::Mixed;

    // 生成网格
    grid_node_arr_ptr = new Array2D<GridNode>(grid_size,grid_size);

    // 网格初始化：类型设置网格点的类型，根据网格边界条件在边界点上赋初始值。 
    grid_initialization(); //移交工作
}



bool Grid::is_regular()
{
    // TODO 
    // 我们只考虑圆是子集的情形，不是子集就exit
    return true;
}

void Grid::grid_initialization()  //TOCHECK
{
    // 设置正方形边界
    for(int idx = 0;idx < grid_size;++idx)
    {
        at(idx,0).type           = GridNodeType::boundary;
        at(idx,grid_size-1).type = GridNodeType::boundary;
        at(0,idx).type           = GridNodeType::boundary;
        at(grid_size-1,idx).type = GridNodeType::boundary;
    }

    // 设置圆的内部和边界 //TODO 
    //考虑从圆的底部开始中心开花，循环到顶部，只算圆的范围多一圈，节省计算量
    int circ_bottom   = ceil ((circ_center.second - circ_radius) / grid_length);
    int circ_top      = floor((circ_center.second + circ_radius) / grid_length) + 1;
    int circ_center_x = floor(circ_center.first / grid_length);
    for(int y = circ_bottom; y<circ_top; y++)
    {
        int x = 0;
        bool this_row_not_finished = true;
        while(this_row_not_finished)
        {
            // left point on grid  (circ_center_x - x,y)
            // right point on grid (circ_center_x + x,y)
            double r_left = std::sqrt(
                std::pow( (circ_center_x - x)*grid_length - circ_center.first , 2) 
              + std::pow( y*grid_length                   - circ_center.second, 2));
            double r_right = std::sqrt(
                std::pow( (circ_center_x + x)*grid_length - circ_center.first , 2) 
              + std::pow( y*grid_length                   - circ_center.second, 2));
            if(r_left < circ_radius)
            {
                at(circ_center_x - x,y).type = GridNodeType::exterier;
            }
            else if(r_left < circ_radius + grid_length)
            {
                at(circ_center_x - x,y).type = GridNodeType::boundary;
            }
            if(r_right < circ_radius)
            {
                at(circ_center_x + x,y).type = GridNodeType::exterier;
            }
            else if(r_right < circ_radius + grid_length)
            {
                at(circ_center_x + x,y).type = GridNodeType::boundary;
                this_row_not_finished = false;
            }
            ++x;
        }
    }
    generate_bit_graph();
}

Grid::~Grid() 
{
    safe_delete();
}

void Grid::safe_delete()
{ 
    delete grid_node_arr_ptr; 
    grid_node_arr_ptr = nullptr;
}

void Grid::generate_bit_graph(int flag)/*default flag = 0*/
{
    std::ofstream f("bitgraph.ppm");
    f << "P3\n"<< grid_size<<' '<<grid_size<<"\n255\n";
    for(int i=0;i<grid_size;++i)
    {
        for(int j=0;j<grid_size;++j)
        {
            GridNodeType temp = at(i,j).type;
            f<<(temp==GridNodeType::boundary)<<' ';
        }
        f <<'\n';
    }  
}