#include "Grid.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

float distance(float x1, float x2, float y1, float y2){
    x1 -= x2;
    y1 -= y2;
    return std::sqrt(x1 * x1 + y1 * y1);
}

int clamp(int i, int min, int max){
    return std::max(std::min(i, max), min);
}

Grid::Grid(const Json::Value &root)
: grid_node_arr_ptr(nullptr)
{
    // 计算网格的大小，步长
    try{
        grid_length = root["grid_h"].asDouble();
        if(grid_length <= 0)
            throw 1;
        grid_size   = floor(1.0 / grid_length) + 1;
        grid_length = 1.0 / (grid_size-1); //重新计算，保证等分
        grid_node_num = grid_size * grid_size;
    }
    catch(...){
        std::cerr << "Invalid parameter grid_h." << std::endl;
        exit(1);
    }

    // 读取圆属性
    if(root["circ_r"].isNull() && root["circ_c"].isNull()){
        circ_center = std::make_pair(0.0, 0.0);
        circ_radius = 0.0;
    }
    else{
        try{
            circ_radius = root["circ_r"].asDouble();
            if(circ_radius < 0)
                throw 1;
            circ_center = std::move(
                std::pair<double, double>{
                    root["circ_c"][0].asDouble(),
                    root["circ_c"][1].asDouble() });
        }
        catch(...){
            std::cerr << "Invalid parameter circ_r or circ_c." << std::endl;
            exit(1);
        }
        
        // 区域是否连通？
        int nearest_x = clamp((int)std::roundf(circ_center.first / grid_length), 0, grid_size - 1),
            nearest_y = clamp((int)std::roundf(circ_center.second / grid_length), 0, grid_size - 1);
        if(((distance(circ_center.first, circ_center.second, 0,               nearest_y      ) <= circ_radius) +
            (distance(circ_center.first, circ_center.second, grid_length - 1, nearest_y      ) <= circ_radius) +
            (distance(circ_center.first, circ_center.second, nearest_x,       0              ) <= circ_radius) +
            (distance(circ_center.first, circ_center.second, nearest_x,       grid_length - 1) <= circ_radius) -
            (distance(circ_center.first, circ_center.second, 0,               0              ) <= circ_radius) -
            (distance(circ_center.first, circ_center.second, grid_length - 1, 0              ) <= circ_radius) -
            (distance(circ_center.first, circ_center.second, 0,               grid_length - 1) <= circ_radius) -
            (distance(circ_center.first, circ_center.second, grid_length - 1, grid_length - 1) <= circ_radius)) >= 2)
        {
            std::cerr << "Domain not connected." << std::endl;
            exit(1);
        }
    }

    // 设置网格的边界条件
    try{
        auto bdry_type_str = root["bdry_type"].asString();
        if(bdry_type_str == std::string("Dirichlet"))
            this->bdryType = GridBdryType::Dirichlet;
        else if(bdry_type_str == std::string("Neumann"))
            this->bdryType = GridBdryType::Neumann;
        else if(bdry_type_str == std::string("DirichletNeumann"))
            this->bdryType = GridBdryType::DirichletNeumann;
        else if(bdry_type_str == std::string("NeumannDirichlet"))
            this->bdryType = GridBdryType::NeumannDirichlet;
        else
            throw 1;
    }
    catch(...){
        std::cerr << "Invalid parameter bdry_type." << std::endl;
        exit(1);
    }
    
    // 生成网格
    grid_node_arr_ptr = new Array2D<GridNode>(grid_size,grid_size);

    // 网格初始化：类型设置网格点的类型，根据网格边界条件在边界点上赋初始值。 
    grid_initialization(); //移交工作

    if(grid_node_num < 4){
        std::cerr << "Point Number no more than 4." << std::endl;
        exit(1);
    }
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

    // 设置圆的内部和边界  
    //考虑从圆的底部开始中心开花，循环到顶部，只算圆的范围多一圈，节省计算量
    int circ_bottom   = floor ((circ_center.second - circ_radius) / grid_length) ;
    int circ_top      = ceil((circ_center.second + circ_radius) / grid_length);
    int circ_center_x = floor(circ_center.first / grid_length);
    for(int y = circ_bottom ; y <= circ_top; y++)
    {
        int x = 0;
        bool this_row_not_finished = true;
        while(this_row_not_finished)
        {
            // left point on grid  (circ_center_x - x,y)
            // right point on grid (circ_center_x + x,y)
           
            double r_current = std::sqrt(
                std::pow( (circ_center_x + x)*grid_length - circ_center.first , 2) 
              + std::pow( y*grid_length                   - circ_center.second, 2));
            
            if(r_current < circ_radius)
            {
                at(circ_center_x + x,y).type = GridNodeType::exterier;
                at(circ_center_x - x,y).type = GridNodeType::exterier;
                grid_node_num -= (1 + (x!=0));
            }
            else if(r_current < circ_radius + grid_length)
            {
                at(circ_center_x + x,y).type = GridNodeType::boundary;
                at(circ_center_x - x,y).type = GridNodeType::boundary;
            }
            else
                this_row_not_finished = false;
            ++x;
        }
    }
    generate_bit_graph(); // TODELETE
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
    for(int i = grid_size - 1;i>= 0;--i)
    {
        for(int j=0;j<grid_size;++j)
        {
            GridNodeType temp = at(j,i).type;
            f<<(temp==GridNodeType::boundary)+0*(temp == GridNodeType::exterier)<<' ';
        }
        f << std::endl;
    } 
    std::cout<< grid_node_num << std::endl; 
}


void Grid::grid_solve(std::function<double(double,double)>f,
                      std::function<double(double,double)>df_dx,
                      std::function<double(double,double)>df_dy)
{
    // 创建矩阵，求解Ax = b
    Eigen::MatrixXd A;
    Eigen::MatrixXd x;
    Eigen::MatrixXd b;

    A.resize(grid_node_num, grid_node_num);
    x.resize(grid_node_num, 1);
    b = x; // lazy boy.....

    // 从grid中写入矩阵 TODO

    // 求解 TODO

    // 从矩阵写入grid
    int       place        = 0; // stand for the current node of the grid with disk
    GridNode *current_node = &at(0,0);
    for(int j = 0; j<grid_node_num; j++)
    {
        current_node->val = x(j);
        current_node = next_node(place);
    }
}

void Grid::grid_output()
{
    std::fstream f("grid_val_output.txt"); //TODO
}

GridNode *Grid::next_node(int &x, int &y)
{
    GridNode *rst = nullptr;
    if(x<0 || x>=grid_size || y<0 || y>=grid_size)
        return rst;  // 一切尽在不言中
    while(rst == nullptr)
    {
        x = (x + 1) % grid_size;
        y = ((x==0) + y) % grid_size;
        if(y == 0)
            return rst;
        if(at(x,y).type != GridNodeType::exterier)
            rst = &at(x,y);
    }
    return rst;
}

GridNode *Grid::next_node(int &idx)
{
    int x = idx % grid_size;
    int y = (idx - x)/grid_size;
    GridNode *rst = next_node(x,y);
    idx = x + grid_size * y;
    return rst;
}