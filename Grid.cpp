#include "Grid.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

double distance(double x1, double y1, double x2, double y2){
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
        inv_grid_length = 1.0 / grid_length;
        grid_size   = floor(inv_grid_length) + 1;
        grid_length = 1.0 / (grid_size-1); //重新计算，保证等分
        grid_node_num = grid_size * grid_size;
    }
    catch(...){
        std::cerr << "Invalid parameter grid_h." << std::endl;
        exit(1);
    }

    // 读取圆属性
    if(root["circ_r"].isNull() && root["circ_c"].isNull()){
        circ_exist = false;
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

    // 是否少于4个点？
    if(grid_node_num < 4){
        std::cerr << "Point Number no more than 4." << std::endl;
        exit(1);
    }

    // 区域是否连通？
    if(circ_exist){
        int nearest_x = clamp((int)std::roundf(circ_center.first* inv_grid_length), 0, grid_size - 1),
            nearest_y = clamp((int)std::roundf(circ_center.second* inv_grid_length), 0, grid_size - 1);
        if(((at(0,               nearest_y      ).type == GridNodeType::exterier) +
            (at(grid_length - 1, nearest_y      ).type == GridNodeType::exterier) +
            (at(nearest_x,       0              ).type == GridNodeType::exterier) +
            (at(nearest_x,       grid_length - 1).type == GridNodeType::exterier) -
            (at(0,               0              ).type == GridNodeType::exterier) -
            (at(grid_length - 1, 0              ).type == GridNodeType::exterier) -
            (at(0,               grid_length - 1).type == GridNodeType::exterier) -
            (at(grid_length - 1, grid_length - 1).type == GridNodeType::exterier)) >= 2)
        {
            std::cerr << "Domain not connected." << std::endl;
            exit(1);
        }
    }
}

void Grid::grid_initialization(){
    grid_node_num = grid_size * grid_size;
    for (size_t i = 0; i < grid_size; i++)
    {
        at(i,0).type           = GridNodeType::squareBoundary;
        at(i,grid_size-1).type = GridNodeType::squareBoundary;
        at(0,i).type           = GridNodeType::squareBoundary;
        at(grid_size-1,i).type = GridNodeType::squareBoundary;
    }
    int min_x = clamp((int)((circ_center.first - circ_radius)           * inv_grid_length), 0, grid_size - 1),
        max_x = clamp((int)std::ceil((circ_center.first + circ_radius)  * inv_grid_length), 0, grid_size - 1),
        min_y = clamp((int)((circ_center.second - circ_radius)          * inv_grid_length), 0, grid_size - 1),
        max_y = clamp((int)std::ceil((circ_center.second + circ_radius) * inv_grid_length), 0, grid_size - 1);
    std::cout << min_x << ' ' << min_y << ',' << max_x << ' ' << max_y << std::endl;
    for (size_t i = min_x; i <= max_x; i++){
        for (size_t j = min_y; j <= max_y; j++){
            if(distance(i * grid_length, j * grid_length,
                        circ_center.first, circ_center.second) <= circ_radius){
                at(i, j).type = GridNodeType::exterier;
                if(i != 0 && at(i - 1, j).type != GridNodeType::exterier)
                    at(i - 1, j).type = GridNodeType::circleBoundary;
                if(j != 0 && at(i, j - 1).type != GridNodeType::exterier)
                    at(i, j - 1).type = GridNodeType::circleBoundary;
                if(i != grid_size - 1 && at(i + 1, j).type != GridNodeType::exterier)
                    at(i + 1, j).type = GridNodeType::circleBoundary;
                if(j != grid_size - 1 && at(i, j + 1).type != GridNodeType::exterier)
                    at(i, j + 1).type = GridNodeType::circleBoundary;
                --grid_node_num;
            }
        }
    }
    generate_bit_graph();
}

/*void Grid::grid_initialization()  //TOCHECK
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
    int circ_bottom   = floor ((circ_center.second - circ_radius)* inv_grid_length) ;
    int circ_top      = ceil((circ_center.second + circ_radius)* inv_grid_length);
    int circ_center_x = floor(circ_center.first* inv_grid_length);
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
}*/

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
            f<<(int)at(j,i).type<<' ';
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