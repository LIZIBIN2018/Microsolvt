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
        if(((at(0,               nearest_y      ).type == GridNodeType::exterior) +
            (at(grid_length - 1, nearest_y      ).type == GridNodeType::exterior) +
            (at(nearest_x,       0              ).type == GridNodeType::exterior) +
            (at(nearest_x,       grid_length - 1).type == GridNodeType::exterior) -
            (at(0,               0              ).type == GridNodeType::exterior) -
            (at(grid_length - 1, 0              ).type == GridNodeType::exterior) -
            (at(0,               grid_length - 1).type == GridNodeType::exterior) -
            (at(grid_length - 1, grid_length - 1).type == GridNodeType::exterior)) >= 2)
        {
            std::cerr << "Domain not connected." << std::endl;
            exit(1);
        }
    }
}

void Grid::grid_initialization(){
    // 边界标记
    grid_node_num = grid_size * grid_size;
    for (size_t i = 0; i < grid_size; i++)
    {
        at(i,0).type           = GridNodeType::squareBoundary;
        at(i,grid_size-1).type = GridNodeType::squareBoundary;
        at(0,i).type           = GridNodeType::squareBoundary;
        at(grid_size-1,i).type = GridNodeType::squareBoundary;
    }
    if(circ_exist){
        int min_x = clamp((int)((circ_center.first - circ_radius)           * inv_grid_length), 0, grid_size - 1),
            max_x = clamp((int)std::ceil((circ_center.first + circ_radius)  * inv_grid_length), 0, grid_size - 1),
            min_y = clamp((int)((circ_center.second - circ_radius)          * inv_grid_length), 0, grid_size - 1),
            max_y = clamp((int)std::ceil((circ_center.second + circ_radius) * inv_grid_length), 0, grid_size - 1);
        // std::cout << min_x << ' ' << min_y << ',' << max_x << ' ' << max_y << std::endl;
        for (size_t i = min_x; i <= max_x; i++){
            for (size_t j = min_y; j <= max_y; j++){
                if(distance(i * grid_length, j * grid_length,
                            circ_center.first, circ_center.second) <= circ_radius){
                    at(i, j).type = GridNodeType::exterior;
                    if(i != 0 && at(i - 1, j).type != GridNodeType::exterior)
                        at(i - 1, j).type = GridNodeType::circleBoundary;
                    if(j != 0 && at(i, j - 1).type != GridNodeType::exterior)
                        at(i, j - 1).type = GridNodeType::circleBoundary;
                    if(i != grid_size - 1 && at(i + 1, j).type != GridNodeType::exterior)
                        at(i + 1, j).type = GridNodeType::circleBoundary;
                    if(j != grid_size - 1 && at(i, j + 1).type != GridNodeType::exterior)
                        at(i, j + 1).type = GridNodeType::circleBoundary;
                    --grid_node_num;
                }
            }
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
    for(int i = grid_size - 1;i>= 0;--i)
    {
        for(int j=0;j<grid_size;++j)
        {
            f<<(int)at(j,i).type<<' ';
        }
        f << std::endl;
    } 
    // std::cout<< grid_node_num << std::endl; 
}


void Grid::grid_solve(std::function<double(double,double)>f,
                      std::function<double(double,double)>df_dx,
                      std::function<double(double,double)>df_dy,
                      std::function<double(double,double)>neg_laplacian_f)
{
    // 创建矩阵，求解Ax = b
    Eigen::MatrixXd A;
    Eigen::MatrixXd x;
    Eigen::MatrixXd b;

    A.resize(grid_node_num, grid_node_num);
    x.resize(grid_node_num, 1);
    b = x; // lazy boy.....

    // 从grid中写入矩阵
    Array2D<int> matrix_idx(grid_size, grid_size);
    for (size_t j = 0, node_idx = 0; j < grid_size; j++)
    {
        for (size_t i = 0; i < grid_size; i++)
        {
            matrix_idx.at(i, j) = at(i, j).type != GridNodeType::exterior ? node_idx++ : -1;
        }
    }
    double invsqr_h = inv_grid_length * inv_grid_length;
    for (size_t j = 0; j < grid_size; j++)
    {
        for (size_t i = 0; i < grid_size; i++)
        {
            int idx = matrix_idx.at(i, j);
            if(at(i, j).type == GridNodeType::interior){
                A(idx, idx) = 4 * invsqr_h;
                A(idx, matrix_idx.at(i - 1, j)) = -invsqr_h;
                A(idx, matrix_idx.at(i + 1, j)) = -invsqr_h;
                A(idx, matrix_idx.at(i, j - 1)) = -invsqr_h;
                A(idx, matrix_idx.at(i, j + 1)) = -invsqr_h;
                b(idx) = neg_laplacian_f(i * grid_length, j * grid_length);
            }
            else if(at(i, j).type == GridNodeType::squareBoundary){
                if(((int)bdryType & 2) == 2) // Neumann
                {
                    b(idx) = f(i * grid_length, j * grid_length) * grid_length * 0.5f;
                    if(j == 0){
                        if(j + 1 >= grid_size || matrix_idx.at(i, j + 1) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = -inv_grid_length;
                        A(idx, matrix_idx.at(i, j + 1)) = inv_grid_length;
                        b(idx) += df_dy(i * grid_length, j * grid_length);
                    }
                    else if(j == grid_size - 1){
                        if(j - 1 < 0 || matrix_idx.at(i, j - 1) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = -inv_grid_length;
                        A(idx, matrix_idx.at(i, j - 1)) = inv_grid_length;
                        b(idx) += -df_dy(i * grid_length, j * grid_length);
                    }
                    else if(i == 0){
                        if(i + 1 >= grid_size || matrix_idx.at(i + 1, j) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = -inv_grid_length;
                        A(idx, matrix_idx.at(i + 1, j)) = inv_grid_length;
                        b(idx) += df_dx(i * grid_length, j * grid_length);
                    }
                    else{
                        if(i - 1 < 0 || matrix_idx.at(i - 1, j) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = -inv_grid_length;
                        A(idx, matrix_idx.at(i - 1, j)) = inv_grid_length;
                        b(idx) += -df_dx(i * grid_length, j * grid_length);
                    }
                }
                else // Dirichlet
                {
                    A(idx, idx) = 1;
                    b(idx) = f(i * grid_length, j * grid_length);
                }
            }
            else if(at(i, j).type == GridNodeType::circleBoundary){
                if(((int)bdryType & 1) == 1) // Neumann
                {
                    double distance_to_center = distance(i * grid_length, j * grid_length, circ_center.first, circ_center.second);
                    double normalX = (i * grid_length - circ_center.first) / distance_to_center,
                            normalY = (j * grid_length - circ_center.second) / distance_to_center;
                    b(idx) = f(i * grid_length, j * grid_length) * grid_length * 0.5f;

                    if(normalY >= 0){
                        if(j + 1 >= grid_size || matrix_idx.at(i, j + 1) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = -inv_grid_length * normalY;
                        A(idx, matrix_idx.at(i, j + 1)) = inv_grid_length * normalY;
                        b(idx) += df_dy(i * grid_length, j * grid_length) * normalY;
                    }
                    else{
                        if(j - 1 < 0 || matrix_idx.at(i, j - 1) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) = inv_grid_length * normalY;
                        A(idx, matrix_idx.at(i, j - 1)) = -inv_grid_length * normalY;
                        b(idx) += df_dy(i * grid_length, j * grid_length) * normalY;
                    }

                    if(normalX >= 0){
                        if(i + 1 >= grid_size || matrix_idx.at(i + 1, j) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) += -inv_grid_length * normalX;
                        A(idx, matrix_idx.at(i + 1, j)) = inv_grid_length * normalX;
                        b(idx) += df_dx(i * grid_length, j * grid_length) * normalX;
                    }
                    else{
                        if(i - 1 < 0 || matrix_idx.at(i - 1, j) == -1){
                            std::cerr << "Invalid mesh." << std::endl;
                            exit(1);
                        }
                        A(idx, idx) += inv_grid_length * normalX;
                        A(idx, matrix_idx.at(i - 1, j)) = -inv_grid_length * normalX;
                        b(idx) += df_dx(i * grid_length, j * grid_length) * normalX;
                    }
                }
                else // Dirichlet
                {
                    A(idx, idx) = 1;
                    b(idx) = f(i * grid_length, j * grid_length);
                }
            }
        }
    }

    // 求解 TODO
    x = A.lu().solve(b);

    // 从矩阵写入grid
    /*int       place        = 0; // stand for the current node of the grid with disk
    GridNode *current_node = &at(0,0);
    for(int j = 0; j<grid_node_num; j++)
    {
        current_node->val = x(j);
        current_node = next_node(place);
    }*/
    for (size_t j = 0, node_idx = 0; j < grid_size; j++)
    {
        for (size_t i = 0; i < grid_size; i++)
        {
            if(matrix_idx.at(i, j) != -1){
                at(i, j).val = x(matrix_idx.at(i, j));
            }
        }
    }
}

void Grid::grid_output(std::string path)
{
    std::ofstream f(path);
    f.precision(17);
    for (size_t j = 0, node_idx = 0; j < grid_size; j++)
    {
        for (size_t i = 0; i < grid_size; i++)
        {
            if(at(i, j).type != GridNodeType::exterior){
                f << at(i, j).val << ' ';
            }
            else{
                f << "* ";
            }
        }
        f << std::endl;
    }
    f.close();
    std::cout << "Job completes with results written to " << path << std::endl;
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
        if(at(x,y).type != GridNodeType::exterior)
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