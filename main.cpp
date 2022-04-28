#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <jsoncpp/json/json.h>
#include <eigen3/Eigen/Dense>
#include "Multigrid/Multigrid.h"
#include "Multigrid/CycleSolver.h"
#include "Multigrid/TransferOperator.h"


// test tool
#include <unistd.h>

using namespace std;
bool read_json(string file, Json::Value &root);
void generate_test_fun(std::function<double(double,double)> testfuns[3],
                       std::function<double(double,double)> d_testfuns[3][3]);
void generate_test_fun_1d(std::function<double(double)> &, std::function<double(double)> &);
int run_multigrid();


int run_test()
{   
    // 文件内容全部读取到root变量中
    Json::Value root;
    if(!read_json("./Multigrid/input.json", root))   // TODO: 如何实现边界猜测的读取？
    {
        if(!read_json("../Multigrid/input.json", root))
        {   
            std::cout << "读取json错误" << std::endl;
            return 1;
        }
    }
        
    // 根据json生成网格
    Multigrid<1> grid1(root);  
    Multigrid<2> grid2(root);
    
    // 根据json生成求解器
    CycleSolver<1> solver1(root);
    CycleSolver<2> solver2(root);
    solver2.test();

     // 生成函数
    std::function<double(double)> f1;
    std::function<double(double)> ddf1;
    std::function<double(double,double)> f2[3];   // 暂时没用
    std::function<double(double,double)> df2[3][3];
    generate_test_fun(f2, df2);
    generate_test_fun_1d(f1,ddf1); 
    
    // 网格求解
    //solver1.solve(grid1,&ddf1);
    solver2.solve(grid2,&(df2[0][2]));

    // 读取输出路径
    std::string output_path;
    try{
        output_path = root["output_path"].asString();
    }catch(...){
        std::cerr << "Invalid parameter output_path." << std::endl;
        exit(1);
    }

    // 输出
    //grid1.grid_output(output_path);
    grid2.grid_output(output_path);

    return 0;
}



int main()
{
    run_test();
}


int run_multigrid()
{
    // 定义网格维度
    constexpr int dim = 1;

    
    // 文件内容全部读取到root变量中
    Json::Value root;
    if(!read_json("./Multigrid/input.json", root))   // TODO: 如何实现边界猜测的读取？
        return 1;
    
    
    // 根据json生成网格
    Multigrid<dim> grid(root);  // TODO: 如何摆脱dim的硬编码？

    // 生成求解器 
    CycleSolver<dim> solver(root);


    // 读取函数序号
    int fun_idx;
    try{ 
        fun_idx = root["func_idx"].asInt();
        if(fun_idx < 0 || fun_idx >= 3)
            throw 1;
    }
    catch(...){
        std::cerr << "Invalid parameter func_idx." << std::endl;
        exit(1);
    }

    // 生成函数
    std::function<double(double)> f1;
    std::function<double(double)> ddf1;
    std::function<double(double,double)> f2[3];   // 暂时没用
    std::function<double(double,double)> df2[3][3];
    generate_test_fun(f2, df2);
    generate_test_fun_1d(f1,ddf1); 
    
    // 网格求解
    if(dim == 1)
        solver.solve(grid,&ddf1);
    else if (dim == 2)
        solver.solve(grid, &df2[fun_idx][2]);
    
    

    std::string output_path;
    try{
        output_path = root["output_path"].asString();
    }catch(...){
        std::cerr << "Invalid parameter output_path." << std::endl;
        exit(1);
    }

    // 输出
    grid.grid_output(output_path);

    return 0;
}

bool read_json(string file, Json::Value &root)
{
    // 打开文件
    ifstream input;
    input.open(file);
    if (!input.is_open())
    {
        cerr << "Error: file is not opened" << endl;
    }

    // 创建Json解析器
    Json::Reader reader;

    // 读取数据，记录在root中
    if (!reader.parse(input, root))
    {
        cout << "reader parse error: " << strerror(errno) << endl;
        return false;
    }
    return true;
}

inline double pow2(double d) { return d * d; }

#define PI 3.14159265358979323846
void generate_test_fun(std::function<double(double,double)> testfuns[3],
                       std::function<double(double,double)> d_testfuns[3][3])
{
    testfuns[0] = [](double x, double y){return exp(y  + sin(x));};
    testfuns[1] = [](double x, double y){return pow2(x - 0.5) + pow2(y - 0.5);};
    testfuns[2] = [](double x, double y){return sin(2 * PI * x) * sin(2 * PI * y);};

    d_testfuns[0][0] = [](double x, double y){return exp(y  + sin(x)) * cos(x);};
    d_testfuns[0][1] = [](double x, double y){return exp(y  + sin(x));};
    d_testfuns[0][2] = [](double x, double y){return exp(y  + sin(x)) * (sin(x) - pow2(cos(x)) - 1);};
    
    d_testfuns[1][0] = [](double x, double y){return 2 * x - 1;};
    d_testfuns[1][1] = [](double x, double y){return 2 * y - 1;};
    d_testfuns[1][2] = [](double x, double y){return -4;};
    
    d_testfuns[2][0] = [](double x, double y){return 2 * PI * cos(2 * PI * x) * sin(2 * PI * y);};
    d_testfuns[2][1] = [](double x, double y){return 2 * PI * sin(2 * PI * x) * cos(2 * PI * y);};
    d_testfuns[2][2] = [](double x, double y){return 4 * PI * PI * sin(2 * PI * x) * sin(2 * PI * y);};
}

void generate_test_fun_1d(std::function<double(double)> &testfuns,
                          std::function<double(double)> &neg_laplacian)
{
    testfuns = [](double x){return exp(sin(x));};
    neg_laplacian = [](double x){return exp(sin(x))*(sin(x)-pow2(cos(x)));};
    //d_testfuns[1] = [](double x){return exp(sin(x)) * (cos(x) * cos(x) - sin(x));};
}
