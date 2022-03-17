#include <iostream>
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>
#include <eigen3/Eigen/Dense>
#include "Grid.h"
#include "Array2D.h"
#include <functional>
using namespace std;


bool read_json(string file, Json::Value &root);
void generate_test_fun(std::function<double(double,double)> testfuns[3],
                       std::function<double(double,double)> d_testfuns[3][3]);

int main()
{
    // 生成测试函数
    std::function<double(double,double)> testfuns[3];
    std::function<double(double,double)> d_testfuns[3][3];
    generate_test_fun(testfuns,d_testfuns);

    // 文件内容全部读取到root变量中
    Json::Value root;
    if(!read_json("input.json", root)) 
        return 1;
    
    // 生成网格
    Grid grid(root); 

    // 求解线性方程组,把结果直接写在网格里（在把网格上的未知数向量化时，我们按字典序排列）
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
    grid.grid_solve(testfuns[fun_idx],
                    d_testfuns[fun_idx][0],
                    d_testfuns[fun_idx][1],
                    d_testfuns[fun_idx][2]);

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
    if(!input.is_open()) { cerr<<"Error: file is not opened"<<endl; }

    // 创建Json解析器
    Json::Reader reader;
    
    // 读取数据，记录在root中
    if(!reader.parse(input,root))
    {
        cout <<"reader parse error: "<<strerror(errno)<<endl;
        return false;
    }
    return true;
}

inline double pow2(double d) { return d * d; }

void generate_test_fun(std::function<double(double,double)> testfuns[3],
                       std::function<double(double,double)> d_testfuns[3][3])
{
    testfuns[0] = [](double x, double y){return exp(y  + sin(x));};
    testfuns[1] = [](double x, double y){return exp(y) * cos(x);};
    testfuns[2] = [](double x, double y){return sin(x) * cos(y) ;};

    d_testfuns[0][0] = [](double x, double y){return  exp(y  + sin(x))*cos(x);};
    d_testfuns[0][1] = [](double x, double y){return  exp(y  + sin(x));};
    d_testfuns[0][2] = [](double x, double y){return  exp(y  + sin(x)) * (sin(x) - pow2(cos(x)) - 1);};
    
    d_testfuns[1][0] = [](double x, double y){return -exp(y) * sin(x);};
    d_testfuns[1][1] = [](double x, double y){return  exp(y) * cos(x);};
    d_testfuns[1][2] = [](double x, double y){return  0;};
    
    d_testfuns[2][0] = [](double x, double y){return  cos(x) * cos(y);};
    d_testfuns[2][1] = [](double x, double y){return -sin(x) * sin(y);};
    d_testfuns[2][2] = [](double x, double y){return  2 * sin(x) * cos(y);};
}