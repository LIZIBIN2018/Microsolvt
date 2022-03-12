#include <iostream>
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>
#include <eigen3/Eigen/Dense>
#include "Grid.h"
#include "Array2D.h"
using namespace std;


bool read_json(string file, Json::Value &root);
void trans_output(Grid &grid);


void testfun() // TODELETE
{
    Array2D<int> arr(2,3);
    int a = 0;
    arr.at(0,1) = 2;
    cout << arr.at(a,1) << endl;
    
    //constexpr int b = arr.at(a,0);
}

int main()
{
    testfun();
    
    // 文件内容全部读取到root变量中
    Json::Value root;
    if(!read_json("input.json", root)) 
        return 1;

    // 生成网格
    Grid grid(root); 

    // 求解线性方程组 TODO（在把网格上的未知数向量化时，我们按字典序排列）

    // 转换输出格式 TODO
    
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

void trans_output(Grid &grid) //TODO 
{

}