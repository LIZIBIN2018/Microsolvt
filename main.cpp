#include <iostream>
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>
#include <eigen3/Eigen/Dense>
using namespace std;

bool read_json(string file, Json::Value &root);


int main()
{
    // 文件内容全部读取到root变量中
    Json::Value root;
    if(!read_json("input.json", root)) 
        return 1;

    // 解释root中信息

    // 生成网格

    // 
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
