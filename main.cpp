#include <iostream>
#include <jsoncpp/json/json.h>
#include <eigen3/Eigen/Dense>
using namespace std;


int main()
{
    Json::Value root;    
    root = 1;
    cout << root << endl;

    Eigen::Matrix2f mtx;
    mtx << 1, 2, 3, 4;
    cout << mtx << endl;

    return 0;
}
