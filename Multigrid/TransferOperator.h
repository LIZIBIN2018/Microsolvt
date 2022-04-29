#pragma once 
#include <iostream>
#include <eigen3/Eigen/Dense>

// 1st layer 
template<int dim>
class RestrictionOperator
{
private:
    virtual double restrict1D(const Eigen::MatrixXd& mat, size_t idx) = 0;
    virtual double restrict2D(const Eigen::MatrixXd& mat, size_t idx, size_t grid_num) = 0;
public:
    Eigen::MatrixXd operator()(const Eigen::MatrixXd &v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        
        size_t halfLength = grid_num >> 1;
        if(dim == 1){
            rst.resize(halfLength, 1);
            for(size_t i = 0; i < halfLength; i++){
                rst(i) = restrict1D(v, i << 1 | 1);  //*2+1
            }
        }
        else if(dim == 2){
            rst.resize(halfLength * halfLength, 1);
            for (size_t row = 0; row < halfLength; row++)
            {
                for (size_t col = 0; col < halfLength; col++)
                {
                    rst(row * halfLength + col) = restrict2D(v, (row << 1 | 1) * grid_num + (col << 1 | 1), grid_num);
                }
            }
        }
        else{
            throw "Invalid dim";
        }
        return rst;
    }
};


template<int dim>
class InterpolationOperator
{
private:
    virtual double interpolate1D(double left, double right) = 0;
    virtual double interpolate2D(double left_up, double right_up, double left_down, double right_down) = 0;
public:
    Eigen::MatrixXd operator()(const Eigen::MatrixXd &v, size_t grid_num)
    {
        Eigen::MatrixXd itp;
        
        size_t doubleLength = grid_num << 1 | 1;
        if(dim == 1){
            itp.resize(doubleLength, 1);
            for (size_t i = 0; i < doubleLength; i++){
                if( i%2==0 )
                    itp(i) = interpolate1D(i == 0 ? v(0) : v((i - 1) >> 1), i == doubleLength - 1 ? v(grid_num - 1) : v(i >> 1));        
                else
                    itp(i) = v(i>>1);
            }
        }
        else if(dim == 2){
            itp.resize(doubleLength * doubleLength, 1);
            for (size_t row = 0; row < doubleLength; row++)
            {
                for (size_t col = 0; col < doubleLength; col++)
                {
                    int index = row * doubleLength + col;
                    double temp = 0;
                    if(row % 2 == 1 && col % 2 == 1)
                    {
                        temp = v((row >> 1) * grid_num + (col >> 1));
                    }
                    else if(row % 2 == 0 && col % 2 == 1)
                    {
                        // row -> row>>1 - 1 , row>>1
                        // col -> col>>1
                        temp = interpolate1D(row==0?0:v(((row>>1)-1)*grid_num + (col>>1)), 
                                                   row==doubleLength-1?0:v((row>>1)*grid_num + (col>>1)));
                    }
                    else if(row % 2 == 1 && col % 2 == 0)
                    {
                        temp = interpolate1D(col==0?0:v((row>>1)*grid_num + (col>>1) - 1), 
                                             col==doubleLength-1?0:v((row>>1)*grid_num + (col>>1)));
                    }
                    else
                    {
                        temp = interpolate2D(
                        row == 0                || col == 0                ? 0 : v(((row>>1) - 1) * grid_num + ((col >> 1) - 1)),
                        row == 0                || col == doubleLength - 1 ? 0 : v(((row>>1) - 1) * grid_num + ((col >> 1))),
                        row == doubleLength - 1 || col == 0                ? 0 : v(((row>>1)) * grid_num + ((col >> 1) - 1)),
                        row == doubleLength - 1 || col == doubleLength - 1 ? 0 : v(((row>>1)) * grid_num + ((col >> 1))));
                    }
                    itp(index) = temp;
                }
            }
        }
        else{
            throw "Invalid dim";
        }
        return itp;
    }
};


// 2nd layer
template<int dim>
class InjectionOperator: public RestrictionOperator<dim>
{
    virtual double restrict1D(const Eigen::MatrixXd& mat, size_t idx){
        return mat(idx, 0);
    }
    virtual double restrict2D(const Eigen::MatrixXd& mat, size_t idx, size_t grid_num)
    {
        return mat(idx, 0);
    }
};


template<int dim>
class FullWeightOperator: public RestrictionOperator<dim>
{
    virtual double restrict1D(const Eigen::MatrixXd& mat, size_t idx){
        return (mat(idx - 1, 0) + mat(idx, 0) * 2 + mat(idx + 1, 0)) * 0.25;
    }
    virtual double restrict2D(const Eigen::MatrixXd& mat, size_t idx, size_t grid_num){
        return (mat(idx - grid_num - 1, 0)     +
                mat(idx - grid_num,     0) * 2 +
                mat(idx - grid_num + 1, 0)     +
                mat(idx - 1       ,     0) * 2 +
                mat(idx           ,     0) * 4 +
                mat(idx + 1       ,     0) * 2 +
                mat(idx + grid_num - 1   )     +
                mat(idx + grid_num,     0) * 2 +
                mat(idx + grid_num + 1, 0)
                ) * 0.0625;
    }
};


template<int dim>
class LinearInterpolationOperator: public InterpolationOperator<dim>
{
    virtual double interpolate1D(double left, double right){
        return (left + right) * 0.5;
    }
    virtual double interpolate2D(double left_up, double right_up, double left_down, double right_down){
        return (left_up + right_up + left_down + right_down) * 0.25;
    }
};


template<int dim>
class QuadraticInterpolationOperator: public InterpolationOperator<dim>
{
    virtual double interpolate1D(double left, double right){
        return 0; // todo
    }
    virtual double interpolate2D(double left_up, double right_up, double left_down, double right_down){
        return 0; // todo
    }
};