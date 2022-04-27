#pragma once 
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <tbb/tbb.h>


// 1st layer 
template<int dim>
class RestrictionOperator
{
private:
    virtual double restrict1D(Eigen::MatrixXd& mat, size_t idx) = 0;
    virtual double restrict2D(Eigen::MatrixXd& mat, size_t idx, size_t grid_num) = 0;
public:
    Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        const size_t threadNum = 16;
        size_t halfLength = grid_num >> 1;
        if(dim == 1){
            if(grid_num & (grid_num + 1) != 0 || v.size() != grid_num)
                throw std::exception("Invalid vector size");
            rst.resize(halfLength, 1);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, threadNum), [&](size_t i){
                size_t bound = std::min(halfLength / threadNum * (i + 1), halfLength);
                for (size_t j = rst.size() / threadNum * i; j < bound; j++){
                    rst(j, 0) = restrict1D(v, j << 1 | 1);
                }
            });
        }
        else if(dim == 2){
            if(grid_num & (grid_num + 1) != 0 || v.size() != grid_num * grid_num)
                throw std::exception("Invalid vector size");
            rst.resize(halfLength * halfLength, 1);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, threadNum), [&](size_t i){
                size_t bound = std::min(halfLength / threadNum * (i + 1), halfLength);
                for (size_t row = 0; row < bound; row++)
                {
                    for (size_t col = 0; col < halfLength; col++)
                    {
                        rst(row * halfLength + col, 0) = restrict2D(v, (row << 1 | 1) * grid_num + (col << 1 | 1), grid_num);
                    }
                }
            });
        }
        else{
            throw std::exception("Invalid dim");
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
    Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd itp;
        const size_t threadNum = 16;
        size_t doubleLength = grid_num << 1 | 1;
        if(dim == 1){
            if(grid_num & (grid_num + 1) != 0 || v.size() != grid_num)
                throw std::exception("Invalid vector size");
            itp.resize(doubleLength, 1);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, threadNum), [&](size_t i){
                size_t bound = std::min(doubleLength / threadNum * (i + 1), doubleLength);
                for (size_t j = doubleLength / threadNum * i; j < bound; j++){
                    itp(j) = interpolate1D(j == 0 ? 0 : v((j - 1) >> 1), j == doubleLength - 1 ? 0 : v(j >> 1, 0));
                }
            });
        }
        else if(dim == 2){
            if(grid_num & (grid_num + 1) != 0 || v.size() != grid_num * grid_num)
                throw std::exception("Invalid vector size");
            itp.resize(doubleLength * doubleLength, 1);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, threadNum), [&](size_t i){
                size_t bound = std::min(doubleLength / threadNum * (i + 1), doubleLength);
                for (size_t row = 0; row < bound; row++)
                {
                    for (size_t col = 0; col < doubleLength; col++)
                    {
                        itp(row * doubleLength + col) = interpolate2D(
                            row == 0            || col == 0            ? 0 : v((row - 1) >> 1) * grid_num + ((col - 1) >> 1),
                            row == 0            || col == grid_num - 1 ? 0 : v((row - 1) >> 1) * grid_num + ((col    ) >> 1),
                            row == grid_num - 1 || col == 0            ? 0 : v((row    ) >> 1) * grid_num + ((col - 1) >> 1),
                            row == grid_num - 1 || col == grid_num - 1 ? 0 : v((row    ) >> 1) * grid_num + ((col    ) >> 1)
                        );
                    }
                }
            });
        }
        else{
            throw std::exception("Invalid dim");
        }
        return itp;
    }
};


// 2nd layer
template<int dim>
class InjectionOperator: public RestrictionOperator<dim>
{
    virtual double restrict1D(Eigen::MatrixXd& mat, size_t idx){
        return mat(idx, 0);
    }
    virtual double restrict2D(Eigen::MatrixXd& mat, size_t idx, size_t grid_num){
        return mat(idx, 0);
    }
};


template<int dim>
class FullWeightOperator: public RestrictionOperator<dim>
{
    virtual double restrict1D(Eigen::MatrixXd& mat, size_t idx){
        return (mat(idx - 1, 0) + mat(idx - 1, 0) * 2 + mat(idx - 1, 0)) * 0.25;
    }
    virtual double restrict2D(Eigen::MatrixXd& mat, size_t idx, size_t grid_num){
        return (mat(idx - grid_num - 1, 0)    +
                mat(idx - grid_num, 0)    * 2 +
                mat(idx - grid_num + 1, 0)    +
                mat(idx - 1, 0)           * 2 +
                mat(idx, 0)               * 4 +
                mat(idx + 1, 0)           * 2 +
                mat(idx + grid_num - 1)       +
                mat(idx + grid_num, 0)    * 2 +
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