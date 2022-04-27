#pragma once 
#include <iostream>
#include <eigen3/Eigen/Dense>



// 1st layer 
template<int dim>
class RestrictionOperator
{
public:
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num) = 0;
};


template<int dim>
class InterpolationOperator
{
public:
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num) = 0;
};


// 2nd layer
template<int dim>
class InjectionOperator: public RestrictionOperator<dim>
{
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        return rst;
    }
};


template<int dim>
class FullWeightOperator: public RestrictionOperator<dim>
{
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        return rst;
    }
};


template<int dim>
class LinearInterpolationOperator: public InterpolationOperator<dim>
{
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        return rst;
    }
};


template<int dim>
class QuadraticInterpolationOperator: public InterpolationOperator<dim>
{
    virtual Eigen::MatrixXd operator()(Eigen::MatrixXd v, size_t grid_num)
    {
        Eigen::MatrixXd rst;
        return rst;
    }
};