#pragma once 
#include <iostream>


template<int dim>
class CycleSolver;


// 1st layer
template<int dim>
class TransferOperator
{
protected:
    CycleSolver<dim> *solver_ptr;
};

// 2nd layer
template<int dim>
class RestrictionOperator: public TransferOperator<dim>
{
    
};


template<int dim>
class InterpolationOperator: public TransferOperator<dim>
{

};


// 3rd layer
template<int dim>
class InjectionOperator: public RestrictionOperator<dim>
{

};


template<int dim>
class FullWeightOperator: public RestrictionOperator<dim>
{

};


template<int dim>
class LinearInterpolationOperator: public InterpolationOperator<dim>
{

};


template<int dim>
class QuadraticInterpolationOperator: public InterpolationOperator<dim>
{

};