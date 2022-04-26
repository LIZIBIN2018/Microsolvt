#pragma once 
#include <iostream>


enum class RstOptType // Restriction Operator
{
    injection,
    fullWeighting
};

enum class ItpOptType // Interpolation Operator
{
    linear,
    quadratic
};

template<int dim>
class CycleSolver;


// 1st layer
template<int dim>
class TransferOperator
{
public:
    TransferOperator(CycleSolver<dim> *solver_ptr_) : solver_ptr(solver_ptr_){}
    virtual ~TransferOperator() {solver_ptr = nullptr;}
protected:
    CycleSolver<dim> *solver_ptr;
};

// 2nd layer
template<int dim>
class RestrictionOperator: public TransferOperator<dim>
{
public:
    RestrictionOperator(CycleSolver<dim> *solver_ptr_) : TransferOperator<dim>(solver_ptr_){}
};


template<int dim>
class InterpolationOperator: public TransferOperator<dim>
{
public:
    InterpolationOperator(CycleSolver<dim> *solver_ptr_) : TransferOperator<dim>(solver_ptr_){}
};


// 3rd layer
template<int dim>
class InjectionOperator: public RestrictionOperator<dim>
{
public:
    InjectionOperator(CycleSolver<dim> *solver_ptr_): RestrictionOperator<dim>(solver_ptr_) { }
};


template<int dim>
class FullWeightOperator: public RestrictionOperator<dim>
{
public:
    FullWeightOperator(CycleSolver<dim> *solver_ptr_): RestrictionOperator<dim>(solver_ptr_) { }
};


template<int dim>
class LinearInterpolationOperator: public InterpolationOperator<dim>
{
public:
    LinearInterpolationOperator(CycleSolver<dim> *solver_ptr_): InterpolationOperator<dim>(solver_ptr_) { }
};


template<int dim>
class QuadraticInterpolationOperator: public InterpolationOperator<dim>
{
public:
    QuadraticInterpolationOperator(CycleSolver<dim> *solver_ptr_): InterpolationOperator<dim>(solver_ptr_) { }
};