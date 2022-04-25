#pragma once
#include <iostream>
#include <utility>
#include <functional>
#include <jsoncpp/json/json.h>
#include "../Array2D.h" //TODO Cmake里解决这个问题

template<int dim>
class Multigrid
{
private:



public:
    Multigrid();
    Multigrid(const Json::Value &root);
    Multigrid(const Multigrid &) = delete;
};

template<int dim>
Multigrid<dim>::Multigrid() { }

template<int dim>
Multigrid<dim>::Multigrid(const Json::Value &root)
{

}