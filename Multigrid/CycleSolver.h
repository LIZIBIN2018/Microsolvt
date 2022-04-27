#pragma once
#include "TransferOperator.h"
#include <functional>
#include "Multigrid.h"
#include <vector>

enum class SolverType
{
    VCycle,
    FullMultigridVCycle
};

template <int dim>
class CycleSolver
{
public: // ctor & dtor
    CycleSolver(Json::Value &root)
    {
        try
        {
            auto rstOptStr = root["restriction_operator"].asString();
            auto itpOptStr = root["interpolation_operator"].asString();
            if (rstOptStr == "injection")
                rst_opt = new InjectionOperator<dim>();
            else if (rstOptStr == "fullWeighting")
                rst_opt = new FullWeightOperator<dim>();
            else
                throw 1;
            if (itpOptStr == "linear")
                itp_opt = new LinearInterpolationOperator<dim>();
            else if (itpOptStr == "quadratic")
                itp_opt = new QuadraticInterpolationOperator<dim>();
            else
                throw 1;
        }
        catch (...)
        {
            std::cerr << "Invalid transfer opetator" << '\n';
        }

        try
        {
            auto cycle_type = root["cycle"].asString();
            if (cycle_type == "V-cycle")
                solver_type = SolverType::VCycle;
            else if (cycle_type == "FullMultigridVCycle")
                solver_type = SolverType::FullMultigridVCycle;
            else
                throw 1;
        }
        catch (...)
        {
            std::cerr << "Invalid cycle type" << std::endl;
            exit(1);
        }

        max_iteration = root["max_iteration"].asint();
        if(max_iteration <= 0)
            throw std::exception("Invalid max_iteration");

        rel_accuracy = root["rel_accuracy"].asDouble();
        if(rel_accuracy <= 0)
            throw std::exception("Invalid rel_accuracy");
    }
    virtual ~CycleSolver() { 
        if(rst_opt) delete rst_opt; 
        if(itp_opt) delete itp_opt; 
    }

    virtual void solve(const Multigrid<dim> &grid,std::function<double(double)> f)
    {
        // 根据v(i) = 1/2*(v(i-1) + v(i+1) + h*h*fj) 构造Av = f, 此处f(j) = h*h*f(v(j))

        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_size,1);
        Eigen::MatrixXd fh = Eigen::MatrixXd;
        fh.resize(grid.grid_size,1);
        for (size_t iter = 0; iter < max_iteration; iter++)
        {
            // TODO 并行
            for(int i = 0;i < grid.grid_size;i++)
            {
                fh(i) = f(i*grid.grid_length)*grid.grid_length*grid.grid_length;
            }
            if(solver_type == SolverType::VCycle)
            {
                vh = VCycle(vh, f, 3, 3, grid.grid_size); //nu1 ,nu2 = ? TODO
            }
            else if(solver_type == SolverType::FullMultigridVCycle)
            {
                vh = FullMultigridVCycle(); //TODO
            }
            for(int i = 0; i < grid_size; i++) //TOCHECK +-1
            {
                grid.data(i) = vh(i);
            } 

            // Calculate Relative Error
            double max_norm = (vh - As[0] * fh).maxCoeff();
            std::cout << "Iteration " << iter << ":    Error = " << max_norm << std::endl;
            if(max_norm / vh.maxCoeff() < rel_accuracy)
                break;
        }
    }

    virtual void solve(const Multigrid<dim> &grid,std::function<double(double, double)> f)
    {
        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_size*grid.grid_size,1);
        Eigen::MatrixXd fh = Eigen::MatrixXd;
        fh.resize(grid.grid_size*grid.grid_size,1);
        // TODO dim==2的情况是如何？
    }

    // TODO 既然传引用，那么还要返回么？
    Eigen::MatrixXd VCycle(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t nu1, size_t nv2, size_t grid_size_cur)
    {
        ////迭代法：v=Rw v + wD^-1f
        // 迭代矩阵装配 TODO
        Eigen::MatrixXd A;
        Eigen::MatrixXd Rw;
        Eigen::MatrixXd wD_inv;
        if(dim == 1)
        {

        }
        else if(dim == 2)
        {

        }
        // 对方程A^h u^h = f^h 迭代nu1次
        for(int i = 0;i<nu1;i++)
        {
            v = Rw*v + f;
        }
        if(grid_size_cur + 1 > coarest)
        {
            auto f_new = (*rst_opt)(f - A*v, grid_size_cur);
            auto v_new = Eigen::MatrixXd::zero((grid_size_cur+1)/2 -1, 1);
            v_new = VCycle(v_new, f_new, nu1, nu2, (grid_size_cur+1)/2 -1);
            v += (*itp_opt)(v_new, grid_size_cur);
        }
        for(int i = 0; i<nu2; i++)
        {
            v = Rw*v + wD_inv*f;
        }
    }
    
    Eigen::MatrixXd FullMultigridVCycle(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t nu0, size_t grid_size_cur)
    {
        Eigen::MatrixXd Rw;
        Eigen::MatrixXd wD_inv;
        // 矩阵装配
        if(dim == 1)
        {

        }
        if(dim ==2 )
        {

        }

        if(grid_size_cur + 1> coarest) 
        {
            auto f_new = (*rst_opt)(f, grid_size_cur);
            auto v_new = FullMultigridVCycle(f_new, nu0, (grid_size_cur+1)/2 - 1);
            auto v = (*itp_opt)(v, grid_size_cur);
        }
        else
        {
            v = Eigen::MatrixXd::Zero(v.rows(), 1);
        }
        for(int i = 0; i<v0; i++)
        {
            v = Rw*v + wD_inv*f;
        }
        return v;
    }

private:                                 // data
    RestrictionOperator<dim>   *rst_opt = nullptr;   // restriction
    InterpolationOperator<dim> *itp_opt = nullptr; // interpolation
    SolverType                  solver_type; 
    size_t                      coarest = 4;
    std::vector<Eigen::MatrixXd>As;
    size_t                      max_iteration;
    double                      rel_accuracy;
};
