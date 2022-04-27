#pragma once
#include "TransferOperator.h"
#include <functional>
#include "Multigrid.h"
#include <utility>
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
        // TODO 并行
        for(int i = 0;i < grid.grid_size;i++)
        {
            fh(i) = f((i+1)*grid.grid_length)*grid.grid_length*grid.grid_length;
        }
        if(solver_type == SolverType::VCycle)
        {
            vh = VCycle(vh, f, 3, 3, grid.grid_size); //nu1 ,nu2 = ? TODO
        }
        else if(solver_type == SolverType::FullMultigridVCycle)
        {
            vh = FullMultigridVCycle(); //TODO
        }
        for(int i = 0; i < grid_size; i++) 
        {
            grid.data(i) = vh(i);
        } 
    }

    virtual void solve(const Multigrid<dim> &grid,std::function<double(double, double)> f)
    {
        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_size*grid.grid_size,1);
        Eigen::MatrixXd fh = Eigen::MatrixXd;
        fh.resize(grid.grid_size*grid.grid_size,1);
        // TODO 并行
        for(int i = 0; i < grid.grid_size; i++)
        {
            for(int j = 0;j < grid.grid_size; j++)
            {
                int index = place2index(i,j);
                fh(index) = f((i+1)*grid.grid_length, (j+1)*grid.grid_length);
            }
        }
        if(solver_type == SolverType::VCycle)
        {
            vh = VCycle(vh, f, 3, 3, grid.grid_size); //nu1 ,nu2 = ? TODO
        }
        else if(solver_type == SolverType::FullMultigridVCycle)
        {
            vh = FullMultigridVCycle(vh,fh,3,grid.grid_size); 
        }
        for(int i = 0; i < grid_size*grid_size; i++) 
        {
            grid.data(i) = vh(i);
        } 
    }

    Eigen::MatrixXd VCycle(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t nu1, size_t nv2, size_t grid_size_cur)
    {
        ////迭代法：v=Rw v + wD^-1f
        // 迭代矩阵装配 TODO
        Eigen::MatrixXd A; //[-1 2 -1]/hh

        // A finished
        A = 2*Eigen::MatrixXd::Identity(grid_size_cur, grid_size_cur);
        for(int i = 0;i < grid_size_cur - 1;i++)
        {
            A(i,i+1) = 1;
            A(i+1,i) = 1;
        }
        A = A*pow(grid_size_cur+1,2);

        if (dim == 1)
        {
            Eigen::MatrixXd Rw;
            double wD_inv;
            Rw = Eigen::MatrixXd::Identity(grid_size_cur, grid_size_cur) - A / pow(grid_size_cur, 2);
            wD_inv = wD_inv = 1 / (3 * pow(1 + grid_size_cur, 2));
            // 对方程A^h u^h = f^h 迭代nu1次
            for (int i = 0; i < nu1; i++)
            {
                v = Rw * v + wD_inv * f;
            }
            if (grid_size_cur + 1 > coarest)
            {
                auto f_new = (*rst_opt)(f - A * v, grid_size_cur);
                auto v_new = Eigen::MatrixXd::zero((grid_size_cur + 1) / 2 - 1, 1);
                v_new = VCycle(v_new, f_new, nu1, nu2, (grid_size_cur + 1) / 2 - 1);
                v += (*itp_opt)(v_new, grid_size_cur);
            }
            for (int i = 0; i < nu2; i++)
            {
                v = Rw * v + wD_inv * f;
            }
        }
        else if (dim == 2)
        {

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


private: // TOOLS
    std::pair<int,int> index2place(int index, double grid_length) 
    { //TODO
        int n = round(1.0/grid_length);
        return std::pair<int, int>{index%};
    }
    int place2index(int i,int j)
    {//TODO
        return 0;
    }

};
