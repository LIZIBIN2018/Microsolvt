#pragma once
#include "TransferOperator.h"
#include "Multigrid.h"
#include <functional>
#include <utility>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <tbb/tbb.h>

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
        auto rstOptStr = root["restriction_operator"].asString();
        auto itpOptStr = root["interpolation_operator"].asString();
        if (rstOptStr == "injection")
            rst_opt = new InjectionOperator<dim>();
        else if (rstOptStr == "fullWeighting")
            rst_opt = new FullWeightOperator<dim>();
        else
            throw "Invalid transfer opetator";
        if (itpOptStr == "linear")
            itp_opt = new LinearInterpolationOperator<dim>();
        else if (itpOptStr == "quadratic")
            itp_opt = new QuadraticInterpolationOperator<dim>();
        else
            throw "Invalid transfer opetator";

        auto cycle_type = root["cycle"].asString();
        if (cycle_type == "V-cycle")
            solver_type = SolverType::VCycle;
        else if (cycle_type == "FullMultigridVCycle")
            solver_type = SolverType::FullMultigridVCycle;
        else
            throw "Invalid cycle type";

        max_iteration = root["max_iteration"].asInt();
        if(max_iteration <= 0)
            throw "Invalid max_iteration";

        rel_accuracy = root["rel_accuracy"].asDouble();
        if(rel_accuracy <= 0)
            throw "Invalid rel_accuracy";
    }
    virtual ~CycleSolver() { 
        if(rst_opt) delete rst_opt; 
        if(itp_opt) delete itp_opt; 
    }

    double maxNorm(Eigen::MatrixXd& v){
        double m = 0.0;
        for (size_t i = 0; i < v.size(); i++)
        {
            if (std::abs(v(i)) > m)
                m = std::abs(v(i));
        }
        return m;
    }

    void solve(Multigrid<dim> &grid, std::function<double(double)> f)
    {
        // 根据v(i) = 1/2*(v(i-1) + v(i+1) + h*h*fj) 构造Av = f, 此处f(j) = h*h*f(v(j))

        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_size, 1);
        Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(grid.grid_size, 1);
        fh.resize(grid.grid_size, 1);
        
        for (int i = 0; i < grid.grid_size; i++)
        {
            fh(i) = f((i + 1) * grid.grid_length);
        }

        for (size_t iter = 0; iter < max_iteration; iter++)
        {
            if(solver_type == SolverType::VCycle)
            {
                VCycle(vh, fh, 3, 3, grid.grid_size); //nu1 ,nu2 = ? TODO
            }
            else if(solver_type == SolverType::FullMultigridVCycle)
            {
                FullMultigridVCycle(vh, fh, 3, grid.grid_size); //TODO
            }
            for(int i = 0; i < grid.grid_size; i++) 
            {
                grid.data(i) = vh(i);
            }

            // Calculate Relative Error
            Eigen::MatrixXd errorVector = getResidual(vh, fh, grid.grid_size);
            double residue = maxNorm(errorVector);
            double rel_error = residue / maxNorm(vh);
            std::cout << "Iteration " << iter << ":    AbsError = " << residue << ", RelError = " << rel_error << std::endl;
            if(rel_error < rel_accuracy)
                break;
        }
    }
    

    void solve(Multigrid<dim> &grid,std::function<double(double, double)> f)
    {
        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_size * grid.grid_size, 1);
        Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(grid.grid_size * grid.grid_size, 1);
        fh.resize(grid.grid_size * grid.grid_size, 1);
        for (size_t iter = 0; iter < max_iteration; iter++)
        {
            for (int i = 0; i < grid.grid_size; i++)
            {
                for (int j = 0; j < grid.grid_size; j++)
                {
                    fh(place2index(i, j, grid.grid_size)) = f((i + 1) * grid.grid_length, (j + 1) * grid.grid_length);
                }
            }
            if(solver_type == SolverType::VCycle)
            {
                VCycle(vh, fh, 3, 3, grid.grid_size); //nu1 ,nu2 = ? TODO
            }
            else if(solver_type == SolverType::FullMultigridVCycle)
            {
                FullMultigridVCycle(vh, fh, 3, grid.grid_size);
            }
            for (int i = 0; i < grid.grid_size * grid.grid_size; i++)
            {
                grid.data(i) = vh(i);
            } 

            // Calculate Relative Error
            Eigen::MatrixXd errorVector = getResidual(vh, fh, grid.grid_size);
            double max_norm = maxNorm(errorVector);
            double rel_error = max_norm / maxNorm(vh);
            std::cout << "Iteration " << iter << ":    Error = " << rel_error << std::endl;
            if(rel_error < rel_accuracy)
                break;
        }
    }

    void Relax(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t grid_size_cur){
        Eigen::MatrixXd v0 = v;
        const size_t threadNum = 16;
        const double omega = 2.0 / 3;



        //TODO parallel for有大问题
        if(dim == 1){
            tbb::parallel_for((size_t)0, threadNum, (size_t)1,   [&](size_t i){
                size_t bound = std::min((grid_size_cur / threadNum + 1)* (i + 1), grid_size_cur);
                for (size_t j = grid_size_cur / threadNum * i; j < bound; j++){
                    v(j) = v0(j) - omega * 0.5 * (2 * v0(j) - (j == 0 ? 0 : v0(j - 1)) - (j == grid_size_cur - 1 ? 0 : v0(j + 1)))
                            + omega / 2 / ((grid_size_cur + 1) * (grid_size_cur + 1)) * f(j);  //TOCHECK 
                }
            });
        }
        else if(dim == 2){
            tbb::parallel_for((size_t)0, threadNum, (size_t)1,   [&](size_t i){
                size_t bound = std::min((grid_size_cur / threadNum + 1)* (i + 1), grid_size_cur);
                for (size_t row = 0; row < bound; row++)
                {
                    for (size_t col = 0; col < grid_size_cur; col++)
                    {
                        size_t idx = row * grid_size_cur + col;
                        v(idx) = v0(idx) 
                                - omega * 0.25 * (4 * v0(idx) 
                                    - (col == 0                 || row == 0                 ? 0 : v0(place2index(row - 1, col - 1, grid_size_cur))) 
                                    - (col == 0                 || row == grid_size_cur - 1 ? 0 : v0(place2index(row + 1, col - 1, grid_size_cur))) 
                                    - (col == grid_size_cur - 1 || row == 0                 ? 0 : v0(place2index(row - 1, col + 1, grid_size_cur))) 
                                    - (col == grid_size_cur - 1 || row == grid_size_cur - 1 ? 0 : v0(place2index(row + 1, col + 1, grid_size_cur))))
                                + omega / 4 / ((grid_size_cur + 1) * (grid_size_cur + 1)) * f(idx); //TOCHECK
                    }
                }
            });
        }
        else{
            throw "Invalid dim";
        }
    }

    Eigen::MatrixXd getResidual(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t grid_size){
        Eigen::MatrixXd r = v;
        const size_t threadNum = 16;
        if(dim == 1){
            tbb::parallel_for((size_t)0, threadNum, (size_t)1, [&](size_t i){
                size_t bound = std::min((grid_size / threadNum + 1) * (i + 1), grid_size);
                for (size_t j = grid_size / threadNum * i; j < bound; j++){
                    r(j) = f(j) - pow(grid_size, 2) * (2 * v(j) - (j == 0 ? 0 : v(j - 1)) - (j == grid_size - 1 ? 0 : v(j + 1)));
                }
            });
        }
        else if(dim == 2){
            tbb::parallel_for((size_t)0, threadNum, (size_t)1, [&](size_t i){
                size_t bound = std::min((grid_size / threadNum + 1) * (i + 1), grid_size);
                for (size_t row = 0; row < bound; row++)
                {
                    for (size_t col = 0; col < grid_size; col++)
                    {
                        size_t idx = row * grid_size + col;
                        r(idx) = f(idx) - pow(grid_size, 2) * (4 * v(idx) 
                                    - (col == 0             || row == 0             ? 0 : v(place2index(row - 1, col - 1, grid_size))) 
                                    - (col == 0             || row == grid_size - 1 ? 0 : v(place2index(row + 1, col - 1, grid_size))) 
                                    - (col == grid_size - 1 || row == 0             ? 0 : v(place2index(row - 1, col + 1, grid_size))) 
                                    - (col == grid_size - 1 || row == grid_size - 1 ? 0 : v(place2index(row + 1, col + 1, grid_size))));
                    }
                }
            });
        }
        else{
            throw "Invalid dim";
        }
        return r;
    }

    void VCycle(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t nu1, size_t nu2, size_t grid_size_cur)
    {
        const size_t threadNum = 16;

        // std::cout << "f = " << f << std::endl;
        // 对方程A^h u^h = f^h 迭代nu1次
        for (int i = 0; i < nu1; i++)
        {
            Relax(v, f, grid_size_cur);
            // std::cout << "Size = " << grid_size_cur << "Residue = " << getResidue(v, f, grid_size_cur) << std::endl;
        }

            // std::cout << grid_size_cur << std::endl;
        if (grid_size_cur + 1 > coarest)
        {
            auto f_new = (*rst_opt)(getResidual(v, f, grid_size_cur), grid_size_cur);
            // std::cout << f(0) << ',' << f_new(0) << std::endl;
            Eigen::MatrixXd v_new = Eigen::MatrixXd::Zero(grid_size_cur >> 1, 1);
            VCycle(v_new, f_new, nu1, nu2, grid_size_cur >> 1);
            v += (*itp_opt)(v_new, grid_size_cur >> 1);
        }

        for (int i = 0; i < nu2; i++)
        {
            Relax(v, f, grid_size_cur);
        }
    }
    
    void FullMultigridVCycle(Eigen::MatrixXd &v, Eigen::MatrixXd &f, size_t nu0, size_t grid_size_cur)
    {
        if(grid_size_cur + 1> coarest) 
        {
            auto f_new = (*rst_opt)(f, grid_size_cur);
            FullMultigridVCycle(v, f_new, nu0, grid_size_cur >> 1);
            v = (*itp_opt)(v, grid_size_cur);
        }
        else
        {
            v = Eigen::MatrixXd::Zero(v.rows(), 1);
        }
        for(int i = 0; i<nu0; i++)
        {
            VCycle(v, f, 3, 3, grid_size_cur);
        }
    }

private: // data
    RestrictionOperator<dim>   *rst_opt = nullptr;   // restriction
    InterpolationOperator<dim> *itp_opt = nullptr;   // interpolation
    SolverType                  solver_type; 
    size_t                      coarest = 8;
    size_t                      max_iteration;
    double                      rel_accuracy;


private: // TOOLS
    size_t place2index(size_t row, size_t col, size_t grid_size)
    {
        return row * grid_size + col;
    }
};
