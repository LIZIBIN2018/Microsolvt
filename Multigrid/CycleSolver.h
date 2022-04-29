#pragma once
#include "TransferOperator.h"
#include "Multigrid.h"
#include <functional>
#include <utility>
#include <vector>
#include <eigen3/Eigen/Eigen>

#include "../TimeCounter.h"

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
        try{
            auto rstOptStr = root["restriction_operator"].asString();
            auto itpOptStr = root["interpolation_operator"].asString();

            // 限制算子
            if (rstOptStr == "injection")
                rst_opt = new InjectionOperator<dim>();
            else if (rstOptStr == "fullWeighting")
                rst_opt = new FullWeightOperator<dim>();
            else
                throw "Invalid transfer opetator";

            // 插值算子
            if (itpOptStr == "linear")
                itp_opt = new LinearInterpolationOperator<dim>();
            else if (itpOptStr == "quadratic")
                itp_opt = new QuadraticInterpolationOperator<dim>();
            else
                throw "Invalid transfer opetator";

            // 求解器类型
            auto cycle_type = root["cycle"].asString();
            if (cycle_type == "V-cycle")
                solver_type = SolverType::VCycle;
            else if (cycle_type == "FullMultigridVCycle")
                solver_type = SolverType::FullMultigridVCycle;
            else
                throw "Invalid cycle type";

            // 最大迭代次数
            max_iteration = root["max_iteration"].asInt();
            if(max_iteration <= 0)
                throw "Invalid max_iteration";

            // 停止精度
            rel_accuracy = root["rel_accuracy"].asDouble();
            if(rel_accuracy <= 0)
                throw "Invalid rel_accuracy";

            // 最粗糙网格密度
            coarest = root["coarest"].asInt();
            if(coarest <= 1 || coarest & (coarest - 1) != 0)
                throw "Invalid parameter coarest";
        }
        catch(const char* s){
            std::cerr << s << std::endl;
            exit(1);
        }

        grid_size = root["grid_n"].asInt() - 1; // n-1
        
        // 矩阵装配
        int grid_size_cur = grid_size;
        while(grid_size_cur + 1 >= coarest)
        {
            Eigen::SparseMatrix<double> A;
            generateSparseAh(A,grid_size_cur);
            Ahs.push_back(A);
            grid_size_cur = grid_size_cur >> 1;
        }
    }
    virtual ~CycleSolver() { 
        if(rst_opt) delete rst_opt; 
        if(itp_opt) delete itp_opt; 
    }

    void solve(Multigrid<dim> &grid, void *f) //DONE
    {
        if(dim == 1)
        {
            using Function = std::function<double(double)>;
            auto f_orig = *(static_cast<Function *>(f));
            solve(grid,f_orig);
        }
        else if(dim == 2)
        {
            using Function = std::function<double(double,double)>;
            auto f_orig = *(static_cast<Function *>(f));
            solve(grid,f_orig);
        }
    }


private:
    void solve(Multigrid<dim> &grid, std::function<double(double)> f)
    {
        // 根据v(i) = 1/2*(v(i-1) + v(i+1) + h*h*fj) 构造Av = f, 此处f(j) = f(v(j))*h*h
        const Eigen::SparseMatrix<double> &Ah = Ahs[0]; 
        Eigen::MatrixXd uh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);  //ground truth 
        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);
        Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);
        double h = 1.0/(grid.grid_size + 1);

        for (int i = 0; i < grid.grid_size; i++)
        {
            fh(i) = f((i + 1) * grid.grid_length) * h * h;
            uh(i) = exp(sin((i + 1) * grid.grid_length)) - (1-(i + 1) * grid.grid_length) - ((i + 1) * grid.grid_length)*exp(sin(1));
        }
        for (size_t iter = 0; iter < max_iteration; iter++)
        {
            if(solver_type == SolverType::VCycle)
            {
                VCycle(vh, fh, 2, 1, grid.grid_size); 
            }
            else if(solver_type == SolverType::FullMultigridVCycle)
            {
                vh = FullMultigridVCycle(vh, fh, 3, grid.grid_size); //TODO
            }
            for(int i = 0; i < grid.grid_size; i++) 
            {
                grid.data(i) = vh(i);
            }

            // Calculate Relative Error
            Eigen::MatrixXd errorVector = fh - Ah*vh;
            double error = errorVector.norm()/grid.grid_size;
            double rel_error = error / fh.norm() * grid.grid_size;
            std::cout << "Iteration " << iter << ":    AbsError = " << error << ", RelError = " << rel_error << ", residual" << (fh-Ah*vh).norm() << std::endl;
            if(rel_error < rel_accuracy)
                break;

        }
    }

    void solve(Multigrid<dim> &grid,std::function<double(double, double)> f)
    {
        Eigen::SparseMatrix<double> &Ah = Ahs[0];
        Eigen::MatrixXd uh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);  //ground truth 
        Eigen::MatrixXd vh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);
        Eigen::MatrixXd fh = Eigen::MatrixXd::Zero(grid.grid_node_num, 1);
        double h = 1.0/(grid.grid_size + 1);

        for (int i = 0; i < grid.grid_size; i++)
        {
            for (int j = 0; j < grid.grid_size; j++)
            {
                fh(place2index(i, j, grid.grid_size)) = f((i + 1) * grid.grid_length, (j + 1) * grid.grid_length)*h*h;
            }
        }

        for (size_t iter = 0; iter < max_iteration; iter++)
        {
            if(solver_type == SolverType::VCycle)
            {
                VCycle(vh, fh, 2, 1, grid.grid_size); 
            }
            else if(solver_type == SolverType::FullMultigridVCycle)
            {
                vh = FullMultigridVCycle(vh, fh, 1, grid.grid_size);
            }
            for (int i = 0; i < grid.grid_node_num; i++)
            {
                grid.data(i) = vh(i);
            } 

            // Calculate Relative Error
            Eigen::MatrixXd errorVector = GetResidual(Ah, vh, fh);
            double max_norm = maxNorm(errorVector);
            double rel_error = max_norm / maxNorm(vh);
            std::cout << "Iteration " << iter << ":    Error = " << rel_error << std::endl;
            if(rel_error < rel_accuracy)
                break;
        }
    }

    void VCycle(Eigen::MatrixXd &vh, const Eigen::MatrixXd &fh, size_t nu1, size_t nu2, size_t grid_size_cur)
    {   
        int layer = round(log(double(grid_size + 1)/(grid_size_cur + 1)) / log(2));
        Eigen::SparseMatrix<double> Ah = Ahs[layer];

        for(int i = 0; i < nu1; i++)
        {
            Relax(Ah, vh, fh);
            //std::cout << "当前误差" << GetResidual(Ah,vh,fh).norm() << std::endl;
        }

        if (grid_size_cur + 1 > coarest)
        {
            
            Eigen::MatrixXd f2h = (*rst_opt)(GetResidual(Ah,vh,fh), grid_size_cur);
            Eigen::MatrixXd v2h = Eigen::MatrixXd::Zero(pow(grid_size_cur >> 1,dim), 1);
            VCycle(v2h, f2h, nu1, nu2, grid_size_cur >> 1);
            
            (*itp_opt)(v2h, grid_size_cur >> 1);
            //std::cout <<"插值后的大小" << (*itp_opt)(v2h, grid_size_cur >> 1).size() << std::endl;
            //std::cout <<"vh的大小" << vh.size() << std::endl;
           
            vh += (*itp_opt)(v2h, grid_size_cur >> 1);
            
        }



        for (int i = 0; i < nu2; i++)
        {
            Relax(Ah, vh, fh);
            //std::cout << "当前误差" << GetResidual(Ah,vh,fh).norm() << std::endl;
        }
    }

    void generateSparseAh(Eigen::SparseMatrix<double> &Ah, size_t grid_size)
    {

        size_t mat_size = pow(grid_size, dim);
        Ah.resize(mat_size, mat_size);
        Ah.setIdentity();
        Ah *= 2 * dim;

        if (dim == 1)
        {
            for (int i = 0; i < grid_size - 1; ++i)
            {
                Ah.insert(i, i + 1) = -1;
                Ah.insert(i + 1, i) = -1;
            }
            // std::cout << Ah << std::endl; // TODELETE
        }
        else if (dim == 2)
        {
            for (int i = 0; i < mat_size - grid_size; ++i)
            {
                Ah.insert(i, i + grid_size) = -1;
                Ah.insert(i + grid_size, i) = -1;
            }
            for (int i = 0; i < mat_size; ++i)
            {
                if (i % grid_size == 0)
                    continue;
                Ah.insert(i - 1, i) = -1;
                Ah.insert(i, i - 1) = -1;
            }
            // std::cout << Ah << std::endl; //TODELETE
        }
        else
            throw "dim error";
    }

    void generateAh(Eigen::MatrixXd &Ah, size_t grid_size)
    {
        size_t mat_size = pow(grid_size, dim);
        Ah = Eigen::MatrixXd::Identity(mat_size,mat_size)*2*dim;
        if(dim == 1)
        {
            for(int i = 0; i < grid_size - 1; ++i)
            {
                Ah(i,i+1) = -1;
                Ah(i+1,i) = -1;
            }
            //std::cout << Ah << std::endl; // TODELETE
        }
        else if(dim == 2)
        {
            for(int i = 0;i < mat_size - grid_size; ++i)
            {
                Ah(i, i+grid_size) = -1;
                Ah(i+grid_size, i) = -1;

            }
            for(int i = 0; i < mat_size; ++i)
            {
                if(i%grid_size == 0)
                    continue;
                Ah(i-1, i) = -1;
                Ah(i, i-1) = -1;
            }
            //std::cout << Ah << std::endl; //TODELETE
        }
        else
            throw "dim error";
    }

    void Relax(const Eigen::MatrixXd &A, Eigen::MatrixXd &v,const Eigen::MatrixXd &f)
    {
        v += -omega/(2*dim)*A*v + omega/(2*dim)*f;
    }
    
    Eigen::MatrixXd GetResidual(const Eigen::MatrixXd &A,const Eigen::MatrixXd &v,const Eigen::MatrixXd &f)
    {
        return f - A*v;
    }

    Eigen::MatrixXd FullMultigridVCycle(Eigen::MatrixXd vh, Eigen::MatrixXd &fh, size_t nu0, size_t grid_size_cur)
    {
        if(grid_size_cur + 1> coarest) 
        {
            auto f2h = (*rst_opt)(fh, grid_size_cur);
            auto v2h = (*rst_opt)(vh, grid_size_cur);
            vh = (*itp_opt)(FullMultigridVCycle(v2h, f2h, nu0, grid_size_cur >> 1),grid_size_cur >> 1);
        }
        else
        {
            vh = Eigen::MatrixXd::Zero(grid_size_cur, 1);
        }

        for(int i = 0; i<nu0; i++)
        {
            VCycle(vh, fh, 3, 3, grid_size_cur);
        }
        return vh;
    }

    double maxNorm(Eigen::MatrixXd& v) //DONE
    {
        double m = 0.0;
        for (size_t i = 0; i < v.size(); i++)
        {
            if (std::abs(v(i)) > m)
                m = std::abs(v(i));
        }
        return m;
    }

private: // data
    RestrictionOperator<dim>                   *rst_opt = nullptr;   // restriction
    InterpolationOperator<dim>                 *itp_opt = nullptr;   // interpolation
    SolverType                                  solver_type; 
    size_t                                      coarest;
    size_t                                      max_iteration;
    size_t                                      grid_size;
    double                                      rel_accuracy;
    double                                      omega = 2.0/3.0;
    std::vector<Eigen::SparseMatrix<double>>    Ahs;   //pre calculate


private: // TOOLS
    size_t place2index(size_t row, size_t col, size_t grid_size)
    {
        return row * grid_size + col;
    }
};
