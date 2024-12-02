/*
    Solvers for the Ax=b system
*/
/*
    Solvers for the Ax=b system
*/
#ifndef _HT_SOLVER_H_
#define _HT_SOLVER_H_

#include "KokkosKernels_config.h"
#include "KokkosSparse_pcg.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_default_types.hpp"
#include <iostream>

#include "GwDomain.h"
#include "GwMatrix.h"
#include "RTMatrix.h"
#include "HTMatrix.h"

#include <vector>
#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosKernels_Handle.hpp>
#include <KokkosSparse_spiluk.hpp>
#include <KokkosSparse_sptrsv.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosBlas1_dot.hpp>
#include <cmath>
#include <memory>
#include <fstream>

#include <Eigen/Eigenvalues>
// using namespace std;
// using namespace Eigen;

using lno_view_t = Kokkos::View<int*>;
using lno_nnz_view_t = Kokkos::View<int*>;
using scalar_view_t = Kokkos::View<double*>;

template <typename execution_space>
class HTSolver {

    typedef typename KokkosSparse::CrsMatrix<default_scalar, default_lno_t, execution_space, void, default_size_type> crsMat_t;
    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
    typedef typename scalar_view_t::value_type scalar_t;
    typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type::value_type size_type;
    typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type::value_type lno_t;

public:
    int iter, iter_max, gsteps, nrow, nnz;
    real eps, eps_min;

    lno_view_t ht_ptr;
    lno_nnz_view_t ht_ind;
    scalar_view_t ht_val, ht_rhs, ht_vecx, ht_diag;

    inline void init(HTMatrix &htA, GwDomain &gdom) {
        nrow = htA.ht_nrow; nnz = htA.ht_nnz;
        iter_max = 40;
        gsteps = 20;
        iter_max = iter_max * gsteps;
        eps_min = 1e-7;
#if SERGHEI_KOKKOSKERNELS_SOLVER
        ht_ptr = lno_view_t("ht_ptr", htA.ht_nrow + 1);
        ht_ind = lno_nnz_view_t("ht_ind", htA.ht_nnz);
        ht_val = scalar_view_t("ht_val", htA.ht_nnz);
        ht_rhs = scalar_view_t("ht_rhs", htA.ht_nrow);
        ht_vecx = scalar_view_t("ht_vecx", htA.ht_nrow);
        ht_diag = scalar_view_t("ht_diag", htA.ht_nrow);
        Kokkos::deep_copy(ht_ptr, htA.ht_ptr);
        Kokkos::deep_copy(ht_ind, htA.ht_ind);
        Kokkos::deep_copy(ht_val, htA.ht_val);
        Kokkos::deep_copy(ht_rhs, htA.ht_rhs);
        Kokkos::deep_copy(ht_diag, htA.ht_diag);
#endif
    }

/*------------使用Eigen库求解矩阵------------------*/
//     void Gauss_Seidel(RTMatrix &rtA) {
//         using namespace Eigen;
//         int n = rtA.rt_nrow;
//         int nnz = rtA.rt_nnz;

//         // 使用指针访问Kokkos的数据
//         const int* ptr = rtA.rt_ptr.data();
//         const int* ind = rtA.rt_ind.data();
//         const double* val = rtA.rt_val.data();
//         const double* b_data = rtA.rt_rhs.data();
//         double* x_data = rtA.rt_x.data();

//     // 使用ILU预处理
// //     ilu_decompose(rtA);

// //     // 将ILU分解后的L和U矩阵转换为Eigen格式
// //     Map<const SparseMatrix<double, RowMajor>> eigenL(n, n, rtA.rt_L_values.size(), rtA.rt_L_row_map.data(), rtA.rt_L_entries.data(), rtA.rt_L_values.data());
// //     Map<const SparseMatrix<double, RowMajor>> eigenU(n, n, rtA.rt_U_values.size(), rtA.rt_U_row_map.data(), rtA.rt_U_entries.data(), rtA.rt_U_values.data());
// //        BiCGSTAB< SparseMatrix<double> > solver;
// //         solver.setTolerance(eps_min);
// //         solver.setMaxIterations(iter_max);
// //         solver.compute(eigenMatrix);

//         // 将指针转换为Eigen的数据结构
//         Map<const SparseMatrix<double, RowMajor>> eigenMatrix(n, n, nnz, ptr, ind, val);
//         Map<const VectorXd> b(b_data, n);
//         Map<VectorXd> x(x_data, n);


//      // //    使用BiCGSTAB求解

//      //    BiCGSTAB<SparseMatrix<double>> solver;
//      //    solver.setTolerance(eps_min);
//      //    solver.setMaxIterations(iter_max);
//      //    solver.compute(eigenMatrix);
	


// 	//    lscg最小二乘共轭梯度法
// 		LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
//           solver.setTolerance(eps_min);
//           solver.setMaxIterations(iter_max);
//           solver.compute(eigenMatrix);



//         if (solver.info() != Success) {
//             std::cerr << "Decomposition failed!" << std::endl;
//             return;
//         };

//         x = solver.solve(b);

//         if (solver.info() != Success) {
//             std::cerr << "Solving failed!" << std::endl;
//             return;
//         }

//         // 结果已经存在于x中，不需要额外复制
//     }
/*------------使用Eigen库求解矩阵------------------*/

/*---------ILU预处理--------------*/
// ILU (Incomplete LU) preconditioner
// void ilu_decompose(RTMatrix &rtA) {
//     using namespace KokkosSparse;
//     using namespace KokkosSparse::Experimental;
//     using namespace KokkosKernels;

//     // Types for the ILU preconditioner
//     using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
//         <default_size_type, default_lno_t, default_scalar, execution_space, execution_space, execution_space>;
//     using crsMat_t = CrsMatrix<default_scalar, default_lno_t, execution_space, void, default_size_type>;

//     // Copy data to Kokkos Views
//     crsMat_t rtmatA("rtmatA", rtA.rt_nrow, rtA.rt_nrow, rtA.rt_nnz, rt_val, rt_ptr, rt_ind);

//     // Create kernel handle and set parameters
//     KernelHandle kh;
//     kh.create_spiluk_handle(SparseILU2, 0); // SparseILU2 is a higher level of ILU, use SparseILU for level 0
//     auto spiluk_handle = kh.get_spiluk_handle();

//     // Allocate memory for L and U factors
//     lno_view_t L_row_map("L_row_map", rtA.rt_nrow + 1);
//     lno_nnz_view_t L_entries("L_entries", rtA.rt_nnz);
//     scalar_view_t L_values("L_values", rtA.rt_nnz);

//     lno_view_t U_row_map("U_row_map", rtA.rt_nrow + 1);
//     lno_nnz_view_t U_entries("U_entries", rtA.rt_nnz);
//     scalar_view_t U_values("U_values", rtA.rt_nnz);

//     // Set views in the handle
//     spiluk_handle->set_level_schedule(0); // Level 0 scheduling

//     // Perform ILU factorization
//     Experimental::spiluk_symbolic(&kh, rtmatA.numRows(), rtmatA.numRows(), rtmatA.nnz(),
//                                   spiluk_handle, rtmatA.graph.row_map, rtmatA.graph.entries);
//     Experimental::spiluk_numeric(&kh, rtmatA.numRows(), rtmatA.numRows(), rtmatA.nnz(),
//                                  spiluk_handle, rtmatA.graph.row_map, rtmatA.graph.entries, rtmatA.values);

//     // Copy the results back to the original RTMatrix
//     Kokkos::deep_copy(L_row_map, spiluk_handle->get_L_row_map());
//     Kokkos::deep_copy(L_entries, spiluk_handle->get_L_entries());
//     Kokkos::deep_copy(L_values, spiluk_handle->get_L_values());

//     Kokkos::deep_copy(U_row_map, spiluk_handle->get_U_row_map());
//     Kokkos::deep_copy(U_entries, spiluk_handle->get_U_entries());
//     Kokkos::deep_copy(U_values, spiluk_handle->get_U_values());

//     rtA.rt_L_row_map = L_row_map;
//     rtA.rt_L_entries = L_entries;
//     rtA.rt_L_values = L_values;
//     rtA.rt_U_row_map = U_row_map;
//     rtA.rt_U_entries = U_entries;
//     rtA.rt_U_values = U_values;

//     kh.destroy_spiluk_handle();
// }


/*---------ILU预处理--------------*/











    /*
        Gauss-Seidel solver
    */
    void Gauss_Seidel(HTMatrix &htA) {
        int n = htA.ht_nrow;

        // Convert Kokkos::View to std::vector for easier manipulation
        std::vector<double> x(n);
        std::vector<double> b(n);

        // Manually copy data from Kokkos::View to std::vector
        for (int i = 0; i < n; ++i) {
            b[i] = htA.ht_rhs(i);
            x[i] = htA.ht_x(i);
        }

        for (int iter = 0; iter < iter_max; ++iter) {
            std::vector<double> x_old = x;

            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                double diag = 0.0;

                for (int j = htA.ht_ptr(i); j < htA.ht_ptr(i + 1); ++j) {
                    int col = htA.ht_ind(j);
                    if (col == i) {
                        diag = htA.ht_val(j);
                    } else {
                        sum += htA.ht_val(j) * x[col];
                    }
                }

                if (diag != 0) {
                    x[i] = (b[i] - sum) / diag;
                }
            }

            // Check for convergence
            double error = 0.0;
            for (int i = 0; i < n; ++i) {
                error += std::pow(x[i] - x_old[i], 2);
            }
            error = std::sqrt(error);

            if (error < eps_min) {
                break;
            }
        }

        // Manually copy data back from std::vector to Kokkos::View
        for (int i = 0; i < n; ++i) {
            htA.ht_x(i) = x[i];
        }
    }
//输出矩阵
    /*
        Gauss-Seidel solver
    */

void save_matrix_and_vector_to_files(const HTMatrix &htA, const std::string &matrix_file, const std::string &vector_file) {
    int n = htA.ht_nrow;
    int nnz = htA.ht_nnz;

    // 将矩阵A保存到txt文件中
    std::ofstream matrix_out(matrix_file);
    if (!matrix_out.is_open()) {
        std::cerr << "Unable to open file for matrix output!" << std::endl;
        return;
    }

    matrix_out << n << " " << n << " " << nnz << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = htA.ht_ptr[i]; j < htA.ht_ptr[i + 1]; ++j) {
            matrix_out << i << " " << htA.ht_ind[j] << " " << htA.ht_val[j] << std::endl;
        }
    }

    matrix_out.close();

    // 将向量b保存到txt文件中
    std::ofstream vector_out(vector_file);
    if (!vector_out.is_open()) {
        std::cerr << "Unable to open file for vector output!" << std::endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        vector_out << htA.ht_rhs[i] << std::endl;
    }

    vector_out.close();
}


  
//输出矩阵

/*-------------检查矩阵是否对称正定---------------------*/
// bool is_symmetric(const RTMatrix &rtA)
// {
// 	using namespace std;
// 	using namespace Eigen;
//      int rt_nrow;
//      std::vector<int> rt_ptr;     // Row pointers
//      std::vector<int> rt_ind;     // Column indices
//      std::vector<double> rt_val;  // Non-zero values
//      std::vector<double> rt_diag; // Diagonal values
//      for (int i = 0; i < rtA.rt_nrow; ++i)
//      {
//           for (int j = rtA.rt_ptr[i]; j < rtA.rt_ptr[i + 1]; ++j)
//           {
//                int row = i;
//                int col = rtA.rt_ind[j];
//                double val = rtA.rt_val[j];

//                // Check if the corresponding off-diagonal element is present and equal
//                bool has_corresponding_element = false;
//                for (int k = rtA.rt_ptr[col]; k < rtA.rt_ptr[col + 1]; ++k)
//                {
//                     if (rtA.rt_ind[k] == row)
//                     {
//                          if (rtA.rt_val[k] != val)
//                          {
//                               return false;
//                          }
//                          has_corresponding_element = true;
//                          break;
//                     }
//                }

//                // If the corresponding off-diagonal element is not present, the matrix is not symmetric
//                if (!has_corresponding_element)
//                {
//                     return false;
//                }
//           }
//      }
//      return true;
// }
bool is_positive_definite(const HTMatrix &htA)
{

	using namespace std;
	using namespace Eigen;
     int ht_nrow;
     std::vector<int> ht_ptr;     // Row pointers
     std::vector<int> ht_ind;     // Column indices
     std::vector<double> ht_val;  // Non-zero values
     std::vector<double> ht_diag; // Diagonal values
     Eigen::MatrixXd matrix(htA.ht_nrow, htA.ht_nrow);

     // Construct the matrix from ht_val
     for (int i = 0; i < htA.ht_nrow; ++i)
     {
          for (int j = htA.ht_ptr[i]; j < htA.ht_ptr[i + 1]; ++j)
          {
               int col = htA.ht_ind[j];
               matrix(i, col) = htA.ht_val[j];
               if (i != col)
               {
                    matrix(col, i) = htA.ht_val[j];
               }
          }
     }

     // Compute eigenvalues
     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix);
     if (eigensolver.info() != Eigen::Success)
     {
          std::cerr << "Eigenvalue computation failed." << std::endl;
          return false;
     }
     // Check if all eigenvalues are positive
     Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
     std::cout << "Eigenvalues:" << eigenvalues << std::endl;

     for (int i = 0; i < eigenvalues.size(); ++i)
     {
          if (eigenvalues[i] <= 0)
          {
               return false;
          }
     }

     return true;
}

/*-------------------矩阵正定性检查---------------*/
	/*
		----------------------------------------------------------
		----------------------------------------------------------
							Iterative Solvers
		----------------------------------------------------------
		----------------------------------------------------------
	*/

	/*
		PCG from KokkosKernels
	*/
	// Top-level PCG solver 预处理共轭梯度法
	void kkpcg(HTMatrix &htA)
	{
		bool usePreconditioner = 1;
		//usePreconditioner 变量标志是否使用预条件器，默认为使用
		const unsigned cg_iteration_limit = 1000000;
		//cg_iteration_limit 定义了PCG算法的最大迭代次数。
		const double   cg_iteration_tolerance     = 1e-8 ;
		//cg_iteration_tolerance 定义了PCG算法的收敛容限。

		decompose(htA);
		// decompose(A) 对矩阵进行分解，
		Kokkos::deep_copy (ht_ind, htA.ht_ind);
		Kokkos::deep_copy (ht_val, htA.ht_val);
		Kokkos::deep_copy (ht_rhs, htA.ht_rhs);
		Kokkos::deep_copy (ht_diag, htA.ht_diag);

		crsMat_t htmatA = crsMat_t("htmatA", htA.ht_nrow, htA.ht_nrow, htA.ht_nnz, ht_val, ht_ptr, ht_ind);

		// initialize KokkosKernels solver
		KokkosKernels::Experimental::Example::CGSolveResult cg_result ;
		typedef KokkosKernels::Experimental::KokkosKernelsHandle
		  < size_type, lno_t, scalar_t, execution_space, execution_space, execution_space > KernelHandle;
		KernelHandle ht_kh;
		ht_kh.create_gs_handle();
		//Kokkos::Impl::Timer timer1;
		KokkosKernels::Experimental::Example::pcgsolve(ht_kh, htmatA, ht_rhs, ht_vecx, ht_diag
		    , cg_iteration_limit, cg_iteration_tolerance, &cg_result, usePreconditioner);
		Kokkos::fence();
		// 调用 pcgsolve 函数进行PCG求解，传入了求解句柄 kh、矩阵 matA、右侧向量 rhs、
		// 解向量 vecx、对角元素 diag、迭代次数限制、迭代收敛容限等参数，
		// 并将求解结果保存在 cg_result 中。
		//solve_time = timer1.seconds();
		//std::cout  << "DEFAULT SOLVE: " << "(P)CG_NUM_ITER = [" << cg_result.iteration << "], " << "RESIDUAL = [" << cg_result.norm_res << "]"<< std::endl ;
		ht_kh.destroy_gs_handle();
		// 调用 destroy_gs_handle() 方法销毁求解句柄

		Kokkos::deep_copy (htA.ht_x, ht_vecx);
		// 将解向量 vecx 复制回 A 中
	}


	/*
    	CG Solver
    */
    int cg(HTMatrix &htA, GwDomain &gdom)	{
		real rho, rhoOld, alpha, beta;
		rhoOld = 0.0;
		decompose(htA);
		// initialize residual
		get_residual(htA);
		iter = 0;	eps = 1.0;
		while (iter < iter_max & eps > eps_min)	{
			precJACO(htA);
			rho = dot(htA.ht_r, htA.ht_z, htA.ht_nrow);
			if (iter == 0)	{
				beta = 0.0;
				Kokkos::deep_copy(htA.ht_p, htA.ht_z);
			}
			else {
				beta = rho / rhoOld;
				kxpy(htA.ht_p, beta, htA.ht_p, htA.ht_z, htA.ht_nrow);
			}
			mul_MV(htA.ht_q, htA, htA.ht_p);
			alpha = rho / dot(htA.ht_q, htA.ht_p, htA.ht_nrow);
			update_X(htA, alpha);
			kxpy(htA.ht_r, -alpha, htA.ht_q, htA.ht_r, htA.ht_nrow);
			rhoOld = rho;
			eps = pow(dot(htA.ht_r, htA.ht_r, htA.ht_nrow),0.5);
			iter += 1;
    	}
		// std::cerr << "      >> CG solver converges in " << iter << " iterations! eps = " << eps <<"\n";
    	return iter;
    }

	/*
    	Jacobi Preconditioner
    */
    void precJACO(HTMatrix htA)	{
		Kokkos::parallel_for( htA.ht_nrow , KOKKOS_LAMBDA(int idom) {
			htA.ht_z(idom) = htA.ht_r(idom) * htA.ht_diag(idom);
		});
    }

	/*
		Reset A.x
	*/
	void reset(realArr x, int n)	{
		Kokkos::parallel_for(n , KOKKOS_LAMBDA(int idom) {x(idom) = 0.0;});
	}

	/*
		Get residual
	*/
	void get_residual(HTMatrix &htA)	{
		Kokkos::parallel_for( htA.ht_nrow , KOKKOS_LAMBDA(int idom) {
			int icol;
			htA.ht_r(idom) = htA.ht_rhs(idom);
			for (icol = htA.ht_ptr(idom); icol < htA.ht_ptr(idom+1); icol++)	{
				htA.ht_r(idom) -= htA.ht_val(icol) * htA.ht_x(htA.ht_ind(icol));
			}
		});
	}

	/*
		Dot product of two vectors
	*/
	real dot(realArr v1, realArr v2, int n)	{
		real out;
		Kokkos::parallel_reduce( n , KOKKOS_LAMBDA (int idx, real &prod) {
			prod += v1(idx) * v2(idx);
		} , out);
		return out;
	}

	/*
		kX+Y
	*/
	void kxpy(realArr out, real k, realArr x, realArr y, int n)	{
		Kokkos::parallel_for( n , KOKKOS_LAMBDA(int idx) {
			out(idx) = k * x(idx) + y(idx);
		});
	}

	/*
		Matrix - Vector Multiplication
	*/
	void mul_MV(realArr out, HTMatrix htA, realArr x)	{
		Kokkos::parallel_for( htA.ht_nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			out(idx) = 0.0;
			for (icol = htA.ht_ptr(idx); icol < htA.ht_ptr(idx+1); icol++)	{
				out(idx) += htA.ht_val(icol) * x(htA.ht_ind(icol));
			}
		});
	}

	/*
		Update solution
	*/
	void update_X(HTMatrix htA, real alpha)	{
		Kokkos::parallel_for( htA.ht_nrow , KOKKOS_LAMBDA(int idx) {
			htA.ht_x(idx) += alpha * htA.ht_p(idx);
		});
	}

	/*
		Get diagonal of Matrix
	*/
	void decompose(HTMatrix htA)	{
		Kokkos::parallel_for(htA.ht_nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			htA.ht_diag(idx) = 0.0;
			for (icol = htA.ht_ptr(idx); icol < htA.ht_ptr(idx+1); icol++)	{
				if (htA.ht_ind(icol) == idx)	{htA.ht_diag(idx) = 1.0 / htA.ht_val(icol);}
			}
		});
	}

};
#endif
