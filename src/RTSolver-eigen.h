/*
    Solvers for the Ax=b system
*/
/*
    Solvers for the Ax=b system
*/
#ifndef _RT_SOLVER_H_
#define _RT_SOLVER_H_

#include "KokkosKernels_config.h"
#include "KokkosSparse_pcg.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_default_types.hpp"
#include <iostream>
#include "GwDomain.h"
#include "GwMatrix.h"
#include "RTMatrix.h"
#include <vector>
#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <Eigen/SVD>
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
class RTSolver {

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

    lno_view_t rt_ptr;
    lno_nnz_view_t rt_ind;
    scalar_view_t rt_val, rt_rhs, rt_vecx, rt_diag;

    inline void init(RTMatrix &rtA, GwDomain &gdom) {
        nrow = rtA.rt_nrow; nnz = rtA.rt_nnz;
        iter_max = 400;
        gsteps = 20;
        iter_max = iter_max * gsteps;
        eps_min = 1e-7;
#if SERGHEI_KOKKOSKERNELS_SOLVER
        rt_ptr = lno_view_t("rt_ptr", rtA.rt_nrow + 1);// 每一行非零元素的起始索引
        rt_ind = lno_nnz_view_t("rt_ind", rtA.rt_nnz);// 非零元素的列索引
        rt_val = scalar_view_t("rt_val", rtA.rt_nnz);// 非零元素的值
        rt_rhs = scalar_view_t("rt_rhs", rtA.rt_nrow); // 右侧向量b
        rt_vecx = scalar_view_t("rt_vecx", rtA.rt_nrow);// 近似解向量x
        rt_diag = scalar_view_t("rt_diag", rtA.rt_nrow);
        Kokkos::deep_copy(rt_ptr, rtA.rt_ptr);
        Kokkos::deep_copy(rt_ind, rtA.rt_ind);
        Kokkos::deep_copy(rt_val, rtA.rt_val);
        Kokkos::deep_copy(rt_rhs, rtA.rt_rhs);
        Kokkos::deep_copy(rt_diag, rtA.rt_diag);
#endif
    }

//!
	//!使用Eigen库求解矩阵
//!

//计算矩阵秩和奇异值
void computeRank(const Eigen::MatrixXd& matrix) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    double tolerance = matrix.norm() * matrix.size() * std::numeric_limits<double>::epsilon();
    int rank = 0;
    for (int i = 0; i < svd.singularValues().size(); ++i) {
        if (svd.singularValues()(i) > tolerance) {
            ++rank;
        }
    }
    std::cout << "Matrix rank: " << rank << std::endl;//打印矩阵秩
    std::cout << "Singular values:\n" << svd.singularValues() << std::endl;//打印奇异值
}
//!
    void Eigen(RTMatrix &rtA) {
        using namespace Eigen;
        int n = rtA.rt_nrow;
        int nnz = rtA.rt_nnz;

        // 使用指针访问Kokkos的数据
        const int* ptr = rtA.rt_ptr.data();
        const int* ind = rtA.rt_ind.data();
        const double* val = rtA.rt_val.data();
        const double* b_data = rtA.rt_rhs.data();
        double* x_data = rtA.rt_x.data();

    // 使用ILU预处理
    // ilu_decompose(rtA);

    // 将ILU分解后的L和U矩阵转换为Eigen格式
    // Map<const SparseMatrix<double, RowMajor>> eigenL(n, n, rtA.rt_L_values.size(), rtA.rt_L_row_map.data(), rtA.rt_L_entries.data(), rtA.rt_L_values.data());
    // Map<const SparseMatrix<double, RowMajor>> eigenU(n, n, rtA.rt_U_values.size(), rtA.rt_U_row_map.data(), rtA.rt_U_entries.data(), rtA.rt_U_values.data());
    //    BiCGSTAB< SparseMatrix<double> > solver;
    //     solver.setTolerance(eps_min);
    //     solver.setMaxIterations(iter_max);
    //     solver.compute(eigenMatrix);

        // 将指针转换为Eigen的数据结构
        Map<const SparseMatrix<double, RowMajor>> eigenMatrix(n, n, nnz, ptr, ind, val);
        Map<const VectorXd> b(b_data, n);
        Map<VectorXd> x(x_data, n);


    //  //!    使用BiCGSTAB求解

        BiCGSTAB<SparseMatrix<double>> solver;
        solver.setTolerance(eps_min);
        solver.setMaxIterations(iter_max);
        solver.compute(eigenMatrix);

	// //!    lscg最小二乘共轭梯度法
	// 	LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
    //       solver.setTolerance(eps_min);
    //       solver.setMaxIterations(iter_max);
    //       solver.compute(eigenMatrix);


// computeRank(eigenMatrix);//计算矩阵的秩和奇异值


        if (solver.info() != Success) {
            std::cerr << "Decomposition failed!" << std::endl;
            return;
        };

        x = solver.solve(b);
        

        if (solver.info() != Success) {
            std::cerr << "Solving failed!" << std::endl;
            return;
        }

        // 结果已经存在于x中，不需要额外复制
    }

//!


//!/*---------ILU预处理--------------*/
//ILU (Incomplete LU) preconditioner
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

//! 打印RTMatrix矩阵的函数
// void printRTMatrix(const RTMatrix &rtA) {
//     int n = rtA.rt_nrow;
//     for (int i = 0; i < n; ++i) {
//         for (int j = rtA.rt_ptr[i]; j < rtA.rt_ptr[i + 1]; ++j) {
//             int col = rtA.rt_ind[j];
//             double val = rtA.rt_val[j];
//             std::cout << "A[" << i << "][" << col << "] = " << val << " ";
//         }
//         std::cout << " | b[" << i << "] = " << rtA.rt_rhs[i] << " | x[" << i << "] = " << rtA.rt_x[i] << std::endl;
//     }
// }
void printRTMatrix(const RTMatrix &rtA) {
    int n = rtA.rt_nrow; // 获取矩阵的行数
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { // 遍历每一列
            double val = 0.0; // 默认值为0
            for (int k = rtA.rt_ptr[i]; k < rtA.rt_ptr[i + 1]; ++k) {
                if (rtA.rt_ind(k) == j) {
                    val = rtA.rt_val(k); // 如果找到非零值，则更新val
                    break; // 找到后退出循环
                }
            }
            std::cout << "A[" << i << "][" << j << "] = " << val << " ";
        }
        std::cout << " | b[" << i << "] = " << rtA.rt_rhs[i] << " | x[" << i << "] = " << rtA.rt_x[i] << std::endl;
    }
}
//!
   //! Gauss-Seidel solver
//!
    void Gauss_Seidel(RTMatrix &rtA) {
        int n = rtA.rt_nrow;

        // Convert Kokkos::View to std::vector for easier manipulation
        std::vector<double> x(n);
        std::vector<double> b(n);

        // Manually copy data from Kokkos::View to std::vector
        for (int i = 0; i < n; ++i) {
            b[i] = rtA.rt_rhs(i);
            x[i] = rtA.rt_x(i);
        }

        for (int iter = 0; iter < iter_max; ++iter) {
            std::vector<double> x_old = x;

            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                double diag = 0.0;

                for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i + 1); ++j) {
                    int col = rtA.rt_ind(j);
                    if (col == i) {
                        diag = rtA.rt_val(j);
                    } else {
                        sum += rtA.rt_val(j) * x[col];
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
            rtA.rt_x(i) = x[i];
        }


    // 打印矩阵
    // printRTMatrix(rtA);
    }


//!kokkos solver 2 Jacobi-BICGSTAB
/*
    Jacobi Preconditioner
*/
void Jacobi_Precondition(RTMatrix &rtA) {
    int n = rtA.rt_nrow;
    realArr &rt_diag = rtA.rt_diag;
    realArr &rt_val = rtA.rt_val;
    intArr &rt_ptr = rtA.rt_ptr;
    intArr &rt_ind = rtA.rt_ind;
    realArr &rt_z = rtA.rt_z;

    // Jacobi Preconditioning: z = M^-1 * r, where M is the diagonal of A
    Kokkos::parallel_for("Jacobi_Precondition", n, KOKKOS_LAMBDA(int i) {
        double diag_val = 0.0;
        // 对角线元素直接使用
        for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i + 1); ++j) {
            int col = rtA.rt_ind(j);
            if (col == i) {
                diag_val = rtA.rt_val(j);  // 仅保存对角元素
            }
        }
        
        if (diag_val != 0.0) {
            rtA.rt_z(i) = rtA.rt_r(i) / diag_val;
        } else {
            rtA.rt_z(i) = 0.0;  // 如果对角元素为零，直接设为零
        }
    });
}

/*
    BICGSTAB Solver using Jacobi Preconditioner
*/
int Jacobi_BICGSTAB_Solve(RTMatrix &rtA) {
    int n = rtA.rt_nrow;
    int iter = 0;
    real eps = 1.0, eps_min = 1e-6, rho = 0.0, rhoOld = 0.0;
    real alpha, beta;
    realArr &rt_r = rtA.rt_r;
    realArr &rt_z = rtA.rt_z;
    realArr &rt_p = rtA.rt_p;
    realArr &rt_q = rtA.rt_q;
    realArr &rt_x = rtA.rt_x;

    // Initialize residual
    get_residual(rtA);
    while (iter < iter_max && eps > eps_min) {
        // Jacobi Preconditioning
        Jacobi_Precondition(rtA);
        
        rho = dot(rt_r, rt_z, n);  // rho = r' * z
        if (iter == 0) {
            beta = 0.0;
            Kokkos::deep_copy(rt_p, rt_z);  // p = z
        } else {
            beta = rho / rhoOld;
            kxpy(rt_p, beta, rt_p, rt_z, n);  // p = z + beta * p
        }

        mul_MV(rt_q, rtA, rt_p);  // q = A * p
        alpha = rho / dot(rt_q, rt_p, n);  // alpha = rho / (p' * q)
        
        update_X(rtA, alpha);  // x = x + alpha * p
        kxpy(rt_r, -alpha, rt_q, rt_r, n);  // r = r - alpha * q

        Jacobi_Precondition(rtA);  // Apply Jacobi Preconditioning to new residual

        rhoOld = rho;
        eps = std::sqrt(dot(rt_r, rt_r, n));  // Error check: eps = sqrt(r' * r)
        iter += 1;
    }

    return iter;
}

/*
    Dot product of two vectors
*/
real dot(realArr v1, realArr v2, int n) {
    real out;
    Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(int idx, real &prod) {
        prod += v1(idx) * v2(idx);
    }, out);
    return out;
}

/*
    kX + Y
*/
void kxpy(realArr out, real k, realArr x, realArr y, int n) {
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idx) {
        out(idx) = k * x(idx) + y(idx);
    });
}

/*
    Matrix - Vector Multiplication
*/
void mul_MV(realArr out, RTMatrix &rtA, realArr x) {
    Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
        int icol;
        out(idx) = 0.0;
        for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); icol++) {
            out(idx) += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
        }
    });
}

/*
    Update solution
*/
void update_X(RTMatrix &rtA, real alpha) {
    Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
        rtA.rt_x(idx) += alpha * rtA.rt_p(idx);
    });
}

/*
    Get residual
*/
void get_residual(RTMatrix &rtA) {
    Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
        int icol;
        rtA.rt_r(idx) = rtA.rt_rhs(idx);
        for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); icol++) {
            rtA.rt_r(idx) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
        }
    });
}


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

bool is_positive_definite(const RTMatrix &rtA)
{

	using namespace std;
	using namespace Eigen;
     int rt_nrow;
     std::vector<int> rt_ptr;     // Row pointers
     std::vector<int> rt_ind;     // Column indices
     std::vector<double> rt_val;  // Non-zero values
     std::vector<double> rt_diag; // Diagonal values
     Eigen::MatrixXd matrix(rtA.rt_nrow, rtA.rt_nrow);

     // Construct the matrix from rt_val
     for (int i = 0; i < rtA.rt_nrow; ++i)
     {
          for (int j = rtA.rt_ptr[i]; j < rtA.rt_ptr[i + 1]; ++j)
          {
               int col = rtA.rt_ind[j];
               matrix(i, col) = rtA.rt_val[j];
               if (i != col)
               {
                    matrix(col, i) = rtA.rt_val[j];
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
	// void kkpcg(RTMatrix &rtA)
	// {
	// 	bool usePreconditioner = 1;
	// 	//usePreconditioner 变量标志是否使用预条件器，默认为使用
	// 	const unsigned cg_iteration_limit = 1000000;
	// 	//cg_iteration_limit 定义了PCG算法的最大迭代次数。
	// 	const double   cg_iteration_tolerance     = 1e-8 ;
	// 	//cg_iteration_tolerance 定义了PCG算法的收敛容限。

	// 	decompose(rtA);
	// 	// decompose(A) 对矩阵进行分解，
	// 	Kokkos::deep_copy (rt_ind, rtA.rt_ind);
	// 	Kokkos::deep_copy (rt_val, rtA.rt_val);
	// 	Kokkos::deep_copy (rt_rhs, rtA.rt_rhs);
	// 	Kokkos::deep_copy (rt_diag, rtA.rt_diag);

	// 	crsMat_t rtmatA = crsMat_t("rtmatA", rtA.rt_nrow, rtA.rt_nrow, rtA.rt_nnz, rt_val, rt_ptr, rt_ind);

	// 	// initialize KokkosKernels solver
	// 	KokkosKernels::Experimental::Example::CGSolveResult cg_result ;
	// 	typedef KokkosKernels::Experimental::KokkosKernelsHandle
	// 	  < size_type, lno_t, scalar_t, execution_space, execution_space, execution_space > KernelHandle;
	// 	KernelHandle rt_kh;
	// 	rt_kh.create_gs_handle();
	// 	//Kokkos::Impl::Timer timer1;
	// 	KokkosKernels::Experimental::Example::pcgsolve(rt_kh, rtmatA, rt_rhs, rt_vecx, rt_diag
	// 	    , cg_iteration_limit, cg_iteration_tolerance, &cg_result, usePreconditioner);
	// 	Kokkos::fence();
	// 	// 调用 pcgsolve 函数进行PCG求解，传入了求解句柄 kh、矩阵 matA、右侧向量 rhs、
	// 	// 解向量 vecx、对角元素 diag、迭代次数限制、迭代收敛容限等参数，
	// 	// 并将求解结果保存在 cg_result 中。
	// 	//solve_time = timer1.seconds();
	// 	//std::cout  << "DEFAULT SOLVE: " << "(P)CG_NUM_ITER = [" << cg_result.iteration << "], " << "RESIDUAL = [" << cg_result.norm_res << "]"<< std::endl ;
	// 	rt_kh.destroy_gs_handle();
	// 	// 调用 destroy_gs_handle() 方法销毁求解句柄

	// 	Kokkos::deep_copy (rtA.rt_x, rt_vecx);
	// 	// 将解向量 vecx 复制回 A 中
	// }


	// /*
    // 	CG Solver
    // */
    // int cg(RTMatrix &rtA, GwDomain &gdom)	{
	// 	real rho, rhoOld, alpha, beta;
	// 	rhoOld = 0.0;
	// 	decompose(rtA);
	// 	// initialize residual
	// 	get_residual(rtA);
	// 	iter = 0;	eps = 1.0;
	// 	while (iter < iter_max & eps > eps_min)	{
	// 		precJACO(rtA);
	// 		rho = dot(rtA.rt_r, rtA.rt_z, rtA.rt_nrow);
	// 		if (iter == 0)	{
	// 			beta = 0.0;
	// 			Kokkos::deep_copy(rtA.rt_p, rtA.rt_z);
	// 		}
	// 		else {
	// 			beta = rho / rhoOld;
	// 			kxpy(rtA.rt_p, beta, rtA.rt_p, rtA.rt_z, rtA.rt_nrow);
	// 		}
	// 		mul_MV(rtA.rt_q, rtA, rtA.rt_p);
	// 		alpha = rho / dot(rtA.rt_q, rtA.rt_p, rtA.rt_nrow);
	// 		update_X(rtA, alpha);
	// 		kxpy(rtA.rt_r, -alpha, rtA.rt_q, rtA.rt_r, rtA.rt_nrow);
	// 		rhoOld = rho;
	// 		eps = pow(dot(rtA.rt_r, rtA.rt_r, rtA.rt_nrow),0.5);
	// 		iter += 1;
    // 	}
	// 	// std::cerr << "      >> CG solver converges in " << iter << " iterations! eps = " << eps <<"\n";
    // 	return iter;
    // }

	// /*
    // 	Jacobi Preconditioner
    // */
    // void precJACO(RTMatrix rtA)	{
	// 	Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idom) {
	// 		rtA.rt_z(idom) = rtA.rt_r(idom) * rtA.rt_diag(idom);
	// 	});
    // }

	// /*
	// 	Reset A.x
	// */
	// void reset(realArr x, int n)	{
	// 	Kokkos::parallel_for(n , KOKKOS_LAMBDA(int idom) {x(idom) = 0.0;});
	// }

	// /*
	// 	Get residual
	// */
	// void get_residual(RTMatrix &rtA)	{
	// 	Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idom) {
	// 		int icol;
	// 		rtA.rt_r(idom) = rtA.rt_rhs(idom);
	// 		for (icol = rtA.rt_ptr(idom); icol < rtA.rt_ptr(idom+1); icol++)	{
	// 			rtA.rt_r(idom) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
	// 		}
	// 	});
	// }

	// /*
	// 	Dot product of two vectors
	// */
	// real dot(realArr v1, realArr v2, int n)	{
	// 	real out;
	// 	Kokkos::parallel_reduce( n , KOKKOS_LAMBDA (int idx, real &prod) {
	// 		prod += v1(idx) * v2(idx);
	// 	} , out);
	// 	return out;
	// }

	// /*
	// 	kX+Y
	// */
	// void kxpy(realArr out, real k, realArr x, realArr y, int n)	{
	// 	Kokkos::parallel_for( n , KOKKOS_LAMBDA(int idx) {
	// 		out(idx) = k * x(idx) + y(idx);
	// 	});
	// }

	// /*
	// 	Matrix - Vector Multiplication
	// */
	// void mul_MV(realArr out, RTMatrix rtA, realArr x)	{
	// 	Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
	// 		int icol;
	// 		out(idx) = 0.0;
	// 		for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++)	{
	// 			out(idx) += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
	// 		}
	// 	});
	// }

	// /*
	// 	Update solution
	// */
	// void update_X(RTMatrix rtA, real alpha)	{
	// 	Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
	// 		rtA.rt_x(idx) += alpha * rtA.rt_p(idx);
	// 	});
	// }

	// /*
	// 	Get diagonal of Matrix
	// */
	// void decompose(RTMatrix rtA)	{
	// 	Kokkos::parallel_for(rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
	// 		int icol;
	// 		rtA.rt_diag(idx) = 0.0;
	// 		for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++)	{
	// 			if (rtA.rt_ind(icol) == idx)	{rtA.rt_diag(idx) = 1.0 / rtA.rt_val(icol);}
	// 		}
	// 	});
	// }

};
#endif
