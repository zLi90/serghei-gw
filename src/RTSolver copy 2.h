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


#include "KokkosSparse_bicgstab.hpp"



using lno_view_t = Kokkos::View<int*>;
using lno_nnz_view_t = Kokkos::View<int*>;
using scalar_view_t = Kokkos::View<double*>;


using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosKernels::Experimental::Example;



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
        // Kokkos::deep_copy(rt_ind, rtA.rt_ind);
        // Kokkos::deep_copy(rt_val, rtA.rt_val);
        // Kokkos::deep_copy(rt_rhs, rtA.rt_rhs);
        // Kokkos::deep_copy(rt_diag, rtA.rt_diag);
#endif
    }

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
   //! BICGSTAB solver with ILU preconditioner kokkos-kernels
//!

void kkpbicgstab(RTMatrix &rtA) {
    using namespace KokkosKernels::Experimental;
    using namespace Example;

    // 参数配置
    constexpr int fill_lev = 1;       // ILU填充级别
    const unsigned iter_max = 8000; // 从rtA获取迭代限制
    const double eps_min = 1e-7;     // 收敛容差

    // 确保初始化已完成
    assert(rtA.rt_ptr.data() && "Call init() first!");

    // 数据同步到设备端（如果init中注释了deep_copy）
    Kokkos::deep_copy(rt_ind, rtA.rt_ind);    // 列索引
    Kokkos::deep_copy(rt_val, rtA.rt_val);    // 非零值
    Kokkos::deep_copy(rt_rhs, rtA.rt_rhs);    // 右侧向量b
    Kokkos::deep_copy(rt_diag, rtA.rt_diag);  // 初始解

    // 构建CRS矩阵（使用init初始化的视图）
    crsMat_t matA("matA", nrow, nrow, nnz, 
                  rt_val, rt_ptr, rt_ind);

    // 创建求解器和预条件子句柄
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
        typename crsMat_t::size_type,
        typename crsMat_t::ordinal_type,
        typename crsMat_t::value_type,
        typename crsMat_t::execution_space,
        typename crsMat_t::memory_space,
        typename crsMat_t::memory_space>;

    KernelHandle solver_handle;
    KernelHandle prec_handle;

    // 配置句柄参数（示例配置）
    solver_handle.set_team_work_size(16);
    prec_handle.set_team_work_size(16);


    // 调用BICGSTAB求解器
    BICGSTABResult result;

    bicgsolve(solver_handle, prec_handle, matA, rt_rhs, rt_vecx, 
                iter_max, eps_min, &result, fill_lev);   
                
                
                Kokkos::fence();
                Kokkos::deep_copy(rtA.rt_x, rt_vecx);



    // 输出诊断信息
    if (result.norm_res > eps_min) {
        std::cerr << "WARNING: BICGSTAB未收敛! 残差: " 
                  << result.norm_res << std::endl;
    }
    printf("BICGSTAB迭代 %lu次, 残差: %.2e, 总时间: %.2fs\n",
           result.iteration, result.norm_res, result.iter_time);
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
                // printf("Gauss-Seidel converged in %d iterations, error value: %e\n", iter, error);
                
                break;
            }

            // else{
            //     printf("Gauss-Seidel can't converged in %d iterations, error value: %e\n", iter, error);
            //     continue;
            // }

            // if (iter == iter_max - 1) {
            //     printf("!!!!!!Gauss-Seidel did not converge after %d iterations, error value: %e !!!!\n", iter, error);
            // }
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
    // printRTMatrix(rtA);
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



};
#endif
