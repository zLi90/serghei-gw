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
#include "RTMatrix.h"
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


//#include "KokkosSparse_bicgstab.hpp"



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

// void printRTMatrix(const RTMatrix &rtA) {
//     int n = rtA.rt_nrow; // 获取矩阵的行数
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) { // 遍历每一列
//             double val = 0.0; // 默认值为0
//             for (int k = rtA.rt_ptr[i]; k < rtA.rt_ptr[i + 1]; ++k) {
//                 if (rtA.rt_ind(k) == j) {
//                     val = rtA.rt_val(k); // 如果找到非零值，则更新val
//                     break; // 找到后退出循环
//                 }
//             }
//             std::cout << "A[" << i << "][" << j << "] = " << val << " ";
//         }
//         std::cout << " | b[" << i << "] = " << rtA.rt_rhs[i] << " | x[" << i << "] = " << rtA.rt_x[i] << std::endl;
//     }
// }



void printRTMatrix(const RTMatrix &rtA) {
    const double THRESHOLD = 1e10;  // 定义异常大值的阈值
    int n = rtA.rt_nrow;
    bool foundAbnormal = false;  // 标记是否发现异常值

    for (int i = 0; i < n; ++i) {
        // 检查右端向量b
        if (std::isnan(rtA.rt_rhs[i]) || std::isinf(rtA.rt_rhs[i]) || std::abs(rtA.rt_rhs[i]) > THRESHOLD) {
            std::cout << "异常值 - b[" << i << "] = " << rtA.rt_rhs[i] << std::endl;
            foundAbnormal = true;
        }

        // 检查解向量x
        if (std::isnan(rtA.rt_x[i]) || std::isinf(rtA.rt_x[i]) || std::abs(rtA.rt_x[i]) > THRESHOLD) {
            std::cout << "异常值 - x[" << i << "] = " << rtA.rt_x[i] << std::endl;
            foundAbnormal = true;
        }

        for (int j = 0; j < n; ++j) {
            double val = 0.0;
            for (int k = rtA.rt_ptr[i]; k < rtA.rt_ptr[i + 1]; ++k) {
                if (rtA.rt_ind[k] == j) {
                    val = rtA.rt_val[k];
                    break;
                }
            }

            // 检查矩阵元素
            if (std::isnan(val)) {
                std::cout << "NaN异常 - A[" << i << "][" << j << "]" << std::endl;
                foundAbnormal = true;
            } else if (std::isinf(val)) {
                std::cout << "无穷大异常 - A[" << i << "][" << j << "] = " << val << std::endl;
                foundAbnormal = true;
            } else if (std::abs(val) > THRESHOLD) {
                std::cout << "过大值异常 - A[" << i << "][" << j << "] = " << val << std::endl;
                foundAbnormal = true;
            } 
            // else if (val != 0.0 && std::abs(val) < 1e-10) {
            //     std::cout << "过小值异常 - A[" << i << "][" << j << "] = " << val << std::endl;
            //     foundAbnormal = true;
            // }
        }
    }

    if (!foundAbnormal) {
        std::cout << "未检测到异常数值" << std::endl;
    } else {
        std::cout << "\n检测到以上异常值，可能导致求解结果为NaN" << std::endl;
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

// printf("rtA.rt_x(%d)=%f\n", i, rtA.rt_x(i));
        if (isnan(rtA.rt_x(i))) {

// printf("!!!!!rtA.rt_x(%d)=%f\n", i, rtA.rt_x(i));
            printRTMatrix(rtA);
        }
	}


// 打印矩阵
// printRTMatrix(rtA);
}





// Helper functions (dot, kxpy, mul_MV, get_residual remain mostly the same)
real dot(realArr v1, realArr v2, int n) {
    real out;
    Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(int idx, real &prod) {
        prod += v1(idx) * v2(idx);
    }, out);
    return out;
}

void kxpy(realArr out, real k, realArr x, realArr y, int n) {
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idx) {
        out(idx) = k * x(idx) + y(idx);
    });
}

void kxpy2(realArr out, real k, realArr x, double y, int n) {
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idx) {
        out(idx) = k * x(idx) + y;
    });
}

void mul_MV(realArr out, RTMatrix &rtA, realArr x) {
    Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
        out(idx) = 0.0;
        for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); ++icol) {
            out(idx) += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
        }
    });
}

void get_residual(RTMatrix &rtA) {
    Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
        rtA.rt_r(idx) = rtA.rt_rhs(idx);
        for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); ++icol) {
            rtA.rt_r(idx) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
        }
    });
}

// Modified Jacobi Preconditioner with input/output parameters
void Jacobi_Precondition(RTMatrix &rtA, realArr r, realArr z) {
    int n = rtA.rt_nrow;
    Kokkos::parallel_for("Jacobi_Precondition", n, KOKKOS_LAMBDA(int i) {
        double diag_val = 0.0;
        for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i+1); ++j) {
            if (rtA.rt_ind(j) == i) {
                diag_val = rtA.rt_val(j);
                break;
            }
        }
        z(i) = (diag_val != 0.0) ? r(i)/diag_val : 0.0;
    });
}

// GMRES with Jacobi Preconditioning
int Jacobi_GMRES_Solve(RTMatrix &rtA ) {

   

	
	int m_restart = 30;
    const int n = rtA.rt_nrow;
    const real eps_min = 1e-6;
    int iter = 0;
    real eps = 1.0;

    // realArr rt_z("z", n), rt_w("w", n);
	realArr &rt_z = rtA.rt_z;
	realArr &rt_w = rtA.rt_w;

    Kokkos::View<real**> v("v", m_restart+1, n);
    Kokkos::View<real*> c("c", m_restart), s("s", m_restart), g("g", m_restart+1);
    Kokkos::View<real**> H("H", m_restart+1, m_restart);

    while (iter < iter_max && eps > eps_min) {
        // Initial residual and preconditioning
        get_residual(rtA);
        Jacobi_Precondition(rtA, rtA.rt_r, rt_z);
        real beta = sqrt(dot(rt_z, rt_z, n));
        if (beta < eps_min) break;

		    // 使用 Kokkos::parallel_for 进行元素-wise 除法
			Kokkos::parallel_for(rt_z.extent(0), KOKKOS_LAMBDA(const int i) {
				rt_z(i) /= beta;
			});

        // Initialize v[0] and g
        Kokkos::deep_copy(Kokkos::subview(v, 0, Kokkos::ALL()), rt_z);
        g(0) = beta;
        for (int i=1; i<=m_restart; ++i) g(i) = 0.0;

        int k;
        for (k=0; k<m_restart; ++k) {
            // Compute w = M^{-1}*A*v[k]
            mul_MV(rt_w, rtA, Kokkos::subview(v, k, Kokkos::ALL()));
            Jacobi_Precondition(rtA, rt_w, rt_w);

            // Modified Gram-Schmidt
            for (int i=0; i<=k; ++i) {
                H(i,k) = dot(Kokkos::subview(v, i, Kokkos::ALL()), rt_w, n);
                kxpy(rt_w, -H(i,k), Kokkos::subview(v, i, Kokkos::ALL()), rt_w, n);
            }
            H(k+1,k) = sqrt(dot(rt_w, rt_w, n));
            if (H(k+1,k) != 0.0)
                kxpy2(Kokkos::subview(v, k+1, Kokkos::ALL()), 1.0/H(k+1,k), rt_w, 0.0, n);

            // Apply previous Givens rotations
            for (int i=0; i<k; ++i) {
                real temp =  c(i)*H(i,k) + s(i)*H(i+1,k);
                H(i+1,k) = -s(i)*H(i,k) + c(i)*H(i+1,k);
                H(i,k) = temp;
            }

            // Compute new Givens rotation
            real gamma = sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k));
            if (gamma == 0) gamma = 1e-16;
            c(k) = H(k,k)/gamma;
            s(k) = H(k+1,k)/gamma;
            H(k,k) = gamma;
            
            // Update g
            real gk = c(k)*g(k) + s(k)*g(k+1);
            g(k+1) = -s(k)*g(k) + c(k)*g(k+1);
            g(k) = gk;

            eps = fabs(g(k+1));
            if (eps < eps_min) { k++; break; }
        }

        // Solve upper triangular system H(0:k,0:k) y = g(0:k)
        Kokkos::View<real*> y("y", k);
        for (int i=k-1; i>=0; --i) {
            y(i) = g(i);
            for (int j=i+1; j<k; ++j) y(i) -= H(i,j)*y(j);
            y(i) /= H(i,i);
        }

        // Update solution
        Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idx) {
            real delta = 0.0;
            for (int j=0; j<k; ++j) delta += v(j, idx)*y(j);
            rtA.rt_x(idx) += delta;
        });

        iter += k;
        eps = fabs(g(k));
    }
    return iter;
}



};
#endif
