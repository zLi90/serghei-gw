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
template <typename execution_space>
class RTSolver {

	typedef typename KokkosSparse::CrsMatrix<default_scalar, default_lno_t, execution_space, void, default_size_type> crsMat_t;
	typedef typename crsMat_t::StaticCrsGraphType graph_t;
	typedef typename graph_t::row_map_type::non_const_type lno_view_t;
	typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
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

	inline void init(RTMatrix &rtA, GwDomain &gdom)	{
		nrow = rtA.rt_nrow;	nnz = rtA.rt_nnz;
		iter_max = 40;
		gsteps = 20;
		iter_max = iter_max * gsteps;
		eps_min = 1e-7;
		#if SERGHEI_KOKKOSKERNELS_SOLVER
		rt_ptr = lno_view_t("rt_ptr", rtA.rt_nrow+1);
		rt_ind = lno_nnz_view_t("rt_ind", rtA.rt_nnz);
		rt_val = scalar_view_t("rt_val", rtA.rt_nnz);
		rt_rhs = scalar_view_t("rt_rhs", rtA.rt_nrow);
		rt_vecx = scalar_view_t("rt_vecx", rtA.rt_nrow);
		rt_diag = scalar_view_t("rt_diag", rtA.rt_nrow);
		Kokkos::deep_copy (rt_ptr, rtA.rt_ptr);
		#endif
	}

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
	void kkpcg(RTMatrix &rtA)
	{
		bool usePreconditioner = 1;
		//usePreconditioner 变量标志是否使用预条件器，默认为使用
		const unsigned cg_iteration_limit = 1000000;
		//cg_iteration_limit 定义了PCG算法的最大迭代次数。
		const double   cg_iteration_tolerance     = 1e-8 ;
		//cg_iteration_tolerance 定义了PCG算法的收敛容限。

		decompose(rtA);
		// decompose(A) 对矩阵进行分解，
		Kokkos::deep_copy (rt_ind, rtA.rt_ind);
		Kokkos::deep_copy (rt_val, rtA.rt_val);
		Kokkos::deep_copy (rt_rhs, rtA.rt_rhs);
		Kokkos::deep_copy (rt_diag, rtA.rt_diag);

		crsMat_t rtmatA = crsMat_t("rtmatA", rtA.rt_nrow, rtA.rt_nrow, rtA.rt_nnz, rt_val, rt_ptr, rt_ind);

		// initialize KokkosKernels solver
		KokkosKernels::Experimental::Example::CGSolveResult cg_result ;
		typedef KokkosKernels::Experimental::KokkosKernelsHandle
		  < size_type, lno_t, scalar_t, execution_space, execution_space, execution_space > KernelHandle;
		KernelHandle rt_kh;
		rt_kh.create_gs_handle();
		//Kokkos::Impl::Timer timer1;
		KokkosKernels::Experimental::Example::pcgsolve(rt_kh, rtmatA, rt_rhs, rt_vecx, rt_diag
		    , cg_iteration_limit, cg_iteration_tolerance, &cg_result, usePreconditioner);
		Kokkos::fence();
		// 调用 pcgsolve 函数进行PCG求解，传入了求解句柄 kh、矩阵 matA、右侧向量 rhs、
		// 解向量 vecx、对角元素 diag、迭代次数限制、迭代收敛容限等参数，
		// 并将求解结果保存在 cg_result 中。
		//solve_time = timer1.seconds();
		//std::cout  << "DEFAULT SOLVE: " << "(P)CG_NUM_ITER = [" << cg_result.iteration << "], " << "RESIDUAL = [" << cg_result.norm_res << "]"<< std::endl ;
		rt_kh.destroy_gs_handle();
		// 调用 destroy_gs_handle() 方法销毁求解句柄

		Kokkos::deep_copy (rtA.rt_x, rt_vecx);
		// 将解向量 vecx 复制回 A 中
	}

	/*
    	CG Solver
    */
    int cg(RTMatrix &rtA, GwDomain &gdom)	{
		real rho, rhoOld, alpha, beta;
		rhoOld = 0.0;
		decompose(rtA);
		// initialize residual
		get_residual(rtA);
		iter = 0;	eps = 1.0;
		while (iter < iter_max & eps > eps_min)	{
			precJACO(rtA);
			rho = dot(rtA.rt_r, rtA.rt_z, rtA.rt_nrow);
			if (iter == 0)	{
				beta = 0.0;
				Kokkos::deep_copy(rtA.rt_p, rtA.rt_z);
			}
			else {
				beta = rho / rhoOld;
				kxpy(rtA.rt_p, beta, rtA.rt_p, rtA.rt_z, rtA.rt_nrow);
			}
			mul_MV(rtA.rt_q, rtA, rtA.rt_p);
			alpha = rho / dot(rtA.rt_q, rtA.rt_p, rtA.rt_nrow);
			update_X(rtA, alpha);
			kxpy(rtA.rt_r, -alpha, rtA.rt_q, rtA.rt_r, rtA.rt_nrow);
			rhoOld = rho;
			eps = pow(dot(rtA.rt_r, rtA.rt_r, rtA.rt_nrow),0.5);
			iter += 1;
    	}
		// std::cerr << "      >> CG solver converges in " << iter << " iterations! eps = " << eps <<"\n";
    	return iter;
    }

	/*
    	Jacobi Preconditioner
    */
    void precJACO(RTMatrix rtA)	{
		Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idom) {
			rtA.rt_z(idom) = rtA.rt_r(idom) * rtA.rt_diag(idom);
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
	void get_residual(RTMatrix &rtA)	{
		Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idom) {
			int icol;
			rtA.rt_r(idom) = rtA.rt_rhs(idom);
			for (icol = rtA.rt_ptr(idom); icol < rtA.rt_ptr(idom+1); icol++)	{
				rtA.rt_r(idom) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
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
	void mul_MV(realArr out, RTMatrix rtA, realArr x)	{
		Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			out(idx) = 0.0;
			for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++)	{
				out(idx) += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
			}
		});
	}

	/*
		Update solution
	*/
	void update_X(RTMatrix rtA, real alpha)	{
		Kokkos::parallel_for( rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
			rtA.rt_x(idx) += alpha * rtA.rt_p(idx);
		});
	}

	/*
		Get diagonal of Matrix
	*/
	void decompose(RTMatrix rtA)	{
		Kokkos::parallel_for(rtA.rt_nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			rtA.rt_diag(idx) = 0.0;
			for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++)	{
				if (rtA.rt_ind(icol) == idx)	{rtA.rt_diag(idx) = 1.0 / rtA.rt_val(icol);}
			}
		});
	}

};
#endif
