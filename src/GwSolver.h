/*
	Solvers for the Ax=b system
*/
#ifndef _GW_SOLVER_H_
#define _GW_SOLVER_H_

#include "GwDomain.h"
#include "GwMatrix.h"

class GwSolver {

public:
	int iter, iter_max, gsteps, nrow, nnz;
	real eps, eps_min;

	inline void init(GwMatrix &A, GwDomain &gdom)	{
		nrow = A.nrow;	nnz = A.nnz;
		iter_max = 40;
		gsteps = 20;
		iter_max = iter_max * gsteps;
		eps_min = 1e-7;
	}

	/*
		----------------------------------------------------------
		----------------------------------------------------------
							Iterative Solvers
		----------------------------------------------------------
		----------------------------------------------------------
	*/

	/*
    	CG Solver
    */
    int cg(GwMatrix &A, GwDomain &gdom)	{
		real rho, rhoOld, alpha, beta;
		rhoOld = 0.0;
		decompose(A);
		// initialize residual
		get_residual(A);
		iter = 0;	eps = 1.0;
		while (iter < iter_max & eps > eps_min)	{
			precJACO(A);
			rho = dot(A.r, A.z, A.nrow);
			if (iter == 0)	{
				beta = 0.0;
				Kokkos::deep_copy(A.p, A.z);
			}
			else {
				beta = rho / rhoOld;
				kxpy(A.p, beta, A.p, A.z, A.nrow);
			}
			mul_MV(A.q, A, A.p);
			alpha = rho / dot(A.q, A.p, A.nrow);
			update_X(A, alpha);
			kxpy(A.r, -alpha, A.q, A.r, A.nrow);
			rhoOld = rho;
			eps = pow(dot(A.r, A.r, A.nrow),0.5);
			iter += 1;
    	}
		// std::cerr << "      >> CG solver converges in " << iter << " iterations! eps = " << eps <<"\n";
    	return iter;
    }

	/*
    	Jacobi Preconditioner
    */
    void precJACO(GwMatrix A)	{
		Kokkos::parallel_for( A.nrow , KOKKOS_LAMBDA(int idom) {
			A.z(idom) = A.r(idom) / A.diag(idom);
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
	void get_residual(GwMatrix &A)	{
		Kokkos::parallel_for( A.nrow , KOKKOS_LAMBDA(int idom) {
			int icol;
			A.r(idom) = A.rhs(idom);
			for (icol = A.ptr(idom); icol < A.ptr(idom+1); icol++)	{
				A.r(idom) -= A.val(icol) * A.x(A.ind(icol));
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
	void mul_MV(realArr out, GwMatrix A, realArr x)	{
		Kokkos::parallel_for( A.nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			out(idx) = 0.0;
			for (icol = A.ptr(idx); icol < A.ptr(idx+1); icol++)	{
				out(idx) += A.val(icol) * x(A.ind(icol));
			}
		});
	}

	/*
		Update solution
	*/
	void update_X(GwMatrix A, real alpha)	{
		Kokkos::parallel_for( A.nrow , KOKKOS_LAMBDA(int idx) {
			A.x(idx) += alpha * A.p(idx);
		});
	}

	/*
		Get diagonal of Matrix
	*/
	void decompose(GwMatrix A)	{
		Kokkos::parallel_for(A.nrow , KOKKOS_LAMBDA(int idx) {
			int icol;
			A.diag(idx) = 0.0;
			for (icol = A.ptr(idx); icol < A.ptr(idx+1); icol++)	{
				if (A.ind(icol) == idx)	{A.diag(idx) = A.val(icol);}
			}
		});
	}

};
#endif
