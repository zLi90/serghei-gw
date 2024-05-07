/*
	The sparse matrix
*/
#ifndef _RT_MATRIX_H_
#define _RT_MATRIX_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"

class RTMatrix {

public:
	int rt_nrow, rt_ncol, rt_nnz, rt_nx, rt_ny, rt_nz;
	intArr rt_ptr, rt_ind, rt_ptrT, rt_indT;
	realArr rt_val, rt_diag, rt_rhs, rt_x, rt_valT, rt_lt, rt_ut;
	// views for cg solver
	realArr rt_r, rt_z, rt_p, rt_q;

	// initialize
	// GwMatrix(GwDomain &gdom)	{
	void init(GwDomain &gdom)	{
		int ii, jj, kk, ndom = gdom.nCell;
		rt_nx = gdom.nx;	rt_ny = gdom.ny;	rt_nz = gdom.nz;
		rt_nrow = ndom;//矩阵行数是nxhc*nyhc*nzhc
		rt_ncol = ndom;//矩阵列数是nxhc*nyhc*nzhc
		// Number of non-zeros for 1D-z, 2D-xz, 3D-xyz simulations
		if (rt_nx == 1 & rt_ny == 1 & rt_nz > 1)	{
			rt_nnz = rt_nrow + (rt_nz-2)*2 + 2;
		}
		else if (rt_nx > 1 & rt_ny == 1 & rt_nz > 1)	{
			rt_nnz = rt_nrow + 4*2 + ((rt_nx-2)*2 + (rt_nz-2)*2)*3 + (rt_nx-2)*(rt_nz-2)*4;
		}
		else if (rt_nx > 1 & rt_ny > 1 & rt_nz > 1)	{
			rt_nnz = rt_nrow + 3*8 + 4*((rt_nx-2)*4 + (rt_ny-2)*4 + (rt_nz-2)*4) +
		    	5*((rt_nx-2)*(rt_nz-2)*2 + (rt_nx-2)*(rt_ny-2)*2 + (rt_ny-2)*(rt_nz-2)*2) +
		    	6*((rt_nx-2)*(rt_ny-2)*(rt_nz-2));
		}
		else	{printf("ERROR : RTDomain must be 1D-z, 2D-xz or 3D-xyz!\n");}

		rt_ptr = intArr("rt_ptr", rt_nrow+1);		rt_ind = intArr("rt_ind", rt_nnz);
		rt_val = realArr("rt_val", rt_nnz);			rt_diag = realArr("rt_diag", rt_nrow);
		rt_rhs = realArr("rt_rhs", rt_nrow);			rt_x = realArr("rt_x", rt_nrow);
		rt_r = realArr("rt_r", rt_nrow);				rt_z = realArr("rt_z", rt_nrow);
		rt_p = realArr("rt_p", rt_nrow);				rt_q = realArr("rt_q", rt_nrow);
		// Get ptr for the CRS matrix
		for (int idx = 0; idx < rt_nrow+1; idx++)	{
			gdom.unpackIndices(idx, kk, jj, ii);
			// gdom.unpackIndicesGw(idx, nz, ny, nx, kk, jj, ii);
			if (idx == rt_nrow) {rt_ptr(idx) = rt_nnz;}
			else	{rt_ptr(idx) = get_irow(ii, jj, kk, rt_nx, rt_ny, rt_nz);}
		}

	}

	// get starting index
	inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
		int rt_nrowi = 0, rt_nlayer1 = 0, rt_nlayerk = 0, rt_ncol1 = 0, rt_ncolj = 0, rt_irow;
	    if (k == 0)   {
	        if (j == 0)   {
	            // 3D
	            if (ny > 1) {if (i > 0)   {rt_nrowi = 4 + (i-1)*5;}}
	            // 2D x-z
	            else {if (i > 0)   {rt_nrowi = 3 + (i-1)*4;}}
	        }
	        else {
	            // 3D
	            rt_ncol1 = 4*2 + 5*(nx-2);
	            rt_ncolj = (j-1)*5*2 + (j-1)*6*(nx-2);
	            if (j == ny-1)   {if (i > 0)   {rt_nrowi = 4 + (i-1)*5;}}
	            else {if (i > 0)   {rt_nrowi = 5 + (i-1)*6;}}
	        }
	    }
	    else {
	        if (nx > 1) {
	            // 3D
	            if (ny > 1) {
	                rt_nlayer1 = 4*4 + 5*2*(nx-2+ny-2) + 6*(nx-2)*(ny-2);
	                rt_nlayerk = (k-1)*(5*4 + 6*2*(nx-2+ny-2) + 7*(nx-2)*(ny-2));
	            }
	            // 2D
	            else {
	                rt_nlayer1 = 3*2 + 4*(nx-2);
	                rt_nlayerk = (k-1)*(4*2 + 5*(nx-2));
	            }
	        }
	        else {
	            // 1D
	            rt_nlayer1 = 2;    rt_nlayerk = (k-1)*3;
	        }
	        if (k == nz-1)    {
	            if (j == 0)   {
	                if (ny > 1) {if (i > 0)   {rt_nrowi = 4 + (i-1)*5;}}
	                else {if (i > 0)   {rt_nrowi = 3 + (i-1)*4;}}
	            }
	            else {
	                rt_ncol1 = 4*2 + 5*(nx-2);
	                rt_ncolj = (j-1)*5*2 + (j-1)*6*(nx-2);
	                if (j == ny-1)   {if (i > 0)   {rt_nrowi = 4 + (i-1)*5;}}
	                else {if (i > 0)   {rt_nrowi = 5 + (i-1)*6;}}
	            }
	        }
	        else {
	            if (j == 0)   {
	                if (ny > 1) {if (i > 0)   {rt_nrowi = 5 + (i-1)*6;}}
	                else {if (i > 0)   {rt_nrowi = 4 + (i-1)*5;}}
	            }
	            else {
	                rt_ncol1 = 5*2 + 6*(nx-2);
	                rt_ncolj = (j-1)*6*2 + (j-1)*7*(nx-2);
	                if (j == ny-1)   {if (i > 0)   {rt_nrowi = 5 + (i-1)*6;}}
	                else {if (i > 0)   {rt_nrowi = 6 + (i-1)*7;}}
	            }
	        }

	    }
	    rt_irow = rt_nlayer1 + rt_nlayerk + rt_ncol1 + rt_ncolj + rt_nrowi;
	    return rt_irow;

	}


};

#endif
