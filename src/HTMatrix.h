/*
	The sparse matrix
*/
#ifndef _HT_MATRIX_H_
#define _HT_MATRIX_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"

class HTMatrix {

public:
	int ht_nrow, ht_ncol, ht_nnz, ht_nx, ht_ny, ht_nz;
	intArr ht_ptr, ht_ind, ht_ptrT, ht_indT;
	realArr ht_val, ht_diag, ht_rhs, ht_x, ht_valT, ht_lt, ht_ut;
	// views for cg solver
	realArr ht_r, ht_z, ht_p, ht_q;

	// initialize
	// GwMatrix(GwDomain &gdom)	{
	void init(GwDomain &gdom)	{
		int ii, jj, kk, ndom = gdom.nCell;
		ht_nx = gdom.nx;	ht_ny = gdom.ny;	ht_nz = gdom.nz;
		ht_nrow = ndom;//矩阵行数是nxhc*nyhc*nzhc
		ht_ncol = ndom;//矩阵列数是nxhc*nyhc*nzhc
		// Number of non-zeros for 1D-z, 2D-xz, 3D-xyz simulations
		if (ht_nx == 1 & ht_ny == 1 & ht_nz > 1)	{
			ht_nnz = ht_nrow + (ht_nz-2)*2 + 2;
		}
		else if (ht_nx > 1 & ht_ny == 1 & ht_nz > 1)	{
			ht_nnz = ht_nrow + 4*2 + ((ht_nx-2)*2 + (ht_nz-2)*2)*3 + (ht_nx-2)*(ht_nz-2)*4;
		}
		else if (ht_nx > 1 & ht_ny > 1 & ht_nz > 1)	{
			ht_nnz = ht_nrow + 3*8 + 4*((ht_nx-2)*4 + (ht_ny-2)*4 + (ht_nz-2)*4) +
		    	5*((ht_nx-2)*(ht_nz-2)*2 + (ht_nx-2)*(ht_ny-2)*2 + (ht_ny-2)*(ht_nz-2)*2) +
		    	6*((ht_nx-2)*(ht_ny-2)*(ht_nz-2));
		}
		else	{printf("ERROR : HTDomain must be 1D-z, 2D-xz or 3D-xyz!\n");}
		ht_ptr = intArr("ht_ptr", ht_nrow+1);		ht_ind = intArr("ht_ind", ht_nnz);
		ht_val = realArr("ht_val", ht_nnz);			ht_diag = realArr("ht_diag", ht_nrow);
		ht_rhs = realArr("ht_rhs", ht_nrow);			ht_x = realArr("ht_x", ht_nrow);
		ht_r = realArr("ht_r", ht_nrow);				ht_z = realArr("ht_z", ht_nrow);
		ht_p = realArr("ht_p", ht_nrow);				ht_q = realArr("ht_q", ht_nrow);
		// Get ptr for the CRS matrix
		for (int idx = 0; idx < ht_nrow+1; idx++)	{
			gdom.unpackIndices(idx, kk, jj, ii);
			// gdom.unpackIndicesGw(idx, nz, ny, nx, kk, jj, ii);
			if (idx == ht_nrow) {ht_ptr(idx) = ht_nnz;}
			else	{ht_ptr(idx) = get_irow(ii, jj, kk, ht_nx, ht_ny, ht_nz);}
		}

	}

	// get starting index
	inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
		int ht_nrowi = 0, ht_nlayer1 = 0, ht_nlayerk = 0, ht_ncol1 = 0, ht_ncolj = 0, ht_irow;
	    if (k == 0)   {
	        if (j == 0)   {
	            // 3D
	            if (ny > 1) {if (i > 0)   {ht_nrowi = 4 + (i-1)*5;}}
	            // 2D x-z
	            else {if (i > 0)   {ht_nrowi = 3 + (i-1)*4;}}
	        }
	        else {
	            // 3D
	            ht_ncol1 = 4*2 + 5*(nx-2);
	            ht_ncolj = (j-1)*5*2 + (j-1)*6*(nx-2);
	            if (j == ny-1)   {if (i > 0)   {ht_nrowi = 4 + (i-1)*5;}}
	            else {if (i > 0)   {ht_nrowi = 5 + (i-1)*6;}}
	        }
	    }
	    else {
	        if (nx > 1) {
	            // 3D
	            if (ny > 1) {
	                ht_nlayer1 = 4*4 + 5*2*(nx-2+ny-2) + 6*(nx-2)*(ny-2);
	                ht_nlayerk = (k-1)*(5*4 + 6*2*(nx-2+ny-2) + 7*(nx-2)*(ny-2));
	            }
	            // 2D
	            else {
	                ht_nlayer1 = 3*2 + 4*(nx-2);
	                ht_nlayerk = (k-1)*(4*2 + 5*(nx-2));
	            }
	        }
	        else {
	            // 1D
	            ht_nlayer1 = 2;    ht_nlayerk = (k-1)*3;
	        }
	        if (k == nz-1)    {
	            if (j == 0)   {
	                if (ny > 1) {if (i > 0)   {ht_nrowi = 4 + (i-1)*5;}}
	                else {if (i > 0)   {ht_nrowi = 3 + (i-1)*4;}}
	            }
	            else {
	                ht_ncol1 = 4*2 + 5*(nx-2);
	                ht_ncolj = (j-1)*5*2 + (j-1)*6*(nx-2);
	                if (j == ny-1)   {if (i > 0)   {ht_nrowi = 4 + (i-1)*5;}}
	                else {if (i > 0)   {ht_nrowi = 5 + (i-1)*6;}}
	            }
	        }
	        else {
	            if (j == 0)   {
	                if (ny > 1) {if (i > 0)   {ht_nrowi = 5 + (i-1)*6;}}
	                else {if (i > 0)   {ht_nrowi = 4 + (i-1)*5;}}
	            }
	            else {
	                ht_ncol1 = 5*2 + 6*(nx-2);
	                ht_ncolj = (j-1)*6*2 + (j-1)*7*(nx-2);
	                if (j == ny-1)   {if (i > 0)   {ht_nrowi = 5 + (i-1)*6;}}
	                else {if (i > 0)   {ht_nrowi = 6 + (i-1)*7;}}
	            }
	        }

	    }
	    ht_irow = ht_nlayer1 + ht_nlayerk + ht_ncol1 + ht_ncolj + ht_nrowi;
	    return ht_irow;

	}


};

#endif
