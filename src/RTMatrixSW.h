/*
	The sparse matrix
*/
#ifndef _RT_MATRIX_SW_H_
#define _RT_MATRIX_SW_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "Domain.h"


class RTMatrixSW {

public:
	int rt_nrow, rt_ncol, rt_nnz, rt_nx, rt_ny, rt_cg_iter;
	intArr rt_ptr, rt_ind, rt_ptrT, rt_indT;
	realArr rt_val, rt_diag, rt_rhs, rt_x, rt_valT, rt_lt, rt_ut;
	// views for cg solver
	realArr rt_r, rt_z, rt_p, rt_q,rt_w;


	// initialize
	
	void init(Domain &dom)	{
		int ii, jj, ndom = dom.nCell;
		rt_nx = dom.nx;	rt_ny = dom.ny;
		//! std::cout << "rt_nx = " << rt_nx << " rt_ny = " << rt_ny << " rt_nz = " << rt_nz << std::endl;
		rt_nrow = ndom;//矩阵行数是nxhc*nyhc*nzhc
		rt_ncol = ndom;//矩阵列数是nxhc*nyhc*nzhc
		// Number of non-zeros for 1D-x, 1D-y,2D-xy simulations

        // 1D-x
        if (rt_nx > 1 & rt_ny == 1 )	{
			rt_nnz = rt_nrow + (rt_nx-2)*2 + 2;
		}
        // 1D-y
        else if (rt_nx == 1 & rt_ny > 1 )	{
			rt_nnz = rt_nrow + (rt_ny-2)*2 + 2;
		}
        // 2D-xy
        else if (rt_nx > 1 & rt_ny > 1 )	{
			rt_nnz = rt_nrow + 4*3 + ((rt_nx-2)*2 + (rt_ny-2)*2)*5 + (rt_nx-2)*(rt_ny-2)*8;
		}
		
		
		else	{printf("ERROR : RTDomainSW must be 1D-x, 1D-y or 2D-xy!\n");}
//rtA.rt_ptr数组，存储了稀疏矩阵中每一行非零元素的起始索引
		rt_ptr = intArr("rt_ptr", rt_nrow+1);		rt_ind = intArr("rt_ind", rt_nnz);
		rt_val = realArr("rt_val", rt_nnz);			rt_diag = realArr("rt_diag", rt_nrow);
		rt_rhs = realArr("rt_rhs", rt_nrow);			rt_x = realArr("rt_x", rt_nrow);
		rt_r = realArr("rt_r", rt_nrow);				rt_z = realArr("rt_z", rt_nrow);
		rt_p = realArr("rt_p", rt_nrow);				rt_q = realArr("rt_q", rt_nrow);
        rt_w = realArr("rt_w", rt_nrow);
		// Get ptr for the CRS matrix
		for (int idx = 0; idx < rt_nrow+1; idx++)	{
			dom.unpackIndices(idx, jj, ii);
			
			if (idx == rt_nrow) 
			{
				rt_ptr(idx) = rt_nnz;
				}
			else	
			{
				rt_ptr(idx) = get_irow(ii, jj, rt_nx, rt_ny);
				}
// std::cout<<"------rt_ptr(idx)------\n"<<rt_ptr(idx)<<std::endl;
// std::cout<<"------rt_nnz------\n"<<rt_nnz<<std::endl;			
		}
	}


	// get starting index
	//计算在三维网格中，给定位置 (i, j) 对应的行在ptr数组（用于存储非零元素的索引）中的起始位置
// inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
	inline int get_irow(int i, int j,  int nx, int ny) {
    int rt_nrowi = 0, rt_nlayer1 = 0, rt_nlayerk = 0, rt_ncol1 = 0, rt_ncolj = 0, rt_irow;

    int JDY = 3;
    int BDY = 4;
    int MDY = 5;
    int NDY = 5;


   


        //!第一列
        if (j == 0)
        {
            
            if ( nx > 1)
            {
                if (ny > 1) {//! 2D xy
                        if (i > 0)
                            {
                                rt_nrowi = JDY + (i - 1) * BDY;
                            }
                }
                else { //!1D-y
                rt_nrowi = (JDY-1) + (nx -1)*(BDY-1);

                }
            }
            else {//!1D-x
                rt_nrowi = JDY -1;

            }
        }
        //! 非第一列
        else  
            {
                if (nx>1){
                    rt_nlayer1 = JDY * 2 + BDY * (nx - 2);
                    rt_nlayerk = (i - 1) * BDY * 2 + (nx - 2)*(j - 1) * NDY;


                    if (j == ny-1){//!最后一列
                        if (i>0){
                            rt_nrowi = JDY + (i - 1) * BDY;
                        }

                    }
                    else {//!中间列
                        rt_nrowi = BDY + (i - 1)*NDY;

                    }
                
                
                }
                else {//! 1Dx,非第一列
                    rt_nrowi = (JDY - 1)+(j-1)*(BDY-1);


                }


            }
 
   
    
    rt_irow = rt_nlayer1 + rt_nlayerk + rt_ncol1 + rt_ncolj + rt_nrowi;
    return rt_irow;

	}


};

#endif
