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
	int rt_nrow, rt_ncol, rt_nnz, rt_nx, rt_ny, rt_nz,  rt_cg_iter;
	intArr rt_ptr, rt_ind, rt_ptrT, rt_indT;
	realArr rt_val, rt_diag, rt_rhs, rt_x, rt_valT, rt_lt, rt_ut;
	// views for cg solver
	realArr rt_r, rt_z, rt_p, rt_q,rt_w;


	// initialize
	// GwMatrix(GwDomain &gdom)	{
	void init(GwDomain &gdom)	{
		int ii, jj, kk, ndom = gdom.nCell;
		rt_nx = gdom.nx;	rt_ny = gdom.ny;	rt_nz = gdom.nz;
		//! std::cout << "rt_nx = " << rt_nx << " rt_ny = " << rt_ny << " rt_nz = " << rt_nz << std::endl;
		rt_nrow = ndom;//矩阵行数是nxhc*nyhc*nzhc
		rt_ncol = ndom;//矩阵列数是nxhc*nyhc*nzhc
		// Number of non-zeros for 1D-z, 2D-xz, 3D-xyz simulations
		if (rt_nx == 1 & rt_ny == 1 & rt_nz > 1)	{
			rt_nnz = rt_nrow + (rt_nz-2)*2 + 2;
		}
		else if (rt_nx > 1 & rt_ny == 1 & rt_nz > 1)	{
			// rt_nnz = rt_nrow + 4*2 + ((rt_nx-2)*2 + (rt_nz-2)*2)*3 + (rt_nx-2)*(rt_nz-2)*4;
			//todo 修改
			rt_nnz = rt_nrow + 4*3 + ((rt_nx-2)*2 + (rt_nz-2)*2)*5 + (rt_nx-2)*(rt_nz-2)*8;
			// rt_nnz = 165;
		}
		else if (rt_nx > 1 & rt_ny > 1 & rt_nz > 1)	{
			// rt_nnz = rt_nrow + 3*8 + 4*((rt_nx-2)*4 + (rt_ny-2)*4 + (rt_nz-2)*4) +
		    // 	5*((rt_nx-2)*(rt_nz-2)*2 + (rt_nx-2)*(rt_ny-2)*2 + (rt_ny-2)*(rt_nz-2)*2) +
		    // 	6*((rt_nx-2)*(rt_ny-2)*(rt_nz-2));
            // rt_nnz = rt_nrow + 7*8 + 
            //     11*((rt_nx-2)*4 + (rt_ny-2)*4 + (rt_nz-2)*4) +
		    // 	17*((rt_nx-2)*(rt_nz-2)*2 + (rt_nx-2)*(rt_ny-2)*2 + (rt_ny-2)*(rt_nz-2)*2) +
		    // 	26*((rt_nx-2)*(rt_ny-2)*(rt_nz-2));

            rt_nnz = rt_nrow 
             + 8*6  // 8个角点，每个角点有7个非零元素
             + ((rt_nx-2)*4 + (rt_ny-2)*4 + (rt_nz-2)*4)*9  // 边上的点，每条边有 (rt_nx-2) 个点，每个点有10个非零元素
             + ((rt_nx-2)*(rt_ny-2)*2 + (rt_nx-2)*(rt_nz-2)*2 + (rt_ny-2)*(rt_nz-2)*2)*13  // 面上的点，每个面有 (rt_nx-2)*(rt_ny-2) 个点，每个点有14个非零元素
             + (rt_nx-2)*(rt_ny-2)*(rt_nz-2)*18;  // 内部的点，每个点有19个非零元素（19点差分格式）
		}
		else	{printf("ERROR : RTDomain must be 1D-z, 2D-xz or 3D-xyz!\n");}
//rtA.rt_ptr数组，存储了稀疏矩阵中每一行非零元素的起始索引
		rt_ptr = intArr("rt_ptr", rt_nrow+1);		rt_ind = intArr("rt_ind", rt_nnz);
		rt_val = realArr("rt_val", rt_nnz);			rt_diag = realArr("rt_diag", rt_nrow);
		rt_rhs = realArr("rt_rhs", rt_nrow);			rt_x = realArr("rt_x", rt_nrow);
		rt_r = realArr("rt_r", rt_nrow);				rt_z = realArr("rt_z", rt_nrow);
		rt_p = realArr("rt_p", rt_nrow);				rt_q = realArr("rt_q", rt_nrow);
        rt_w = realArr("rt_w", rt_nrow);
		// Get ptr for the CRS matrix
		for (int idx = 0; idx < rt_nrow+1; idx++)	{
			gdom.unpackIndices(idx, kk, jj, ii);
			// gdom.unpackIndicesGw(idx, nz, ny, nx, kk, jj, ii);
			if (idx == rt_nrow) 
			{
				rt_ptr(idx) = rt_nnz;
				}
			else	
			{
				rt_ptr(idx) = get_irow(ii, jj, kk, rt_nx, rt_ny, rt_nz);
				}
// std::cout<<"------rt_ptr(idx)------\n"<<rt_ptr(idx)<<std::endl;
// std::cout<<"------rt_nnz------\n"<<rt_nnz<<std::endl;			
		}
	}


	// get starting index
	//计算在三维网格中，给定位置 (i, j, k) 对应的行在ptr数组（用于存储非零元素的索引）中的起始位置
// inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
	inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
    int rt_nrowi = 0, rt_nlayer1 = 0, rt_nlayerk = 0, rt_ncol1 = 0, rt_ncolj = 0, rt_irow;

    int JDY = 7;
    int BDY = 10;
    int MDY = 14;
    int NDY = 19;

    //!矩阵第一层
    if (k == 0)
    { 
        //!第一层第一列
        if (j == 0)
        {
            //! 3D
            if (ny > 1)
            {
                if (i > 0)
                {
                    //rt_nrowi = 4 + (i - 1) * 5;
                    rt_nrowi = 7 + (i - 1) * 10;//?
                    rt_irow = JDY + (i - 1) * BDY;
                }
            }
            //! 2D x-z
            else
            {
                //!第一层第一列非第一行
                if (i > 0)
                {
                    // rt_nrowi = 3 + (i - 1) * 4;
                    // todo 修改
                    rt_nrowi = 4 + (i - 1) * 6;
                }
            }
        }
        //!第一层非第一列
        else
        {
            // 3D
            // rt_ncol1 = 4 * 2 + 5 * (nx - 2);
            // rt_ncolj = (j - 1) * 5 * 2 + (j - 1) * 6 * (nx - 2);

            // rt_ncol1 = 7 * 4 + 10 * (nx - 2);    //? 第一层第一列的总非零元素数量
            rt_ncol1 = JDY*2 + BDY*(nx-2);
            // rt_ncolj = (j - 1) * 10 * 4 + (j - 1) * 14 * (nx - 2);     //? 第一层非第一列的总非零元素数量 
            rt_ncolj = (j - 1) * BDY*2 + (j - 1) * MDY*(nx-2);
            if (j == ny - 1)
            {
                if (i > 0)
                {
                    // rt_nrowi = 4 + (i - 1) * 5;
                    // rt_nrowi = 7 + (i - 1) * 10;    //? 最后一列的非零元素数量
                    rt_nrowi = JDY + (i - 1) * BDY;
                }
            }
            else
            {
                if (i > 0)
                {
                    // rt_nrowi = 5 + (i - 1) * 6;
                    // rt_nrowi = 10 + (i - 1) * 14;    //? 中间列的非零元素数量
                    rt_nrowi = BDY + (i - 1) * MDY;
                }
            }
        }
    }
    //!矩阵非第一层
    else
    {
        if (nx > 1)
        {
            // 3D
            if (ny > 1)
            {
                // rt_nlayer1 = 4 * 4 + 5 * 2 * (nx - 2 + ny - 2) + 6 * (nx - 2) * (ny - 2);
                // rt_nlayerk = (k - 1) * (5 * 4 + 6 * 2 * (nx - 2 + ny - 2) + 7 * (nx - 2) * (ny - 2));
                // rt_nlayer1 = 7 * 8 + 10 * 4 * (nx - 2 + ny - 2) + 14 * 2 * ((nx - 2) * (ny - 2));   //? 第一层的总非零元素数量
                rt_nlayer1 = 4*JDY + BDY*2*(nx-2+ny-2) + MDY*(nx-2)*(ny-2);
                // rt_nlayerk = (k - 1) * (10 * 8 + 14 * 4 * (nx - 2 + ny - 2) + 19 * 2 * ((nx - 2) * (ny - 2)));   //? 中间层的总非零元素数量            
                rt_nlayerk = (k - 1) * (JDY*4 + MDY*2*(nx-2+ny-2) + NDY*(nx-2)*(ny-2));
            }
            // 2D
            else
            {
                // rt_nlayer1 = 3 * 2 + 4 * (nx - 2);
                // rt_nlayerk = (k - 1) * (4 * 2 + 5 * (nx - 2));
                // todo 修改
                rt_nlayer1 = 4 * 2 + 6 * (nx - 2);
                rt_nlayerk = (k - 1) * (6 * 2 + 9 * (nx - 2));
            }
        }
        //! nx=1
        else
        {
            // 1D
            rt_nlayer1 = 2;
            rt_nlayerk = (k - 1) * 3;
        }

        //!矩阵最后一层
        if (k == nz - 1)
        {
            //!最后一层第一列
            if (j == 0)
            {
                //! 3D
                if (ny > 1)
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 4 + (i - 1) * 5;
                        // rt_nrowi = 7 + (i - 1) * 10;   //? 最后一层第一列的非零元素数量
                        rt_nrowi = 4 + (i - 1) * 5;
                    }
                }
                //! 2D x-z
                else
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 3 + (i - 1) * 4;
						//!todo DXY修改
                        rt_nrowi = 4 + (i - 1) * 6;
                    }
                }
            }
            //!最后一层非第一列
            else
            {
                // rt_ncol1 = 4 * 2 + 5 * (nx - 2);
                // rt_ncolj = (j - 1) * 5 * 2 + (j - 1) * 6 * (nx - 2);

                // rt_ncol1 = 7 * 4 + 10 * (nx - 2);  //? 最后一层第一列的总非零元素数量
                rt_ncol1 = JDY*2 + BDY*(nx-2);
                // rt_ncolj = (j - 1) * 10 * 4 + (j - 1) * 14 * (nx - 2);  //? 最后一层非第一列的总非零元素数量               
                rt_ncolj = (j - 1) * BDY*2 + (j - 1) * MDY*(nx-2);
                if (j == ny - 1)
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 4 + (i - 1) * 5;
                        // rt_nrowi = 7 + (i - 1) * 10;    //? 最后一列的非零元素数量
                        rt_nrowi = JDY + (i - 1) * BDY;
                    }
                }
                else
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 5 + (i - 1) * 6;
                        rt_nrowi = 10 + (i - 1) * 14;    //? 中间列的非零元素数量
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
            }
        }
        //!非第一层非最后一层，中间层
        else
        {
            //!中间层第一列
            if (j == 0)
            {
                // 3D
                if (ny > 1)
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 5 + (i - 1) * 6;
                        // rt_nrowi = 10 + (i - 1) * 14;   //? 中间层第一列的非零元素数量
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
                // 2D x-z
                else
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 4 + (i - 1) * 5;
						//!todo DXY修改
                        rt_nrowi = 6 + (i - 1) * 9;
                    }
                }
            }
            //!中间层非第一列
            else
            {
                // rt_ncol1 = 5 * 2 + 6 * (nx - 2);
                // rt_ncolj = (j - 1) * 6 * 2 + (j - 1) * 7 * (nx - 2);

                // rt_ncol1 = 10 * 4 + 14 * (nx - 2);   // ?中间层第一列的总非零元素数量
                rt_ncol1 = BDY*2 + MDY*(nx-2);
                // rt_ncolj = (j - 1) * 14 * 4 + (j - 1) * 19 * (nx - 2);   //? 中间层非第一列的总非零元素数量                
                rt_ncolj = (j - 1) * MDY*2 + (j - 1) * NDY*(nx-2);
                if (j == ny - 1)
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 5 + (i - 1) * 6;
                        // rt_nrowi = 10 + (i - 1) * 14;    //? 最后一列的非零元素数量
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
                else
                {
                    if (i > 0)
                    {
                        // rt_nrowi = 6 + (i - 1) * 7;
                        // rt_nrowi = 14 + (i - 1) * 19;     //? 中间列的非零元素数量
                        rt_nrowi = MDY + (i - 1) * NDY;
                    }
                }
            }
        }
    }
    rt_irow = rt_nlayer1 + rt_nlayerk + rt_ncol1 + rt_ncolj + rt_nrowi;
    return rt_irow;

	}


};

#endif
