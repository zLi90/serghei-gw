/*
    Surface Water Transport Module
*/
#ifndef _RT_FUNCTION_SW_H_
#define _RT_FUNCTION_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include <Kokkos_Core.hpp>
#include "define.h"

#include "const.h"
#include "BC.h"
#include "Domain.h"
#include "State.h"
#include "SArray.h"
#include "Indexing.h"
#include "RTStateSW.h"
#include <set>
#include <math.h>

class RTFunctionSW {
public:
    // 显式对流-弥散求解器（二维地表水版）
    KOKKOS_INLINE_FUNCTION 
    void transport(RTStateSW &rtsw, State &state, Domain &dom, ExternalBoundaries &ebc) {
        const real dx = dom.dx();
        const real dy = dom.dx();
        rtsw.dt = dom.dt;
        
        Kokkos::parallel_for("sw_transport", dom.nCellMem, 
            KOKKOS_LAMBDA(const int iGlob) {

            int i, j, id0, id1, id2, nxhc, nyhc;
            nxhc = dom.nx+2*hc;
            nyhc = dom.ny+2*hc;   
            unpackIndicesUniformGrid(iGlob, nyhc, nxhc, j, i);  
            

            // real h = 0.1; // 固定水深
            real h = state.h(iGlob);

            //复制上一时间步浓度
            // rtsw.c(iGlob, 0) = rtsw.c(iGlob, 1);

            //!上一时刻质量=上一时刻浓度*上一时刻水深*面积
            // rtsw.solute_mass(iGlob) = rtsw.c(iGlob, 0) * state.h(iGlob) * dx * dy; // 质量守恒
            rtsw.solute_mass(iGlob) = rtsw.c(iGlob, 0) * h * dx * dy; // 质量守恒


       
            
            // 跳过无效单元
            // if (state.isnodata(iGlob)) return;
            
            real adv_x = 0.0, adv_y = 0.0;
            real disp_x = 0.0, disp_y = 0.0;
            
            
            
            // 仅在水深大于阈值时计算
            if (h > state.hmin) {

                // printf("state.hmin = %f\n", state.hmin);
                // 计算流速分量
                real u = -state.hu(iGlob);
                real v = -state.hv(iGlob);

                //!debug
                // real u = 0; // 固定流速
                // real v = -0.01; // 固定流速
          
                // 计算弥散系数（简化版）
                real ALPHAL = 0.1;
                real ALPHAT = 0.01;
                real DMOLECULAR = 0; // 分子扩散系数
                
    
                real Dx, Dy;
                real vel = sqrt(u*u + v*v);
                if (vel < 1e-6) {
                    Dx = DMOLECULAR;
                    Dy = DMOLECULAR;
                } 
                else {
                    Dx = ALPHAT * vel + (ALPHAL - ALPHAT) * u*u / vel + DMOLECULAR;
                    Dy = ALPHAT * vel + (ALPHAL - ALPHAT) * v*v / vel + DMOLECULAR;
                }
                // 防止NaN
                if (isnan(Dx)) Dx = 0.0;
                if (isnan(Dy)) Dy = 0.0;
            


                // int east = -1, west = -1, north = -1, south = -1;
                // 获取相邻单元索引
                int east = (i < nxhc-1) ? iGlob + 1 : -1; // 东边单元
                int west = (i > 0) ? iGlob - 1 : -1; // 西边单元
                int north = (j > 0) ? iGlob - nxhc : -1; // 北边单元
                int south = (j < nyhc-1) ? iGlob + nxhc : -1; // 南边单元
                
               
            // c3: c_i-1, c_i, c_i+1
            SArray<real,3> c3;
            // c4: cxp, cxm, cyp, cym
            SArray<real,4> c4;
            
                
                // X方向通量
                if (east != -1 && state.h(east) > state.hmin) {
                    real c_east = rtsw.c(east, 1);
                    // real c_edge = (u > 0) ? rtsw.c(iGlob, 1) : c_east;//上游迎风格式
                    // real c_edge = (rtsw.c(iGlob, 1) + c_east) / 2.0;//中心差分格式
                    
                    // TVD scheme
                    if (state.h(iGlob)>0){
                        if (u > 0){
                            id0 = iGlob - 1;
                            id1 = iGlob;
                            id2 = iGlob + 1;
                        }
                        else {
                            id0 = iGlob + 2;
                            id1 = iGlob + 1;
                            id2 = iGlob;
                        }
                        c3(0) = rtsw.c(id0,1);   c3(1) = rtsw.c(id1,1);   c3(2) = rtsw.c(id2,1);
                                   if (state.h(id0) <= 0.0 || state.isnodata(id0) == 1)    {c3(0) = rtsw.c(id1,1);}
                        // get edge concentration with tvd scheme
                        c4(0) = tvd(c3, u, rtsw.dt, dom.dx());
                    }
                    else {
                        u = 0.0; c4(0) = 0.0;
                    }
                    // adv_x -= u * dy * state.h(east) * c_edge; // 考虑水深
                    adv_x -= u * dy * state.h(east) * c4(0); // TVD
                    disp_x += Dx * dy * (c_east - rtsw.c(iGlob, 1)) / dx;
                }
     
                
                if (west != -1 && state.h(west) > state.hmin) {
                    real c_west = rtsw.c(west, 1);
                    
                    // real c_edge = (u < 0) ? rtsw.c(iGlob, 1) : c_west;//上游迎风格式
                    // real c_edge = (rtsw.c(iGlob, 1) + c_west) / 2.0;//中心差分格式
                    // TVD scheme
                    if (state.h(iGlob)>0){
                        if (u > 0){
                            id0 = iGlob - 2;
                            id1 = iGlob - 1;
                            id2 = iGlob;
                        }
                        else {
                            id0 = iGlob + 1;
                            id1 = iGlob;
                            id2 = iGlob - 1;
                        }
                        c3(0) = rtsw.c(id0,1);   c3(1) = rtsw.c(id1,1);   c3(2) = rtsw.c(id2,1);
                                   if (state.h(id0) <= 0.0 || state.isnodata(id0) == 1)    {c3(0) = rtsw.c(id1,1);}
                        // get edge concentration with tvd scheme
                        c4(1) = tvd(c3, u, rtsw.dt, dom.dx());
                    }
                    else {
                        u = 0.0; c4(1) = 0.0;
                    }
                    // adv_x += u * dy * state.h(west) * c_edge; // 考虑水深
                    adv_x += u * dy * state.h(west) * c4(1); // TVD
                    disp_x += Dx * dy * (c_west - rtsw.c(iGlob, 1)) / dx;
                }
 
                
                // Y方向通量
                if (north != -1 && state.h(north) > state.hmin) {
                    real c_north = rtsw.c(north, 1);
                    
                    // real c_edge = (v < 0) ?  rtsw.c(iGlob, 1) : c_north;
                    // real c_edge = (rtsw.c(iGlob, 1) + c_north) / 2.0;
                    
                    
                    
                    // TVD scheme
                    if (state.h(iGlob)>0){
                        // North界面 (当前单元与上方单元之间)
                        if (v > 0) { // 向下流动
                            id0 = iGlob - 2*nxhc; // 更上方单元
                            id1 = iGlob - nxhc;   // 上方单元
                            id2 = iGlob;          // 当前单元
                        } else { // 向上流动
                            id0 = iGlob + nxhc;   // 下方单元
                            id1 = iGlob;          // 当前单元
                            id2 = iGlob - nxhc;   // 上方单元
                        }
                        c3(0) = rtsw.c(id0,1);   c3(1) = rtsw.c(id1,1);   c3(2) = rtsw.c(id2,1);
                        if (state.h(id0) <= 0.0 || state.isnodata(id0) == 1)    {c3(0) = rtsw.c(id1,1);}
                        // get edge concentration with tvd scheme
                        c4(2) = tvd(c3, v, rtsw.dt, dom.dx());     
                    }
                    else {
                        v = 0.0; c4(2) = 0.0;
                    }
                    // adv_y += v * dx * state.h(north) * c_edge; // 考虑水深
                    adv_y += v * dx * state.h(north) * c4(2); //with TVD
                    disp_y += Dy * dx * state.h(north)  * (c_north - rtsw.c(iGlob, 1)) / dy;
                }
          
                
                if (south != -1 && state.h(south) > state.hmin) {
                    real c_south = rtsw.c(south, 1);
                    real c_edge = (v > 0) ? rtsw.c(iGlob, 1) : c_south;
                    // real c_edge = (rtsw.c(iGlob, 1) + c_south) / 2.0;
                    
                    // TVD scheme
                if (state.h(iGlob)>0){
                    // South界面 (当前单元与下方单元之间)
                    if (v > 0) { // 向下流动
                        id0 = iGlob - nxhc;   // 上方单元
                        id1 = iGlob;          // 当前单元
                        id2 = iGlob + nxhc;   // 下方单元
                    } else { // 向上流动
                        id0 = iGlob + 2*nxhc; // 更下方单元
                        id1 = iGlob + nxhc;   // 下方单元
                        id2 = iGlob;          // 当前单元
                    }
                    c3(0) = rtsw.c(id0,1);   c3(1) = rtsw.c(id1,1);   c3(2) = rtsw.c(id2,1);
                    if (state.h(id0) <= 0.0 || state.isnodata(id0) == 1)    {c3(0) = rtsw.c(id1,1);}
                    // get edge concentration with tvd scheme
                    c4(3) = tvd(c3, v, rtsw.dt, dom.dx());
                }
                else {
                    v = 0.0; c4(3) = 0.0;
                }
                    
                    // adv_y -= v * dx * state.h(south) * c_edge; // 考虑水深
                    adv_y -= v * dx * state.h(south) * c4(3);  // tvd

                    disp_y += Dy * dx * state.h(south) * (c_south - rtsw.c(iGlob, 1)) / dy;
                }
            
                
                if (state.h(iGlob) > state.hmin) {

                
                // 更新溶质质量 (m = c * h * dx * dy)
                rtsw.solute_mass(iGlob) += rtsw.dt * (adv_x + disp_x + adv_y + disp_y);
                

                // 从质量转换为浓度: c = m / (h * area)
                real area = dx * dy; // 单元面积
                rtsw.c(iGlob, 1) = rtsw.solute_mass(iGlob) / (h * area);
                }
                // printf("iGlob=%d, c(iGlob,0)=%f, c(iGlob,1)=%f,  adv_y=%f, disp_y=%f\n", 
                //     iGlob, rtsw.c(iGlob, 0), rtsw.c(iGlob, 1),  adv_y, disp_y);
               
                }
        
                });       

                // apply boundary conditions
                for (int k = 0; k < ebc.rtbc.size(); k ++) {
                    ebc.rtbc[k].applyrtbc(rtsw, dom);
                }
               

                // 更新浓度
                Kokkos::parallel_for("update_concentration", dom.nCellMem,
                            KOKKOS_LAMBDA (const int iGlob) {
                                rtsw.c(iGlob, 0) = rtsw.c(iGlob,1);
                            }); 
        
    }

        // 2nd-order TVD scheme with superbee limiter
    KOKKOS_INLINE_FUNCTION real tvd(const SArray<real,3> &c3, real u, real dt, real dx) {
        real r, phi = 0.0, r1 = 1.0, r2 = 2.0;
        real co = myfabs(u) * dt / dx;//Courant number
        real c0 = c3(0), c1 = c3(1), c2 = c3(2);
        if (c2 != c1)   {
            r = (c1 - c0) / (c2 - c1);
            if (2.0*r < 1.0)    {r1 = 2.0*r;}
            if (r < 2.0)    {r2 = r;}
            if (r1 > 0.0 || r2 > 0.0)   {phi = max(r1, r2);}
        }
        return c1 + 0.5*phi*(1.0-co)*(c2-c1);
    }





 
};

#endif
#endif