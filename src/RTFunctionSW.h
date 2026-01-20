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




real calculate_dt(State &state, Domain &dom) {
    real max_vel = 0;
    Kokkos::parallel_reduce("max_vel", dom.nCellMem,
        KOKKOS_LAMBDA(const int iGlob, real& maxv){
            real u = state.hu4rtsw(iGlob,0)/state.h4rtsw(iGlob,0);
            
            real v = state.hv4rtsw(iGlob,0)/state.h4rtsw(iGlob,0);
            real vel = sqrt(u*u + v*v);
            if (vel > maxv) maxv = vel;
        }, Kokkos::Max<real>(max_vel));
    
    const real CFL = 0.5;
    real dt = CFL * dom.dx() / (max_vel + 1e-6);
    return dt;
}






    // 显式对流-弥散求解器（二维地表水版）
    KOKKOS_INLINE_FUNCTION 
    void transport(RTStateSW &rtsw, State &state, Domain &dom, ExternalBoundaries &ebc,  const SourceSinkData &ss) {
        const real dx = dom.dx();
        const real dy = dom.dx();


        // 时间步长控制
        // 在transport函数开头
// rtsw.dt = calculate_dt(state, dom);
// dom.dt = rtsw.dt;
// printf("transport: dt = %f\n", rtsw.dt);


        rtsw.dt = dom.dt;
        // rtsw.dt = 0.1;

        

        Kokkos::parallel_for("const_h", dom.nCellMem, 
        KOKKOS_LAMBDA(const int iGlob){
            int i, j, id0, id1, id2, nxhc, nyhc;
            nxhc = dom.nx+2*hc;
            nyhc = dom.ny+2*hc;   
            unpackIndicesUniformGrid(iGlob, nyhc, nxhc, j, i); 
                   
            // state.h4rtsw(iGlob, 1) = 0.001;
            // state.h4rtsw(iGlob, 0) = 0.001; 

        }); 
        
        Kokkos::parallel_for("sw_transport", dom.nCellMem, 
            KOKKOS_LAMBDA(const int iGlob) {



            int i, j, id0, id1, id2, nxhc, nyhc;
            nxhc = dom.nx+2*hc;
            nyhc = dom.ny+2*hc;   
            unpackIndicesUniformGrid(iGlob, nyhc, nxhc, j, i);  
            
            
real hmin = 0; // 水深阈值
          

            //复制上一时间步浓度
            // rtsw.c(iGlob, 0) = rtsw.c(iGlob, 1);

            
            //!初始化：初始质量 m^n = C^n·h^n·dx·dy
            
            // rtsw.solute_mass(iGlob) = rtsw.c(iGlob, 0) * rtsw.h(iGlob) * dx * dy; // 质量守恒
            // rtsw.solute_mass(iGlob, 0) = rtsw.c(iGlob, 0) * state.h4rtsw(iGlob, 0) * dx * dy; // 质量守恒
rtsw.solute_mass(iGlob, 0) = rtsw.c(iGlob, 0) ;

       
                // 获取相邻单元索引
                int east = (i < nxhc-1) ? iGlob + 1 : -1; // 东边单元
                int west = (i > 0) ? iGlob - 1 : -1; // 西边单元
                int north = (j > 0) ? iGlob - nxhc : -1; // 北边单元
                int south = (j < nyhc-1) ? iGlob + nxhc : -1; // 南边单元            
            
            
                // 跳过无效单元
            // if (state.isnodata(iGlob)) return;
            
            real adv_x = 0.0, adv_y = 0.0;
            real disp_x = 0.0, disp_y = 0.0;
            
            
            
            // 仅在水深大于阈值时计算
            // if (state.h4rtsw(iGlob, 0) > hmin) {

                // printf("hmin = %f\n", hmin);
                // 计算流速分量
                // 计算流速分量 - 关键修改：右边界使用左侧单元流速
                real u, v;
                u = state.hu4rtsw(iGlob, 0) / max(state.h4rtsw(iGlob, 0), 1e-6);
             
                v = -state.hv4rtsw(iGlob, 0) / max(state.h4rtsw(iGlob, 0), 1e-6);

                //!debug
                // real u = 0.01; // 固定流速
                // real v = 0; // 固定流速
          
                // 计算弥散系数（简化版）
                real ALPHAL = 0;
                real ALPHAT = 0;
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
            



                
               
            // c3: c_i-1, c_i, c_i+1
            SArray<real,3> c3;
            // c4: cxp, cxm, cyp, cym
            SArray<real,4> c4;
            
                
                // X方向通量
                if (east != -1 && state.h4rtsw(east, 0) > hmin) {//非esat边界，即可能是west边界
                    real c_east = rtsw.c(east, 0);
                    // real c_edge = (u > 0) ? rtsw.c(iGlob, 1) : c_east;//上游迎风格式
                    // real c_edge = (rtsw.c(iGlob, 1) + c_east) / 2.0;//中心差分格式
                    
                    // TVD scheme
                    if (state.h4rtsw(iGlob, 0) > hmin){
                        if (u >= 0){
                            id0 = iGlob - 1;
                            id1 = iGlob;
                            id2 = iGlob + 1;
                        }
                        else {
                            id0 = iGlob + 2;
                            id1 = iGlob + 1;
                            id2 = iGlob;
                        }
                        c3(0) = rtsw.c(id0,0);   c3(1) = rtsw.c(id1,0);   c3(2) = rtsw.c(id2,0);
                        if (state.h4rtsw(id0, 0) <= 0.0 || state.isnodata(id0) == 1)    {c3(0) = rtsw.c(id1,0);}
                        // get edge concentration with tvd scheme
                        c4(0) = tvd(c3, u, rtsw.dt, dom.dx());
                    }
                    else {
                        u = 0.0; c4(0) = 0.0;
                    }
                    // adv_x -= u * dy * 0.5 * (state.h4rtsw(iGlob, 0) + state.h4rtsw(east, 0)) * c4(0); // TVD
                    adv_x -= u * c4(0) / dx; // TVD
                                      
                    disp_x += Dx * (rtsw.c(east, 0) - rtsw.c(iGlob, 0)) / dx/dx;
                   
                    
                }
     
                
                if (west != -1 && state.h4rtsw(west, 0) > hmin) {//非west边界，即可能是east边界
                    real c_west = rtsw.c(west, 0);
                    
                    // real c_edge = (u < 0) ? rtsw.c(iGlob, 1) : c_west;//上游迎风格式
                    // real c_edge = (rtsw.c(iGlob, 1) + c_west) / 2.0;//中心差分格式
                    // TVD scheme
                    if (state.h4rtsw(iGlob, 0)>hmin){
                        if (u >= 0){
                            id0 = iGlob - 2;
                            id1 = iGlob - 1;
                            id2 = iGlob;
                        }
                        else {
                            id0 = iGlob + 1;//!
                            id1 = iGlob;
                            id2 = iGlob - 1;
                        }
                        c3(0) = rtsw.c(id0,0);   c3(1) = rtsw.c(id1,0);   c3(2) = rtsw.c(id2,0);
                                   if (state.h4rtsw(id0, 0) <= 0.0 || state.isnodata(id0) == 1)    {
                                    // printf("state.isnodata(id0)=%d,state.h4rtsw(id0, 0) =%e, u=%e\n",state.isnodata(id0),state.h4rtsw(id0, 0),u);
                                   
                                    // printf("111111\n");
                                    c3(0) = rtsw.c(id1,0);
                                }
                        // get edge concentration with tvd scheme
                        c4(1) = tvd(c3, u, rtsw.dt, dom.dx());
                    }
                    else {
                        u = 0.0; 
                        c4(1) = 0.0;
                        //todo 添加当水深小于最小水深时，disp_x=0
                    }
                    // adv_x += u * dy * rtsw.h(west) * c_edge; // 考虑水深
                    // adv_x += u * dy * rtsw.h(west) * c4(1); // TVD
                    // adv_x += u * dy * 0.5 * (state.h4rtsw(iGlob, 0) + state.h4rtsw(west, 0)) * c4(1); // TVD
                    
                    
                    adv_x += u  * c4(1) /dx; // TVD
                
                    disp_x += Dx * (rtsw.c(west, 0) - rtsw.c(iGlob, 0)) / dx /dx;
                }
 
                
// Y方向通量
// 非上边界处理
if (north != -1 && state.h4rtsw(north, 0) > hmin) {
                        real c_north = rtsw.c(north, 0);
                    
                    real c_edge = (v < 0) ?  rtsw.c(iGlob, 0) : c_north;
    if (state.h4rtsw(iGlob, 0) > hmin) {
        // 修正索引顺序
        if (v > 0) { // 向下流动
                            id0 = iGlob - 2*nxhc; // 更上方单元
                            id1 = iGlob - nxhc;   // 上方单元
                            id2 = iGlob;          // 当前单元
        } else { // 向上流动
            id0 = iGlob + nxhc;   // !南边单元
            id1 = iGlob;          // 当前单元
            id2 = iGlob - nxhc;   // 北边单元
        }
        c3(0) = rtsw.c(id0,0); c3(1) = rtsw.c(id1,0); c3(2) = rtsw.c(id2,0);
        if (state.h4rtsw(id0, 0) <= 0.0 || state.isnodata(id0) == 1) {c3(0) = c3(1); }
        c4(2) = tvd(c3, v, rtsw.dt, dom.dx()); // 使用 dy
    }
    else {
        v = 0.0; c4(2) = 0.0;
    }
    // 使用界面平均水深
    // adv_y += v * dx* 0.5 * (state.h4rtsw(iGlob, 0) + state.h4rtsw(north, 0)) * c4(2);
    adv_y += v  * c4(2) /dx;
   
    // adv_y += v * dx * rtsw.h(north) * c_edge; // 考虑水深
    // disp_y += Dy * dx * state.h(north)  * (c_north - rtsw.c(iGlob, 0)) / dy;
    disp_y += Dy   * (c_north - rtsw.c(iGlob, 0)) / dy /dy;

}
          
//非下边界处理                
    if (south != -1 && state.h4rtsw(south, 0) > hmin) { 
    real c_south = rtsw.c(south, 0); // 南边单元的浓度
    real c_edge = (v > 0) ? rtsw.c(iGlob, 0) : c_south;
    

    if (state.h4rtsw(iGlob, 0) > hmin) { // 当前单元有有效水深
        // 南边界面处理
        if (v > 0) { // 向下流动（当前单元 → 南边单元）
                        id0 = iGlob - nxhc;   // 上方单元
                        id1 = iGlob;          // 当前单元
                        id2 = iGlob + nxhc;   // 下方单元
        } else { // 向上流动（南边单元 → 当前单元）
            id0 = iGlob + 2 * nxhc;   // 更南边的单元
            id1 = iGlob + nxhc;       // 南边单元
            id2 = iGlob;              // 当前单元
        }

        // 提取相邻单元的浓度值
        c3(0) = rtsw.c(id0, 0);
        c3(1) = rtsw.c(id1, 0);
        c3(2) = rtsw.c(id2, 0);

        // 处理无效单元（如边界外或干单元）
        if (state.h4rtsw(id0, 0) <= 0.0 || state.isnodata(id0) == 1) { c3(0) = c3(1); }

        // 使用TVD格式计算南边界面浓度
        c4(3) = tvd(c3, v, rtsw.dt, dy); // 使用 dy 替代 dx
    } else {
        v = 0.0; // 当前单元无有效水深时，速度置零
        c4(3) = 0.0;
    }


    // adv_y -= v * dx * rtsw.h(south) * c_edge; // 考虑水深
    // adv_y -= v * dom.dx() * 0.5 * (state.h4rtsw(iGlob, 0) + state.h4rtsw(south, 0)) * c4(3);
    adv_y -= v * c4(3) / dx;
    // disp_y += Dy * dx * state.h(south) * (c_south - rtsw.c(iGlob, 0)) / dy;
    disp_y += Dy *  (c_south - rtsw.c(iGlob, 0)) / dy /dy;

}
            
                
                if (state.h4rtsw(iGlob, 0) > state.hmin) {
                   

                
                // 更新溶质质量 (m = c * h * dx * dy)
                rtsw.solute_mass(iGlob, 1) = rtsw.solute_mass(iGlob, 0) + rtsw.dt * (adv_x + disp_x + adv_y + disp_y);
                // rtsw.solute_mass(iGlob) = rtsw.c(iGlob, 0) * state.h4rtsw(iGlob, 0) * dx * dy + rtsw.dt * (adv_x + disp_x + adv_y + disp_y);
      
                // 从质量转换为浓度: c = m / (h * area)
                real area = dx * dy; // 单元面积
                // rtsw.c(iGlob, 1) = rtsw.solute_mass(iGlob, 1) / (state.h4rtsw(iGlob, 1) * area);
                rtsw.c(iGlob, 1) = rtsw.solute_mass(iGlob, 1);

//! 添加地表降雨中溶质浓度含量
if (dom.isRain) {
    // int ii = dom.getIndex(iGlob);
    real RainCon = 1e-15;
    // if (i< nxhc-1){
    rtsw.c(iGlob, 1) += ss.rainRate(iGlob) * RainCon * rtsw.dt / state.h4rtsw(iGlob, 0);
    // }
    
}


                } 
// else{
//     printf("lower than state.hmin=%e\n",state.hmin);
// }              
                // }        
                });    



                

                // apply boundary conditions
                for (int k = 0; k < ebc.rtbc.size(); k ++) {
                    ebc.rtbc[k].applyrtbc(rtsw, dom);
                }
               

                // 更新浓度
                Kokkos::parallel_for("update_concentration", dom.nCellMem,
                            KOKKOS_LAMBDA (const int iGlob) {
                                rtsw.c(iGlob, 0) = rtsw.c(iGlob,1);
                                rtsw.c(iGlob, 2) = rtsw.c(iGlob, 1);
                                rtsw.csw4gw(iGlob) = rtsw.c(iGlob, 1);
                                // if (iGlob == 74){
                                // printf("iGlob=%d , 222rtsw.c(iGlob, 1)=%e, rtsw.c(iGlob, 0)=%e\n",iGlob,rtsw.c(iGlob, 1),rtsw.c(iGlob, 0));
                                // printf("222rtsw.csw4gw(iGlob)=%f, rtsw.c(iGlob, 0)=%f, state.hmin=%e\n",rtsw.csw4gw(iGlob), rtsw.c(iGlob, 0), state.hmin);
                            // }                              
                            }); 

        
    }
    
    //!计算地表地下溶质交换量
    inline void computeGwExchange_For_RT(State &state, RTStateSW &rtsw , const Domain &dom) {
        Kokkos::parallel_for( dom.nCell , KOKKOS_LAMBDA (int idom) {
            int ii = dom.getIndex(idom);//idom without halo, ii with halo

// printf("idom=%d, ii=%d,rtsw.c(ii, 1)=%e \n",idom,ii, rtsw.c(ii, 1));
            
            rtsw.c(ii, 1) += rtsw.dt * rtsw.ConQss(idom)  / state.h4rtsw(ii, 0);
            if (rtsw.c(ii, 1)<0) {rtsw.c(ii, 1)=0.0;}
            if (isnan(rtsw.c(ii, 1))) {rtsw.c(ii, 1)=0.0;}
            rtsw.c(ii, 0) = rtsw.c(ii,1);
            // if (idom = 11){
// printf("idom=%d, ii=%d, 333rtsw.c(ii, 1)=%e,rtsw.ConQss(idom)=%e, rtsw.dt=%e, state.h4rtsw(ii, 0)=%e\n",idom,ii, rtsw.c(ii, 1),rtsw.ConQss(idom), rtsw.dt, state.h4rtsw(ii, 0));
            //     }
            
            
            // if(rtsw.c(ii, 1)<TOL12) {rtsw.c(ii, 1)=0.0;}
        });
    }



    // 修改后的TVD函数
//     KOKKOS_INLINE_FUNCTION real tvd(const SArray<real,3> &c3, real u, real dt, real dy) {
//         real r, phi = 0.0;
//         real co = myfabs(u) * dt / dy; // 使用 dy 替代 dx
//         real c0 = c3(0), c1 = c3(1), c2 = c3(2);
//         if (c2 != c1) {
//             r = (c1 - c0) / (c2 - c1);
//             real r1 = (2.0 * r < 1.0) ? 2.0 * r : 1.0; // min(2r, 1)
//             real r2 = (r < 2.0) ? r : 2.0;            // min(r, 2)
//             phi = max(r1, r2);
//             phi = min(phi, 1.0); // 确保 phi 不超过 1.0
//             if (r < 0.0) phi = 0.0; // 确保非负
//         }
//         if (u >= 0) {
//     return c1 + 0.5 * phi * (c2 - c1);
// } else {
//     return c1 - 0.5 * phi * (c2 - c1);
// }
//     }


       // 修改后的TVD函数
    KOKKOS_INLINE_FUNCTION real tvd(const SArray<real,3> &c3, real u, real dt, real dy) {
        real r, phi = 0.0;
        real co = myfabs(u) * dt / dy; // 使用 dy 替代 dx
        real c0 = c3(0), c1 = c3(1), c2 = c3(2);
        if (c2 != c1) {
            r = (c1 - c0) / (c2 - c1);
            real r1 = (2.0 * r < 1.0) ? 2.0 * r : 1.0; // min(2r, 1)
            real r2 = (r < 2.0) ? r : 2.0;            // min(r, 2)
            phi = max(r1, r2);
            phi = min(phi, 1.0); // 确保 phi 不超过 1.0
            if (r < 0.0) phi = 0.0; // 确保非负
        }
        return c1 + 0.5 * phi * (1.0 - co) * (c2 - c1);
    }






 
};

#endif
#endif