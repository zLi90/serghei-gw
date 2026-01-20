/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _METEO_STATE_H_
#define _METEO_STATE_H_

#include "../define.h"
// #include "SArray.h" // 如果 realArr 定义在 define.h 中，通常不需要显式包含 SArray.h，视具体 define.h 内容而定

class MeteoState {

public:
    // Meta data
    real lat, lon, elev;
    real angstA, angstB; // Angstrom parameters

    int start_year, start_month, start_day, start_doy; // Simulation start date

    // Data Size
    size_t num_records; // 修改为 size_t

    // Time series data (Kokkos Views)
    // 1D Arrays for each variable
    realArr time;   // seconds
    realArr irrad;  // J/m2/d
    realArr tmin;   // C
    realArr tmax;   // C
    realArr vap;    // hPa
    realArr wind;   // m/s
    realArr rain;   // mm/d

    // ADDED for Evapotranspiration module
    realArr e0;     // Potential evaporation from open water [cm/day]
    realArr es0;    // Potential evaporation from bare soil [cm/day]
    realArr et0;    // Potential evapotranspiration from reference crop [cm/day]

    // Helper to get data at a specific simulation time (Step Function)
    // 根据模拟时间(秒)获取气象数据索引
    KOKKOS_INLINE_FUNCTION
    int get_index_at_time(real current_time) const {
        // Assuming time is sorted and uniform (daily, 86400s)
        // index = floor(current_time / 86400.0)
        int idx = static_cast<int>(current_time / 86400.0);
        
        // Handle boundaries
        if (idx < 0) idx = 0;
        // 注意：num_records 是 size_t，需要转换比较
        if (idx >= static_cast<int>(num_records)) idx = static_cast<int>(num_records) - 1;
        
        return idx;
    }

    // Allocate memory
    // 参数 n 改为 size_t
    void allocate(size_t n) {
        num_records = n;
        
        // 直接使用 std::string label 构造，第二个参数必须是 size_t
        time  = realArr("meteo_time", n);
        irrad = realArr("meteo_irrad", n);
        tmin  = realArr("meteo_tmin", n);
        tmax  = realArr("meteo_tmax", n);
        vap   = realArr("meteo_vap", n);
        wind  = realArr("meteo_wind", n);
        rain  = realArr("meteo_rain", n);
        
        // Allocate ADDED arrays
        e0    = realArr("meteo_e0", n);
        es0   = realArr("meteo_es0", n);
        et0   = realArr("meteo_et0", n);
    }
};

#endif