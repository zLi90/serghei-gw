#ifndef _AFGEN_H_
#define _AFGEN_H_

#include "../define.h"

class Afgen {
public:
    // Linear interpolation function
    // table: 1D array [x1, y1, x2, y2, ...]
    KOKKOS_INLINE_FUNCTION
    static real lookup(const realArr& table, real x) {
        int n = table.extent(0);
        if (n < 2) return 0.0;

        // Check boundaries
        if (x <= table(0)) return table(1);
        if (x >= table(n-2)) return table(n-1);

        // Linear search (for small tables this is fine on GPU)
        for (int i = 0; i < n - 2; i += 2) {
            if (x >= table(i) && x < table(i+2)) {
                real x1 = table(i);
                real y1 = table(i+1);
                real x2 = table(i+2);
                real y2 = table(i+3);
                
                real slope = (y2 - y1) / (x2 - x1);
                return y1 + slope * (x - x1);
            }
        }
        return table(n-1); // Should not reach here
    }
};

#endif