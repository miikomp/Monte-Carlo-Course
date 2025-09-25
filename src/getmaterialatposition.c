#include "header.h"

int getMaterialAtPosition(double x, double y, double z) {

    /* For now there is no geometry so we can just return the first material idx */
    x = x; y = y; z = z; // suppress unused variable warnings

    if (DATA.n_mats > 0)
        return 0;
    else
        return -1;
}