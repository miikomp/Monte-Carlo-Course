#include "header.h"

long getMaterialAtPosition(double x, double y, double z) {
    /* Find current cell */

    int err;
    long cell_idx = cellSearch(x, y, z, &err);
    Cell *cell;

    if (cell_idx >= 0 && err == CELL_ERR_OK)
        cell = &DATA.cells[cell_idx];
    else
        return -1;

    /* Get material idx */
    
    return cell->mat_idx;
}