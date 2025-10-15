#include "header.h"

long getMaterialAtPosition(double x, double y, double z, int* err) {

    /* Find current cell */

    cellSearchRes res = cellSearch(x, y, z, 0.0, 0.0, 0.0);
    long cell_idx = res.cell_idx;
    *err = res.err;
    Cell *cell;

    if (cell_idx >= 0 && err && *err == CELL_ERR_OK)
        cell = &DATA.cells[cell_idx];
    else
        return -1;

    /* Get material idx */
    
    return cell->mat_idx;
}