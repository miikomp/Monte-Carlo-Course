#include "header.h"

long getMaterialAtPosition(double x, double y, double z, int* err) {

    /* Find current cell */

    long cell_idx = cellSearch(x, y, z, err, NULL, NULL, NULL, NULL);
    Cell *cell;

    if (cell_idx >= 0 && err && *err == CELL_ERR_OK)
        cell = &DATA.cells[cell_idx];
    else
        return -1;

    /* Get material idx */
    
    return cell->mat_idx;
}