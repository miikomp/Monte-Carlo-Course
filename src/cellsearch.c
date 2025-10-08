#include "header.h"

int cellSearch(double x, double y, double z, int *err) {

    const double EPS = 1e-12;

    if (err)
        *err = CELL_ERR_OK;

    int ret = -1;
    double ret_min_abs = INFINITY;

    /* Loop over all cells */

    for (size_t c = 0; c < DATA.n_cells; c++)
    {
        Cell *cell = &DATA.cells[c];

        /* Loop over all surfaces defining this cell */
        double min_abs = INFINITY;
        size_t s;
        for (s = 0; s < cell->n_surfs; s++)
        {
            int side = cell->sides[s];
            Surface *surf = &DATA.surfs[cell->surf_idxs[s]];

            /* Do surface test */

            double res = surfaceTest(surf->type, surf->params, surf->n_params, x, y, z);
            double abs_res = fabs(res);
            if (abs_res < min_abs)
                min_abs = abs_res;

            /* Check result. Sign convention means that product of surfaceTest and
               side flag is positive when the point fulfills the criteria. */

            if (side * res < -EPS)
                break;
        }

        /* Check if all surfaces passed the test */

        if (s == cell->n_surfs)
        {
            if (ret < 0)
            {
                ret = (int)c;
                ret_min_abs = min_abs;
            }
            else
            {
                /* More than one cell satisfied the criteria. Need to check for overlap */

                if (!(ret_min_abs <= EPS && min_abs <= EPS))
                {
                    /* True overlap */
                    if (err && *err != CELL_ERR_OVERLAP)
                        *err = CELL_ERR_OVERLAP;
                }
                else if (min_abs < ret_min_abs)
                {
                    /* On boundary of multiple surfaces */
                    
                    ret = (int)c;
                    ret_min_abs = min_abs;
                }
            }
        }
    }

    if (ret < 0)
    {
        if (err)
            *err = CELL_ERR_UNDEFINED;
    }

    return ret;
}
