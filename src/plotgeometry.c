#include "header.h"

long plotGeometry() {
    fprintf(stdout, "\nPlotting geometry...\n");

    // Test surface test function
    double x = 0.0, y = 0.0, z = 0.0;
    for (size_t i = 0; i < DATA.n_surf; ++i) {
        Surface *s = &DATA.surfs[i];
        double res = surfaceTest(s->type, s->params, s->n_params, x, y, z);
        fprintf(stdout, "Surface %s test at (%.2f, %.2f, %.2f): %.2f\n", s->name, x, y, z, res);
    }

    fprintf(stdout, "DONE.\n");
    return 0;
}