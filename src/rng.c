#include "header.h"


void sampleDirection(double *u, double *v, double *w, xoshiro256ss_state *rng) {
    double a, b, s;
    
    do
    {
        /* Two random numbers between [-1, 1] */
        
        a = 2.0*randd(rng) - 1.0;
        b = 2.0*randd(rng) - 1.0;

        /* Determine if the point falls inside the unit circle */

        s = a*a + b*b;
    } 
    while (s >= 1.0);

    /* Compute factor to normalize the vector */
    
    double factor = 2.0 * sqrt(1.0 - s);

    /* Put direction vector */

    *u = factor * a;
    *v = factor * b;
    *w = 1.0 - 2.0 * s;
}