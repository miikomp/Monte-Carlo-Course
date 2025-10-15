#include "header.h"

double sampleMaxwellianEnergy(Neutron *n, double T)
{
    /* Input is the distribution parameter, i.e, nuclear temperature in MeV */
    /* Algorithm used from MCNP4C manual page 2-44 */
    double a, b, c, d, R;
    do
    {
        a = randd(&n->state);
        b = randd(&n->state);

        R = a*a + b*b;

    } while (R <= 1.0);

    c = randd(&n->state);
    d = randd(&n->state);
    
    return -T * (a*a * log(c)/R + log(d));
}

void sampleNeutronDirection(Neutron *n) {
    double a, b, s;

    do
    {
        /* Two random numbers between [-1, 1] */
        
        a = 2.0*randd(&n->state) - 1.0;
        b = 2.0*randd(&n->state) - 1.0;

        /* Determine if the point falls inside the unit circle */

        s = a*a + b*b;
    } 
    while (s >= 1.0);

    /* Compute factor to normalize the vector */
    
    double factor = 2.0 * sqrt(1.0 - s);

    /* Put direction vector */

    n->u = factor * a;
    n->v = factor * b;
    n->w = 1.0 - 2.0 * s;
}

double sampleDistanceToCollision(Neutron *n)
{
    /* First we need to get the material we are in */
    int err;
    int mat_idx = getMaterialAtPosition(n->x, n->y, n->z, &err);

    if (mat_idx < 0 || mat_idx >= (int)DATA.n_mats)
    {
        /* Neutron is outside the geometry */

        return -1.0;
    }

    /* Get material */

    Material *mat = &DATA.mats[mat_idx];

    /* Get total macroscopic cross section */

    double sigma_t = getTotalMacroscopicXS(n->E, mat);
    if (sigma_t <= 0.0)
        return -1.0;

    /* Sample distance to collision */

    double xi = randd(&n->state);
    if (xi <= 0.0)
        xi = 1.0 - 1e-12;
    double d = -log(1.0 - xi) / sigma_t;
    return d;
}

int sampleCollisionNuclide(Neutron *n)
{
    /* Get active material */

    Material *mat = &DATA.mats[n->mat_idx];

    /* Sample nuclide in material based on the macroscopic cross sections */

    double total_macro = getTotalMacroscopicXS(n->E, mat);

    if (total_macro <= 0.0)
        return -1;

    /* Scale random number with total macroscopic cross section */

    double xi = randd(&n->state) * total_macro;
    double sum = 0.0;

    /* Loop over nuclides until cumulated cross section exceeds scaled random number */
    
    for (size_t i = 0; i < mat->n_nucs; ++i) 
    {
        MaterialNuclide *mnuc = &mat->nucs[i];
        double xs = getTotalMicroscopicXS(n->E, &mnuc->nuc_data);
        if (xs < 0.0)
            continue;

        sum += mnuc->N_i * xs * BARN_TO_CM2;

        if (sum >= xi) 
            return (int)i;
    }

    /* If we reach here something went wrong */

    return (mat->n_nucs > 0) ? (int)(mat->n_nucs - 1) : -1;
}

int sampleInteractionType(Neutron *n, Nuclide *nuc)
{
    /* Get total microscopic cross section */

    double sigma_t = getTotalMicroscopicXS(n->E, nuc);
    if (sigma_t <= 0.0)
        return -1;

    /* Sample interaction type based on microscopic cross sections */

    double xi = randd(&n->state) * sigma_t;
    double sum = 0.0;

    for (size_t i = 1; i < nuc->n_xs; ++i) 
    {
        XsTable *xs_table = &nuc->xs[i];
        if (xs_table->mt < 1)
            continue;
        double xs = getMicroscopicXS(n->E, xs_table);
        if (xs < 0.0)
            continue;
        sum += xs;
        if (sum >= xi) 
            return xs_table->mt;
    }

    /* If we reach here something went wrong */

    return -1;
}

void sampleIsotropicDirection(xoshiro256ss_state *state, double *u, double *v, double *w)
{
    double a, b, s;
    do
    {
        a = 2.0 * randd(state) - 1.0;
        b = 2.0 * randd(state) - 1.0;
        s = a * a + b * b;
    }
    while (s >= 1.0 || s == 0.0);

    double factor = 2.0 * sqrt(1.0 - s);
    *u = factor * a;
    *v = factor * b;
    *w = 1.0 - 2.0 * s;
}
