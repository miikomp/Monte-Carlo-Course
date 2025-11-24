#include "header.h"

void handleElasticScatter(Neutron *n, Nuclide *nuc)
{
    double Vt, ut, vt, wt;

    /* Depending on neutron energy the target velocity may be ignored or sampled */
    if (n->E < E_FG_LIMIT)
    {
        /* --- Free gas model --- */

        /* Sample target energy from Maxwellian distribution */

        double Et = sampleMaxwellianEnergy(n, DATA.mats[n->mat_idx].kT);

        /* Calculate speed (cm/s) */

        Vt = sqrt(fmax(0.0, 2.0 * Et / (AMU_TO_MEV_C2 * nuc->AW)));

        /* Sample target direction randomly */

        sampleIsotropicDirection(&n->state, &ut, &vt, &wt);
    }
    else
    {
        /* --- Stationary target ---*/

        Vt = 0.0;
        ut = 0.0;
        vt = 0.0;
        wt = 1.0;
    }

    /* Neutron speed (cm/s) */

    double Vn = sqrt(fmax(0.0, 2.0 * n->E / AMU_TO_MEV_C2));

    /* Mass ratio A = m_target / m_neutron */

    double A = nuc->AW * INV_MASS_NEUTRON;
    double denom = A + 1.0;

    /* Laboratory velocities */

    double vn_x = Vn * n->u;
    double vn_y = Vn * n->v;
    double vn_z = Vn * n->w;

    double vt_x = Vt * ut;
    double vt_y = Vt * vt;
    double vt_z = Vt * wt;

    /* Centre-of-mass velocity */

    double vcm_x = (vn_x + A * vt_x) / denom;
    double vcm_y = (vn_y + A * vt_y) / denom;
    double vcm_z = (vn_z + A * vt_z) / denom;

    /* Neutron velocity in COM frame */

    double vn_cx = vn_x - vcm_x;
    double vn_cy = vn_y - vcm_y;
    double vn_cz = vn_z - vcm_z;
    double Vc = sqrt(vn_cx * vn_cx + vn_cy * vn_cy + vn_cz * vn_cz);

    if (!(Vc > 0.0))
    {
        /* Degenerate case: retain direction but ensure it is normalised */
        double norm = sqrt(n->u * n->u + n->v * n->v + n->w * n->w);
        if (norm > 0.0)
        {
            n->u /= norm;
            n->v /= norm;
            n->w /= norm;
        }
        return;
    }

    /* Sample new COM direction isotropically */
    double uc, vc, wc;
    sampleIsotropicDirection(&n->state, &uc, &vc, &wc);

    /* New neutron velocity in COM frame */

    double vn_cx_p = Vc * uc;
    double vn_cy_p = Vc * vc;
    double vn_cz_p = Vc * wc;

    /* Transform back to laboratory frame */

    double vn_x_p = vn_cx_p + vcm_x;
    double vn_y_p = vn_cy_p + vcm_y;
    double vn_z_p = vn_cz_p + vcm_z;

    double Vn_p = sqrt(vn_x_p * vn_x_p + vn_y_p * vn_y_p + vn_z_p * vn_z_p);

    if (Vn_p <= 0.0)
    {
        sampleIsotropicDirection(&n->state, &n->u, &n->v, &n->w);
        n->E = 0.0;
        return;
    }

    /* Put new neutron direction and energy */

    n->u = vn_x_p / Vn_p;
    n->v = vn_y_p / Vn_p;
    n->w = vn_z_p / Vn_p;
    n->E = 0.5 * AMU_TO_MEV_C2 * Vn_p * Vn_p;
}

void handleInelasticScatter(Neutron *n, Nuclide *nuc, int mt)
{
    /* Get xsdata for inelastic level scattering */

    int xs_idx = nuc->mt_idx[mt];
    XsTable xstable = nuc->xs[xs_idx];

    /* Check for discrete vs. continuum */

    if (mt == MT_INELASTIC_CONTINUUM)
    {
        fprintf(stderr, "[ERROR] Continuum inelastic scattering not implemented.\n");
        n->status = NEUTRON_DEAD_TERMINATED;
        return;
    }
    
    /* --- Discrete level inelastic scattering --- */
    
    /* Mass ratio */
    
    double A = nuc->AW * INV_MASS_NEUTRON;
    double denom = A + 1.0;

    /* Neutron lab speed */

    double Vn = sqrt(2.0 * n->E / AMU_TO_MEV_C2);

    /* Neutron lab direction vector */

    double vn_x = Vn * n->u;
    double vn_y = Vn * n->v;
    double vn_z = Vn * n->w;

    /* CM velocity (assume target stationary) */

    double vcm_x = vn_x / denom;
    double vcm_y = vn_y / denom;
    double vcm_z = vn_z / denom;

    /* Neutron velocity in CM */

    double vn_cx = vn_x - vcm_x;
    double vn_cy = vn_y - vcm_y;
    double vn_cz = vn_z - vcm_z;

    double Vc = sqrt(vn_cx*vn_cx + vn_cy*vn_cy + vn_cz*vn_cz);

    if (!(Vc > 0.0)) 
    {
        /* Invalid, keep direction normalized and zero energy */
        double norm = sqrt(n->u*n->u + n->v*n->v + n->w*n->w);
        if (norm > 0.0) 
        { 
            n->u/=norm; 
            n->v/=norm; 
            n->w/=norm; 
        }
        n->E = 0.0;
        return;
    }

    /* CM relative kinetic energy */

    double Erel = 0.5 * AMU_TO_MEV_C2 * Vc * Vc * (1.0 + 1.0 / A);

    /* Apply MT Q value to relative energy */

    double Q = xstable.Q;
    double Erel_p = Erel + Q;

    if (Erel <= 0.0)
    {
        /* Invalid, terminate history */
        n->E = 0.0;
        n->status = NEUTRON_DEAD_TERMINATED;
        return;
    }

    /* Scale neutron speed */

    double Vc_p  = Vc * sqrt(Erel_p / Erel);


    /* Sample new direction in CM isotropically */
    
    double uc, vc, wc;
    sampleIsotropicDirection(&n->state, &uc, &vc, &wc);

    /* New neutron velocity in CM */
    
    double vn_cx_p = Vc_p * uc;
    double vn_cy_p = Vc_p * vc;
    double vn_cz_p = Vc_p * wc;
    
    /* Back to lab */
    
    double vn_x_p = vn_cx_p + vcm_x;
    double vn_y_p = vn_cy_p + vcm_y;
    double vn_z_p = vn_cz_p + vcm_z;

    const double Vn_p = sqrt(vn_x_p*vn_x_p + vn_y_p*vn_y_p + vn_z_p*vn_z_p);

    if (Vn_p <= 0.0) 
    {
        n->E = 0.0;
        sampleIsotropicDirection(&n->state, &n->u, &n->v, &n->w);
        return;
    }

    /* Put data in neutron */

    n->u = vn_x_p / Vn_p;
    n->v = vn_y_p / Vn_p;
    n->w = vn_z_p / Vn_p;
    n->E = 0.5 * AMU_TO_MEV_C2 * Vn_p * Vn_p;

    return;
}

void handleFission(Neutron *n, Nuclide *nuc)
{
    if (!nuc->has_nubar)
    {
        n->fission_yield = 0;
        return;
    }
    
    /* Get average fission yield for the nuclide at specific neutron energy */

    NubarTable nubar_tab = nuc->nubar;

    double nubar;

    /* If below lower boundary use lower boundary value */

    if (n->E <= nubar_tab.E[0]) 
    {
        nubar = nubar_tab.nu[0];
    } 

    /* If above upper boundary use upper boundary value */

    else if (n->E >= nubar_tab.E[nubar_tab.n - 1]) 
    {
        nubar = nubar_tab.nu[nubar_tab.n - 1];
    }

    /* Find the gap where current neutron energy is in using binary search */

    else 
    {
        size_t lo = 0, hi = nubar_tab.n - 1;
        while (hi - lo > 1) {
            size_t mid = (lo + hi) >> 1;
            if (nubar_tab.E[mid] <= n->E)
                lo = mid;
            else
                hi = mid;
        }
        /* Use linear interpolation */

        nubar = linlin(nubar_tab.E[lo], nubar_tab.nu[lo],
                       nubar_tab.E[lo + 1], nubar_tab.nu[lo + 1], n->E);
    }

    /* Sample number of neutrons emitted */

    int k = (int)floor(nubar);
    double xi = randd(&n->state);

    if (xi < (nubar - (double)k))
        k++;

    n->fission_yield = k;

    /* Kill neutron by fission */

    n->status = NEUTRON_DEAD_FISSION;
}

double getMajorantXS(const double E)
{
    const MajorantXsTable *maj = &DATA.majorant_xs;

    if (maj->n == 0 || !maj->E || !maj->xs)
        return -1.0;

    double sigma_M;

    if (E <= maj->E[0])
    {
        sigma_M = maj->xs[0];
    }
    else if (E >= maj->E[maj->n - 1])
    {
        sigma_M = maj->xs[maj->n - 1];
    }
    else
    {
        size_t lo = 0, hi = maj->n - 1;
        while (hi - lo > 1)
        {
            size_t mid = (lo + hi) >> 1;
            if (maj->E[mid] <= E)
                lo = mid;
            else
                hi = mid;
        }

        sigma_M = linlin(maj->E[lo], maj->xs[lo],
                         maj->E[lo + 1], maj->xs[lo + 1], E);
    }

    return sigma_M;
}

double getTotalMacroscopicXS(const double E, Material* mat)
{
    /* Get material total macroscopic cross section table */

    MacroXsTable total = mat->macro_xs[0];

    if (total.n == 0 || !total.E || !total.xs)
        return -1.0;

    double sigma_t;

    /* If below lower boundary use lower boundary value */

    if (E <= total.E[0]) 
    {
        sigma_t = total.xs[0];
    } 

    /* If above upper boundary use upper boundary value */

    else if (E >= total.E[total.n - 1]) 
    {
        sigma_t = total.xs[total.n - 1];
    } 

    /* Find the gap where current neutron energy is in using binary search */

    else 
    {
        size_t lo = 0, hi = total.n - 1;
        while (hi - lo > 1) {
            size_t mid = (lo + hi) >> 1;
            if (total.E[mid] <= E)
                lo = mid;
            else
                hi = mid;
        }
        /* Use linear interpolation */

        sigma_t = linlin(total.E[lo], total.xs[lo],
                         total.E[lo + 1], total.xs[lo + 1], E);
    }

    return sigma_t;
}

double getTotalMicroscopicXS(const double E, Nuclide* nuc)
{
    /* Get nuclide total microscopic cross section table */

    if (nuc->n_xs == 0 || !nuc->xs || !nuc->xs[0].E || !nuc->xs[0].xs)
        return -1.0;

    XsTable total = nuc->xs[0];

    double sigma_t;

    /* If below lower boundary use lower boundary value */

    if (E <= total.E[0]) 
    {
        sigma_t = total.xs[0];
    } 

    /* If above upper boundary use upper boundary value */

    else if (E >= total.E[total.n - 1]) 
    {
        sigma_t = total.xs[total.n - 1];
    } 

    /* Find the gap where current neutron energy is in using binary search */

    else 
    {
        size_t lo = 0, hi = total.n - 1;
        while (hi - lo > 1) {
            size_t mid = (lo + hi) >> 1;
            if (total.E[mid] <= E)
                lo = mid;
            else
                hi = mid;
        }
        /* Use linear interpolation */

        sigma_t = linlin(total.E[lo], total.xs[lo],
                         total.E[lo + 1], total.xs[lo + 1], E);
    }

    return sigma_t;
}

double getMicroscopicXS(const double E, XsTable* xs_table)
{
    /* Get nuclide microscopic cross section table for specific MT */

    if (xs_table->n == 0 || !xs_table->E || !xs_table->xs)
        return -1.0;

    double sigma;

    /* If below lower boundary use lower boundary value */

    if (E <= xs_table->E[0]) 
    {
        sigma = xs_table->xs[0];
    } 

    /* If above upper boundary use upper boundary value */

    else if (E >= xs_table->E[xs_table->n - 1]) 
    {
        sigma = xs_table->xs[xs_table->n - 1];
    } 

    /* Find the gap where current neutron energy is in using binary search */

    else 
    {
        size_t lo = 0, hi = xs_table->n - 1;
        while (hi - lo > 1) {
            size_t mid = (lo + hi) >> 1;
            if (xs_table->E[mid] <= E)
                lo = mid;
            else
                hi = mid;
        }
        /* Use linear interpolation */

        sigma = linlin(xs_table->E[lo], xs_table->xs[lo],
                       xs_table->E[lo + 1], xs_table->xs[lo + 1], E);
    }

    return sigma;
}

double getVelocityCmPerS(const double E)
{
    if (E < 0.0)
        return -1.0;
    return C_LIGHT * sqrt(2.0 * E / (AMU_TO_MEV_C2 * MASS_NEUTRON));
}

int checkNeutronCutoff(Neutron *n)
{
    if (n->E < GLOB.energy_cutoff)
    {
        n->status = NEUTRON_DEAD_TERMINATED;
        return 1;
    }
    else if (n->time >= GLOB.time_cutoff)
    {
        n->status = NEUTRON_DEAD_TERMINATED;
        return 1;
    }
    else if (n->genc >= GLOB.generation_cutoff)
    {
        n->status = NEUTRON_DEAD_LEAKAGE;
        return 1;
    }
    return 0;
}
