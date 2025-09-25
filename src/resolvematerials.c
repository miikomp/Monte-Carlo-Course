

#include "header.h"

int resolveMaterials(TempNucDataLib *lib, size_t nlib) {

    fprintf(stdout, "\nProcessing materials...\n");

    for (size_t im = 0; im < DATA.n_mats; ++im) 
    {
        Material *M = &DATA.mats[im];
        fprintf(stdout, "\nAdding %zu nuclide(s) to material \"%s\":\n",
                (size_t)M->n_nucs, M->name);

        /* First pass: map ZA to lib entry, copy identity (name, AW), pick temperature */

        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            MaterialNuclide *mc = &M->nucs[c];
            const int ZA = mc->nuc_data.ZA;
            if (ZA <= 0) 
            {
                fprintf(stderr, "[ERROR] Material \"%s\" contains invalid ZA=%d.\n", M->name, ZA);
                exit(EXIT_FAILURE);
            }
            const int Z  = ZA / 1000;
            const int A  = ZA % 1000;

            /* find matching nuclide in library */

            const TempNucDataLib *N = NULL;
            for (size_t i = 0; i < nlib; ++i) 
            {
                if (lib[i].Z == Z && lib[i].A == A) 
                { 
                    N = &lib[i]; 
                    break; 
                }
            }
            if (!N) 
            {
                fprintf(stderr, "[ERROR] Material \"%s\" references ZA=%d which is not in the XS library.\n",
                        M->name, ZA);
                exit(EXIT_FAILURE);
            }

            snprintf(mc->nuc_data.name, sizeof mc->nuc_data.name, "%s", N->name);
            mc->nuc_data.Z  = Z;
            mc->nuc_data.A  = A;
            mc->nuc_data.ZA = ZA;
            mc->nuc_data.AW = N->AW;

            if (N->n_var == 0) 
            {
                fprintf(stderr, "[ERROR] No temperature variants found for ZA=%d.\n", ZA);
                exit(EXIT_FAILURE);
            }

            /* Choose closest temperature variant (tie choose higher T) */
            /* No temperature adjustments are made */

            size_t vbest = 0;
            double dbest = fabs(N->var[0].temp - M->temp);
            for (size_t v = 1; v < N->n_var; ++v) 
            {
                double d = fabs(N->var[v].temp - M->temp);
                if (d < dbest || (d == dbest && N->var[v].temp > N->var[vbest].temp)) {
                    dbest = d; vbest = v;
                }
            }
            const Nuclide *src = &N->var[vbest];

            /* Copy the chosen variant into material nuclide data */

            mc->nuc_data.temp       = src->temp;
            if (src->AW != 0.0)
                mc->nuc_data.AW     = src->AW;
            mc->nuc_data.has_nubar  = src->has_nubar;
            memcpy(mc->nuc_data.mt_idx, src->mt_idx, sizeof mc->nuc_data.mt_idx);

            mc->nuc_data.n_xs = src->n_xs;
            mc->nuc_data.xs   = (XsTable*)malloc(src->n_xs * sizeof *mc->nuc_data.xs);
            mc->nuc_data.nubar.n  = 0;
            mc->nuc_data.nubar.E  = NULL;
            mc->nuc_data.nubar.nu = NULL;
            if (!mc->nuc_data.xs) 
            { 
                fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                exit(EXIT_FAILURE); 
            }

            for (size_t k = 0; k < src->n_xs; ++k) 
            {
                const XsTable *ts = &src->xs[k];
                XsTable *td = &mc->nuc_data.xs[k];
                td->mt = ts->mt;
                td->Q  = ts->Q;
                td->n  = ts->n;
                td->E  = (double*)malloc(ts->n * sizeof(double));
                td->xs = (double*)malloc(ts->n * sizeof(double));
                if (!td->E || !td->xs) 
                { 
                    fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                    exit(EXIT_FAILURE); 
                }
                memcpy(td->E,  ts->E,  ts->n * sizeof(double));
                memcpy(td->xs, ts->xs, ts->n * sizeof(double));
            }
            if (src->has_nubar && src->nubar.n > 0) 
            {
                mc->nuc_data.has_nubar = 1;
                mc->nuc_data.nubar.n   = src->nubar.n;
                mc->nuc_data.nubar.E   = (double*)malloc(src->nubar.n * sizeof(double));
                mc->nuc_data.nubar.nu  = (double*)malloc(src->nubar.n * sizeof(double));
                if (!mc->nuc_data.nubar.E || !mc->nuc_data.nubar.nu) 
                { 
                    fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                    exit(EXIT_FAILURE); 
                }
                memcpy(mc->nuc_data.nubar.E,  src->nubar.E,  src->nubar.n * sizeof(double));
                memcpy(mc->nuc_data.nubar.nu, src->nubar.nu, src->nubar.n * sizeof(double));
            }

            /* Print per-nuclide summary */
            
            size_t n_modes = (mc->nuc_data.n_xs > 0 ? mc->nuc_data.n_xs - 1 : 0);
            fprintf(stdout, "  %5d - %s with %zu reaction modes (using %.0fK data).\n",
                    mc->nuc_data.ZA, mc->nuc_data.name, n_modes, mc->nuc_data.temp);
        }

        /* Second pass: compute fractions and number densities using AW and mdens */
        
        fprintf(stdout, "\n  Calculating densities and fractions:\n");

        double sum_atoms = 0.0, sum_w = 0.0;
        for (size_t c = 0; c < M->n_nucs; ++c) 
        { 
            sum_atoms += M->nucs[c].atom_frac; 
            sum_w += M->nucs[c].mass_frac; 
        }

        const double Ntot_target = (M->adens > 0.0) ? (M->adens * 1.0e24) : 0.0;
        

        if (sum_atoms > 0.0) 
        {
            /* Atomic counts provided */

            for (size_t c = 0; c < M->n_nucs; ++c) 
                M->nucs[c].atom_frac /= sum_atoms;
        

            double Abar = 0.0;
            for (size_t c = 0; c < M->n_nucs; ++c) 
                Abar += M->nucs[c].atom_frac * M->nucs[c].nuc_data.AW;

            if (M->mdens > 0.0 && Abar > 0.0) 
            {
                /* Use mass density to set absolute scale */

                const double Ntot = (M->mdens * NA) / Abar;  /* atoms/cm^3 */

                for (size_t c = 0; c < M->n_nucs; ++c) 
                {
                    M->nucs[c].N_i = M->nucs[c].atom_frac * Ntot;
                    M->nucs[c].mass_frac = (M->nucs[c].atom_frac * M->nucs[c].nuc_data.AW) / Abar;
                }
            } 
            else if (M->adens > 0.0) 
            {
                /* Use atomic density atoms/b·cm */

                const double Ntot = Ntot_target;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    M->nucs[c].N_i = M->nucs[c].atom_frac * Ntot;
        
                /* Derive mass fractions from x_i and Abar */

                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].mass_frac = (M->nucs[c].atom_frac * M->nucs[c].nuc_data.AW) / Abar;
            } 
            else
            {
                /* No absolute density use relative scale */

                for (size_t c = 0; c < M->n_nucs; ++c) M->nucs[c].N_i = M->nucs[c].atom_frac;

                /* mass_frac from relative N_i */

                double Wsum = 0.0;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    Wsum += M->nucs[c].N_i * M->nucs[c].nuc_data.AW;
                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].mass_frac = (Wsum > 0.0) ? (M->nucs[c].N_i * M->nucs[c].nuc_data.AW / Wsum) : 0.0;
            }
        } 
        else 
        {
            /* Mass fractions provided: w_i already normalized */

            if (M->mdens > 0.0) 
            {
                /* Convert w_i + rho to absolute N_i */

                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].N_i = (M->mdens * M->nucs[c].mass_frac / M->nucs[c].nuc_data.AW) * NA;

                /* Derive atomic fractions */

                double Nsum = 0.0;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    Nsum += M->nucs[c].N_i;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    M->nucs[c].atom_frac = (Nsum > 0.0) ? (M->nucs[c].N_i / Nsum) : 0.0;

            } 
            else if (M->adens > 0.0) 
            {
                /* Use mdens to set total atoms */
            
                double denom = 0.0;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    denom += (M->nucs[c].nuc_data.AW > 0.0) ? (M->nucs[c].mass_frac / M->nucs[c].nuc_data.AW) : 0.0;
                if (denom <= 0.0) 
                {
                    for (size_t c = 0; c < M->n_nucs; ++c) 
                        M->nucs[c].N_i = (M->nucs[c].mass_frac / fmax(M->nucs[c].nuc_data.AW, 1e-300));
                } 
                else 
                {
                    for (size_t c = 0; c < M->n_nucs; ++c) 
                    {
                        double xi = (M->nucs[c].mass_frac / M->nucs[c].nuc_data.AW) / denom;
                        M->nucs[c].atom_frac = xi;
                        M->nucs[c].N_i = xi * Ntot_target;
                    }
                }
            } 
            else 
            {
                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].N_i = (M->nucs[c].mass_frac / fmax(M->nucs[c].nuc_data.AW, 1e-300));

                /* Derive atomic fractions */

                double Nsum = 0.0; 
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    Nsum += M->nucs[c].N_i;

                for (size_t c = 0; c < M->n_nucs; ++c) 
                    M->nucs[c].atom_frac = (Nsum > 0.0) ? (M->nucs[c].N_i / Nsum) : 0.0;
            }
        }

        /* Print summary of material */

        double xsum = 0.0, wsum = 0.0, Nsum = 0.0, mass_sum = 0.0, Abar_sum = 0.0;

        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            xsum     += M->nucs[c].atom_frac;
            wsum     += M->nucs[c].mass_frac;
            Nsum     += M->nucs[c].N_i;
            mass_sum += M->nucs[c].N_i * M->nucs[c].nuc_data.AW;
            Abar_sum += M->nucs[c].atom_frac * M->nucs[c].nuc_data.AW;
        }

        /* Compute derived densities for verification */

        const double rho_calc   = mass_sum / NA;        /* g/cm^3 */
        const double ndens_calc = Nsum * 1.0e-24;       /* atoms/b*cm */

        fprintf(stdout,
            "  dens=%.6E %s, T=%.1fK\n"
            "  components=%zu, sum(x)=%.6f, sum(w)=%.6f\n"
            "  Abar=%.6f g/mol/atom, N_tot=%.6e 1/cm^3\n"
            "  mdens_calc=%.6E g/cm^3, adens_calc=%.6E atoms/b*cm\n\n",
            (M->adens > 0.0) ? M->adens : M->mdens,
            (M->adens > 0.0) ? "atoms/b*cm" : "g/cm3",
            M->temp,
            M->n_nucs, 
            xsum, 
            wsum,
            Abar_sum, 
            Nsum, 
            rho_calc, 
            ndens_calc
        );

        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            const MaterialNuclide *mc = &M->nucs[c];
            size_t n_modes = (mc->nuc_data.n_xs > 0 ? mc->nuc_data.n_xs - 1 : 0);
            fprintf(stdout,
                "  %5d - %5s : T=%.0fK, AW=%.6E, afrac=%.6E, mfrac=%.6E, N_i=%.6E atoms/b*cm\n",
                mc->nuc_data.ZA, 
                mc->nuc_data.name,
                mc->nuc_data.temp,
                mc->nuc_data.AW, 
                mc->atom_frac, 
                mc->mass_frac, 
                mc->N_i * 1e-24
                );
            if (VERBOSITY >= 1) 
            {
                fprintf(stdout, "  %5zu Reaction modes:\n", n_modes);
                for (size_t k = 0; k < mc->nuc_data.n_xs; ++k) 
                {
                    XsTable *tab = &mc->nuc_data.xs[k];
                    fprintf(stdout, "     MT = %3d, Q = %11.4E, points = %zu\n", tab->mt, tab->Q, tab->n);
                }
                if (mc->nuc_data.has_nubar)
                    fprintf(stdout, "    Nubar grid, points = %zu\n", mc->nuc_data.nubar.n);
                fprintf(stdout, "\n");
            }
        }
    } 

    fprintf(stdout, "\nDONE.\n");

    /* ########################################################################################## */
    /* Free the parsed temporary library */

    for (size_t i = 0; i < nlib; ++i) {
        for (size_t v = 0; v < lib[i].n_var; ++v) 
        {
            Nuclide *vv = &lib[i].var[v];
            for (size_t k = 0; k < vv->n_xs; ++k) 
            {
                free(vv->xs[k].E);
                free(vv->xs[k].xs);
            }
            free(vv->xs);
            if (vv->has_nubar) 
            { 
                free(vv->nubar.E); 
                free(vv->nubar.nu); 
            }
        }
        free(lib[i].var);
    }
    free(lib);

    return EXIT_SUCCESS;
}