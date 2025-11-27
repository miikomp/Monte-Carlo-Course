#include "header.h"

static const XsTable *findXsTable(const Nuclide *nu, int mt);

int computeMacroXs(void) 
{
    fprintf(stdout, "\nComputing macroscopic cross sections...\n");

    /* Loop over materials*/

    for (size_t im = 0; im < DATA.n_mats; ++im) 
    {
        Material *M = &DATA.mats[im];

        /* Reset material macro tables */

        if (!M->macro_xs || M->n_macro_xs == 0) 
        {
            M->macro_xs = NULL;
            M->n_macro_xs = 0;
        }

        for (size_t i = 0; i < M->n_macro_xs; ++i) 
        {
            free(M->macro_xs[i].E);
            free(M->macro_xs[i].xs);
        }
        free(M->macro_xs);
        M->macro_xs = NULL;
        M->n_macro_xs = 0;

        if (M->n_nucs == 0) 
        {
            fprintf(stdout, "\nMaterial \"%s\" has no nuclides. Skipping...\n", M->name);
            continue;
        }

        /* Collect all unique MT identifiers present in material nuclides */

        int *mt_list = NULL;
        size_t mt_count = 0, mt_cap = 0;

        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            const Nuclide *nu = &M->nucs[c].nuc_data;

            for (size_t k = 0; k < nu->n_xs; ++k) 
            {
                const int mt = nu->xs[k].mt;
                if (mt < 0)
                    continue;

                bool present = false;
                for (size_t j = 0; j < mt_count; ++j) 
                {
                    if (mt_list[j] == mt) 
                    {
                        present = true;
                        break;
                    }
                }
                if (present)
                    continue;

                if (mt_count == mt_cap) 
                {
                    size_t new_cap = mt_cap ? mt_cap * 2 : 8;
                    int *tmp = (int*)realloc(mt_list, new_cap * sizeof(int));
                    if (!tmp) 
                    {
                        free(mt_list);
                        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                        exit(EXIT_FAILURE);
                    }
                    mt_list = tmp;
                    mt_cap = new_cap;
                }
                mt_list[mt_count++] = mt;
            }
        }

        if (mt_count == 0) 
        {
            fprintf(stdout, "\nMaterial \"%s\" has no microscopic cross sections. Skipping...\n", M->name);
            free(mt_list);
            continue;
        }

        MacroXsTable *macro = (MacroXsTable*)calloc(mt_count, sizeof(MacroXsTable));
        if (!macro) 
        {
            free(mt_list);
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            exit(EXIT_FAILURE);
        }

        size_t produced = 0;

        /* Loop over all unique MTs found */

        for (size_t m_idx = 0; m_idx < mt_count; ++m_idx) 
        {
            const int mt = mt_list[m_idx];

            double *E_union = NULL;
            size_t n_union = 0, cap_union = 0;

            /* Gather union energy grid from nuclides that contain this MT */

            for (size_t c = 0; c < M->n_nucs; ++c) 
            {
                const Nuclide *nu = &M->nucs[c].nuc_data;
                const XsTable *tab = findXsTable(nu, mt);
                if (!tab || tab->n == 0)
                    continue;

                if (cap_union < n_union + tab->n) 
                {
                    size_t new_cap = cap_union ? cap_union : tab->n;
                    while (new_cap < n_union + tab->n)
                        new_cap *= 2;
                    double *tmp = (double*)realloc(E_union, new_cap * sizeof(double));
                    if (!tmp) {
                        free(E_union);
                        free(mt_list);
                        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                        exit(EXIT_FAILURE);
                    }
                    E_union = tmp;
                    cap_union = new_cap;
                }

                memcpy(&E_union[n_union], tab->E, tab->n * sizeof(double));
                n_union += tab->n;
            }

            if (n_union == 0) {
                free(E_union);
                continue;
            }

            /* Sort energies and make union grid */

            qsort(E_union, n_union, sizeof(double), cmpDouble);

            size_t unique_count = 0;
            for (size_t i = 0; i < n_union; ++i) 
            {
                if (i == 0) 
                {
                    E_union[unique_count++] = E_union[i];
                    continue;
                }
                double a = E_union[i];
                double b = E_union[unique_count - 1];
                double tol = 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));

                /* Check for duplicate grid points */

                if (fabs(a - b) > tol)
                    E_union[unique_count++] = a;
            }

            double *E_final = (double*)realloc(E_union, unique_count * sizeof(double));
            if (!E_final) 
            {
                free(E_union);
                free(mt_list);
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
            E_union = E_final;

            double *macro_values = (double*)malloc(unique_count * sizeof(double));
            if (!macro_values) 
            {
                free(E_union);
                free(mt_list);
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }

            for (size_t j = 0; j < unique_count; ++j) 
            {
                double E = E_union[j];
                double sum = 0.0;

                for (size_t c = 0; c < M->n_nucs; ++c) 
                {
                    const MaterialNuclide *mn = &M->nucs[c];
                    const Nuclide *nu = &mn->nuc_data;
                    const XsTable *tab = findXsTable(nu, mt);
                    if (!tab || tab->n == 0)
                        continue;

                    double micro;
                    if (E <= tab->E[0])
                    {
                        micro = tab->xs[0];
                    }
                    else if (E >= tab->E[tab->n - 1])
                    {
                        micro = tab->xs[tab->n - 1];
                    }
                    else
                    {
                        size_t lo = 0, hi = tab->n - 1;
                        while (hi - lo > 1)
                        {
                            size_t mid = (lo + hi) >> 1;
                            if (tab->E[mid] <= E)
                                lo = mid;
                            else
                                hi = mid;
                        }
                        micro = linlin(tab->E[lo], tab->xs[lo], tab->E[lo + 1], tab->xs[lo + 1], E);
                    }

                    sum += mn->N_i * micro * BARN_TO_CM2;
                }

                macro_values[j] = sum;
            }

            /* Drop leading zero-valued points to avoid degenerate tables */
            
            size_t start = 0;
            while (start + 1 < unique_count && macro_values[start] == 0.0)
                start++;
            if (start > 0)
            {
                size_t new_n = unique_count - start;
                memmove(E_union, E_union + start, new_n * sizeof(double));
                memmove(macro_values, macro_values + start, new_n * sizeof(double));
                unique_count = new_n;
            }

            MacroXsTable *dest = &macro[produced++];
            dest->mt = mt;
            dest->n = unique_count;
            dest->E = E_union;
            dest->xs = macro_values;
        }

        free(mt_list);

        if (produced == 0) 
        {
            free(macro);
            fprintf(stdout, "\nMaterial \"%s\": no macroscopic cross sections generated.\n", M->name);
            M->macro_xs = NULL;
            M->n_macro_xs = 0;
            continue;
        }

        if (produced < mt_count) 
        {
            MacroXsTable *tmp = (MacroXsTable*)realloc(macro, produced * sizeof(MacroXsTable));
            if (!tmp) {
                free(macro);
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
            macro = tmp;
        }

        M->macro_xs = macro;
        M->n_macro_xs = produced;

        if (VERBOSITY >= 2)
        {
            fprintf(stdout, "\n  Material \"%s\" macroscopic tables:\n", M->name);
            for (size_t i = 0; i < M->n_macro_xs; ++i) 
            {
                const MacroXsTable *tab = &M->macro_xs[i];
                fprintf(stdout, "    MT = %4d, points = %zu\n", tab->mt, tab->n);
            }
        }
        else
        {
            fprintf(stdout, "  Material \"%s\": %zu macroscopic cross section groups.\n", M->name, M->n_macro_xs);
        }
    }

    
    fprintf(stdout, "\nDONE.\n");

    /* Calculate memory usage for xsdata in DATA */

    size_t mem_bytes = sizeof(DATA);

    /* Loop over everything to count size */

    for (size_t im = 0; im < DATA.n_mats; im++) 
    {
        Material *M = &DATA.mats[im];
        mem_bytes += sizeof(*M);
        mem_bytes += M->n_nucs * sizeof(MaterialNuclide);
        for (size_t c = 0; c < M->n_nucs; c++) 
        {
            MaterialNuclide *mc = &M->nucs[c];
            mem_bytes += sizeof(*mc);

            /* Add memory for cross section tables */

            mem_bytes += mc->nuc_data.n_xs * sizeof(XsTable);
            for (size_t k = 0; k < mc->nuc_data.n_xs; ++k) 
            {
                XsTable *tab = &mc->nuc_data.xs[k];
                mem_bytes += tab->n * sizeof(double);
                mem_bytes += tab->n * sizeof(double);
            }

            /* Add memory for nubar if present */

            if (mc->nuc_data.has_nubar && mc->nuc_data.nubar.n > 0) 
            {
                mem_bytes += mc->nuc_data.nubar.n * sizeof(double);
                mem_bytes += mc->nuc_data.nubar.n * sizeof(double);
            }
        }
        for (size_t k = 0; k < M->n_macro_xs; k++) 
        {
            MacroXsTable *tab = &M->macro_xs[k];
            mem_bytes += sizeof(*tab);
            mem_bytes += tab->n * sizeof(double);
            mem_bytes += tab->n * sizeof(double);
        }
    }

    /* Print memory footprint */

    GLOB.mem_xsdata = mem_bytes;
    fprintf(stdout, "\nMemory allocated for XS-data: %.2f MB\n", mem_bytes / (1024.0 * 1024.0));

    /* ########################################################################################## */
    /* Write MATLAB output for verification */

    if (DATA.n_mats > 0) 
    {
        char mpath[MAX_PATH];
        snprintf(mpath, sizeof mpath, "%s_macroxs.m", GLOB.inputfname);

        FILE *m = fopen(mpath, "w");
        if (!m) 
        {
            fprintf(stderr, "[ERROR] Could not open MATLAB macro output file \"%s\" for writing.\n", mpath);
        } 
        else 
        {
            fprintf(m, "%% Auto-generated macroscopic cross sections\n");
            fprintf(m, "%% Units: E in MeV, Sigma in 1/cm\n\n");

            for (size_t im = 0; im < DATA.n_mats; ++im) 
            {
                const Material *M = &DATA.mats[im];

                for (size_t k = 0; k < M->n_macro_xs; ++k) 
                {
                    const MacroXsTable *tab = &M->macro_xs[k];
                    char varname[256];
                    snprintf(varname, sizeof varname, "macro_mat_%s_MT%d", M->name, tab->mt);

                    fprintf(m, "%% Material: %s (index=%zu), MT=%d, points=%zu\n", M->name, im, tab->mt, tab->n);
                    fprintf(m, "%s = [\n", varname);
                    for (size_t j = 0; j < tab->n; ++j)
                        fprintf(m, "  %.16e  %.16e\n", tab->E[j], tab->xs[j]);
                    fprintf(m, "];\n\n");
                }
            }

            fclose(m);
            fprintf(stdout, "\nWrote all macroscopic XS-tables to \"%s\".\n", mpath);
        }
    }
    return 0;
}

/* ############################################################################################## */

static const XsTable *findXsTable(const Nuclide *nu, int mt) {
    size_t map_size = sizeof(nu->mt_idx) / sizeof(nu->mt_idx[0]);
    if (mt >= 0 && (size_t)mt < map_size) {
        int idx = nu->mt_idx[mt];
        if (idx >= 0 && (size_t)idx < nu->n_xs)
            return &nu->xs[idx];
    }

    for (size_t k = 0; k < nu->n_xs; ++k) {
        if (nu->xs[k].mt == mt)
            return &nu->xs[k];
    }
    return NULL;
}
