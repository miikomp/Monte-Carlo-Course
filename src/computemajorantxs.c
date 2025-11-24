#include "header.h"

static double evalMacroXsAtEnergy(const MacroXsTable *tab, double E)
{
    if (!tab || tab->n == 0 || !tab->E || !tab->xs)
        return -1.0;

    if (E <= tab->E[0])
        return tab->xs[0];

    if (E >= tab->E[tab->n - 1])
        return tab->xs[tab->n - 1];

    size_t lo = 0, hi = tab->n - 1;
    while (hi - lo > 1)
    {
        size_t mid = (lo + hi) >> 1;
        if (tab->E[mid] <= E)
            lo = mid;
        else
            hi = mid;
    }

    return linlin(tab->E[lo], tab->xs[lo],
                  tab->E[lo + 1], tab->xs[lo + 1], E);
}

int computeMajorantXS(void)
{
    fprintf(stdout, "\nComputing majorant cross section...\n");

    DATA.majorant_xs.E  = NULL;
    DATA.majorant_xs.xs = NULL;
    DATA.majorant_xs.n  = 0;

    if (DATA.n_mats == 0)
    {
        fprintf(stderr, "[ERROR] No materials available to construct a majorant.\n");
        return EXIT_FAILURE;
    }

    double *E_union   = NULL;
    size_t  n_union   = 0;
    size_t  cap_union = 0;
    size_t  mats_used = 0;

    /* Build the union energy grid from all material total macroscopic tables */
    for (size_t im = 0; im < DATA.n_mats; ++im)
    {
        const Material *mat = &DATA.mats[im];
        if (!mat->macro_xs || mat->n_macro_xs == 0)
        {
            if (VERBOSITY >= 1)
                fprintf(stdout, "  Material \"%s\" has no macroscopic XS. Skipping.\n", mat->name);
            continue;
        }

        const MacroXsTable *total = &mat->macro_xs[0]; /* total is always index 0 */
        if (total->n == 0 || !total->E)
        {
            if (VERBOSITY >= 1)
                fprintf(stdout, "  Material \"%s\" has no total macroscopic XS. Skipping.\n", mat->name);
            continue;
        }

        mats_used++;

        if (cap_union < n_union + total->n)
        {
            size_t new_cap = cap_union ? cap_union : total->n;
            while (new_cap < n_union + total->n)
                new_cap *= 2;

            double *tmp = (double*)realloc(E_union, new_cap * sizeof(double));
            if (!tmp)
            {
                free(E_union);
                fprintf(stderr, "[ERROR] Memory allocation failed while building majorant grid.\n");
                exit(EXIT_FAILURE);
            }
            E_union   = tmp;
            cap_union = new_cap;
        }

        memcpy(&E_union[n_union], total->E, total->n * sizeof(double));
        n_union += total->n;
    }

    if (mats_used == 0 || n_union == 0)
    {
        free(E_union);
        fprintf(stderr, "[ERROR] Failed to assemble union energy grid for majorant cross section.\n");
        return EXIT_FAILURE;
    }

    /* Sort and uniquify the union grid */
    qsort(E_union, n_union, sizeof(double), cmpDouble);

    size_t unique_count = 0;
    for (size_t i = 0; i < n_union; ++i)
    {
        if (i == 0)
        {
            E_union[unique_count++] = E_union[i];
            continue;
        }

        double a   = E_union[i];
        double b   = E_union[unique_count - 1];
        double tol = 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));

        if (fabs(a - b) > tol)
            E_union[unique_count++] = a;
    }

    if (unique_count == 0)
    {
        free(E_union);
        fprintf(stderr, "[ERROR] Union grid for majorant is empty after deduplication.\n");
        return EXIT_FAILURE;
    }

    double *E_final = (double*)realloc(E_union, unique_count * sizeof(double));
    if (!E_final)
    {
        free(E_union);
        fprintf(stderr, "[ERROR] Memory allocation failed while finalising majorant grid.\n");
        exit(EXIT_FAILURE);
    }
    E_union = E_final;

    double *sigma_M = (double*)malloc(unique_count * sizeof(double));
    if (!sigma_M)
    {
        free(E_union);
        fprintf(stderr, "[ERROR] Memory allocation failed for majorant cross section values.\n");
        exit(EXIT_FAILURE);
    }

    for (size_t j = 0; j < unique_count; ++j)
    {
        double E = E_union[j];
        double sigma_max = 0.0;
        bool   have_val  = false;

        for (size_t im = 0; im < DATA.n_mats; ++im)
        {
            const Material *mat = &DATA.mats[im];
            if (!mat->macro_xs || mat->n_macro_xs == 0)
                continue;

            const MacroXsTable *total = &mat->macro_xs[0];
            if (total->n == 0 || !total->E)
                continue;

            double sigma_t = evalMacroXsAtEnergy(total, E);
            if (sigma_t < 0.0)
                continue;

            if (!have_val || sigma_t > sigma_max)
                sigma_max = sigma_t;

            have_val = true;
        }

        if (!have_val)
        {
            free(E_union);
            free(sigma_M);
            fprintf(stderr, "[ERROR] Majorant evaluation failed at energy %.6e MeV.\n", E);
            return EXIT_FAILURE;
        }

        sigma_M[j] = sigma_max;
    }

    DATA.majorant_xs.n  = unique_count;
    DATA.majorant_xs.E  = E_union;
    DATA.majorant_xs.xs = sigma_M;

    fprintf(stdout, "  Union grid points: %zu (materials contributing: %zu)\n", unique_count, mats_used);

    /* Append majorant to MATLAB macroxs output for verification */
    if (GLOB.inputfname && GLOB.inputfname[0] != '\0')
    {
        char mpath[MAX_PATH];
        snprintf(mpath, sizeof mpath, "%s_macroxs.m", GLOB.inputfname);

        FILE *m = fopen(mpath, "a");
        if (!m)
        {
            fprintf(stderr, "[WARNING] Could not open \"%s\" to append majorant XS.\n", mpath);
        }
        else
        {
            fprintf(m, "%% Majorant cross section (union grid across all materials)\n");
            fprintf(m, "majorant_sigma = [\n");
            for (size_t j = 0; j < unique_count; ++j)
                fprintf(m, "  %.16e  %.16e\n", E_union[j], sigma_M[j]);
            fprintf(m, "];\n\n");
            fclose(m);
        }
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}
