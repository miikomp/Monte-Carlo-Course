#include "header.h"

void *xrealloc(void *p, size_t nbytes);

int parseXsDataFile(const char *abspath, double Tlist, Nuclide *out,
                    int *headerZ, int *headerA, double *headerAW,
                    int *has_nubar, int *n_react_modes);

int processXsData(void) {

    char oldcwd[MAX_PATH];
    if (!getcwd(oldcwd, sizeof oldcwd)) 
    {
        fprintf(stderr, "[ERROR] getcwd failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }

    char list_dir[MAX_PATH]  = {0};
    char list_base[MAX_PATH] = {0};

    size_t len = strlen(GLOB.xslibpath);
    size_t slash_pos = len;
    for (size_t i = 0; i < len; ++i) 
    {
        char c = GLOB.xslibpath[i];
        if (c == '/' || c == '\\') 
            slash_pos = i; 
    }

    if (slash_pos == len) 
        snprintf(list_base, sizeof list_base, "%s", GLOB.xslibpath);
    else 
    {
        size_t dlen = slash_pos;
        if (dlen >= sizeof list_dir) dlen = sizeof list_dir - 1;
        memcpy(list_dir, GLOB.xslibpath, dlen);
        list_dir[dlen] = '\0';
        snprintf(list_base, sizeof list_base, "%s", GLOB.xslibpath + slash_pos + 1);
    }

    int changed_dir = 0;
    if (list_dir[0] != '\0') 
    {
        if (chdir(list_dir) != 0) 
        {
            fprintf(stderr, "[ERROR] chdir(\"%s\") failed: %s\n", list_dir, strerror(errno));
            return EXIT_FAILURE;
        }
        changed_dir = 1;
    }

    /* Open the xslist file */

    FILE *list = fopen(list_base, "r");
    if (!list) 
    {
        fprintf(stderr, "[ERROR] Could not open cross section library file at \"%s\".\n",
                GLOB.xslibpath);
        if (changed_dir) 
            if (chdir(oldcwd) != 0)
                fprintf(stderr, "[ERROR] chdir(\"%s\") failed: %s\n", oldcwd, strerror(errno));
        return EXIT_FAILURE;
    }

    /* ########################################################################################## */
    /* First parse all the xsdata into a temporary library */

    fprintf(stdout, "\nProcessing cross section library \"%s\"...\n\n", GLOB.xslibpath);

    NuclideData *lib = NULL;
    size_t nlib = 0;

    char line[4096];
    long lnum = 0;

    /* Skip first line */
    if (!fgets(line, sizeof line, list)) 
    {
        fprintf(stderr, "[ERROR] Empty cross section library file \"%s\".\n", GLOB.xslibpath);
        fclose(list);
        return EXIT_FAILURE;
    }

    ++lnum;

    /* Read lines */

    while (fgets(line, sizeof line, list)) 
    {
        ++lnum;

        /* skip empty/whitespace-only lines */

        char *p = line;
        while (isspace((unsigned char)*p)) ++p;
        if (*p == '\0') continue;

        /* ZA  Name  T[K]  path */
        
        char *tokZA  = strtok(p,     DELIMS);
        char *tokNm  = strtok(NULL,  DELIMS);
        char *tokT   = strtok(NULL,  DELIMS);
        char *tokPth = strtok(NULL,  DELIMS);

        if (!tokZA || !tokNm || !tokT || !tokPth) 
        {
            fprintf(stderr, "[ERROR] xslist line %ld: expected \"ZA Name T[K] path\".\n", lnum);
            fclose(list);
            return 1;
        }

        char *end = NULL;
        long ZA = strtol(tokZA, &end, 10);
        if (!end || *end) 
        {
            fprintf(stderr, "[ERROR] xslist line %ld: bad ZA \"%s\".\n", lnum, tokZA);
            fclose(list);
            return 1;
        }
        double Tlist = strtod(tokT, &end);
        if (!end || *end) 
        {
            fprintf(stderr, "[ERROR] xslist line %ld: bad temperature \"%s\".\n", lnum, tokT);
            fclose(list);
            return 1;
        }
        const char *name = tokNm;
        const char *path = tokPth;

        /* Find or create the species record by ZA names assumed same, first is kept */

        int Z = (int)ZA / 1000;
        int A = (int)ZA % 1000;
        NuclideData *N = NULL;
        for (size_t i = 0; i < nlib; ++i) {
            if (lib[i].Z == Z && lib[i].A == A) {
                N = &lib[i];
                break;
            }
        }
        if (!N) {
            lib = (NuclideData*)xrealloc(lib, (nlib + 1) * sizeof *lib);
            N = &lib[nlib++];
            memset(N, 0, sizeof *N);
            snprintf(N->name, sizeof N->name, "%s", name ? name : "");
            N->Z  = Z;
            N->A  = A;
            N->AW = 0.0;
            N->n_var = 0;
            N->var   = NULL;
        }

        /* Parse the file into a nuclide with xsdata and calculate MT=1 at xs[0] */

        Nuclide var;
        int hZ=0, hA=0, has_nubar=0, nmods=0;
        double hAW=0.0;

        if (parseXsDataFile(path, Tlist, &var, &hZ, &hA, &hAW, &has_nubar, &nmods) != 0) {
            fprintf(stderr, "[ERROR] Failed parsing \"%s\" (xslist line %ld).\n", path, lnum);
            fclose(list);

            /* Free everything */
            
            for (size_t i = 0; i < nlib; ++i) 
            {
                for (size_t v = 0; v < lib[i].n_var; ++v) 
                {
                    Nuclide *vv = &lib[i].var[v];
                    for (size_t k = 0; k < vv->n_xs; ++k) 
                    {
                        free(vv->xs[k].E);
                        free(vv->xs[k].xs);
                    }
                    free(vv->xs);
                    if (vv->has_nubar) { free(vv->nubar.E); free(vv->nubar.nu); }
                }
                free(lib[i].var);
            }
            free(lib);
            return 1;
        }

        if (N->AW == 0.0) 
            N->AW = hAW;
        else if (fabs(N->AW - hAW) > 1e-8) 
        {
            fprintf(stderr, "[WARN] AW mismatch for ZA=%ld: existing=%.8f, file=%.8f (%s)\n",
                    ZA, N->AW, hAW, path);
        }

        /* Append this temperature variant */

        N->var = (Nuclide*)xrealloc(N->var, (N->n_var + 1) * sizeof *N->var);
        N->var[N->n_var++] = var;

        /* Print summary line */

        fprintf(stdout,
            "[NOTE] Parsed xsdata: ZA=%ld (%s), T=%.1f K, AW=%.3f, nubar=%s, MTs=%d, MT=1 with %zu points.\n",
            ZA, N->name,
            var.temp,
            N->AW,
            has_nubar ? "yes" : "no",
            nmods,
            var.xs[0].n);
    }

    /* Close and restore wrkdir*/

    fclose(list);
    if (changed_dir) 
        if (chdir(oldcwd) != 0) 
            fprintf(stderr, "[ERROR] Failed to restore working directory to \"%s\": %s\n", oldcwd, strerror(errno));
    
    /* ########################################################################################## */
    /* Write all MTs computed/parsed into file for verification */

    char mpath[MAX_PATH];
    snprintf(mpath, sizeof mpath, "%s_xs.m", GLOB.fname);

    FILE *m = fopen(mpath, "w");
    if (!m) 
    {
        fprintf(stderr, "[ERROR] Could not open MATLAB output file \"%s\" for writing.\n", mpath);
    } 
    else 
    {
        fprintf(m, "%% Auto-generated cross sections\n");
        fprintf(m, "%% Units: E in MeV, sigma in barns; nubar is dimensionless\n\n");

        for (size_t i = 0; i < nlib; ++i) 
        {
            const NuclideData *N = &lib[i];
            const int ZA = N->Z * 1000 + N->A;

            for (size_t v = 0; v < N->n_var; ++v) 
            {
                const Nuclide *var = &N->var[v];
                const int Tk = (int)lrint(var->temp);

                /* All MTs, including computed total */
                for (size_t k = 0; k < var->n_xs; ++k) 
                {
                    const XsTable *tab = &var->xs[k];
                    char varname[128];
                    snprintf(varname, sizeof varname, "xs_%d_MT%d_%dK", ZA, tab->mt, Tk);

                    fprintf(m, "%% ZA=%d (%s), T=%g K, MT=%d, Q=%.8g, points=%zu\n",
                            ZA, N->name, (double)var->temp, tab->mt, tab->Q, tab->n);
                    fprintf(m, "%s = [\n", varname);
                    for (size_t j = 0; j < tab->n; ++j)
                        fprintf(m, "  %.16e  %.16e\n", tab->E[j], tab->xs[j]);
                    fprintf(m, "];\n\n");
                }

                /* NUBAR if present */
                if (var->has_nubar && var->nubar.n > 0 && var->nubar.E && var->nubar.nu) 
                {
                    char vnn[128];
                    snprintf(vnn, sizeof vnn, "nubar_%d_%dK", ZA, Tk);
                    fprintf(m, "%% ZA=%d (%s), T=%g K, nubar points=%zu\n",
                            ZA, N->name, (double)var->temp, var->nubar.n);
                    fprintf(m, "%s = [\n", vnn);
                    for (size_t j = 0; j < var->nubar.n; ++j)
                        fprintf(m, "  %.16e  %.16e\n", var->nubar.E[j], var->nubar.nu[j]);
                    fprintf(m, "];\n\n");
                }
            }
        }
        fclose(m);
        fprintf(stdout, "[NOTE] Wrote all MT tables (and nubar) to \"%s\".\n", mpath);
    }

    /* ########################################################################################## */
    /* Resolve all materials */    

    fprintf(stdout, "\nResolving materials...\n\n");

    /* Free everything */

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

/**
 * @brief Parse one xs data file into a Nuclide variant. Computes MT=1 on a union grid.
 * Returns number of reaction modes parsed and whether nubar was present.
 * 
 * @param abspath path to file
 * @param Tlist temperature listed in the file
 * @param out Nuclide to fill
 * @param headerZ output Z
 * @param headerA output A
 * @param headerAW output atomic weight
 * @param has_nubar output whether nubar was present
 * @param n_react_modes output number of reaction modes parsed
 * @return int 
 */
int parseXsDataFile(const char *abspath, double Tlist,
                             Nuclide *out,
                             int *headerZ, int *headerA, double *headerAW,
                             int *has_nubar, int *n_react_modes)
{
    FILE *f = fopen(abspath, "r");
    if (!f) 
    {
        fprintf(stderr, "[ERROR] Cannot open xs data file: %s\n", abspath);
        return -1;
    }

    /* Header line */

    char sym[64] = {0};
    int Z=0, A=0; double AW=0.0, Tfile=0.0;
    if (fscanf(f, "%63s %d %d %lf %lf", sym, &Z, &A, &AW, &Tfile) != 5) 
    {
        fprintf(stderr, "[ERROR] Bad header in %s\n", abspath);
        fclose(f);
        return -1;
    }
    *headerZ = Z; *headerA = A; *headerAW = AW;

    /* NNU line */
    int NNU = 0;
    if (fscanf(f, "%d", &NNU) != 1) 
    {
        fprintf(stderr, "[ERROR] Bad NNU line in %s\n", abspath);
        fclose(f);
        return -1;
    }

    out->temp = Tlist;
    out->n_xs   = 0;
    out->xs     = NULL;
    out->has_nubar = (NNU > 0) ? 1 : 0;
    out->nubar.n = 0; out->nubar.E = out->nubar.nu = NULL;

    if (NNU > 0) 
    {
        out->nubar.n = (size_t)NNU;
        out->nubar.E  = (double*)xrealloc(NULL, (size_t)NNU * sizeof(double));
        out->nubar.nu = (double*)xrealloc(NULL, (size_t)NNU * sizeof(double));
        for (int i = 0; i < NNU; ++i) 
        {
            if (fscanf(f, "%lf %lf", &out->nubar.E[i], &out->nubar.nu[i]) != 2) 
            {
                fprintf(stderr, "[ERROR] Bad NNU pair in %s\n", abspath);
                fclose(f);
                return -1;
            }
        }
    }

    /* Read all reaction blocks (MT Q NE + NE pairs) */

    int parsed_modes = 0;

    for (;;) 
    {
        int mt = 0, ne = 0;
        double Q = 0.0;

        /* Peek next token if EOF break */

        int r = fscanf(f, "%d %lf %d", &mt, &Q, &ne);
        if (r == EOF || r == 0) 
            break;      
        if (r != 3) 
        {
            fprintf(stderr, "[ERROR] Bad MT block header in %s (expected: MT Q NE)\n", abspath);
            fclose(f);
            return -1;
        }
        if (ne < 2) 
        {
            fprintf(stderr, "[ERROR] NE < 2 for MT=%d in %s\n", mt, abspath);
            fclose(f);
            return -1;
        }

        XsTable xt;
        xt.mt = mt; xt.Q = Q; xt.n = (size_t)ne;
        xt.E  = (double*)xrealloc(NULL, (size_t)ne * sizeof(double));
        xt.xs = (double*)xrealloc(NULL, (size_t)ne * sizeof(double));

        for (int i = 0; i < ne; ++i) 
        {
            if (fscanf(f, "%lf %lf", &xt.E[i], &xt.xs[i]) != 2) 
            {
                fprintf(stderr, "[ERROR] Bad E-xs pair for MT=%d in %s\n", mt, abspath);
                fclose(f);
                return -1;
            }
            if (i > 0 && !(xt.E[i] > xt.E[i-1])) 
            {
                fprintf(stderr, "[ERROR] Non-ascending energies for MT=%d in %s\n", mt, abspath);
                fclose(f);
                return -1;
            }
            if (xt.xs[i] < 0.0) 
            {
                fprintf(stderr, "[ERROR] Negative cross section for MT=%d in %s\n", mt, abspath);
                fclose(f);
                return -1;
            }
        }

        /* Append */

        out->xs = (XsTable*)xrealloc(out->xs, (out->n_xs + 1) * sizeof *out->xs);
        out->xs[out->n_xs++] = xt;
        parsed_modes++;
    }

    fclose(f);

    /* ---- Compute MT=1 store it at xs[0] ---- */

    /* Figure out the size of the union grid needed */
    /* Count all breakspoints in reaction modes */
    size_t total_pts = 0;
    for (size_t k = 0; k < out->n_xs; ++k) 
    {
        if (out->xs[k].mt == 1) 
            continue;
        total_pts += out->xs[k].n;
    }

    /* If there are no reaction modes do nothing */
    if (total_pts > 0) {
        double *Eunion = (double*)xrealloc(NULL, total_pts * sizeof(double));
        size_t m = 0;
        for (size_t k = 0; k < out->n_xs; ++k) 
        {
            if (out->xs[k].mt == 1) 
                continue;
            memcpy(&Eunion[m], out->xs[k].E, out->xs[k].n * sizeof(double));
            m += out->xs[k].n;
        }

        /* Sort and build */

        qsort(Eunion, m, sizeof(double), cmpDouble);
        size_t u = 0;
        for (size_t i = 0; i < m; ++i) 
        {
            if (i == 0) 
            { 
                Eunion[u++] = Eunion[i]; 
                continue; 
            }

            double a = Eunion[i], b = Eunion[u-1];
            double tol = 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));
            if (fabs(a - b) > tol) Eunion[u++] = a;
        }

        /* Sum partials at each union energy using linear interpolation */

        double *xs_tot = (double*)xrealloc(NULL, u * sizeof(double));
        for (size_t j = 0; j < u; ++j) 
        {
            double Ej = Eunion[j];
            double sum = 0.0;

            for (size_t k = 0; k < out->n_xs; ++k) 
            {
                if (out->xs[k].mt == 1) 
                    continue;

                const double *x = out->xs[k].E;
                const double *y = out->xs[k].xs;
                size_t n = out->xs[k].n;

                double val;
                if (Ej <= x[0])       
                    val = y[0];
                else if (Ej >= x[n-1]) 
                    val = y[n-1];
                else 
                {
                    /* Search to find bounds to interpolate within */

                    size_t lo = 0, hi = n - 1;
                    while (hi - lo > 1) 
                    {
                        size_t mid = (lo + hi) >> 1;
                        if (x[mid] <= Ej) 
                            lo = mid; 
                        else 
                            hi = mid;
                    }

                    /* Interpolate */

                    val = linlin(x[lo], y[lo], x[lo+1], y[lo+1], Ej);
                }
                sum += val;
            }

            /* Clamp */

            if (sum >= 0.0 || sum < -1e-12)
                xs_tot[j] = sum;
            else
                xs_tot[j] = 0.0;
        }

        /* Prepend MT=1 table at index 0 (keep partials after it) */

        XsTable *newxs = (XsTable*)xrealloc(NULL, (out->n_xs + 1) * sizeof *newxs);

        newxs[0].mt = 1;
        newxs[0].Q  = 0.0;
        newxs[0].n  = u;
        newxs[0].E  = Eunion;
        newxs[0].xs = xs_tot;

        for (size_t k = 0; k < out->n_xs; ++k) newxs[k+1] = out->xs[k];

        free(out->xs);
        out->xs   = newxs;
        out->n_xs = out->n_xs + 1;
    }

    *has_nubar      = (NNU > 0) ? 1 : 0;
    *n_react_modes  = parsed_modes;

    return 0;
}

/**
 * @brief Helper realloc that exits on failure.
 * 
 * @param p ptr to realloc 
 * @param nbytes size in bytes
 * @return void* ptr to reallocated memory 
 */
void *xrealloc(void *p, size_t nbytes) {
    void *q = realloc(p, nbytes);
    if (!q) { fprintf(stderr, "[ERROR] Memory allocation failed.\n"); exit(EXIT_FAILURE); }
    return q;
}