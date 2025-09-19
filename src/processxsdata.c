#include "header.h"

void *xrealloc(void *p, size_t nbytes);

int parseXsDataFile(const char *abspath, double Tlist, Nuclide *out,
                    int *headerZ, int *headerA, double *headerAW,
                    int *has_nubar, int *n_react_modes);

int processXsData(void) {

    /* Save current working directory */

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

    /* Open the xslist file and change working directory */

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

    /* Skip first line in xslist */

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
        /* This way different temperatures of the same nuclide are grouped */

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
            "Parsed xsdata: ZA = %ld (%s), T = %.1fK, AW = %.3f, nubar = %s, MTs = %d, MT = 1 calculated on union grid.\n",
            ZA, N->name,
            var.temp,
            N->AW,
            has_nubar ? "yes" : "no",
            nmods
        );

        for (size_t k = 0; k < var.n_xs; ++k) 
        {
            XsTable *tab = &var.xs[k];
            fprintf(stdout, "  MT = %3d, Q = %13.4E, points = %zu\n", tab->mt, tab->Q, tab->n);
        }
        if (var.has_nubar)
            fprintf(stdout, "  Nubar grid, points = %zu\n", var.nubar.n);
        fprintf(stdout, "\n");
    }

    /* Close and restore wrkdir */

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
        fprintf(stdout, "\nWrote all MT tables (and nubar) to \"%s\".\n", mpath);
    }

    /* ########################################################################################## */
    /* Resolve all materials */
    /* This means copying the required xsdata, by picking closest temperature variant */
    /* Computes missing fractions/densities */

    fprintf(stdout, "\nProcessing materials...\n");

    for (size_t im = 0; im < DATA.n_mats; ++im) 
    {
        Material *M = &DATA.mats[im];
        fprintf(stdout, "\n\nResolving material \"%s\" with %zu nuclide(s):\n",
                M->name, (size_t)M->n_nucs);

        /* First pass: map ZA to lib entry, copy identity (name, AW), pick temperature */

        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            MaterialNuclide *mc = &M->nucs[c];
            int Z = mc->ZA / 1000;
            int A = mc->ZA % 1000;

            /* find matching nuclide in library */

            const NuclideData *N = NULL;
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
                        M->name, mc->ZA);
                exit(EXIT_FAILURE);
            }

            snprintf(mc->name, sizeof mc->name, "%s", N->name);
            mc->AW = N->AW;

            if (N->n_var == 0) 
            {
                fprintf(stderr, "[ERROR] No temperature variants found for ZA=%d.\n", mc->ZA);
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

            mc->data.temp   = src->temp;
            mc->data.n_xs     = src->n_xs;
            mc->data.xs       = (XsTable*)malloc(src->n_xs * sizeof *mc->data.xs);
            mc->data.has_nubar = src->has_nubar;
            mc->data.nubar.n   = 0;
            mc->data.nubar.E   = NULL;
            mc->data.nubar.nu  = NULL;
            if (!mc->data.xs) 
            { 
                fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                exit(EXIT_FAILURE); 
            }

            for (size_t k = 0; k < src->n_xs; ++k) 
            {
                const XsTable *ts = &src->xs[k];
                XsTable *td = &mc->data.xs[k];
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
                mc->data.has_nubar = 1;
                mc->data.nubar.n   = src->nubar.n;
                mc->data.nubar.E   = (double*)malloc(src->nubar.n * sizeof(double));
                mc->data.nubar.nu  = (double*)malloc(src->nubar.n * sizeof(double));
                if (!mc->data.nubar.E || !mc->data.nubar.nu) 
                { 
                    fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                    exit(EXIT_FAILURE); 
                }
                memcpy(mc->data.nubar.E,  src->nubar.E,  src->nubar.n * sizeof(double));
                memcpy(mc->data.nubar.nu, src->nubar.nu, src->nubar.n * sizeof(double));
            }

            /* Print per-nuclide summary */
            
            size_t n_modes = (mc->data.n_xs > 0 ? mc->data.n_xs - 1 : 0);
            fprintf(stdout, "  %d - %s with %zu reaction modes (%.1fK data linked).\n",
                    mc->ZA, mc->name, n_modes, mc->data.temp);
        }

        /* Second pass: compute fractions and number densities using AW and mdens */

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
                Abar += M->nucs[c].atom_frac * M->nucs[c].AW;

            if (M->mdens > 0.0 && Abar > 0.0) 
            {
                /* Use mass density to set absolute scale */

                const double Ntot = (M->mdens * NA) / Abar;  /* atoms/cm^3 */

                for (size_t c = 0; c < M->n_nucs; ++c) 
                {
                    M->nucs[c].N_i = M->nucs[c].atom_frac * Ntot;
                    M->nucs[c].mass_frac = (M->nucs[c].atom_frac * M->nucs[c].AW) / Abar;
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
                    M->nucs[c].mass_frac = (M->nucs[c].atom_frac * M->nucs[c].AW) / Abar;
            } 
            else
            {
                /* No absolute density use relative scale */

                for (size_t c = 0; c < M->n_nucs; ++c) M->nucs[c].N_i = M->nucs[c].atom_frac;

                /* mass_frac from relative N_i */

                double Wsum = 0.0;
                for (size_t c = 0; c < M->n_nucs; ++c) 
                    Wsum += M->nucs[c].N_i * M->nucs[c].AW;
                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].mass_frac = (Wsum > 0.0) ? (M->nucs[c].N_i * M->nucs[c].AW / Wsum) : 0.0;
            }
        } 
        else 
        {
            /* Mass fractions provided: w_i already normalized */

            if (M->mdens > 0.0) 
            {
                /* Convert w_i + rho to absolute N_i */

                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].N_i = (M->mdens * M->nucs[c].mass_frac / M->nucs[c].AW) * NA;

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
                    denom += (M->nucs[c].AW > 0.0) ? (M->nucs[c].mass_frac / M->nucs[c].AW) : 0.0;
                if (denom <= 0.0) 
                {
                    for (size_t c = 0; c < M->n_nucs; ++c) 
                        M->nucs[c].N_i = (M->nucs[c].mass_frac / fmax(M->nucs[c].AW, 1e-300));
                } 
                else 
                {
                    for (size_t c = 0; c < M->n_nucs; ++c) 
                    {
                        double xi = (M->nucs[c].mass_frac / M->nucs[c].AW) / denom;
                        M->nucs[c].atom_frac = xi;
                        M->nucs[c].N_i = xi * Ntot_target;
                    }
                }
            } 
            else 
            {
                for (size_t c = 0; c < M->n_nucs; ++c)
                    M->nucs[c].N_i = (M->nucs[c].mass_frac / fmax(M->nucs[c].AW, 1e-300));

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
            mass_sum += M->nucs[c].N_i * M->nucs[c].AW;
            Abar_sum += M->nucs[c].atom_frac * M->nucs[c].AW;
        }

        /* Compute derived densities for verification */

        const double rho_calc   = mass_sum / NA;        /* g/cm^3 */
        const double ndens_calc = Nsum * 1.0e-24;       /* atoms/b*cm */

        fprintf(stdout,
            "\n  Summary for material \"%s\":\n"
            "  dens=%.6E %s, T=%.1fK\n"
            "  components=%zu, sum(x)=%.6f, sum(w)=%.6f\n"
            "  Abar=%.6f g/mol/atom, N_tot=%.6e 1/cm^3\n"
            "  mdens_calc=%.6E g/cm^3, adens_calc=%.6E atoms/b*cm\n\n"
            "  Nuclides:\n",
            
            M->name,
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
            size_t n_modes = (mc->data.n_xs > 0 ? mc->data.n_xs - 1 : 0);
            fprintf(stdout,
                "  %5d - %5s : AW=%10.6f  afrac=%.6f  mfrac=%.6f  N_i=%.6e atoms/b*cm  (MTs=%2zu, T=%.1fK)\n",
                mc->ZA, 
                mc->name, 
                mc->AW, 
                mc->atom_frac, 
                mc->mass_frac, 
                mc->N_i * 1e-24, 
                n_modes, 
                mc->data.temp);
        }
    }

    /* Calculate memory usage for DATA */

    size_t mem_bytes = sizeof(DATA);

    /* Loop over everything to count size */

    for (size_t im = 0; im < DATA.n_mats; ++im) 
    {
        Material *M = &DATA.mats[im];
        mem_bytes += sizeof(*M);
        mem_bytes += M->n_nucs * sizeof(MaterialNuclide);
        for (size_t c = 0; c < M->n_nucs; ++c) 
        {
            MaterialNuclide *mc = &M->nucs[c];
            mem_bytes += sizeof(*mc);

            /* Add memory for cross section tables */

            mem_bytes += mc->data.n_xs * sizeof(XsTable);
            for (size_t k = 0; k < mc->data.n_xs; ++k) 
            {
                XsTable *tab = &mc->data.xs[k];
                mem_bytes += tab->n * sizeof(double);
                mem_bytes += tab->n * sizeof(double);
            }
            /* Add memory for nubar if present */

            if (mc->data.has_nubar && mc->data.nubar.n > 0) 
            {
                mem_bytes += mc->data.nubar.n * sizeof(double);
                mem_bytes += mc->data.nubar.n * sizeof(double);
            }
        }
    }

    /* Print memory footprint */

    fprintf(stdout, "\nMemory allocated for XS data: %.2f MB\n", mem_bytes / (1024.0 * 1024.0));

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

/* ############################################################################################## */
/* Functions */

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

    while(1) 
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