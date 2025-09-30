#include "header.h"

static int parseLong(const char *s, long *out);

static int parseULong(const char *s, uint64_t *out);

static int parseUInt(const char *s, uint32_t *out);

static int parseDouble(const char *s, double *out);

/* small growable array for components while parsing */
typedef struct {
    MaterialNuclide *v;
    size_t n;
    size_t cap;
} CompVec;

long readInput() {
    long np = 0l, lnum = 0l;

    /* Open file for reading */
    
    FILE *fp = fopen(GLOB.inputf, "r");
    if (fp == NULL) 
    {
        fprintf(stderr, "[ERROR] Cannot open file \"%s\".\n", GLOB.inputf);
        exit(EXIT_FAILURE);
    }

    /* ########################################################################################## */

    /* Try to parse keyword arguments and related options */

    fprintf(stdout, "\nReading input file \"%s\"...\n", GLOB.inputf);

    char line[4096];

    /* Read lines until end-of-file */

    while (fgets(line, sizeof(line), fp)) 
    {
        ++lnum;

        /* Strip comments */
        
        char *cmnt = strchr(line, '#');
        if (cmnt)
            *cmnt = '\0';

        /* Tokenize with whitespaces */

        char *tok = strtok(line, DELIMS);

        /* Skip empty lines */

        if (!tok)
            continue;
            
    /* ########################################################################################## */
    /* --- Handle keywords --- */

        /* ###################################################################################### */
        /* --- mat --- (Material definition)*/
        if (!strcmp(tok, "mat"))
        {
            /* Try to get the material parameters from the same line */

            const char *nameTok = strtok(NULL, DELIMS);
            const char *densTok = strtok(NULL, DELIMS);
            const char *tempTok = strtok(NULL, DELIMS);

            /* Parse material block */

            Material M;
            memset(&M, 0, sizeof M);

            /* Check for valid parameters and put into Material */

            snprintf(M.name, sizeof(M.name), "%s", nameTok);

            double density;
            if (!parseDouble(densTok, &density) || density == 0.0)
            {
                fprintf(stderr, "[ERROR] Bad density \"%s\" for material \"%s\" (line %ld).\n", densTok, M.name, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Negative value is mass density g/cm^3 and positive is atomic density atoms/b*cm*/

            if (density < 0.0)
            {
                M.mdens = -density;
                M.adens = 0.0;
            }
            else
            {
                M.mdens = 0.0;
                M.adens = density;
            }

             /* Temperature in Kelvin */

            if (!parseDouble(tempTok, &M.T) || M.T <= 0.0) 
            {
                fprintf(stderr, "[ERROR] Bad temperature for material \"%s\" (line %ld).\n", M.name, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Temperature in MeV */

            M.kT = M.T * BOLTZMANN;

            /* Store header line number for errors */

            long mat_lnum = lnum;

            /* Read ZA + fraction pairs until empty line or comment put into T struct */

            CompVec cv = {0};

            while (1) {
                
                /* If EOF reached break loop */

                if (!fgets(line, sizeof line, fp)) 
                    break;
                
                /* Increment line number */

                lnum++;

                /* Remove leading whitespace */

                char *s = line; 
                while (isspace((unsigned char)*s)) 
                    s++;
                
                /* Stop at empty line or comment */

                if (*s == '\0' || *s == '#') 
                    break;

                /* Parse ZA and fraction */

                char *tokZA = strtok(s, DELIMS);
                char *tokFr = strtok(NULL, DELIMS);

                if (!tokZA || !tokFr) 
                {
                    fprintf(stderr, "[ERROR] Expected \"ZA fraction\" (line %ld) in material \"%s\".\n", lnum, M.name);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                long ZA;
                if (!parseLong(tokZA, &ZA)) 
                {
                    fprintf(stderr, "[ERROR] Bad ZA \"%s\" on line %ld.\n", tokZA, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double frac;
                if (!parseDouble(tokFr, &frac)) 
                {
                    fprintf(stderr, "[ERROR] Bad fraction \"%s\" on line %ld.\n", tokFr, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                if (frac == 0.0) 
                {
                    fprintf(stderr, "[ERROR] Zero fraction is not allowed on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Initialize a new material nuclide */

                MaterialNuclide mnuc;
                memset(&mnuc, 0, sizeof mnuc);
                mnuc.N_i   = 0.0;
                mnuc.N_tot = 0.0;
                mnuc.nuc_data.ZA = (int)ZA;
                mnuc.nuc_data.Z  = (int)ZA / 1000;
                mnuc.nuc_data.A  = (int)ZA % 1000;

                /* Negative value is mass fractions and positive is atomic fractions */

                if (frac > 0.0) 
                {
                    mnuc.atom_frac = frac;
                    mnuc.mass_frac = 0.0;
                }
                else 
                {
                    mnuc.atom_frac = 0.0;
                    mnuc.mass_frac = -frac;
                }

                /* Expand vector */

                if (cv.n == cv.cap) 
                {
                    cv.cap = cv.cap ? cv.cap * 2 : 4;

                    cv.v = (MaterialNuclide*)realloc(cv.v, cv.cap * sizeof *cv.v);
                    if (!cv.v) 
                    { 
                        fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                        fclose(fp);
                        exit(EXIT_FAILURE); 
                    }
                }

                /* Append to vector and increment counter */

                cv.v[cv.n++] = mnuc;
            }

            /* If no pairs parsed exit */

            if (cv.n == 0) 
            {
                fprintf(stderr, "[ERROR] Material \"%s\" has no components (line %ld).\n", M.name, mat_lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Normalize fractions */
            /*the sum that ends up above 0 is the type in which values were given */

            double sum_atoms = 0.0, sum_w = 0.0;

            for (size_t i = 0; i < cv.n; ++i) 
            { 
                sum_atoms += cv.v[i].atom_frac; 
                sum_w += cv.v[i].mass_frac; 
            }

            /* Check for mismatch in fraction types */

            if (sum_atoms > 0.0 && sum_w > 0.0) 
            {
                fprintf(stderr, "[ERROR] Mixed fraction types in material \"%s\" (line %ld).\n", M.name, mat_lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Normalize */

            if (sum_atoms > 0.0) 
            {
                for (size_t i = 0; i < cv.n; ++i) 
                    cv.v[i].atom_frac /= sum_atoms;
            } 
            else if (sum_w > 0.0)
            {
                for (size_t i = 0; i < cv.n; ++i) 
                    cv.v[i].mass_frac /= sum_w;    
            }
            else
            {
                /* WTF? should not be here */

                fprintf(stderr, "[ERROR] Zero fractions in material \"%s\" (line %ld).\n", M.name, mat_lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Append T vector into Material */

            M.n_nucs = cv.n;
            M.nucs   = cv.v;

            /* Append to DATA.mats by reallocating the array */

            DATA.mats = (Material*)realloc(DATA.mats, (DATA.n_mats + 1) * sizeof *DATA.mats);
            if (!DATA.mats) 
            { 
                fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                exit(EXIT_FAILURE); 
            }

            /* Put in array and increment counter */

            DATA.mats[DATA.n_mats++] = M;

            if (VERBOSITY >= 2) 
            {
                /* Print summary for successfully parsed material*/

                fprintf(stdout, "Parsed material \"%s\" at %.6E %s at %.1fK with %zu nuclide(s).\n", 
                        M.name, 
                        (M.adens > 0.0) ? M.adens : M.mdens,
                        (M.adens > 0.0) ? "atoms/b*cm" : "g/cm3",
                        M.T, 
                        M.n_nucs
                    );
            }
            np++;
        }

        /* ###################################################################################### */
        /* --- src --- (External source definition)*/
        else if (!strcmp(tok, "src"))
        {
            /* Read input parameters */

            char *a1 = strtok(NULL, DELIMS);
            char *a2 = strtok(NULL, DELIMS);
            char *a3 = strtok(NULL, DELIMS);
            char *a4 = strtok(NULL, DELIMS);
            char *a5 = strtok(NULL, DELIMS);
            if (!a1 || !a2 || !a3 || !a4 || !a5) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            uint32_t src_type;
            double x, y, z, E;
            if (!parseUInt(a1, &src_type) ||!parseDouble(a2, &x) ||!parseDouble(a3, &y) ||!parseDouble(a4, &z) ||!parseDouble(a5, &E)) 
            {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Check for valid values. For now only type 1, monoenergetic point source is supported */

            if (src_type < 1 || E <= 0.0) 
            {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            if (DATA.src == NULL) 
            {
                DATA.src = (SourceDefinition*)malloc(sizeof(SourceDefinition));
                if (!DATA.src) 
                { 
                    fprintf(stderr,"[ERROR] Memory allocation failed.\n"); 
                    fclose(fp);
                    exit(EXIT_FAILURE); 
                }
            }

            /* Put source data */
            DATA.src_type = (int)src_type;

            /* Monoenergetic point source */
            if (src_type == 1) 
            {
                DATA.src->mono.E = E;
                DATA.src->mono.x = x;
                DATA.src->mono.y = y;
                DATA.src->mono.z = z;
            } 
            else 
            {
                fprintf(stderr, "[ERROR] Source type %d not implemented (line %ld).\n", src_type, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            if (VERBOSITY >= 2)
                fprintf(stdout, "Parsed type %d source at (%.2f, %.2f, %.2f) with E=%.2f MeV.\n", src_type, x, y, z, E);

            np++;
        }

        /* ###################################################################################### */        
        /* --- det --- (Detector definition for scoring reaction channels)*/
        else if (!strcmp(tok, "det"))
        {
            if (DATA.n_detectors >= (size_t)MAX_NUM_DETECTORS) 
            {
                fprintf(stderr, "[ERROR] Only %d detectors supported (line %ld).\n", MAX_NUM_DETECTORS, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            Detector *det = (Detector*)malloc(sizeof(Detector));
            if (!det)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            memset(det, 0, sizeof(Detector));

            /* Read detector name */

            char *name_tok = strtok(NULL, DELIMS);
            snprintf(det->name, sizeof(det->name), "%s", name_tok);

            /* Read options */

            bool has_dr = false;
            bool dr_all = false;
            bool has_dm = false;
            char mat_name_buf[128] = {0};
            bool has_de = false;
            double e_min = 0.0, e_max = 0.0;
            uint64_t n_bins = 0;
            EnergyBinSpacing energy_spacing = ENERGY_BIN_SPACING_LINEAR;

            char *opt = NULL;
            while ((opt = strtok(NULL, DELIMS)) != NULL)
            {
                if (!strcmp(opt, "dr"))
                {
                    char *resp_tok = strtok(NULL, DELIMS);
                    if (!resp_tok)
                    {
                        fprintf(stderr, "[ERROR] Missing detector response after 'dr' (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }
                    has_dr = true;
                    if (!strcmp(resp_tok, "all"))
                    {
                        dr_all = true;
                    }
                    else
                    {
                        fprintf(stderr, "[ERROR] Unsupported detector response '%s' (line %ld).\n", resp_tok, lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }
                }
                else if (!strcmp(opt, "dm"))
                {
                    char *mat_tok = strtok(NULL, DELIMS);
                    if (!mat_tok)
                    {
                        fprintf(stderr, "[ERROR] Missing material name after 'dm' (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }
                    snprintf(mat_name_buf, sizeof(mat_name_buf), "%s", mat_tok);
                    has_dm = true;
                }
                else if (!strcmp(opt, "de"))
                {
                    char *emin_tok = strtok(NULL, DELIMS);
                    char *emax_tok = strtok(NULL, DELIMS);
                    char *nbins_tok = strtok(NULL, DELIMS);
                    char *spacing_tok = strtok(NULL, DELIMS);
                    if (!emin_tok || !emax_tok || !nbins_tok)
                    {
                        fprintf(stderr, "[ERROR] Incomplete energy grid specification on line %ld.\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (!spacing_tok)
                        spacing_tok = "1";

                    if (!parseDouble(emin_tok, &e_min) || !parseDouble(emax_tok, &e_max))
                    {
                        fprintf(stderr, "[ERROR] Invalid energy bounds on line %ld.\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (!parseULong(nbins_tok, &n_bins) || n_bins == 0)
                    {
                        fprintf(stderr, "[ERROR] Invalid bin count in detector energy grid (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    long spacing_val = 0;
                    if (!parseLong(spacing_tok, &spacing_val))
                    {
                        fprintf(stderr, "[ERROR] Invalid grid spacing flag in detector energy grid (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (spacing_val == 0)
                        energy_spacing = ENERGY_BIN_SPACING_LOG;
                    else if (spacing_val == 1)
                        energy_spacing = ENERGY_BIN_SPACING_LINEAR;
                    else
                    {
                        fprintf(stderr, "[ERROR] Grid spacing flag must be 0 (log) or 1 (lin) on line %ld.\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (e_max <= e_min)
                    {
                        fprintf(stderr, "[ERROR] Detector energy max must exceed min (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }
                    if (energy_spacing == ENERGY_BIN_SPACING_LOG && e_min <= 0.0)
                    {
                        fprintf(stderr, "[ERROR] Logarithmic energy bins require Emin > 0 (line %ld).\n", lnum);
                        free(det);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    has_de = true;
                }
                else
                {
                    fprintf(stderr, "[ERROR] Unknown detector option '%s' (line %ld).\n", opt, lnum);
                    free(det);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }

            /* If detector response defined create a reaction rate detector */

            if (has_dr)
            {
                if (!dr_all)
                {
                    fprintf(stderr, "[ERROR] Only 'dr all' is supported currently (line %ld).\n", lnum);
                    free(det);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                if (!has_dm)
                {
                    fprintf(stderr, "[ERROR] Reaction rate detectors require a material (use 'dm <name>') (line %ld).\n", lnum);
                    free(det);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                det->type = DETECTOR_TYPE_REACTION_RATE;
                ReactionRateDetector *rr = &det->data.reaction_rate;
                memset(rr, 0, sizeof(*rr));
                snprintf(rr->material_name, sizeof(rr->material_name), "%s", mat_name_buf);
                rr->material_index = -1;
                rr->energy_grid.enabled = has_de;
                if (has_de)
                {
                    rr->energy_grid.enabled = true;
                    rr->energy_grid.n_bins = (size_t)n_bins;
                    rr->energy_grid.E_min = e_min;
                    rr->energy_grid.E_max = e_max;
                    rr->energy_grid.spacing = energy_spacing;
                }
                else
                    rr->energy_grid.enabled = false;

                if (VERBOSITY >= 2)
                {
                    fprintf(stdout, "Parsed reaction-rate detector '%s' (material: %s)%s.\n",
                            det->name, rr->material_name, has_de ? " with energy grid" : "");
                }
            }
            /* Otherwise create an energy spectrum detector */

            else
            {
                det->type = DETECTOR_TYPE_ENERGY_SPECTRUM;
                EnergySpectrumDetector *es = &det->data.energy_spectrum;
                memset(es, 0, sizeof(*es));
                es->material_index = -1;
                es->has_material_filter = has_dm;
                if (has_dm)
                    snprintf(es->material_name, sizeof(es->material_name), "%s", mat_name_buf);

                if (!has_de)
                {
                    fprintf(stderr, "[ERROR] Energy spectrum detectors require an energy grid ('de').\n");
                    free(det);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                es->grid.enabled = true;
                es->grid.n_bins = (size_t)n_bins;
                es->grid.E_min = e_min;
                es->grid.E_max = e_max;
                es->grid.spacing = energy_spacing;

                if (VERBOSITY >= 2)
                {
                    const char *mat_desc = es->has_material_filter ? es->material_name : "all materials";
                    fprintf(stdout, "Parsed energy-spectrum detector '%s' (%s) with %zu bins between %E - %E\n",
                            det->name, mat_desc, es->grid.n_bins, es->grid.E_min, es->grid.E_max);
                }
            }

            DATA.detectors[DATA.n_detectors] = det;
            DATA.n_detectors++;
            np++;
        }


        /* ###################################################################################### */        
        /* --- set --- (Settings keyword followed by subwords and arguments) */

        else if (!strcmp(tok, "set"))
        {
            /* Read the  following keyword */

            char *subkey = strtok(NULL, DELIMS);
            if (!subkey) 
            {
                fprintf(stderr, "[ERROR] Missing specifier keyword for \"set\" on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* ################################################################################## */
            /* --- Handle the specifier keywords --- */

            /* --- Neutron energy cutoff */
            if (!strcmp(subkey, "ecut"))
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double ecut;
                if (!parseDouble(a1, &ecut) || ecut < 0.0) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set ecut\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.energy_cutoff = ecut;
                np++;
            }
            /* --- Collision number cutoff */
            else if (!strcmp(subkey, "ccut"))
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                double ccut;
                if (!parseDouble(a1, &ccut) || ccut < 0) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set ccut\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.collision_cutoff = (long)ccut;
                np++;
            }
            /* -- Generation cut off */
            else if (!strcmp(subkey, "gcut"))
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double gcut;
                if (!parseDouble(a1, &gcut) || gcut < 0) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set gcut\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.generation_cutoff = gcut;
                np++;
            }
            /* -- Time cut off */
            else if (!strcmp(subkey, "tcut"))
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double tcut;
                if (!parseDouble(a1, &tcut) || tcut < 0) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set tcut\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.time_cutoff = tcut;
                np++;
            }
            /* --- Cross section library file path */
            else if (!strcmp(subkey, "xslibpath"))
            {   
                char *path = strtok(NULL, DELIMS);
                if (!path) {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Check for valid length */

                size_t len = strlen(path);
                if (len >= MAX_STR_LEN) {
                    fprintf(stderr, "[ERROR] Cross section library path too long on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                snprintf(GLOB.xslibpath, MAX_STR_LEN, "%s", path);
                np++;
            }
            /* --- Criticality source simulation population */
            else if (!strcmp(subkey, "pop")) 
            {
                char *a1 = strtok(NULL, DELIMS);
                char *a2 = strtok(NULL, DELIMS);
                char *a3 = strtok(NULL, DELIMS);
                if (!a1 || !a2) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                long n_particles, n_generations, n_inactive;
                if (!parseLong(a1, &n_particles) || !parseLong(a2, &n_generations)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                if (!parseLong(a3, &n_inactive))
                    n_inactive = 0;

                /* Check for valid values */

                if (n_generations < 1 || n_particles < 1 || n_inactive < 0) {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.n_generations = n_generations;
                GLOB.n_particles = n_particles;
                GLOB.n_inactive = n_inactive;
                if (GLOB.mode != RUNMODE_CHECK)
                    GLOB.mode = RUNMODE_CRITICALITY;
                np++;
            }
            /* --- RNG seed */
            else if (!strcmp(subkey, "seed")) 
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                uint64_t seed;
                if (!parseULong(a1, &seed)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.seed = seed;
                np++;
            }
            else 
            {
                fprintf(stderr, "[WARNING] Unknown sub-keyword \"%s\" for \"set\" on line %ld. Skipping...\n", subkey, lnum);
            }
        }

        /* ###################################################################################### */
    
        else
        {
            /* Some unknown keyword, print warning and continue to next line */

            fprintf(stderr, "[WARNING] Skipping line %ld with unknown keyword argument: \"%s\"...\n", lnum, tok);
        }
    }

    /* ########################################################################################## */
    /* Return number of arguments succesfully parsed */

    return np;
}

/**
 * @brief Parse a signed 64-bit integer from provided str
 * 
 * @param s str to parse from
 * @param out ptr where parsed result is put
 * @return int 0 on failure 1 on success
 */
static int parseLong(const char *s, long *out) {
    if (!s || !out)
        return 0;
    errno = 0;
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (errno || end == s || *end != '\0') 
        return 0;
    *out = v;
    return 1;
}

/**
 * @brief Parse a unsigned 64-bit integer from provided str
 * 
 * @param s str to parse from
 * @param out ptr where parsed result is put
 * @return int 0 on failure 1 on success
 */
static int parseULong(const char *s, uint64_t *out) {
    if (!s || !out)
        return 0;
    errno = 0;
    char *end = NULL;
    unsigned long v = strtoul(s, &end, 10);
    if (errno || end == s || *end != '\0') return 0;
    *out = (uint64_t)v;
    return 1;
}

/**
 * @brief Parse a unsigned 32-bit integer from provided str
 * 
 * @param s str to parse from
 * @param out ptr where parsed result is put
 * @return int 0 on failure 1 on success
 */
static int parseUInt(const char *s, uint32_t *out) {
    if (!s || !out)
        return 0;
    errno = 0;
    char *end = NULL;
    unsigned long v = strtoul(s, &end, 10);
    if (errno || end == s || *end != '\0') return 0;
    *out = (uint32_t)v;
    return 1;
}

/**
 * @brief Parse a 64-bit floating point number from provided str
 * 
 * @param s str to parse from
 * @param out ptr where parsed result is put
 * @return int 0 on failure 1 on success
 */
static int parseDouble(const char *s, double *out) {
    if (!s || !out)
        return 0;
    errno = 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if (errno || end == s || *end != '\0') 
        return 0;
    *out = v;
    return 1;
}
