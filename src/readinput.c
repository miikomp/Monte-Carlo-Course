#include "header.h"

#define MAX_N_PARAMS 256

static int parseLong(const char *s, long *out);

static int parseULong(const char *s, uint64_t *out);

static int parseUInt(const char *s, uint32_t *out);

static int parseDouble(const char *s, double *out);

static void readLatticeUniverseNames(FILE *fp, long *lnum, Lattice *lat);

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

    /* Standard rand used for colour randomization */
    srand((unsigned int)time(NULL));

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

            if (!nameTok || !densTok || !tempTok) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Parse material block */

            Material M;
            memset(&M, 0, sizeof M);

            /* Check for valid parameters and put into Material */
            if (!strcmp(nameTok, "outside"))
            {
                fprintf(stderr, "[ERROR] Material name cannot be \"outside\" (line %ld).\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

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

            /* See if colour is provided, randomise if not */
            char *rgbtok = strtok(NULL, DELIMS);
            if (rgbtok && !strcmp(rgbtok, "rgb"))
            {
                /* Read three colour channels */
                
                char *r_tok = strtok(NULL, DELIMS);
                char *g_tok = strtok(NULL, DELIMS);
                char *b_tok = strtok(NULL, DELIMS);

                int r, g, b;

                if (!parseUInt(r_tok, (uint32_t*)&r) || r < 0 || r > 255 ||
                    !parseUInt(g_tok, (uint32_t*)&g) || g < 0 || g > 255 ||
                    !parseUInt(b_tok, (uint32_t*)&b) || b < 0 || b > 255)
                {
                    fprintf(stderr, "[ERROR] Bad RGB values for material \"%s\" (line %ld).\n", M.name, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                M.rgb[0] = r;
                M.rgb[1] = g;
                M.rgb[2] = b;

                if (VERBOSITY >= 2)
                    fprintf(stdout, "Material \"%s\" colour set to RGB(%u,%u,%u).\n", M.name, r, g, b);
            }
            else
            {
                /* Randomize all three channels */

                M.rgb[0] = (uint32_t)(rand() % 256);
                M.rgb[1] = (uint32_t)(rand() % 256);
                M.rgb[2] = (uint32_t)(rand() % 256);

                if (VERBOSITY >= 2)
                    fprintf(stdout, "Material \"%s\" colour randomized to RGB(%u,%u,%u).\n", M.name, M.rgb[0], M.rgb[1], M.rgb[2]);
            }

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

            if (!a1 || !a2) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            uint32_t src_type;
            if (!parseUInt(a2, &src_type)) 
            {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Allocate new source */

            if (DATA.src == NULL) 
            {
                DATA.src = (ParticleSource*)malloc(sizeof(ParticleSource));
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

            if (src_type == SRC_MONO_POINT) 
            {
                /* Put source name */
                snprintf(DATA.src->npoint.name, sizeof(DATA.src->npoint.name), "%s", a2);

                char *a2 = strtok(NULL, DELIMS);
                char *a3 = strtok(NULL, DELIMS);
                char *a4 = strtok(NULL, DELIMS);
                char *a5 = strtok(NULL, DELIMS);

                if (!a2 || !a3 || !a4 || !a5)
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                parseDouble(a2, &DATA.src->npoint.E);
                parseDouble(a3, &DATA.src->npoint.x);
                parseDouble(a4, &DATA.src->npoint.y);
                parseDouble(a5, &DATA.src->npoint.z);
            }
            /* Fissile material source */

            else if (src_type == SRC_FISS_MAT)
            {
                /* Put source name */
                snprintf(DATA.src->fission.name, sizeof(DATA.src->fission.name), "%s", a2);

                char *a3 = strtok(NULL, DELIMS);

                if (!a3)
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Put material name */

                snprintf(DATA.src->fission.mat_name, sizeof(DATA.src->fission.mat_name), "%s", a3);
                DATA.src->fission.mat_idx = -1;

                /* Try to read boundaries */

                char *a4 = strtok(NULL, DELIMS);
                char *a5 = strtok(NULL, DELIMS);
                char *a6 = strtok(NULL, DELIMS);
                char *a7 = strtok(NULL, DELIMS);
                char *a8 = strtok(NULL, DELIMS);
                char *a9 = strtok(NULL, DELIMS);

                if (!a4 || !a5 || !a6 || !a7 || !a8 || !a9)
                {
                    /* Default to geometry bounds */

                    DATA.src->fission.xmin = NAN;
                    DATA.src->fission.xmax = NAN;
                    DATA.src->fission.ymin = NAN;
                    DATA.src->fission.ymax = NAN;
                    DATA.src->fission.zmin = NAN;
                    DATA.src->fission.zmax = NAN;
                }
                else
                {
                    parseDouble(a4, &DATA.src->fission.xmin);
                    parseDouble(a5, &DATA.src->fission.xmax);
                    parseDouble(a6, &DATA.src->fission.ymin);
                    parseDouble(a7, &DATA.src->fission.ymax);
                    parseDouble(a8, &DATA.src->fission.zmin);
                    parseDouble(a9, &DATA.src->fission.zmax);
                }

            }
            else 
            {
                fprintf(stderr, "[ERROR] Source type %d not implemented (line %ld).\n", DATA.src_type, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            if (VERBOSITY >= 2)
                fprintf(stdout, "Parsed a source definition of type: %d.\n", DATA.src_type);

            np++;
        }

        /* ###################################################################################### */        
        /* --- det --- (Detector definition for scoring reaction channels)*/
        else if (!strcmp(tok, "det"))
        {
            Detector det;
            memset(&det, 0, sizeof(det));

            char *name_tok = strtok(NULL, DELIMS);
            if (!name_tok)
            {
                fprintf(stderr, "[ERROR] Missing detector name on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            snprintf(det.name, sizeof(det.name), "%s", name_tok);

            DetectorResponseBin *resp_bins = NULL;
            size_t resp_count = 0u;
            size_t resp_cap = 0u;

            bool seen_dx = false;
            bool seen_dy = false;
            bool seen_dz = false;
            bool seen_de = false;
            bool seen_dt = false;
            bool seen_dm = false;
            char mat_filter_name[MAX_STR_LEN] = {0};

            char *opt = NULL;
            while ((opt = strtok(NULL, DELIMS)) != NULL)
            {
                if (!strcmp(opt, "dr"))
                {
                    char *mt_tok = strtok(NULL, DELIMS);
                    char *target_tok = strtok(NULL, DELIMS);
                    if (!mt_tok || !target_tok)
                    {
                        fprintf(stderr, "[ERROR] Incomplete 'dr' specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    long mt_val = 0;
                    if (!parseLong(mt_tok, &mt_val))
                    {
                        fprintf(stderr, "[ERROR] Invalid MT value \"%s\" on line %ld.\n", mt_tok, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (resp_count == resp_cap)
                    {
                        size_t new_cap = resp_cap ? resp_cap * 2u : 4u;
                        DetectorResponseBin *tmp = (DetectorResponseBin*)realloc(resp_bins, new_cap * sizeof(DetectorResponseBin));
                        if (!tmp)
                        {
                            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                            free(resp_bins);
                            fclose(fp);
                            exit(EXIT_FAILURE);
                        }
                        resp_bins = tmp;
                        resp_cap = new_cap;
                    }

                    DetectorResponseBin *bin = &resp_bins[resp_count++];
                    memset(bin, 0, sizeof(*bin));
                    bin->response_id = mt_val;
                    bin->rmat_id = -1;
                    bin->multiplier = 1.0;
                    snprintf(bin->rmat_name, sizeof(bin->rmat_name), "%s", target_tok);
                }
                else if (!strcmp(opt, "dx") || !strcmp(opt, "dy") || !strcmp(opt, "dz"))
                {
                    bool *seen_flag = NULL;
                    DetectorCartesianAxis *axis = NULL;
                    const char *label = opt;

                    if (!strcmp(opt, "dx"))
                    {
                        seen_flag = &seen_dx;
                        axis = &det.mesh_x;
                    }
                    else if (!strcmp(opt, "dy"))
                    {
                        seen_flag = &seen_dy;
                        axis = &det.mesh_y;
                    }
                    else
                    {
                        seen_flag = &seen_dz;
                        axis = &det.mesh_z;
                    }

                    if (*seen_flag)
                    {
                        fprintf(stderr, "[ERROR] Duplicate '%s' specification on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    char *min_tok = strtok(NULL, DELIMS);
                    char *max_tok = strtok(NULL, DELIMS);
                    char *bins_tok = strtok(NULL, DELIMS);
                    if (!min_tok || !max_tok || !bins_tok)
                    {
                        fprintf(stderr, "[ERROR] Incomplete '%s' specification on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    double min_val = 0.0;
                    double max_val = 0.0;
                    if (!parseDouble(min_tok, &min_val) || !parseDouble(max_tok, &max_val))
                    {
                        fprintf(stderr, "[ERROR] Invalid '%s' bounds on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (max_val <= min_val)
                    {
                        fprintf(stderr, "[ERROR] '%s' max must exceed min on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    uint64_t bins_tmp = 0;
                    if (!parseULong(bins_tok, &bins_tmp) || bins_tmp == 0)
                    {
                        fprintf(stderr, "[ERROR] Invalid '%s' bin count on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (bins_tmp > (uint64_t)SIZE_MAX)
                    {
                        fprintf(stderr, "[ERROR] '%s' bin count too large on line %ld.\n", label, lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    axis->active = true;
                    axis->type = DETECTOR_MESH_AXIS_UNIFORM;
                    axis->n_bins = (size_t)bins_tmp;
                    axis->origin = min_val;
                    axis->min = min_val;
                    axis->max = max_val;
                    axis->spacing = (max_val - min_val) / (double)axis->n_bins;
                    axis->edges = NULL;
                    *seen_flag = true;
                }
                else if (!strcmp(opt, "de"))
                {
                    if (seen_de)
                    {
                        fprintf(stderr, "[ERROR] Duplicate 'de' specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    char *type_tok = strtok(NULL, DELIMS);
                    char *emin_tok = strtok(NULL, DELIMS);
                    char *emax_tok = strtok(NULL, DELIMS);
                    char *bins_tok = strtok(NULL, DELIMS);
                    if (!type_tok || !emin_tok || !emax_tok || !bins_tok)
                    {
                        fprintf(stderr, "[ERROR] Incomplete energy bin specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    long type_val = 0;
                    if (!parseLong(type_tok, &type_val) ||
                        (type_val != ENERGY_GRID_LIN && type_val != ENERGY_GRID_LOG))
                    {
                        fprintf(stderr, "[ERROR] Energy grid type must be 1 (lin) or 2 (log) on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    double emin = 0.0;
                    double emax = 0.0;
                    if (!parseDouble(emin_tok, &emin) || !parseDouble(emax_tok, &emax))
                    {
                        fprintf(stderr, "[ERROR] Invalid energy bounds on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (emax <= emin)
                    {
                        fprintf(stderr, "[ERROR] Energy max must exceed min on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (type_val == ENERGY_GRID_LOG && emin <= 0.0)
                    {
                        fprintf(stderr, "[ERROR] Logarithmic energy bins require Emin > 0 on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    uint64_t bins_tmp = 0;
                    if (!parseULong(bins_tok, &bins_tmp) || bins_tmp == 0)
                    {
                        fprintf(stderr, "[ERROR] Invalid number of energy bins on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (bins_tmp > (uint64_t)SIZE_MAX)
                    {
                        fprintf(stderr, "[ERROR] Energy bin count too large on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    det.energy.active = true;
                    det.energy.n_bins = (size_t)bins_tmp;
                    det.energy.type = (EnergyGridType)type_val;
                    det.energy.spacing = (type_val == ENERGY_GRID_LIN) ?
                                          DETECTOR_GRID_SPACING_LINEAR : DETECTOR_GRID_SPACING_LOG;
                    det.energy.min = emin;
                    det.energy.max = emax;
                    det.energy.edges = NULL;
                    seen_de = true;
                }
                else if (!strcmp(opt, "dt"))
                {
                    if (seen_dt)
                    {
                        fprintf(stderr, "[ERROR] Duplicate 'dt' specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    char *tmin_tok = strtok(NULL, DELIMS);
                    char *tmax_tok = strtok(NULL, DELIMS);
                    char *bins_tok = strtok(NULL, DELIMS);
                    if (!tmin_tok || !tmax_tok || !bins_tok)
                    {
                        fprintf(stderr, "[ERROR] Incomplete time bin specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    double tmin = 0.0;
                    double tmax = 0.0;
                    if (!parseDouble(tmin_tok, &tmin) || !parseDouble(tmax_tok, &tmax))
                    {
                        fprintf(stderr, "[ERROR] Invalid time bounds on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (tmax <= tmin)
                    {
                        fprintf(stderr, "[ERROR] Time max must exceed min on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    uint64_t bins_tmp = 0;
                    if (!parseULong(bins_tok, &bins_tmp) || bins_tmp == 0)
                    {
                        fprintf(stderr, "[ERROR] Invalid number of time bins on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (bins_tmp > (uint64_t)SIZE_MAX)
                    {
                        fprintf(stderr, "[ERROR] Time bin count too large on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    det.time.active = true;
                    det.time.n_bins = (size_t)bins_tmp;
                    det.time.spacing = DETECTOR_GRID_SPACING_LINEAR;
                    det.time.min = tmin;
                    det.time.max = tmax;
                    det.time.edges = NULL;
                    seen_dt = true;
                }
                else if (!strcmp(opt, "dm"))
                {
                    if (seen_dm)
                    {
                        fprintf(stderr, "[ERROR] Duplicate 'dm' specification on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    char *mat_tok = strtok(NULL, DELIMS);
                    if (!mat_tok)
                    {
                        fprintf(stderr, "[ERROR] Missing material after 'dm' on line %ld.\n", lnum);
                        free(resp_bins);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    snprintf(mat_filter_name, sizeof(mat_filter_name), "%s", mat_tok);
                    seen_dm = true;
                }
                else
                {
                    fprintf(stderr, "[ERROR] Unknown detector option '%s' on line %ld.\n", opt, lnum);
                    free(resp_bins);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }

            if (resp_count == 0)
            {
                resp_bins = (DetectorResponseBin*)calloc(1, sizeof(DetectorResponseBin));
                if (!resp_bins)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                resp_bins[0].response_id = 0;
                resp_bins[0].rmat_id = -1;
                resp_bins[0].multiplier = 1.0;
                resp_bins[0].rmat_name[0] = '\0';
                resp_count = 1;
            }

            det.responses.active = true;
            det.responses.n_bins = resp_count;
            det.responses.bins = resp_bins;

            if (seen_dm)
            {
                det.has_material_filter = true;
                det.material_filter_index = -1;
                snprintf(det.material_filter_name, sizeof(det.material_filter_name), "%s", mat_filter_name);
            }

            if (VERBOSITY >= 2)
            {
                size_t n_resp = det.responses.n_bins;
                size_t n_time = det.time.active ? det.time.n_bins : 1u;
                size_t n_energy = det.energy.active ? det.energy.n_bins : 1u;
                size_t nx = det.mesh_x.active ? det.mesh_x.n_bins : 1u;
                size_t ny = det.mesh_y.active ? det.mesh_y.n_bins : 1u;
                size_t nz = det.mesh_z.active ? det.mesh_z.n_bins : 1u;
                fprintf(stdout,
                        "Parsed detector '%s' with %zu response bin(s), time=%zu, energy=%zu, mesh=(%zu,%zu,%zu).\n",
                        det.name, n_resp, n_time, n_energy, nx, ny, nz);
            }

            Detector *new_arr = NULL;
            if (DATA.n_detectors == 0)
                new_arr = (Detector*)calloc(1u, sizeof(Detector));
            else
                new_arr = (Detector*)realloc(DATA.detectors, (DATA.n_detectors + 1u) * sizeof(Detector));

            if (!new_arr)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                free(resp_bins);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            DATA.detectors = new_arr;
            DATA.detectors[DATA.n_detectors] = det;
            DATA.n_detectors++;
            np++;
        }
        /* ###################################################################################### */        
        /* --- surf --- (Surface definition) */

        else if (!strcmp(tok, "surf"))
        {
            char *name = strtok(NULL, DELIMS);
            char *type = strtok(NULL, DELIMS);

            if (!name || !type) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Check for duplicate definition */
            for (size_t s = 0; s < DATA.n_surf; s++)
            {
                Surface *S = &DATA.surfs[s];
                if (!strcmp(name, S->name))
                {
                    fprintf(stderr, "[ERROR] Duplicate definition of surface \"%s\" on line %ld\n", name, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }

            /* Read all params into an array */
            
            double params[MAX_N_PARAMS];
            size_t nparams = 0;
            char* paramTok = strtok(NULL, DELIMS);

            while (nparams < MAX_N_PARAMS && paramTok != NULL)
            {
                if (!parseDouble(paramTok, &params[nparams++]))
                {
                    fprintf(stderr, "[ERROR] Invalid surface parameter \"%s\" on line %ld.\n", paramTok, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                paramTok = strtok(NULL, DELIMS);   
            }

            /* Create new surface from input */

            Surface S = {0};
            S.t_idx = -1;

            /* Copy all parameters */

            S.n_params = nparams;
            snprintf(S.name, sizeof(S.name), "%s", name);

            S.params = (double*)calloc(nparams, sizeof(double));
            if (!S.params)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            memcpy(S.params, params, nparams * sizeof(double));

            /* Put correct type */

            if (!strcmp(type, "px"))
                S.type = SURF_PLANEX;
            else if (!strcmp(type, "py"))
                S.type = SURF_PLANEY;
            else if (!strcmp(type, "pz"))
                S.type = SURF_PLANEZ;
            else if (!strcmp(type, "plane"))
                S.type = SURF_PLANE;
            else if (!strcmp(type, "sqr"))
                S.type = SURF_SQR;
            else if (!strcmp(type, "tri"))
                S.type = SURF_TRI;
            else if (!strcmp(type, "hexx"))
                S.type = SURF_HEXX;
            else if (!strcmp(type, "hexy"))
                S.type = SURF_HEXY;
            else if (!strcmp(type, "sph"))
                S.type = SURF_SPH;
            else if (!strcmp(type, "cylx"))
                S.type = SURF_CYLX;
            else if (!strcmp(type, "cyly"))
                S.type = SURF_CYLY;
            else if (!strcmp(type, "cylz") || !strcmp(type, "cyl"))
                S.type = SURF_CYLZ;
            else if (!strcmp(type, "cone"))
                S.type = SURF_CONE;
            else if (!strcmp(type, "cube"))
                S.type = SURF_CUBE;
            else if (!strcmp(type, "cuboid"))
                S.type = SURF_CUBOID;
            else if (!strcmp(type, "torus"))
                S.type = SURF_TORUS;
            else if (!strcmp(type, "inf"))
                S.type = SURF_INF;
            else
            {
                fprintf(stderr, "[ERROR] Unknown surface type \"%s\" on line %ld.\n", type, lnum);
                free(S.params);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Put surface into DATA */

            if (DATA.n_surf == 0)
            {
                DATA.surfs = (Surface*)calloc(1, sizeof(Surface));
            }
            else
            {
                DATA.surfs = (Surface*)realloc(DATA.surfs, (DATA.n_surf + 1) * sizeof(Surface));
            }
            if (!DATA.surfs)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            DATA.surfs[DATA.n_surf++] = S;
            np++;

            if (VERBOSITY >= 1)
                fprintf(stdout, "Parsed surface \"%s\" of type: %s (%d), with %zu params\n", name, type, S.type, nparams);
        }

        /* ###################################################################################### */
        /* --- cell --- (Cell definition) */

        else if (!strcmp(tok, "cell"))
        {
            /* Read required parameters */
            char *name = strtok(NULL, DELIMS);
            char *uni = strtok(NULL, DELIMS);
            char *fillTok = strtok(NULL, DELIMS);

            if (!name || !uni || !fillTok)
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Check for duplicate definition */
            for (size_t c = 0; c < DATA.n_cells; c++)
            {
                Cell *C = &DATA.cells[c];
                if (!strcmp(name, C->name))
                {
                    fprintf(stderr, "[ERROR] Duplicate definition of cell \"%s\" on line %ld\n", name, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }

            /* Check if universe filled or material */

            bool is_filled = false;
            char *fill_uni;
            if (!strcmp(fillTok, "fill"))
            {
                /* Filled with universe*/

                fill_uni = strtok(NULL, DELIMS);
                if (!fill_uni)
                {
                    fprintf(stderr, "[ERROR] Missing universe name after \"fill\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                is_filled = true;
            }
            else 
            {
                /* Filled with material */

                is_filled = false;
            }

            /* Read list of intersections into an array */
            size_t n_surfs = 0;
            char surf_names[MAX_N_PARAMS][MAX_STR_LEN];
            int sides[MAX_N_PARAMS];

            char *surfTok = strtok(NULL, DELIMS);

            while (n_surfs < MAX_N_PARAMS && surfTok != NULL)
            {
                /* Check if complement specifier given before surface name */

                if (!strncmp(surfTok, "-", 1))
                {
                    sides[n_surfs] = -1;
                    surfTok++;
                }
                else
                    sides[n_surfs] = 1;
                
                strncpy(surf_names[n_surfs++], surfTok, MAX_STR_LEN);
                surfTok = strtok(NULL, DELIMS);  
            }

            /* Create new cell from input */

            Cell C = {0};
            C.mat_idx = -1;
            C.filluni_idx = -1;
            C.uni_idx = -1;

            /* Copy all parameters */

            C.n_surfs = n_surfs;
            snprintf(C.name, sizeof(C.name), "%s", name);
            snprintf(C.uni_name, sizeof(C.uni_name), "%s", uni);

            /* Depending on filling type, put given name */

            if (is_filled)
            {
                snprintf(C.filluni_name, sizeof(C.filluni_name), "%s", fill_uni);
                C.unifilled = true;
            }
            else
            {
                snprintf(C.mat_name, sizeof(C.mat_name), "%s", fillTok);
                C.unifilled = false;
            }
            /* Allocate arrays */

            C.surf_names = (char*)calloc(n_surfs, MAX_STR_LEN);
            C.surf_idxs = (int*)calloc(n_surfs, sizeof(int));
            memset(C.surf_idxs, -1, n_surfs * sizeof(int));
            C.sides = (int*)calloc(n_surfs, sizeof(int));
            if (!C.surf_names || !C.surf_idxs || !C.sides)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            memcpy(C.surf_names, surf_names, n_surfs * MAX_STR_LEN);
            memcpy(C.sides, sides, n_surfs * sizeof(int));

            /* Put cell into DATA.cells */

            if (DATA.n_cells == 0)
            {
                DATA.cells = (Cell*)calloc(1, sizeof(Cell));
            }
            else
            {
                DATA.cells = (Cell*)realloc(DATA.cells, (DATA.n_cells + 1) * sizeof(Cell));
            }
            if (!DATA.cells)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            DATA.cells[DATA.n_cells++] = C;
            np++;

            if (VERBOSITY >= 1)
            {
                if (is_filled)
                    fprintf(stdout, "Parsed cell \"%s\", filled with universe \"%s\", in universe \"%s\" with %zu defining surface(s) given.\n", C.name, C.filluni_name, C.uni_name, C.n_surfs);
                else
                    fprintf(stdout, "Parsed cell \"%s\", of material \"%s\", in universe \"%s\" with %zu defining surface(s) given.\n", C.name, C.mat_name, C.uni_name, C.n_surfs);
            }
        }

        /* ###################################################################################### */        
        /* --- lat --- (Lattice definition) */
        else if (!strcmp(tok, "lat"))
        {
            /* Read required arguments */

            char* latUni = strtok(NULL, DELIMS);
            char* typeTok = strtok(NULL, DELIMS);

            if (!latUni || !typeTok)
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Create new lattice item */

            Lattice L = {0};
            L.uni_idx = -1;
            L.bb_surf_idx = -1;
            snprintf(L.name, sizeof(L.name), "%s", latUni);

            /* Parse type */
            long type;
            if (!parseLong(typeTok, &type))
            {
                fprintf(stderr, "[ERROR] Invalid lattice type \"%s\" on line %ld.\n", typeTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            if (type == 1)
                L.type = LAT_SQUARE_FINITE;
            else if (type == 2)
                L.type = LAT_HEXX_FINITE;
            else if (type == 3)
                L.type = LAT_HEXY_FINITE;
            else if (type == -1)
                L.type = LAT_SQUARE_INFINITE;
            else if (type == -2)
                L.type = LAT_HEXX_INFINITE;
            else if (type == -3)
                L.type = LAT_HEXY_INFINITE;

            /* Read parameters according to type */

            switch (L.type)
            {
                case LAT_SQUARE_INFINITE:
                case LAT_HEXX_INFINITE:
                case LAT_HEXY_INFINITE:
                {
                    char* x0 = strtok(NULL, DELIMS);
                    char* y0 = strtok(NULL, DELIMS);
                    char* p = strtok(NULL, DELIMS);
                    char* u = strtok(NULL, DELIMS);

                    if (!x0 || !y0 || !p || !u)
                    {
                        fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (!parseDouble(x0, &L.x0) || !parseDouble(y0, &L.y0) || !parseDouble(p, &L.pitch))
                    {
                        fprintf(stderr, "[ERROR] Invalid value for \"lat\" on line %ld.\n", lnum);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    L.nx = L.ny = L.nz = 1;
                    L.n_unis = 1;

                    L.uni_names = (char*)calloc(L.n_unis, MAX_STR_LEN);
                    if (!L.uni_names)
                    {
                        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    snprintf(L.uni_names, MAX_STR_LEN, "%s", u);

                    break;
                }   
                case LAT_SQUARE_FINITE:
                case LAT_HEXX_FINITE:
                case LAT_HEXY_FINITE:
                {
                    char *x0 = strtok(NULL, DELIMS);
                    char* y0 = strtok(NULL, DELIMS);
                    char* nx = strtok(NULL, DELIMS);
                    char* ny = strtok(NULL, DELIMS);
                    char* p = strtok(NULL, DELIMS);

                    if (!x0 || !y0 || !nx || !ny || !p)
                    {
                        fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (!parseDouble(x0, &L.x0) || !parseDouble(y0, &L.y0) ||
                        !parseLong(nx, &L.nx) || !parseLong(ny, &L.ny) || !parseDouble(p, &L.pitch))
                    {
                        fprintf(stderr, "[ERROR] Invalid value for \"lat\" on line %ld.\n", lnum);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    if (L.nx <= 0 || L.ny <= 0)
                    {
                        fprintf(stderr, "[ERROR] Invalid lattice dimensions \"%ld x %ld\" on line %ld.\n", L.nx, L.ny, lnum);
                        fclose(fp);
                        exit(EXIT_FAILURE);
                    }

                    L.nz = 1;
                    L.n_unis = (size_t)L.nx * (size_t)L.ny * (size_t)L.nz;

                    readLatticeUniverseNames(fp, &lnum, &L);

                    break;
                }
                default:
                {
                    fprintf(stderr, "[ERROR] Unknown lattice type \"%s\" on line %ld.\n", typeTok, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }
            

            /* Put cell into DATA.lats */

            if (DATA.n_lats == 0)
            {
                DATA.lats = (Lattice*)calloc(1, sizeof(Lattice));
            }
            else
            {
                DATA.lats = (Lattice*)realloc(DATA.lats, (DATA.n_lats + 1) * sizeof(Lattice));
            }
            if (!DATA.lats)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            DATA.lats[DATA.n_lats++] = L;
            np++;

            if (VERBOSITY >= 1)
                fprintf(stdout, "Parsed a lattice of type %d, filled with %zu universe(s)\n", (int)L.type, L.n_unis);
            
        }
        
        /* ###################################################################################### */
        /* --- trans --- (Geometry transformation) */

        else if (!strcmp(tok, "trans"))
        {
            /* Read type and target */

            char* typeTok = strtok(NULL, DELIMS);
            char* target = strtok(NULL, DELIMS);

            if (!typeTok || !target)
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Create new transformation and put type and target */

            Transform T = {0};
            snprintf(T.target_name, sizeof(T.target_name), "%s", target);
            
            if (!strcmp(typeTok, "s"))
                T.type = TRA_SURFACE;
            else
            {
                fprintf(stderr, "[ERROR] Unknown type for transformation \"%s\".\n", typeTok);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Read following params into an array */

            double params[MAX_N_PARAMS];
            size_t nparams = 0;
            char* paramTok = strtok(NULL, DELIMS);

            while (nparams < MAX_N_PARAMS && paramTok != NULL)
            {
                if (!parseDouble(paramTok, &params[nparams++]))
                {
                    fprintf(stderr, "[ERROR] Invalid transformation parameter \"%s\" on line %ld.\n", paramTok, lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                paramTok = strtok(NULL, DELIMS);   
            }

            const double eps = 1e-12;

            if (nparams == 3 || nparams == 6)
            {
                T.translation[0] = params[0];
                T.translation[1] = params[1];
                T.translation[2] = params[2];

                if ((fabs(params[0]) > eps) || (fabs(params[1]) > eps) || (fabs(params[2]) > eps))
                    T.action = TRA_TRANSLATION;
            }

            if (nparams == 6)
            {
                double ax = params[3];
                double ay = params[4];
                double az = params[5];

                double a = ax * (M_PI / 180.0);
                double b = ay * (M_PI / 180.0);
                double c = az * (M_PI / 180.0);

                double ca = cos(a), sa = sin(a);
                double cb = cos(b), sb = sin(b);
                double cc = cos(c), sc = sin(c);

                T.rotation[0][0] = cb * cc;
                T.rotation[0][1] = ca * sc + sa * sb * cc;
                T.rotation[0][2] = sa * sc - ca * sb * cc;
                T.rotation[1][0] = -cb * sc;
                T.rotation[1][1] = ca * cc - sa * sb * sc;
                T.rotation[1][2] = sa * cc + ca * sb * sc;
                T.rotation[2][0] = sb;
                T.rotation[2][1] = -sa * cb;
                T.rotation[2][2] = ca * cb;

                if ((fabs(ax) > eps) || (fabs(ay) > eps) || (fabs(az) > eps))
                {
                    if (T.action == TRA_TRANSLATION)
                        T.action = TRA_BOTH;
                    else
                        T.action = TRA_ROTATION;
                }
            }

            if (nparams != 3 && nparams != 6)
            {
                fprintf(stderr, "[ERROR] Invalid number of parameters for \"trans\" was expecting 3 or 6, got %zu\n", nparams);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            if (T.action == TRA_NONE)
            {
                fprintf(stderr, "[ERROR] Empty transformation on line %ld.", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Put cell into DATA.trans */

            if (DATA.n_transforms == 0)
            {
                DATA.transforms = (Transform*)calloc(1, sizeof(Transform));
            }
            else
            {
                DATA.transforms = (Transform*)realloc(DATA.transforms, (DATA.n_transforms + 1) * sizeof(Transform));
            }
            if (!DATA.transforms)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            DATA.transforms[DATA.n_transforms++] = T;
            np++;

            if (VERBOSITY >= 1)
                fprintf(stdout, "Parsed a transformation of type %d, with %s targeting \"%s\"\n", (int)T.type, nparams == 6 ? "rotation and translation" : "translation only", T.target_name);
        }

        /* ###################################################################################### */
        /* --- plot --- (Geometry plot definition) */

        else if (!strcmp(tok, "plot"))
        {
            /* Read plot axis and nominal image size in pixels (pixel count of long side)*/

            char *axTok = strtok(NULL, DELIMS);
            char *boundTok = strtok(NULL, DELIMS);
            char *pixTok = strtok(NULL, DELIMS);
            char *posTok = strtok(NULL, DELIMS);

            if (!axTok || !boundTok || !pixTok || !posTok)
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Create new GeometryPlotter item */

            GeometryPlotter gpl = {0};
            
            /* Parse axis */

            long ax;
            if (!strcmp(axTok, "xy") || !strcmp(axTok, "XY") || (parseLong(axTok, &ax) && ax == 3))
                gpl.axis = 3;
            else if (!strcmp(axTok, "xz") || !strcmp(axTok, "XZ") || (parseLong(axTok, &ax) && ax == 2))
                gpl.axis = 2;
            else if (!strcmp(axTok, "yz") || !strcmp(axTok, "YZ") || (parseLong(axTok, &ax) && ax == 1))
                gpl.axis = 1;
            else
            {
                fprintf(stderr, "[ERROR] Invalid axis specifier \"%s\" on line %ld.\n", axTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Parse boundary type */

            long bound;
            if (!strcmp(boundTok, "c") || !strcmp(boundTok, "C") || (parseLong(boundTok, &bound) && bound == 2))
                gpl.bounds = 2;
            else if (!strcmp(boundTok, "m") || !strcmp(boundTok, "M") || (parseLong(boundTok, &bound) && bound == 1))
                gpl.bounds = 1;
            else if ((parseLong(boundTok, &bound) && bound == 0))
                gpl.bounds = 0;
            else
            {
                fprintf(stderr, "[ERROR] Invalid boundary specifier \"%s\" on line %ld.\n", boundTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Parse and  out pixels to pixx, pixy will be defined later to match the aspect ratio */

            if (!parseLong(pixTok, &gpl.pixx) || gpl.pixx <= 0)
            {
                fprintf(stderr, "[ERROR] Invalid pixel size \"%s\" on line %ld.\n", axTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Parse position along normal axis */

            if (!parseDouble(posTok, &gpl.pos))
                        {
                fprintf(stderr, "[ERROR] Invalid plot plane position \"%s\" on line %ld.\n", posTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Try to read four more parameters for bounds */

            char *min1Tok = strtok(NULL, DELIMS);
            char *max1Tok = strtok(NULL, DELIMS);
            char *min2Tok = strtok(NULL, DELIMS);
            char *max2Tok = strtok(NULL, DELIMS);

            double min1, max1, min2, max2;

            if (!min1Tok || !max1Tok || !min2Tok ||!max2Tok)
                min1 = max1 = min2 = max2 = NAN;
            else
            {
                if (!parseDouble(min1Tok, &min1) || !parseDouble(max1Tok, &max1) ||
                    !parseDouble(min2Tok, &min2) || !parseDouble(max2Tok, &max2))
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                if ((min1 >= max1) || (min2 >= max2))
                {
                    fprintf(stderr, "[ERROR] Invalid boundaries for plot on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);    
                }
            }

            /* Put to gpl */

            if (gpl.axis == 1) // YZ-plane
            {
                gpl.ymin = min1;
                gpl.ymax = max1;
                gpl.zmin = min2;
                gpl.zmax = max2;
            }
            else if (gpl.axis == 2) // XZ-plane
            {
                gpl.xmin = min1;
                gpl.xmax = max1;
                gpl.zmin = min2;
                gpl.zmax = max2;
            }
            else // XY-plane
            {
                gpl.xmin = min1;
                gpl.xmax = max1;
                gpl.ymin = min2;
                gpl.ymax = max2;
            }

            /* Put gpl into DATA.gpls */

            if (DATA.n_gpls == 0)
            {
                DATA.gpls = (GeometryPlotter*)calloc(1, sizeof(GeometryPlotter));
            }
            else
            {
                DATA.gpls = (GeometryPlotter*)realloc(DATA.gpls, (DATA.n_gpls + 1) * sizeof(GeometryPlotter));
            }
            if (!DATA.gpls)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            DATA.gpls[DATA.n_gpls++] = gpl;
            np++;

            if (VERBOSITY >= 1)
                fprintf(stdout, "Parsed a geometry plot on axis \"%ld\", with limits (%lf %lf %lf %lf), boundaries \"%ld\".\n", gpl.axis, min1, max1, min2, max2, gpl.bounds);
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

            /* --- Boundary coefficent */
            else if (!strcmp(subkey, "bc"))
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                long bc;
                if (!parseLong(a1, &bc) || bc < BC_BLACK || bc > BC_PERIODIC) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set bc\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                DATA.boundary_coef = bc;
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

                double n_particles, n_generations, n_inactive;
                if (!parseDouble(a1, &n_particles) || !parseDouble(a2, &n_generations)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                if (!parseDouble(a3, &n_inactive))
                    n_inactive = 0;

                /* Check for valid values */

                if (n_generations < 1.0 || n_particles < 1.0 || n_inactive < 0.0) {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.n_generations = (long)n_generations;
                GLOB.n_particles = (long)n_particles;
                GLOB.n_inactive = (long)n_inactive;
                GLOB.mode = RUNMODE_CRITICALITY;
                np++;
            }
            /* --- External source simulation population */
            else if (!strcmp(subkey, "nps")) 
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

                double n_particles, n_cycles, n_inactive;
                if (!parseDouble(a1, &n_particles) || !parseDouble(a2, &n_cycles)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                if (!parseDouble(a3, &n_inactive))
                    n_inactive = 0;

                /* Check for valid values */

                if (n_cycles < 1.0 || n_particles < 1.0 || n_inactive < 0.0) {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.n_cycles = (long)n_cycles;
                GLOB.n_particles = (long)(n_particles / n_cycles);
                GLOB.n_inactive = (long)n_inactive;
                GLOB.mode = RUNMODE_EXTERNAL_SOURCE;
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
            /* --- Neutron buffer factor */
            else if (!strcmp(subkey, "nbuf")) 
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double nbuf_factor;
                if (!parseDouble(a1, &nbuf_factor)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.nbuf_factor = nbuf_factor;
                np++;
            }
            /* --- Normalisation source rate */
            else if (!strcmp(subkey, "srcrate")) 
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double srcrate;
                if (!parseDouble(a1, &srcrate) || srcrate < 0.0) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.srcrate = srcrate;
                GLOB.norm_mode = NORM_SRCRATE;
                np++;
            }
            /* --- Normalisation power */
            else if (!strcmp(subkey, "power")) 
            {
                char *a1 = strtok(NULL, DELIMS);
                if (!a1) 
                {
                    fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                double power;
                if (!parseDouble(a1, &power) || power < 0.0) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                GLOB.power = power;
                GLOB.norm_mode = NORM_POWER;
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
static void readLatticeUniverseNames(FILE *fp, long *lnum, Lattice *lat)
{
    size_t expected = lat->n_unis;

    if (expected == 0)
        return;

    lat->uni_names = (char*)calloc(expected, MAX_STR_LEN);
    if (!lat->uni_names)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed while reading lattice '%s'.\n",
                lat->name);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    size_t count = 0;
    char *tok = strtok(NULL, DELIMS);
    char line_buf[4096];

    while (1)
    {
        while (tok)
        {
            snprintf(lat->uni_names + count * MAX_STR_LEN, MAX_STR_LEN, "%s", tok);
            if (++count == expected)
                break;
            tok = strtok(NULL, DELIMS);
        }

        if (count == expected)
            break;

        if (!fgets(line_buf, sizeof line_buf, fp))
        {
            fprintf(stderr, "[ERROR] Lattice '%s' expected %zu universe names but reached EOF.\n",
                    lat->name, expected);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        (*lnum)++;

        char *cmnt = strchr(line_buf, '#');
        if (cmnt)
            *cmnt = '\0';

        char *s = line_buf;
        while (isspace((unsigned char)*s))
            s++;

        if (*s == '\0')
        {
            fprintf(stderr, "[ERROR] Empty line encountered while reading lattice '%s' (line %ld).\n",
                    lat->name, *lnum);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        tok = strtok(s, DELIMS);
    }
}

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
