#include "header.h"

#define MAX_N_PARAMS 256

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
        /* --- surf --- (Surface definition) */

        else if (!strcmp(tok, "surf"))
        {
            /* "surf" [NAME] [TYPE] [ARGS] (all in one line) */

            char *name = strtok(NULL, DELIMS);
            char *type = strtok(NULL, DELIMS);

            if (!name || !type) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
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

            Surface S;

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
            else if (!strcmp(type, "cube"))
                S.type = SURF_CUBE;
            else if (!strcmp(type, "cuboid"))
                S.type = SURF_CUBOID;
            else if (!strcmp(type, "torusx"))
                S.type = SURF_TORUSX;
            else if (!strcmp(type, "torusy"))
                S.type = SURF_TORUSY;
            else if (!strcmp(type, "torusz"))
                S.type = SURF_TORUSZ;
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
            snprintf(L.uni_name, sizeof(L.uni_name), "%s", latUni);

            /* Parse type */

            long type = 0;
            if (!parseLong(typeTok, &type))
            {
                fprintf(stderr, "[ERROR] Invalid lattice type \"%s\" on line %ld.\n", typeTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            if (type == 1)
                L.type = LAT_SQUARE_INFINITE;
            else if (type == 2)
                L.type = LAT_SQUARE_FINITE;
            else
            {
                fprintf(stderr, "[ERROR] Unknown lattice type \"%s\" on line %ld.\n", typeTok, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Read parameters according to type */

            if (L.type == LAT_SQUARE_INFINITE)
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

                if (!parseDouble(x0, &L.x0) || !parseDouble(y0, &L.y0) || !parseDouble(p, &L.dx))
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"lat\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Pitch in y is same as x */

                L.dy = L.dx;
                L.n_unis = 1;

                /* Filling universe */
                L.uni_names = (char*)calloc(1, MAX_STR_LEN);
                if (!L.uni_names)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                snprintf(L.uni_names, MAX_STR_LEN, "%s", u);
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
