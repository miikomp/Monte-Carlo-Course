#include "header.h"

static int parseLong(const char *s, long *out);

static int parseULong(const char *s, uint64_t *out);

static int parseUInt(const char *s, uint32_t *out);

static int parseDouble(const char *s, double *out);

/* small growable array for components while parsing */
typedef struct {
    MaterialNuclide *v;
    size_t n, cap;
} CompVec;

long readInput() {
    long np = 0l, lnum = 0l;

    /* Avoid compiler warning */
    parseUInt("\0", NULL);

    /* Open file for reading */
    
    FILE *fp = fopen(GLOB.fname, "r");
    if (fp == NULL) 
    {
        fprintf(stderr, "[ERROR] Cannot open file \"%s\".", GLOB.fname);
        exit(EXIT_FAILURE);
    }

    /* ########################################################################################## */

    /* Try to parse keyword arguments and related options */

    fprintf(stdout, "\nReading input file \"%s\"...\n\n", GLOB.fname);

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
        /* --- seed --- */
    
        if (!strcmp(tok, "seed")) 
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

            /* Put seed */

            GLOB.seed = seed;
            np++;
        }

        /* ###################################################################################### */
        /* --- pop --- */

        else if (!strcmp(tok, "pop")) 
        {
            char *a1 = strtok(NULL, DELIMS);
            char *a2 = strtok(NULL, DELIMS);
            if (!a1 || !a2) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            long n_particles, n_generations;
            if (!parseLong(a1, &n_particles) ||!parseLong(a2, &n_generations)) 
            {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Check for valid values */

            if (n_generations < 1 || n_particles < 1) {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Put outer and inner iterations */

            GLOB.n_generations = n_generations;
            GLOB.n_particles = n_particles;
            np++;
        }

        /* ###################################################################################### */
        /* --- xslibpath --- */
        
        else if (!strcmp(tok, "xslibpath"))
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

            /* Put to GLOB */

            snprintf(GLOB.xslibpath, MAX_STR_LEN, "%s", path);
            np++;
        }

        /* ###################################################################################### */
        /* --- material --- */
        else if (!strcmp(tok, "mat"))
        {
            /* Try to get the material specifications from the same line */

            const char *nameTok = strtok(NULL, DELIMS);
            const char *densTok = strtok(NULL, DELIMS);
            const char *tempTok = strtok(NULL, DELIMS);

            /* Parse material block */

            Material M;
            memset(&M, 0, sizeof M);

            /* Check for valid header */

            snprintf(M.name, sizeof M.name, "%s", nameTok);
            double density;
            if (!parseDouble(densTok, &density) || density == 0.0)
            {
                fprintf(stderr, "[ERROR] Bad density for material \"%s\" (line %ld).\n", M.name, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

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

             /* Parse temperature */

            if (!parseDouble(tempTok, &M.temp)) 
            {
                fprintf(stderr, "[ERROR] Bad temperature for material \"%s\" (line %ld).\n", M.name, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Read ZA + fraction pairs until empty line or comment */

            CompVec cv = {0};

            for (;;) {
                
                /* If EOF reached break loop */

                if (!fgets(line, sizeof line, fp)) 
                    break;
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
                mnuc.ZA = (int)ZA;
                mnuc.N_i = 0.0;

                /* val > 0 -> atomic count, val < 0 -> mass fraction */

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
                cv.v[cv.n++] = mnuc;
            }

            if (cv.n == 0) 
            {
                fprintf(stderr, "[ERROR] Material \"%s\" has no components (line %ld).\n", M.name, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Normalize fractions the sum that ends up above 0 is the type in which values were given */

            double sum_atoms = 0.0, sum_w = 0.0;
            for (size_t i = 0; i < cv.n; ++i) 
            { 
                sum_atoms += cv.v[i].atom_frac; 
                sum_w += cv.v[i].mass_frac; 
            }
            if (sum_atoms > 0.0) 
            {
                for (size_t i = 0; i < cv.n; ++i) 
                    cv.v[i].atom_frac /= sum_atoms;
            } 
            else 
            {
                for (size_t i = 0; i < cv.n; ++i) 
                    cv.v[i].mass_frac /= sum_w;    
            }

            /* Move into global DATA */

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

            /* Print summary for successfully parsed material*/

            fprintf(stdout, "Parsed material \"%s\" at %.6E %s at %.1fK with %zu nuclide(s).\n", 
                    M.name, 
                    (M.adens > 0.0) ? M.adens : M.mdens,
                    (M.adens > 0.0) ? "atoms/b*cm" : "g/cm3",
                    M.temp, 
                    M.n_nucs
                );
            np++;
        }

        /* ###################################################################################### */        
        /* --- Set --- */

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

            /* Read the first value */

            char *value = strtok(NULL, DELIMS);
            if (!value) 
            {
                fprintf(stderr, "[ERROR] Missing value for \"set %s\" on line %ld.\n", subkey, lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* ################################################################################## */
            /* --- Handle the specifier keywords --- */

            if (!strcmp(subkey, "mode")) 
            {
                long mode;
                if (!parseLong(value, &mode)) 
                {
                    fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Check mode */
                
                if (!((mode == RUNMODE_CIRCLE_PI) || 
                    (mode == RUNMODE_BUFFONS_PI) || 
                    (mode == RUNMODE_CHECK))) 
                {
                    fprintf(stderr, "[ERROR] Unknown mode identifier \"%ld\" for \"mode\".\n", mode);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
                
                /* Put mode */
                
                GLOB.mode = mode;
                np++;
            }
            else if (!strcmp(subkey, "needle")) 
            {
                double needle_length;
                if (!parseDouble(value, &needle_length)) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set needle\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Put needle length */

                GLOB.needle_length = needle_length;
                np++;
            }
            else if (!strcmp(subkey, "lines")) 
            {
                double line_spacing;
                if (!parseDouble(value, &line_spacing)) 
                {
                    fprintf(stderr, "[ERROR] Invalid value for \"set lines\" on line %ld.\n", lnum);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }

                /* Put line spacing */

                GLOB.line_spacing = line_spacing;
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
    errno = 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if (errno || end == s || *end != '\0') 
        return 0;
    *out = v;
    return 1;
}
