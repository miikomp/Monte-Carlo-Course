#include "header.h"

static int parseLong(const char *s, long *out);

static int parseULong(const char *s, uint64_t *out);

static int parseUInt(const char *s, uint32_t *out);

static int parseDouble(const char *s, double *out);

long readInput() {
    long np = 0u, lnum = 0u;

    /* Open file for reading */
    
    FILE *fp = fopen(GLOB.fname, "r");
    if (fp == NULL) 
    {
        fprintf(stderr, "[ERROR] Cannot open file \"%s\".", GLOB.fname);
        exit(EXIT_FAILURE);
    }

    /* Try to parse keyword arguments and related options */

    fprintf(stdout, "\nReading input file \"%s\"...\n", GLOB.fname);

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
        
        /* --- Parse line --- */

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

            long n_outer, n_inner;
            if (!parseLong(a1, &n_outer) || !parseLong(a2, &n_inner)) 
            {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Check for valid values */

            if (n_outer <= 1 || n_inner <= 1) {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Put outer and inner iterations */

            GLOB.n_outer = n_outer;
            GLOB.n_inner = n_inner;
            np++;
        }
        else if (!strcmp(tok, "mode")) 
        {
            char *a1 = strtok(NULL, DELIMS);
            if (!a1) 
            {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            long mode;
            if (!parseLong(a1, &mode)) 
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
        else
        {
            /* Some unknown keyword, print warning and continue to next line */

            fprintf(stderr, "[WARNING] Skipping line %ld with unknown keyword argument: \"%s\"...\n", lnum, tok);
        }
    }

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