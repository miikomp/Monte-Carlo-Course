#include "header.h"

static int parseLong(const char *s, long *out);

static int parseULong(const char *s, uint64_t *out);

static int parseDouble(const char *s, double *out);

long readInput() {
    long np = 0;

    /* Open file for reading */
    
    FILE *fp = fopen(GLOB.fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "[ERROR] Cannot open file %s", GLOB.fname);
        exit(EXIT_FAILURE);
    }

    /* Try to parse keyword arguments and related options */

    fprintf(stdout, "Reading input file \"%s\"...\n", GLOB.fname);

    char line[4096];
    long lnum = 0;
    const char DELIMS[] = " \t\r\n";

    /* Read lines until end-of-file */

    while (fgets(line, sizeof(line), fp)) {
        ++lnum;

        /* Strip comments */
        
        char *cmnt = strchr(line, '#');
        if (cmnt)
            *cmnt = '\0';

        /* Tokenize with whitespaces */

        char *tok = strtok(line, DELIMS);
        if (!tok)
            continue;
        
        /* --- Parse line --- */

        if (!strcmp(tok, "seed")) {
            char *a1 = strtok(NULL, DELIMS);
            if (!a1) {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            uint64_t seed;
            if (!parseULong(a1, &seed)) {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Put seed */

            GLOB.seed = seed;
            np++;
        }
        else if (!strcmp(tok, "pop")) {
            char *a1 = strtok(NULL, DELIMS);
            char *a2 = strtok(NULL, DELIMS);
            if (!a1 || !a2) {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            uint64_t n_outer, n_inner;
            if (!parseULong(a1, &n_outer) || !parseULong(a2, &n_inner)) {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Put outer and inner iterations */

            GLOB.n_outer = n_outer;
            GLOB.n_inner = n_inner;
            np++;
        }
        else if (!strcmp(tok, "mode")) {
            char *a1 = strtok(NULL, DELIMS);
            if (!a1) {
                fprintf(stderr, "[ERROR] Incomplete input on line %ld.", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            uint64_t mode;
            if (!parseULong(a1, &mode)) {
                fprintf(stderr, "[ERROR] Invalid input on line %ld.\n", lnum);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            /* Check mode */
            
            if (!(mode == RUNTYPE_CIRCLE_PI) || !(mode == RUNTYPE_CIRCLE_PI)) {
                fprintf(stderr, "[ERROR] Unknown mode identifier \"%ld\" for \"mode\"\n", mode);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
            
            /* Put mode */
            
            GLOB.mode = mode;
            np++;
        }
    }

    fprintf(stdout, "DONE. %ld keyword arguments succesfully parsed\n", np);

    return np;
}


static int parseLong(const char *s, long *out) {
    errno = 0;
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (errno || end == s || *end != '\0') 
        return 0;
    *out = v;
    return 1;
}

static int parseULong(const char *s, uint64_t *out) {
    errno = 0;
    char *end = NULL;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || end == s || *end != '\0') return 0;
    *out = (uint64_t)v;
    return 1;
}

static int parseDouble(const char *s, double *out) {
    errno = 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if (errno || end == s || *end != '\0') 
        return 0;
    *out = v;
    return 1;
}