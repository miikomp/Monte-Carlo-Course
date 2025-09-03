#include "header.h"

int main(int argc, char **argv) {
    
    /* Check for valid usage */

    if (argc < 2) {
        fprintf(stderr, "[ERROR] Usage: %s [OPTS] filename\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Filename is the last argument */

    GLOB.fname = argv[argc - 1];

    /* Try to parse commandline arguments */
    
    for (long n = 1; n < argc - 1; n++) {
        if (!strcmp(argv[n], "-omp")) {
            if (n + 1 >= argc - 1) {
                fprintf(stderr, "[ERROR] Number of threads not given for \"%s\"\n", argv[n]);
                return EXIT_FAILURE;
            }

            /* Put number of threads */

            GLOB.n_threads = labs(strtol(argv[++n], NULL, 10));
        }
    }

    /* Read input file */
    
    readInput();

    /* Exit */

    return EXIT_SUCCESS;
}