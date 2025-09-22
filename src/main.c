#include "header.h"

int main(int argc, char **argv) {
    int val;

    /* ########################################################################################## */

    /* Check for valid usage of executable */

    if (argc < 2) 
    {
        fprintf(stderr, "[ERROR] Usage: %s [OPTS] filename.\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Filename is the last argument */

    GLOB.fname = argv[argc - 1];

    /* Try to parse commandline arguments */
    
    for (long i = 1; i < argc - 1; i++) 
    {
        if (!strcmp(argv[i], "-omp")) 
        {
            if (i + 1 >= argc - 1) 
            {
                fprintf(stderr, "[ERROR] Number of threads not given for \"%s\".\n", argv[i]);
                return EXIT_FAILURE;
            }

            /* Put number of threads */

            GLOB.n_threads = (int)fmin((strtol(argv[++i], NULL, 10)), 16);
        }
        else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) 
        {
            /* Verbosity level defaults to 1 if -v is given and 0 otherwise, 2 is maximum output */

            if (i + 1 >= argc - 1) 
                val = 1;
            else
                val = (int)strtol(argv[++i], NULL, 10);

            /* Set verbosity */

            VERBOSITY = val;
        }
        else 
        {
            fprintf(stderr, "[ERROR] Unknown commandline argument \"%s\".\n", argv[i]);
            return EXIT_FAILURE;
        }
        
    }

    /* ########################################################################################## */

    /* Read input file */
    
    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Processing input\n");
    fprintf(stdout, "------------------------\n");

    long np = readInput();
    GLOB.n_kwargs = np;
    fprintf(stdout, "DONE.\n");
    if (VERBOSITY >= 1)
        fprintf(stdout, "Succesfully parsed %ld keyword arguments.\n", np);
    
    /* Process input data */

    if (processInput() != 0) 
        return EXIT_FAILURE;

    /* ########################################################################################## */

    /* If in transport mode read and process the xsdata and resolve all materials */

    if (GLOB.mode == RUNMODE_TRANSPORT || GLOB.mode == RUNMODE_CHECK)
    {
        fprintf(stdout, "\n------------------------\n");
        fprintf(stdout, "  Processing XS-data\n");
        fprintf(stdout, "------------------------\n");

        /* Read and process the xsdata from the given path and resolve all materials*/

        if (processXsData() != 0) 
        {
            fprintf(stderr, "[ERROR] Could not process cross section library file at \"%s\".\n", GLOB.xslibpath);
            return EXIT_FAILURE;
        }

        /* Compute macroscopic cross sections for all resolved materials and reaction modes */
        
        if (computeMacroXs() != 0)
        {
            fprintf(stderr, "[ERROR] Failed to compute macroscopic cross sections.\n");
            return EXIT_FAILURE;
        }
    }
    /* ########################################################################################## */

    /* Set desired number of threads and get actual provided */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();
    GLOB.n_threads = nt;

    /* Disable dynamic teaming to not mess up thread-private things */

    omp_set_dynamic(0);

    /* Initialize look-up tables */

    initTrigTables();

    /* ########################################################################################## */

    /* Start timer */

    GLOB.t0 = omp_get_wtime();

    /* --- Main loop --- */
    
    /* Dispatch case to correct sub-routine */

    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Starting simulation\n");
    fprintf(stdout, "------------------------\n\n");

    switch (GLOB.mode) {
    case RUNMODE_TRANSPORT:
    {
        /* Fill neutron bank with initial source neutrons */

        if (sampleInitialSource() < 0) {
            fprintf(stderr, "[ERROR] Failed to sample initial source.\n");
            return EXIT_FAILURE;
        }

        /* Run the transport simulation */

        if (runTransport() != 0) {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            return EXIT_FAILURE;
        }

        break;
    }
    case RUNMODE_CIRCLE_PI: 
    {
        /* Get scrambler seed from global seed */

        uint64_t sm = (GLOB.seed << 32) ^ UINT64_C(0x94D049BB133111EB);

        if (runCirclePi(sm) != 0)
        {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            return EXIT_FAILURE;
        }

        break;
    }
    case RUNMODE_BUFFONS_PI: 
    {

        /* Get scrambler seed from global seed */

        uint64_t sm = (GLOB.seed << 32) ^ UINT64_C(0x94D049BB133111EB);

        if (runBuffonsPi(sm) != 0) {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            return EXIT_FAILURE;
        }

        break;
    }
    default: 
    {
        fprintf(stderr, "[ERROR] Mode %ld not implemented.\n", GLOB.mode);
        return EXIT_FAILURE;
    }
    }

    /* ########################################################################################## */

    /* Stop timer */

    GLOB.t1 = omp_get_wtime();

    /* Print runtime information */

    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Runtime: %.4lfs\n", GLOB.t1 - GLOB.t0);
    fprintf(stdout, "------------------------\n\n");

    /* Prepare for termination */

    /* Exit */

    return EXIT_SUCCESS;
}