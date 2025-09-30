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
        else if (!strcmp(argv[i], "-check")) 
        {
            GLOB.mode = RUNMODE_CHECK;
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
        fprintf(stdout, "\nSuccesfully parsed %ld keyword arguments.\n", np);
    
    /* Process input data */
    fprintf(stdout, "\nProcessing input data...\n");

    if (processInput() != 0) 
        return EXIT_FAILURE;

    fprintf(stdout, "DONE.\n");

    /* ########################################################################################## */

    /* If in transport mode read and process the xsdata and resolve all materials */

    if (GLOB.mode == RUNMODE_CRITICALITY || GLOB.mode == RUNMODE_CHECK)
    {
        fprintf(stdout, "\n------------------------\n");
        fprintf(stdout, "  Processing XS-data\n");
        fprintf(stdout, "------------------------\n");

        /* Read and process the xsdata from the given path into a temporary librabry */
        TempNucDataLib *lib = NULL;
        size_t nlib = 0;

        if (processXsData(&lib, &nlib) != 0) 
        {
            fprintf(stderr, "[ERROR] Could not process cross section library file at \"%s\".\n", GLOB.xslibpath);
            return EXIT_FAILURE;
        }

        /* Resolve all materials using the temporary nuclide library */

        if (resolveMaterials(lib, nlib) != 0) 
        {
            fprintf(stderr, "[ERROR] Could not resolve materials.\n");
            return EXIT_FAILURE;
        }

        /* Library freed inside resolveMaterials avoid dangling pointer */

        lib = NULL;
        nlib = 0;

        /* Compute macroscopic cross sections for all resolved materials and reaction modes */
        
        if (computeMacroXs() != 0)
        {
            fprintf(stderr, "[ERROR] Failed to compute macroscopic cross sections.\n");
            return EXIT_FAILURE;
        }

        fprintf(stdout, "\n----------------------------\n");
        fprintf(stdout, "  Preparing simulation\n");
        fprintf(stdout, "----------------------------\n\n");

        /* Initialize results data struct */

        fprintf(stdout, "Clearing results...\n");

        RES.n_generations = GLOB.n_generations;
        RES.avg_scores = (GenerationScores*)calloc((size_t)GLOB.n_generations, sizeof(GenerationScores));
        if (!RES.avg_scores)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            return EXIT_FAILURE;
        }

        /* Calculate memory footprint of results data structure */

        size_t total_bytes = sizeof(ResultsData) + (size_t)GLOB.n_generations * sizeof(GenerationScores);
        fprintf(stdout, "Memory allocated for results: %.2f MB\n", (double)total_bytes / (1024.0 * 1024.0));


        fprintf(stdout, "DONE.\n");

        /* Process detectors */

        if (processDetectors() != 0) 
        {
            fprintf(stderr, "[ERROR] Could not process detectors.\n");
            return EXIT_FAILURE;
        }

        /* Fill neutron bank with initial source neutrons */

        if (sampleInitialSource() < 0) {
            fprintf(stderr, "[ERROR] Failed to sample initial source.\n");
            return EXIT_FAILURE;
        }
    }
    /* ########################################################################################## */

    /* Set desired number of threads and put actual number provided */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();
    GLOB.n_threads = nt;

    /* Disable dynamic teaming to not mess up thread-private things */

    omp_set_dynamic(0);

    /* Initialize look-up tables */

    initTrigTables();

    if (GLOB.mode == RUNMODE_CHECK)
    {
        fprintf(stdout, "\nLaunched in check mode and found no issues. Exiting...\n");
        return EXIT_SUCCESS;
    }

    /* ########################################################################################## */

    /* Start timer */

    GLOB.t0 = omp_get_wtime();

    /* --- Main loop --- */
    
    /* Dispatch case to correct sub-routine */

    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Starting simulation\n");
    fprintf(stdout, "------------------------\n\n");

    switch (GLOB.mode) {
    case RUNMODE_EXTERNAL_SOURCE:
    {
        /* 
           Run the external source simulation:
           Simulates a set number of cycles for a set number of neutrons. 
           Each cycle is independent of the last and last for the entire history of each neutron.
           For super critical systems a cut-off must be in place to avoid infinite loops.
        */

        fprintf(stdout, "Running external source simulation for %ld generations with %ld neutrons each.\n", 
            GLOB.n_generations, GLOB.n_particles);

        if (runExternalSourceSimulation() != 0) {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            return EXIT_FAILURE;
        }

        break;
    }
    case RUNMODE_CRITICALITY:
    {
        /* 
           Run the criticality simulation: 
           Simulates a set number of generations so that each new generation is spawned from
           the fission sites of the last generation. 
        */

        fprintf(stdout, "Running criticality source simulation for %ld generations with %ld neutrons each.\n", 
            GLOB.n_generations, GLOB.n_particles);

        if (runCriticalitySimulation() != 0) {
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

    /* ########################################################################################## */

    /* Process and output transport results */

    fprintf(stdout, "Processing results...\n");

    processTransportResults();
    processDetectorResults();

    fprintf(stdout, "\nDONE.\n");
    

    /* Prepare for termination */

    free(RES.avg_scores);
    RES.avg_scores = NULL;

    /* Exit */

    return EXIT_SUCCESS;
}
