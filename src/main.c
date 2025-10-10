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

    char *filename = argv[argc - 1];
    static char inputfname[MAX_PATH];

    GLOB.inputf = filename;

    /* Strip file extension and path if present */
    const char *dot = strrchr(filename, '.');
    const char *slash = strrchr(filename, '/');

    if (dot && dot != filename && (!slash || dot > slash))
    {
        size_t len = (size_t)(dot - filename);
        if (len >= sizeof(inputfname))
            len = sizeof(inputfname) - 1;
        memcpy(inputfname, filename, len);
        inputfname[len] = '\0';
        GLOB.inputfname = inputfname;
    }
    else
    {
        GLOB.inputfname = filename;
    }

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
        else if (!strcmp(argv[i], "-norun") || !strcmp(argv[i], "--norun")) 
        {
            GLOB.norun = true;
        }
        else if (!strcmp(argv[i], "-checkvolumes") || !strcmp(argv[i], "--checkvolumes")) 
        {
            if (i + 2 >= argc - 1) 
            {
                fprintf(stderr, "[ERROR] Not enough arguments given for \"%s\".\n", argv[i]);
                return EXIT_FAILURE;
            }

            long type = strtol(argv[++i], NULL, 10);
            if (type < 1 || type > 2) 
            {
                fprintf(stderr, "[ERROR] Invalid volume check type %ld given for \"%s\".\n", type, argv[i - 1]);
                return EXIT_FAILURE;
            }
            if (type == 1)
                GLOB.n_points = strtol(argv[++i], NULL, 10);
            else
                GLOB.n_lines = strtol(argv[++i], NULL, 10);

            if (GLOB.n_points <= 0 && type == 1) 
            {
                fprintf(stderr, "[ERROR] Number of random points must be positive.\n");
                return EXIT_FAILURE;
            }
            if (GLOB.n_lines <= 0 && type == 2) 
            {
                fprintf(stderr, "[ERROR] Number of random lines must be positive.\n");
                return EXIT_FAILURE;
            }
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

    /* Set desired number of threads and put actual number provided */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();
    GLOB.n_threads = nt;

    /* Disable dynamic teaming to not mess up thread-private things */

    omp_set_dynamic(0);

    /* ########################################################################################## */

    /* Read and process the xsdata and resolve all materials */

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

    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Processing geometry\n");
    fprintf(stdout, "------------------------\n");

    /* Process geometry universes */

    if (resolveUniverses() != 0)
    {
        fprintf(stderr, "[ERROR] Failed to process universes.\n");
        return EXIT_FAILURE;
    }

    /* Process lattices */

    if (resolveLattices() != 0)
    {
        fprintf(stderr, "[ERROR] Failed to process lattices.\n");
        return EXIT_FAILURE;
    }

    /* Process geometry cells */

    if (resolveCells() != 0)
    {
        fprintf(stderr, "[ERROR] Failed to process geometry cells.\n");
        return EXIT_FAILURE;
    }

    /* Process transformations */

    if (resolveTransformations() != 0)
    {
        fprintf(stderr, "[ERROR] Failed to process transformations.\n");
        return EXIT_FAILURE;
    }

    /* Calculate and put outer bounds */

    if (resolveOuterBounds() != 0)
    {
        fprintf(stderr, "[ERROR] Failed to calculate outer bounds.\n");
        return EXIT_FAILURE;
    }

    /* Plot geometry */
    
    if (plotGeometry() != 0) 
    {
        fprintf(stderr, "[ERROR] Could not plot geometry.\n");
        return EXIT_FAILURE;
    }

    /* Check volumes by sampling random points */

    if (checkVolumes() != 0)
    {
        fprintf(stderr, "[ERROR] Volume checking failed.\n");
        return EXIT_FAILURE;
    }

    /* Check volumes by sampling random lines */

    if (checkVolumes2() != 0)
    {
        fprintf(stderr, "[ERROR] Volume checking failed.\n");
        return EXIT_FAILURE;
    }

    /* ########################################################################################## */

    fprintf(stdout, "\n----------------------------\n");
    fprintf(stdout, "  Preparing simulation\n");
    fprintf(stdout, "----------------------------\n\n");

    /* Initialize results data struct */

    fprintf(stdout, "Clearing results...\n");

    size_t n_res = (size_t)GLOB.n_generations ? GLOB.n_generations > 0 : GLOB.n_cycles;
    RES.n_iterations = n_res;
    RES.avg_scores = (TransportRunScores*)calloc(n_res, sizeof(TransportRunScores));
    if (!RES.avg_scores)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    /* Calculate memory footprint of results data structure */

    size_t total_bytes = sizeof(ResultsData) + n_res * sizeof(TransportRunScores);
    fprintf(stdout, "Memory allocated for results: %.2f MB\n", (double)total_bytes / (1024.0 * 1024.0));


    fprintf(stdout, "DONE.\n");

    /* Process detectors */

    if (processDetectors() != 0) 
    {
        fprintf(stderr, "[ERROR] Could not process detectors.\n");
        return EXIT_FAILURE;
    }

    /* Fill neutron bank with initial source neutrons */

    fprintf(stdout, "\nSampling initial neutron source...\n");

    if (sampleInitialSource() < 0) {
        fprintf(stderr, "[ERROR] Failed to sample initial source.\n");
        return EXIT_FAILURE;
    }
    size_t mem_bytes = DATA.bank_cap * sizeof(Neutron);
    GLOB.mem_nbank = mem_bytes;
    fprintf(stdout, "Memory allocated for neutron bank: %.2f MB\n", (double)mem_bytes / (1024.0 * 1024.0));
    fprintf(stdout, "DONE.\n");

    /* ########################################################################################## */

    /* Initialize look-up tables */

    initTrigTables();
    
    /* If launched in check mode, exit now */

    if (GLOB.norun)
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

        fprintf(stdout, "Running external source simulation for %ld cycles with %ld neutrons each.\n", 
            GLOB.n_cycles, GLOB.n_particles);

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
