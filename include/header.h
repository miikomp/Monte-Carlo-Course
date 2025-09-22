#ifndef HEADER_H
#define HEADER_H

/* --- Standard libraries --- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <omp.h>
#include <errno.h>
#include <time.h>
#include <stdalign.h>
#include <unistd.h>
#include <limits.h>
#include <stdbool.h>

/* Other headers */
#include "rng.h"
#include "data.h"

/* Verbosity level: 0 = standard output, 1 = increased output, 2 = all the output */
extern int VERBOSITY;

/* --- Constants --- */

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define RUNMODE_TRANSPORT   0
#define RUNMODE_CIRCLE_PI   1
#define RUNMODE_BUFFONS_PI  2
#define RUNMODE_CHECK       3

#define PRG_BAR_WIDTH 50
#define TRIG_LOOKUP_TABLE_SIZE 10000

#define UNION_ENERGY_GRID_SIZE 500
#define UNION_ENERGY_MIN 1e-11
#define UNION_ENERGY_MAX 20.0

#define DELIMS " \t\r\n"

#define M_PI 3.14159265358979323846
#define NA 6.02214076e23
#define BOLTZMANN 8.617333262145e-11   /* MeV/K */

#define BARN_TO_CM2 1.0e-24
#define NT_FISSION 1.2895              /* MeV Nuclear temperature of U235 fission */

#define DEFAULT_NEEDLE_LENGTH 0.85
#define DEFAULT_LINE_SPACING  1.0
#define MIN_BANK_SIZE 1000

enum SRC_TYPES{
    SRC_MONO_POINT = 1
};

enum DET_TYPES{
    DET_REACTION_RATE = 1
};

/* inverse of 32-bit max integer */

#define INV_INT32_MAX (1.0/(double)UINT32_MAX)

/* --- Function declaration --- */

/**
 * @brief Reads an input file and setups data and parameters. Checks for valid keywords 
 * at the start of each line. When found tries to parse necessary arguments. 
 * Skips comments from "#" to end of line. Raises errors and warnings when applicable.
 * 
 * @return long – Number of keyword arguments succesfully parsed
 */
long readInput();

/**
 * @brief Initialize trigonometric lookup tables for sin and cos
 * 
 */
void initTrigTables();

/**
 * @brief Single threaded outerloop of the quarter circle pi approximation.
 * 
 * @param sm seed for seed scrambler
 * @param seeds thread private seed storage
 * @return int 0 on success 1 on failure
 */
int runCirclePi(uint64_t sm);

/**
 * @brief Single threaded outerloop of Buffon's needle simulation.
 * 
 * @param sm seed for seed scrambler
 * @param seeds thread private seed storage
 * @return int 0 on success 1 on failure
 */
int runBuffonsPi(uint64_t sm);

/**
 * @brief Summarize an array of results by calculating and printing basic statistics like
 * mean, standard deviation, 95% confidence interval and figure-of-merit (FOM). Used for the two
 * pi estimation methods.
 * 
 * @param results ptr to results array with GLOB.n_generations items
 */
void summarizePiResultsArray(const double *results);

/**
 * @brief Process and validate input data stored in GLOB. Raise errors when applicable.
 * 
 * @return int 0 on success 1 on failure
 */
int processInput();

/**
 * @brief Process cross section library data from given file path and resolve data for all materials
 * 
 * @return int 0 on success 1 on failure
 */
int processXsData(void);

/**
 * @brief Compute macroscopic cross sections for all resolved materials.
 *
 * @return int 0 on success, 1 on failure
 */
int computeMacroXs(void);

/**
 * @brief Samples a source and fills the neutron bank with initial neutrons corresponding to 
 * the number of neutrons per generation.
 * 
 * @return long 0 on succesful source sampling, or -1 on failure
 */
long sampleInitialSource(void);

/**
 * @brief Run the transport simulation.
 * 
 * @return int 0 on success 1 on failure
 */
int runTransport(void);

/* Compare two doubles, used for quicksort */
static inline int cmpDouble(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da < db) ? -1 : (da > db);
}

/* Linear–linear interpolation on one segment [E1, E2].
   Returns xs(E) from (E1, xs1) and (E2, xs2).
   Assumes energies in MeV, cross sections in barns. */
static inline double linlin(double E1, double sigma1,
                            double E2, double sigma2,
                            double E)
{
    /* Handle degenerate/unsorted input */
    if (E2 == E1) 
        return 0.5*(sigma1 + sigma2);
    if (E2 <  E1) 
    {
        double tE=E1; E1=E2; E2=tE;
        double ts=sigma1; sigma1=sigma2; sigma2=ts;
    }
    const double t = (E - E1) / (E2 - E1);
    return sigma1 + t * (sigma2 - sigma1);
}

/* Same as above, but clamps outside the segment:
   E <= E1 -> xs1,  E >= E2 -> xs2. */
static inline double linlin_clamp(double E1, double sigma1,
                                  double E2, double sigma2,
                                  double E)
{
    if (E2 == E1) 
        return 0.5*(sigma1 + sigma2);
    if (E2 <  E1) 
    {
        double tE=E1; E1=E2; E2=tE;
        double ts=sigma1; sigma1=sigma2; sigma2=ts;
    }
    if (E <= E1) 
        return sigma1;
    if (E >= E2) 
        return sigma2;
    const double t = (E - E1) / (E2 - E1);
    return sigma1 + t * (sigma2 - sigma1);
}

/* Look-up tables etc. */

extern double sin_table[TRIG_LOOKUP_TABLE_SIZE];
extern double cos_table[TRIG_LOOKUP_TABLE_SIZE];
extern double tan_table[TRIG_LOOKUP_TABLE_SIZE];

#endif
