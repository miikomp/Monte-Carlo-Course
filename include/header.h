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

/* Other headers */
#include "data.h"

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

#define DEFAULT_NEEDLE_LENGTH 0.85
#define DEFAULT_LINE_SPACING  1.0

#define MAX_STR_LEN 1024
#define MAX_PATH 4096

/* inverse of 32-bit max integer is used when 64-bit random integers are split to two 32-bit ones */

#define INV_INT32_MAX (1.0/(double)UINT32_MAX)

/* --- Data structures --- */

/**
 * @brief Global struct for general information and pointers to other data
 * 
 * @param fname Input filename
 * @param n_threads Number of n_threads
 * @param n_generations Number of outer iterations
 * @param n_particles Number of inner iterations
 * @param n_kwargs Number of keyword arguments succesfully parsed from the input file
 * @param seed Seed for random number generator
 * @param mode Type of calculation to run
 */
typedef struct {
    /* General parameters */
    const char *fname;
    const char *outfname;
    const char *errfname;
    char        xslibpath[MAX_STR_LEN];
    long        n_kwargs;
    uint64_t    seed;     
    long        mode;
    double      t0;
    double      t1;
    int         n_threads;

    /* Iteration parameters */
    long        n_generations;
    long        n_particles;
    long        n_inactive;    

    /* Buffon's needle specific parameters */
    double      needle_length;
    double      line_spacing;
} runInfo;

extern runInfo GLOB;

typedef struct {
    long n_tot;
    long n_hits;
} PiResult;

typedef struct {
    uint64_t s[4];
} xoshiro256ss_state;

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

/* --- Inline function declarations --- */

/* Random number generators for seed scrambling and random double */

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t xoshiro256ss(xoshiro256ss_state *state) {
    const uint64_t result = rotl(state->s[1] * 5, 7) * 9;
    const uint64_t t = state->s[1] << 17;

    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];

    state->s[2] ^= t;
    state->s[3] = rotl(state->s[3], 45);

    return result;
}

/* Uniform double within [0,1] */
static inline double randd(xoshiro256ss_state *state) {
    return (xoshiro256ss(state) >> 11) * (1.0 / 9007199254740992.0);
}


/* splitmix64 for seeding xoshiro256** */
static inline uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

/* Simple seeding for xoshiro256** */
static inline void xoshiro256ss_seed(xoshiro256ss_state *state, uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i)
        state->s[i] = splitmix64(&x);
}

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