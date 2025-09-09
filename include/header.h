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

/* --- Constants --- */

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define RUNMODE_CIRCLE_PI   1
#define RUNMODE_BUFFONS_PI  2
#define RUNMODE_CHECK       3

#define PRG_BAR_WIDTH 50
#define TRIG_LOOKUP_TABLE_SIZE 10000

#define DELIMS " \t\r\n"

#define M_PI 3.14159265358979323846

/* inverse of 32-bit max integer is used when 64-bit random integers are split to two 32-bit ones */

#define INV_INT32_MAX (1.0/(double)UINT32_MAX)

/* --- Data structures --- */

/**
 * @brief Global struct for general information and pointers to other data
 * 
 * @param fname Input filename
 * @param n_threads Number of n_threads
 * @param n_outer Number of outer iterations
 * @param n_inner Number of inner iterations
 * @param n_kwargs Number of keyword arguments succesfully parsed from the input file
 * @param seed Seed for random number generator
 * @param mode Type of calculation to run
 */
typedef struct {
    /* General parameters */
    const char *fname;
    int         n_threads;
    long        n_outer;
    long        n_inner;    
    long        n_kwargs;
    uint64_t    seed;     
    long        mode;
    double      t0;
    double      t1;

    /* Buffon's needle specific parameters */
    double      needle_length;
    double      line_spacing;
} runInfo;

extern runInfo GLOB;

/**
 * @brief Used for storing tallies collected during the simulation
 * 
 * @param n_tot Total number of simulated histories or points
 * @param n_hits Total number of simulated histories that satisfy a condition
 */
typedef struct {
    long   n_tot;
    long   n_hits;  
} Tallies;

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
 * @brief Initialize an empty Tallies struct
 * 
 * @return Tallies 
 */
Tallies initTallies();

/**
 * @brief Initialize trigonometric lookup tables for sin and cos
 * 
 */
void initTrigTables();

/**
 * @brief Runs the quarter circle approximation for pi
 * 
 * @param out ptr to result struct
 * @param seeds array of thread private seeds
 * @return int 0 on success 1 on failure
 */
int runCirclePiInner(PiResult *out, xoshiro256ss_state *state);

/**
 * @brief Runs the specified number of outer iterations for Buffon's Needle
 * 
 * @param out ptr to result struct
 * @param seeds array of thread private seeds
 * @return int 0 on success 1 on failure
 */
int runBuffonsPiInner(PiResult *out, xoshiro256ss_state *state);

/**
 * @brief Single threaded outerloop of the quarter pi approximation.
 * 
 * @param sm seed for seed scrambler
 * @param seeds thread private seed storage
 * @return int 0 on success 1 on failure
 */
int runCirclePi(uint64_t sm, xoshiro256ss_state *state);

/**
 * @brief Single threaded outerloop of Buffon's needle simulation.
 * 
 * @param sm seed for seed scrambler
 * @param seeds thread private seed storage
 * @return int 0 on success 1 on failure
 */
int runBuffonsPi(uint64_t sm, xoshiro256ss_state *state);

/**
 * @brief Summarize an array of results by calculating and printing basic statistics like
 * mean, standard deviation, 95% confidence interval and figure-of-merit (FOM).
 * 
 * @param results ptr to results array with GLOB.n_outer items
 */
void summarizeResultsArray(const double *results);

/**
 * @brief Process and validate input data stored in GLOB. Raise errors when applicable.
 * 
 * @return int 0 on success 1 on failure
 */
int processInput();

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

// Uniform double in [0,1]
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

/* Look-up tables etc. */

extern double sin_table[TRIG_LOOKUP_TABLE_SIZE];
extern double cos_table[TRIG_LOOKUP_TABLE_SIZE];

#endif