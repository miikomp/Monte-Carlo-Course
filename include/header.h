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

#define CACHELINE 64

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
    const char *fname;
    int         n_threads;
    long        n_outer;
    long        n_inner;    
    long        n_kwargs;
    uint64_t    seed;     
    long        mode;
} runInfo;

extern runInfo GLOB;

/**
 * @brief Used for storing tallies collected during the simulation
 * 
 * @param n_tot Total number of simulated histories or points
 * @param n_hits Total number of simulated histories that satisfy a condition
 */
typedef struct {
    long long   n_tot;
    long long   n_hits;  
} Tallies;

typedef struct {
    long long n_tot;
    long long n_hits;
} PiResult;

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
 * @brief Runs the quarter circle approximation for pi
 * 
 * @param out ptr to result struct
 * @param seeds array of thread private seeds
 * @return int 0 on success 1 on failure
 */
int runCirclePi(PiResult *out, uint64_t *seeds);

/**
 * @brief Summarize an array of results by calculating and printing basic statistics
 * 
 * @param results ptr to results array with GLOB.n_outer items
 */
void summarizeResultsArray(const long double *results);

/* --- Inline function declarations --- */

/**
 * @brief Fast thread-safe 64-bit RNG. xorshift* algorithm.
 * 
 * @param s ptr to state (seed)
 * @return uint64_t 
 */
static inline uint64_t xorShift64(uint64_t *s) {
    uint64_t x = *s ? *s : UINT64_C(0x106689D45497FDB5);
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}

/**
 * @brief 64-bit scrambler for seeding other RNGs
 * 
 * @param x ptr to seed
 * @return uint64_t 
 */
static inline uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

/* Random number generators for seed scrambling and random double */
/* All based on literature (Not original work) */

/**
 * @brief xorshift* based 53-bit floating point number RNG. Return between [0, 1]
 * 
 * @param s ptr to state (seed)
 * @return double 
 */
static inline double randd(uint64_t *s) {
    return (xorShift64(s) >> 11) * (1.0/9007199254740992.0);
}

#endif