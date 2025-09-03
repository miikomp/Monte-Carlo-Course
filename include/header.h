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

/* --- Function declaration --- */

/**
 * @brief Reads an input file and setups data and parameters
 * 
 * @return long – Number of keyword arguments succesfully parsed
 */
long readInput();

/* --- Constants --- */

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define RUNTYPE_CIRCLE_PI   1
#define RUNTYPE_BUFFONS_PI  2



/* --- Data structures --- */

/**
 * @brief Global struct for general information and pointers to other data
 * 
 * @param fname Input filename
 * @param n_threads Number of n_threads
 * @param seed Random number generator seed   
 * @param n_outer Number of outer iterations
 * @param n_inner Number of inner iterations
 * @param mode Type of calculation to run
 */
typedef struct {
    const char *fname;
    long        n_threads;
    uint64_t    seed;     
    uint64_t    n_outer;
    uint64_t    n_inner;    
    uint64_t    mode;
} runInfo;

extern runInfo GLOB;
#endif