#ifndef DATA_H
#define DATA_H

#include "header.h"

#define MAX_MT 200
#define MAX_STR_LEN 1024
#define MAX_PATH 4096
#define MAX_COLLISION_BINS 1000

/* --- Cross sectional data structures --- */

// Energy wise microscopic xsdata for a single MT-reaction
typedef struct {
    int     mt;
    double  Q; 
    size_t  n;
    double *E;  // MeV
    double *xs; // barns
} XsTable;

// Precomputed macroscopic cross section table on a material-specific energy grid
typedef struct {
    int     mt;   // reaction identifier
    size_t  n;    // number of energy points
    double *E;    // MeV
    double *xs;   // 1/cm
} MacroXsTable;

// Energy wise nubar data for a single fissile nuclide
typedef struct {
    size_t  n;
    double *E;  // MeV
    double *nu;
} NubarTable;

// xsdata for a single nuclide at a specific temperature
typedef struct {
    char       name[16];
    int        Z;
    int        A;
    int        ZA;      // 1000*Z + A
    double     temp;    // K
    double     AW;      // atomic weight
    int        mt_idx[MAX_MT]; // index for each MT-reaction in xs array, -1 if not present
    size_t     n_xs;
    XsTable   *xs;      // barns
    bool       has_nubar;
    NubarTable nubar;
} Nuclide;

/* --- Material data structures ---*/

// Struct for material nuclide component data
typedef struct {
    double  atom_frac;
    double  mass_frac;
    double  N_i;    // 1/cm^3
    double  N_tot;  // 1/cm^3
    Nuclide nuc_data;
} MaterialNuclide;

// Struct for material data
typedef struct {
    char     name[128];
    double   mdens;   // g/cm3
    double   adens;   // 1/b*cm2
    double   temp;    // K
    double   temp_MeV;
    size_t   n_nucs;
    MaterialNuclide *nucs;
    size_t   n_macro_xs;
    MacroXsTable *macro_xs;
} Material;

/* --- Particle data structures --- */

// Neutron data structure
typedef struct {
    bool alive;
    uint64_t id;
    double x, y, z;       // position cm
    double u, v, w;       // direction cosines
    double E;             // energy MeV
    uint64_t seed;        // RNG seed
    xoshiro256ss_state state;
    double path_length;   // accumulated path length cm
    double fission_yield; // fission neutrons produced
} Neutron;

/* --- Scoring data structures --- */

// Per generation scoring structure. Scalar to allow #pragma omp reduction
typedef struct {
    double   total_path_length;
    double   total_fission_yield;
    uint64_t total_collisions;
    uint64_t total_histories;
    size_t   max_collision_bin;
    double   collision_energy_sum[MAX_COLLISION_BINS];
    uint64_t collision_energy_count[MAX_COLLISION_BINS];
} GenerationScores;

// Reaction rate detector structure
typedef struct {
    char     name[128];         // detector name
    size_t   n_channels;        // number of tracked (nuclide, MT) pairs
    int     *materials;         // material index per channel
    int     *nuclides;          // nuclide index within material
    int     *mt_numbers;        // MT identifiers for each channel
    double  *sum_counts;        // accumulated reaction counts
    double  *sum_sq_counts;     // accumulated squared counts for variance
    double  *sum_weights;       // accumulated particle weights per channel
} ReactionRateDetector;

typedef union {
    ReactionRateDetector rr;
} DetectorDefinition;


/* --- Particle source data structures --- */

typedef struct {
    double E;       // MeV
    double x, y, z; // position cm
} MonoNeutronPointSource;

typedef union {
    MonoNeutronPointSource mono;
} SourceDefinition;

/* --- Global data structures --- */

// Collection data structure for all of the data needed during a run
typedef struct {
    uint64_t  generation;
    size_t    n_mats;
    Material *mats;
    size_t    n_bank;
    size_t    bank_cap;
    Neutron  *bank;
    SourceDefinition *src;
    uint32_t src_type;
    DetectorDefinition *detector;
    uint32_t detector_type;
} runData;

typedef struct {
    int64_t          n_generations;
    GenerationScores avg_scores;
} ResultsData;

// Global struct for general information and pointers to other data
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

    /* Cutoff parameters */
    double      energy_cutoff; // MeV

    /* Allocated memory footprint (bytes)*/
    size_t      mem_xsdata;
    size_t      mem_nbank;

    /* Buffon's needle specific parameters */
    double      needle_length;
    double      line_spacing;
} runInfo;

extern runData DATA;
extern ResultsData RES;
extern runInfo GLOB;

/* --- Misc. data structures --- */

typedef struct {
    long n_tot;
    long n_hits;
} PiResult;

#endif
