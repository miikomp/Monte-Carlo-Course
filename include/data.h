#ifndef DATA_H
#define DATA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define MAX_MT 200
#define MAX_STR_LEN 1024
#define MAX_PATH 4096
#define MAX_COLLISION_BINS 250
#define MAX_NUM_DETECTORS 16

typedef struct {
    uint64_t s[4];
} xoshiro256ss_state;

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
    XsTable   *xs;      // barns idx 0 is total
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
    MacroXsTable *macro_xs; // 1/cm idx 0 is total
} Material;

/* --- Particle data structures --- */

// Neutron data structure
typedef struct {
    int status;           // 0=alive, 1=dead by fission, 2=dead by capture, 3=dead by leakage   
    uint64_t id;
    int mat_idx;          // active material idx in DATA.mats, -1 if not set
    double x, y, z;       // position cm
    double u, v, w;       // direction cosines
    double E;             // energy MeV
    uint64_t seed;        // RNG seed
    xoshiro256ss_state state;
    double path_length;   // accumulated path length cm
    int fission_yield;   // fission neutrons produced
} Neutron;

/* --- Scoring data structures --- */

// Per generation scoring structure. And the accompynying custom reduction clause
typedef struct {
    uint64_t n_histories;
    double   total_path_length;
    long   total_fission_yield;
    uint64_t total_collisions;
    uint64_t total_captures;
    uint64_t total_elastic_scatters;
    uint64_t total_inelastic_scatters;
    uint64_t total_fissions;
    uint64_t total_leakages;
    uint64_t total_unknowns;
    size_t   max_collision_bin;
    double   collision_energy_sum[MAX_COLLISION_BINS];
    uint64_t collision_energy_count[MAX_COLLISION_BINS];
} GenerationScores;

#ifdef _OPENMP
static inline void GenerationScoresReduce(GenerationScores *restrict out,
                                          const GenerationScores *restrict in)
{
    out->total_path_length   += in->total_path_length;
    out->total_fission_yield += in->total_fission_yield;
    out->total_collisions    += in->total_collisions;
    out->n_histories         += in->n_histories;
    out->total_captures      += in->total_captures;
    out->total_elastic_scatters   += in->total_elastic_scatters;
    out->total_inelastic_scatters += in->total_inelastic_scatters;
    out->total_fissions      += in->total_fissions;
    out->total_leakages      += in->total_leakages;
    out->total_unknowns      += in->total_unknowns;
    if (in->max_collision_bin > out->max_collision_bin)
        out->max_collision_bin = in->max_collision_bin;
    for (size_t b = 0; b < MAX_COLLISION_BINS; ++b)
    {
        out->collision_energy_sum[b]   += in->collision_energy_sum[b];
        out->collision_energy_count[b] += in->collision_energy_count[b];
    }
}

#pragma omp declare reduction(+:GenerationScores: GenerationScoresReduce(&omp_out, &omp_in)) \
    initializer(omp_priv = (GenerationScores){0})
#endif

/* --- Detector structures --- */

/* Reaction rate detector structures */
typedef enum {
    RRDET_MODE_ELASTIC = 0,
    RRDET_MODE_FISSION,
    RRDET_MODE_INELASTIC,
    RRDET_MODE_CAPTURE,
    RRDET_MODE_COUNT
} ReactionRateMode;

typedef struct {
    double sum;     // accumulated score
    double sum_sq;  // accumulated score squared
} ReactionRateTally;

typedef struct {
    char  nuclide_name[16];
    int   nuclide_index;   // index inside the Material.nucs array
    size_t n_energy_bins;  // number of energy bins (>=1)
    ReactionRateTally *tallies; // flattened array [n_energy_bins][RRDET_MODE_COUNT]
} ReactionRateNuclide;

typedef enum {
    ENERGY_BIN_SPACING_LOG = 0,
    ENERGY_BIN_SPACING_LINEAR = 1
} EnergyBinSpacing;

typedef struct {
    bool   enabled;
    size_t n_bins;
    double E_min;
    double E_max;
    EnergyBinSpacing spacing;
    double *edges;            // energy bin edges, length n_bins + 1
} EnergyGrid;

typedef struct {
    char   material_name[128];  // tracked material name from input
    int    material_index;      // resolved material index, -1 if unresolved
    size_t n_nuclides;          // length of the nuclides array
    ReactionRateNuclide *nuclides;
    EnergyGrid energy_grid;
    size_t n_energy_bins;       // number of bins used for tallying (1 if grid disabled)
} ReactionRateDetector;

typedef struct {
    bool   has_material_filter;   // true if material filter requested
    char   material_name[128];    // tracked material name from input (empty if global)
    int    material_index;        // resolved material index, -1 if unresolved or unused
    EnergyGrid grid;        // energy bin specification
    double *track_length_sum;     // accumulated track-length tallies per bin
    double *track_length_sum_sq;  // accumulated squared tallies
} EnergySpectrumDetector;

/* General detector structures */
typedef enum {
    DETECTOR_TYPE_REACTION_RATE = 0,
    DETECTOR_TYPE_ENERGY_SPECTRUM,
    DETECTOR_TYPE_COUNT
} DetectorType;

typedef struct {
    DetectorType type;
    char   name[128];           // detector name from input
    uint64_t n_histories;       // number of tallied source histories
    union {
        ReactionRateDetector reaction_rate;
        EnergySpectrumDetector energy_spectrum;
    } data;
} Detector;

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
    size_t n_detectors;
    Detector *detectors[MAX_NUM_DETECTORS];
} runData;

typedef struct {
    int64_t          n_generations;
    GenerationScores *avg_scores;
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
    xoshiro256ss_state rng_state;     
    long        mode;
    double      t0;
    double      t1;
    int         n_threads;

    /* Iteration parameters */
    long        n_generations;
    long        n_particles;
    long        n_inactive;

    /* Cutoff parameters */
    double      energy_cutoff; // MeV
    long        max_collisions; // per neutron

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

/* Temporary data struct for used when parsing nuclide data at multiple temperatures */
/* The reason is that in the future interpolation could be applied to temperature adjust */
/* nuclides but as of now just the closest match is linked to a material */
typedef struct {
    char            name[16];
    int             Z;
    int             A;
    double          AW;
    size_t          n_var;
    Nuclide        *var;
} TempNucDataLib;

#endif
