#ifndef DATA_H
#define DATA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define MAX_MT 200
#define MAX_STR_LEN 256
#define MAX_PATH 4096
#define MAX_COLLISION_BINS 250
#define MAX_TIME_BINS 1000
#define MAX_NUM_DETECTORS 16

/* --- Surface types --- */
typedef enum {
    SURF_UNKNOWN = 0,
    SURF_PLANEX,     // YZ-plane
    SURF_PLANEY,     // XZ-plane
    SURF_PLANEZ,     // XY-plane
    SURF_PLANE,      // Arbitrary plane
    SURF_SPH,        // Sphere
    SURF_CYLX,       // Cylinder along X-axis
    SURF_CYLY,       // Cylinder along Y-axis
    SURF_CYLZ,       // Cylinder along Z-axis
    SURF_SQR,        // Square prism
    SURF_HEXX,       // Hexagonal prism X-type
    SURF_HEXY,       // Hexagonal prism Y-type
    SURF_CUBE,       // Cube
    SURF_CUBOID,     // Cuboid
    SURF_TORUSX,     // Elliptical torus with major radius perpendicular to X-axis
    SURF_TORUSY,     // Elliptical torus with major radius perpendicular to Y-axis
    SURF_TORUSZ      // Elliptical torus with major radius perpendicular to Z-axis
} SurfaceTypes;

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
    double     T;    // K
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
    char     name[MAX_STR_LEN];
    uint32_t rgb[3];  // colour in geometry plots
    double   mdens;   // g/cm3
    double   adens;   // 1/b*cm2
    double   T;       // K
    double   kT;      // MeV
    size_t   n_nucs;
    MaterialNuclide *nucs;
    size_t   n_macro_xs;
    MacroXsTable *macro_xs; // 1/cm idx 0 is total
} Material;

// Struct for a surface
typedef struct {
    char         name[MAX_STR_LEN];
    SurfaceTypes type;
    size_t       n_params;
    double      *params;
} Surface;

// Struct for a cell
typedef struct {
    char    name[MAX_STR_LEN];
    char    mat_name[MAX_STR_LEN]; // Stores material name given as input
    int     mat_idx;               // index in DATA.mats, -1 if outside cell
    char    uni_name[MAX_STR_LEN]; // Stores universe name given as input
    size_t  n_surfs;
    char   *surf_names;            // Defining surface names given as input
    int    *surf_idxs;             // array of surface indices (in DATA.surfs)
    int    *side;                  // > 0 outside, < 0 inside
} Cell;

/* --- Particle data structures --- */

// Neutron data structure
typedef struct {
    int status;           // defined in enum NEUTRON_STATUS
    uint64_t id;
    int mat_idx;          // active material idx in DATA.mats, -1 if not set
    double x, y, z;       // position cm
    double u, v, w;       // direction cosines
    double E;             // energy MeV
    uint64_t seed;        // RNG seed
    xoshiro256ss_state state;
    double path_length;   // accumulated path length cm
    double fast_path_length; // path length in fast region cm
    double time;        // time since birth in seconds
    double slowing_down_time;   // time in fast region seconds
    long genc;      // generation counter
    int fission_yield;   // fission neutrons produced
} Neutron;

/* --- Scoring data structures --- */

// Per generation scoring structure. Only scalars and fixed size arrays to allow use of OpenMP reduction.
typedef struct {
    uint64_t n_histories;
    double   total_path_length;
    double   total_fast_path_length;
    double   total_time;
    double   total_slowing_down_time;
    long   total_fission_yield;
    uint64_t total_collisions;
    uint64_t total_captures;
    uint64_t total_elastic_scatters;
    uint64_t total_inelastic_scatters;
    uint64_t total_fissions;
    uint64_t total_thermal_fissions;
    uint64_t total_fast_fissions;
    uint64_t total_leakages;
    uint64_t total_unknowns;
    uint64_t total_terminated;
    size_t   max_collision_bin;
    double   collision_energy_sum[MAX_COLLISION_BINS];
    uint64_t collision_energy_count[MAX_COLLISION_BINS];
    double   fission_time_yield[MAX_TIME_BINS];
    uint64_t fission_time_events[MAX_TIME_BINS];
} TransportRunScores;

/* Custom OpenMP reduction clause for the TransportRunScores structure */

#ifdef _OPENMP
static inline void TransportRunScoresReduce(TransportRunScores *restrict out,
                                          const TransportRunScores *restrict in)
{
    out->total_path_length   += in->total_path_length;
    out->total_fast_path_length += in->total_fast_path_length;
    out->total_time         += in->total_time;
    out->total_slowing_down_time    += in->total_slowing_down_time;
    out->total_fission_yield += in->total_fission_yield;
    out->total_collisions    += in->total_collisions;
    out->n_histories         += in->n_histories;
    out->total_captures      += in->total_captures;
    out->total_elastic_scatters   += in->total_elastic_scatters;
    out->total_inelastic_scatters += in->total_inelastic_scatters;
    out->total_fissions      += in->total_fissions;
    out->total_thermal_fissions += in->total_thermal_fissions;
    out->total_fast_fissions += in->total_fast_fissions;
    out->total_leakages      += in->total_leakages;
    out->total_unknowns      += in->total_unknowns;
    out->total_terminated      += in->total_terminated;
    if (in->max_collision_bin > out->max_collision_bin)
        out->max_collision_bin = in->max_collision_bin;
    for (size_t b = 0; b < MAX_COLLISION_BINS; ++b)
    {
        out->collision_energy_sum[b]   += in->collision_energy_sum[b];
        out->collision_energy_count[b] += in->collision_energy_count[b];
    }
    for (size_t b = 0; b < MAX_TIME_BINS; ++b)
    {
        out->fission_time_yield[b]   += in->fission_time_yield[b];
        out->fission_time_events[b] += in->fission_time_events[b];
    }
}

#pragma omp declare reduction(+:TransportRunScores: TransportRunScoresReduce(&omp_out, &omp_in)) \
    initializer(omp_priv = (TransportRunScores){0})
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
    char   name[MAX_STR_LEN];           // detector name from input
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

/* ############################################################################################## */
/* --- Global data structures --- */

// Collection data structure for all of the data needed during a run
typedef struct {
    /* System */
    size_t    n_mats;
    Material *mats;
    size_t    n_surf;
    Surface  *surfs;
    size_t    n_cells;
    Cell     *cells;

    /* Bounds */
    double    x_min, x_max;
    double    y_min, y_max;
    double    z_min, z_max;

    /* Runtime */
    uint64_t  cur_gen;
    uint64_t  cur_cycle;
    size_t    n_bank;
    size_t    bank_cap;
    Neutron  *bank;

    /* Sources */
    SourceDefinition *src;
    uint32_t src_type;

    /* Detectors*/
    size_t n_detectors;
    Detector *detectors[MAX_NUM_DETECTORS];
} runData;

typedef struct {
    int64_t          n_iterations;
    TransportRunScores *avg_scores;
} ResultsData;

// Global struct for general information and pointers to other data
typedef struct {
    /* General parameters */
    const char *inputf;
    const char *inputfname;
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
    double      nbuf_factor;
    bool        norun;

    /* Iteration parameters */
    long        n_generations;  // criticality simulation
    long        n_cycles;       // external source simulation
    long        n_particles;    // particles per cycle/generation
    long        n_inactive;     // number of inactive generations/cycles

    /* Cutoff parameters */
    double      energy_cutoff; // MeV
    long        collision_cutoff; // per neutron
    long        generation_cutoff; // in external source mode to limits lengths of fission histories
    double      time_cutoff; // seconds

    /* Allocated memory footprint (bytes)*/
    size_t      mem_xsdata;
    size_t      mem_nbank;
} runInfo;

extern runData DATA;
extern ResultsData RES;
extern runInfo GLOB;

/* --- Misc. data structures --- */

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
