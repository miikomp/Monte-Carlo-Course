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
#include <float.h>
#include <gd.h>

/* Other headers */
#include "rng.h"
#include "endfmts.h"
#include "data.h"

/* Verbosity level: 0 = standard output, 1 = increased output, 2 = all the output */
extern int VERBOSITY;

static inline double getNormalizationFactor(void)
{
    return (GLOB.norm_factor > 0.0) ? GLOB.norm_factor : 1.0;
}

static inline const char *getNormalizationModeDescription(void)
{
    switch (GLOB.norm_mode)
    {
        case NORM_POWER:   return "power";
        case NORM_SRCRATE: return "source rate";
        default:           return "unity";
    }
}

/* --- Constants --- */

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define TRIG_LOOKUP_TABLE_SIZE 10000

#define UNION_ENERGY_GRID_SIZE 500
#define UNION_ENERGY_MIN 1e-11
#define UNION_ENERGY_MAX 20.0

#define DELIMS " \t\r\n"

#define CI95_FACTOR 1.959963984540054
#define M_PI 3.14159265358979323846
#define SQRT3 1.7320508075688772
#define NA 6.02214076e23
#define BOLTZMANN 8.617333262145e-11   /* MeV/K */
#define AMU_TO_MEV_C2 931.49410242     /* MeV/c^2 */
#define MASS_NEUTRON 1.00866491588     /* amu */
#define INV_MASS_NEUTRON (1.0/MASS_NEUTRON)
#define C_LIGHT 2.99792458e10          /* cm/s */

#define BARN_TO_CM2 1.0e-24
#define TNUC_FISSION 1.2895           /* MeV Nuclear temperature of U235 fission */
#define E_FG_LIMIT 0.0002             /* MeV energy cutoff for free gas model at 200eV */
#define E_THERMAL 5.0e-7              /* MeV Thermal threshold energy (Cadmium cutoff) */
#define E_EPITHERMAL 0.1              /* MeV Epithermal threshold energy */
#define STEP_INTPL 1e-6             /* Interpolation distance to cross boundaries and avoid floating point errors */

#define DEFAULT_NEEDLE_LENGTH 0.85
#define DEFAULT_LINE_SPACING  1.0
#define MIN_BANK_SIZE 1000

#define TIME_BIN_WIDTH 1e-9    /* seconds */

/* --- Enums --- */

enum RUN_MODES{
    RUNMODE_EXTERNAL_SOURCE = 0,
    RUNMODE_CRITICALITY,
};

enum NEUTRON_STATUS{
    NEUTRON_ALIVE = 0,
    NEUTRON_DEAD_FISSION,
    NEUTRON_DEAD_CAPTURE,
    NEUTRON_DEAD_LEAKAGE,
    NEUTRON_DEAD_TERMINATED
};

enum SRC_TYPES {
    SRC_MONO_POINT = 1,
    SRC_MONO_REGION,
    SRC_FISS_MAT = 18
};

enum DET_TYPES{
    DET_REACTION_RATE = 1
};

enum CELL_ERRORS {
    CELL_ERR_OK = 0,
    CELL_ERR_UNDEFINED,
    CELL_ERR_OVERLAP
};

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
 * @brief Process and validate input data stored in GLOB. Raise errors when applicable.
 * 
 * @return int 0 on success 1 on failure
 */
int processInput();

/**
 * @brief Process cross section library data from given file path into
 * a temporary nuclide library that is put into the output pointers.
 * 
 * @return int 0 on success 1 on failure
 */
int processXsData(TempNucDataLib **outlib, size_t *noutlib);

/**
 * @brief Resolve material data for all materials using the temporary nuclide library.
 * Handles freeing the temporary library after use.
 * 
 * @param lib Pointer to the temporary nuclide library
 * @param nlib Number of nuclides in the library
 * @return int 0 on success 1 on failure
 */
int resolveMaterials(TempNucDataLib *lib, size_t nlib);

/**
 * @brief Resolve all cell data, populating the surf_idxs arrays and material indices.
 * Checks for valid surface and material names.
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveCells();

/**
 * @brief Search for the cell containing the point (x, y, z).
 * Sets error flag to "undefined" if point is outside all cells, and "overlap" if more than one cell
 * occupy the point. Boundary cases are handled correctly so only true overlaps are flagged.
 * 
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 * @param err Pointer to an error flag
 * @param lx ptr to output resolved local X-coordinate
 * @param ly ptr to output resolved local Y-coordinate
 * @param lz ptr to output resolved local Z-coordinate
 * @param lidx ptr to output lattice element idx
 * @param plx ptr to put parent cell local X-coordinate
 * @param ply ptr to put parent cell local Y-coordinate
 * @param plz ptr to put parent cell local Z-coordinate
 * @param p_cell_idx ptr to put parent cell idx in DATA.cells
 * 
 * @return long Index of the cell containing the point, or -1 on failure
 */
long cellSearch(double x, double y, double z, int *err, double *lx, double *ly, double *lz, long* lidx,
                double *plx, double *ply, double *plz, long *p_cell_idx);


/**
 * @brief Apply boundary condition to the coordinate and direction vector.
 * 
 * @param x 
 * @param y 
 * @param z 
 * @param u 
 * @param v 
 * @param w 
 */
void applyBoundaryConditions(double *x, double *y, double *z, double *u, double *v, double *w);

/**
 * @brief Resolve outer boundary of root universe. Calculates and puts outer bounds into DATA.
 * At the moment the outer boundary must be defined using a single surface. Axial bounds are
 * possible using truncated versions of infinite prisms and cylinders.
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveOuterBounds();

/**
 * @brief Allocates and populates the universe array in DATA. 
 * Universes own indeces to cells held in DATA.cells.
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveUniverses();

/**
 * @brief Resolve all lattice data, populating the uni_idxs arrays and checking for valid universe names.
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveLattices();

/**
 * @brief Resolve all geometric transformations, by linking them to the specified targets.
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveTransformations();

/**
 * @brief Determines material volumes by sampling random points within the outer boundaries.
 * The routine is parallelised using OpenMP. Results are printed to stdout and saved under each material.
 * 
 * @return int 0 on succes, 1 on failure
 */
int checkVolumes();


/**
 * @brief Determines material volumes by sampling random lines between outer bounds of the geometry.
 * The routine is parallelised using OpenMP. Results are printed to stdout and saved under each material.
 * 
 * @return int 0 on succes, 1 on failure
 */
int checkVolumes2();

/**
 * @brief Compute macroscopic cross sections for all resolved materials.
 *
 * @return int 0 on success, 1 on failure
 */
int computeMacroXs(void);

/**
 * @brief Process detector data and put parameters
 * 
 * @return int 0 on success, 1 on failure
 */
int resolveDetectors(void);


/**
 * @brief Compute detector bin idx for a neutron context
 * 
 * @param n ptr to neutron
 * @param det_idx detector idx in DATA.detectors
 * @return long bin idx, or -1 if no valid bin for neutron
 */
long computeDetectorBin(Neutron *n, size_t det_idx);

/**
 * @brief Process scored detector results and output relevant data into file and stdout.
 * 
 */
void processDetectorResults();

/**
 * @brief Samples a source and fills the neutron bank with initial neutrons corresponding to 
 * the number of neutrons per generation.
 * 
 * @return long 0 on succesful source sampling, or -1 on failure
 */
long sampleInitialSource(void);

/**
 * @brief Runs the external source transport simulation. 
 * 
 * @return int 0 on success 1 on failure
 */
int runExternalSourceSimulation(void);

/**
 * @brief Run the criticality source simulation.
 * 
 * @return int 0 on success 1 on failure
 */
int runCriticalitySimulation(void);

/**
 * @brief Process transport results and output relevant data into file and stdout.
 * 
 */
void processTransportResults(void);


/**
 * @brief Return distance to nearest boundary from point (x, y, z) in direction (u, v, w) or INFINITY
 * when no boundary is in line of sight.
 * 
 * @param x 
 * @param y 
 * @param z 
 * @param u 
 * @param v 
 * @param w 
 * @return double 
 */
double distanceToNearestBoundary(double x, double y, double z, double u, double v, double w);

/**
 * @brief Process detector tallies gathered during transport and write outputs.
 *
 */
void processDetectorResults(void);

/**
 * @brief Build a fission neutron bank from fission sites.
 * 
 * @return int 0 on success, 1 on failure
 */
int buildFissionBank(void);

/**
 * @brief Handle an elastic scatter event between a neutron and a nuclide.
 * Depending on neutron energy scattering is handled from either
 *  - Stationary target (E >= 200eV)
 *  - Free gas model (E < 200eV)
 *
 * @param n Pointer to neutron
 * @param nuc Pointer to nuclide
 */
void handleElasticScatter(Neutron *n, Nuclide *nuc);

/**
 * @brief Handle inelastic level scattering.
 * 
 * @param n Pointer to neutron
 * @param nuc Pointer to nuclide
 * @param mt ENDF reaction mode identifier (For getting Q value etc.)
 */
void handleInelasticScatter(Neutron *n, Nuclide *nuc, int mt);

/**
 * @brief Handle a fission event. Samples the number of neutrons produced, creation of fission
 * neutrons is not handled here but when building the fission bank from sampled fission sites.
 *
 * @param n Pointer to neutron
 * @param nuc Pointer to nuclide
 */
void handleFission(Neutron *n, Nuclide *nuc);

/**
 * @brief Get the Total Macroscopic xs for a given material and neutron energy
 * 
 * @param E Energy in MeV
 * @param mat Pointer to material
 * @return double total macroscopic cross section in 1/cm, or -1.0 on failure
 */
double getTotalMacroscopicXS(const double E, Material* mat);

/**
 * @brief Get the Total Microscopic xs for a given nuclide and neutron energy
 * 
 * @param E Energy in MeV
 * @param nuc Pointer to nuclide
 * @return double total microscopic cross section in barns, or -1.0 on failure
 */
double getTotalMicroscopicXS(const double E, Nuclide* nuc);

/**
 * @brief Get the Microscopic xs for a given MT reaction and neutron energy
 * 
 * @param E Energy in MeV
 * @param xs_table Pointer to XS table
 * @return double microscopic cross section in barns, or -1.0 on failure
 */
double getMicroscopicXS(const double E, XsTable* xs_table);

/**
 * @brief Get the Material At Position object
 * 
 * @param x X-coordinate in cm
 * @param y Y-coordinate in cm
 * @param z Z-coordinate in cm
 * @param err ptr to error flag set by cellSearch
 * @return long Material index, or -1 if not found
 */
long getMaterialAtPosition(double x, double y, double z, int* err);

/**
 * @brief Get the Velocity cm/s for a given neutron energy
 * 
 * @param E Neutron energy in MeV
 * @return double Velocity in cm/s, or -1.0 on failure
 */
double getVelocityCmPerS(const double E);

/**
 * @brief Check if any cutt-off has been reached, terminate neutron if yes.
 * 
 * @param n Pointer to Neutron
 * 
 * @return int 0 if no cut-off reached, 1 if cut-off reached and neutron terminated
 */
int checkNeutronCutoff(Neutron *n);

/**
 * @brief Initialize a fission neutron.
 * 
 * @param parent Pointer to the parent neutron
 * @param new_neutron Pointer to the new neutron
 * @param idx Index of the new neutron in the fission bank (used for setting unique IDs)
 */
void initFissionNeutron(Neutron *parent, Neutron *new_neutron, long idx);

/**
 * @brief Test if a point (X, Y, Z) is inside a given surface. Computes surface equation and returns raw value.
 * 
 * @param type Surface type
 * @param params Surface parameters
 * @param nparams Number of parameters
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 * @return double < 0 if inside, > 0 if outside, 0 if on surface
 */
double surfaceTest(SurfaceTypes type, double *params, size_t nparams, double x, double y, double z);

/**
 * @brief Apply transformation to point and direction vector. Direction vector can be NULL.
 * 
 * @param T Pointer to transformation
 * @param x 
 * @param y 
 * @param z 
 * @param u 
 * @param v 
 * @param w 
 */
void applyTransformation(Transform* T, double *x, double *y, double *z, double *u, double *v, double *w);

/**
 * @brief Return distance to surface
 * 
 * @param type Surface type
 * @param params array of surface parameters
 * @param n_params number of surface parameters
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 * @param u X-component of direction
 * @param v Y-component of direction
 * @param w Z-component of direction
 * @return double distance to surface or INFINITY is not in line-of-sight
 */
double surfaceDistance(SurfaceTypes type, double* params, size_t n_params, double x, double y, double z, double u, double v, double w);

/**
 * @brief Plot geometry according to specified plotter instances. 
 * Track plotting is also handled here conditionally
 * 
 * @return long 0 success, 1 on failure
 */
long plotGeometry();

/* --- Inline functions --- */

/* Compare two doubles, used for quicksort */
static inline int cmpDouble(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da < db) ? -1 : (da > db);
}

/* Linear–linear interpolation on one segment [E1, E2].
   Returns xs(E) from (E1, xs1) and (E2, xs2). */
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
