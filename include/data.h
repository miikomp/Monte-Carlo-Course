#ifndef DATA_H
#define DATA_H

#include "header.h"

/* --- Cross sectional data structures --- */

// Energy wise xsdata for a single MT-reaction
typedef struct {
    int     mt;
    double  Q; 
    size_t  n;
    double *E;  // MeV
    double *xs; // barns
} XsTable;

// Energy wise nubar data for a single fissile nuclide
typedef struct {
    size_t  n;
    double *E;  // MeV
    double *nu;
} NubarTable;

// xsdata for a single nuclide at a specific temperature (Nuclide variant)
typedef struct {
    double     temp;    // K
    int        mt_idx[200]; // index for each MT-reaction in xs array, -1 if not present
    size_t     n_xs;
    XsTable   *xs;      // barns
    int        has_nubar;
    NubarTable nubar;
} Nuclide;

// xsdata collection for a single nuclide at all temperatures only used as a temporary 
// data structure prior to resolving materials. The Nuclide struct is the one used during transport
typedef struct {
    char            name[16];
    int             Z;
    int             A;
    double          AW;
    size_t          n_var;
    Nuclide        *var;
} NuclideData;

// Struct for material nuclide component data
typedef struct {
    int     ZA;
    char    name[16];
    double  AW;
    double  atom_frac;
    double  mass_frac;
    double  N_i;
    double  N_tot;
    Nuclide data;
} MaterialNuclide;

// Struct for material data
typedef struct {
    char     name[128];
    double   density;   // g/cm3
    double   temp;      // K
    size_t   n_nucs;
    MaterialNuclide *nucs;
} Material;

// Collection data structure for all of the data
typedef struct {
    size_t    n_mats;
    Material *mats;
} runData;

extern runData DATA;

#endif