#ifndef RNG_H
#define RNG_H

#include <stdint.h>

#include "data.h"

/* --- Function declarations --- */

void sampleNeutronDirection(Neutron *n);

double sampleDistanceToCollision(Neutron *n);

int sampleCollisionNuclide(Neutron *n);

int sampleInteractionType(Neutron *n, Nuclide *nuc);

double sampleMaxwellianEnergy(Neutron *n, double T);

void sampleIsotropicDirection(xoshiro256ss_state *state, double *u, double *v, double *w);

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

#endif
