#include "header.h"

void initFissionNeutron(Neutron *parent, Neutron *new_neutron, long idx)
{
    /* Initialize misc. parameters */

    new_neutron->status = NEUTRON_ALIVE;
    new_neutron->id = (DATA.cur_gen - 1) * GLOB.n_particles + idx;
    new_neutron->mat_idx = -1;
    new_neutron->last_mt = 0;
    new_neutron->last_nuc = 0;
    new_neutron->fission_yield = 0;

    new_neutron->path_length = 0.0;
    new_neutron->fast_path_length = 0.0;
    new_neutron->slowing_down_time = 0.0;

    /* Inherit parent time for time-resolved generational history */

    new_neutron->time = parent->time;

    /* Inherit location of fission site */

    new_neutron->x = parent->x;
    new_neutron->y = parent->y;
    new_neutron->z = parent->z;

    /* Inherit seed */

    new_neutron->seed = parent->seed + idx;

    /* Inherit generation incremented by one */

    new_neutron->genc = parent->genc + 1;

    /* Derive rng state */

    xoshiro256ss_seed(&new_neutron->state, new_neutron->seed);

    /* Isotropic initial direction */

    sampleNeutronDirection(new_neutron);

    /* Energy from Watts distribution */

    new_neutron->E = sampleMaxwellianEnergy(new_neutron, TNUC_FISSION);
}