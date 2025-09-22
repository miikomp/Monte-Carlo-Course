#include "header.h"

int runTransport(void) {

    /*
    1) On all but the first generation either:
        a) Build a fission neutron bank from fission sites of last generation
        b) Re-use last generation's bank
    2) Initialize a GenerationScores struct to hold scores for this generation
    3) Loop over all neutrons in the bank by comb sampling #pragma omp parallel for
        a) While neutron is alive repeatedly:
            i) Sample distance to next collision
            ii) Move neutron to collision site
            iii) Sample reaction type
            iv) Score path length, collisions, fission yield
            v) If fission, sample number of neutrons produced and add to fission bank
            vi) If capture kill neutron
            vii) If scatter sample new energy and direction
                - Target energy is either ignored, or free gas model is used (E < 200eV)
            viii) if maximum number of collisions or energy below cutoff kill neutron
        b) Use #pragma reduction to sum scores into GenerationScores struct
    4) After all neutrons are done, add generation scores to RES.avg_scores in some valid way (single threaded?)
    5) Repeat until all generations are done
    */

    return 0;
}