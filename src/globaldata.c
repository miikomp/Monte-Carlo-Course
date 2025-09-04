#include "header.h"

/* Initialize global data structs */

runInfo GLOB = {
    .fname          = NULL,
    .n_threads      = 1,
    .n_outer        = 1l,
    .n_inner        = 1l,
    .n_kwargs       = 0l,
    .mode           = 0l,
    .seed           = 0u,
};

Tallies initTallies() {
    Tallies ret = {
        .n_tot      = 0ll,
        .n_hits     = 0ll,
        .dis        = 0.0
    };

    return ret;
}