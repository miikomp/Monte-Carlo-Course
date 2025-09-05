#include "header.h"

/* Initialize global data structs */

runInfo GLOB = {
    .fname          = NULL,
    .n_threads      = 1,
    .n_outer        = 0l,
    .n_inner        = 0l,
    .n_kwargs       = 0l,
    .mode           = 0l,
    .seed           = 0u,
    .t0             = 0.0,
    .t1             = 0.0,
};

Tallies initTallies() {
    Tallies ret = {
        .n_tot      = 0ul,
        .n_hits     = 0ul,
    };

    return ret;
}