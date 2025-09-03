#include "header.h"

/* Initialize global data structs */

runInfo GLOB = {
    .fname          = NULL,
    .n_threads    = 1,
    .seed           = 0u,
    .n_outer      = 1u,
    .n_inner      = 1u,
};