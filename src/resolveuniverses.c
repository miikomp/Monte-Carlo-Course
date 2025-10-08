#include "header.h"

int resolveUniverses() {

    fprintf(stdout, "\nProcessing universes...\n");

    if (DATA.n_cells == 0)
    {
        fprintf(stderr, "[ERROR] No cells defined.\n");
        return EXIT_FAILURE;
    }

    char **names = NULL;
    size_t n_unis = 0;

    /* Collect unique universe names from cells */

    for (size_t c = 0; c < DATA.n_cells; ++c)
    {
        Cell *cell = &DATA.cells[c];
        const char *uni_name = cell->uni_name;

        if (!uni_name[0])
        {
            fprintf(stderr, "[ERROR] Cell \"%s\" has an empty universe name.\n", DATA.cells[c].name);
            goto fail;
        }

        bool seen = false;
        for (size_t u = 0; u < n_unis; ++u)
        {
            if (!strcmp(uni_name, names[u]))
            {
                seen = true;
                cell->uni_idx = (int)u;
                break;
            }
        }
        if (seen)
            continue;

        char **tmp = (char**)realloc(names, (n_unis + 1) * sizeof(char*));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            goto fail;
        }
        names = tmp;

        names[n_unis] = (char*)calloc(MAX_STR_LEN, sizeof(char));
        if (!names[n_unis])
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            goto fail;
        }
        snprintf(names[n_unis], MAX_STR_LEN, "%s", uni_name);
        cell->uni_idx = (int)n_unis;
        ++n_unis;
    }

    if (n_unis == 0)
    {
        fprintf(stderr, "[ERROR] Geometry includes no universes.\n");
        goto fail;
    }

    /* Ensure root universe "0" exists */

    size_t root_idx = n_unis;
    for (size_t u = 0; u < n_unis; ++u)
    {
        if (!strcmp(names[u], "0"))
        {
            root_idx = u;
            break;
        }
    }

    if (root_idx == n_unis)
    {
        fprintf(stderr, "[ERROR] Root universe \"0\" is not defined.\n");
        goto fail;
    }

    /* Place root universe at index 0 */

    if (root_idx != 0)
    {
        char *tmp_name = names[0];
        names[0] = names[root_idx];
        names[root_idx] = tmp_name;

        for (size_t c = 0; c < DATA.n_cells; ++c)
        {
            int idx = DATA.cells[c].uni_idx;
            if (idx == (int)root_idx)
                DATA.cells[c].uni_idx = 0;
            else if (idx == 0)
                DATA.cells[c].uni_idx = (int)root_idx;
        }
    }

    /* Put to DATA */

    DATA.unis = (Universe*)calloc(n_unis, sizeof(Universe));
    if (!DATA.unis)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        goto fail;
    }
    DATA.n_unis = n_unis;

    /* Put name and type */

    for (size_t u = 0; u < n_unis; ++u)
    {
        Universe *uni = &DATA.unis[u];
        snprintf(uni->name, sizeof(uni->name), "%s", names[u]);
        uni->type = (u == 0) ? UNI_ROOT : UNI_NORMAL;
        uni->lat_idx = -1;
        uni->n_cells = 0;
        uni->cell_idxs = NULL;
    }

    /* Free name buffer */

    for (size_t u = 0; u < n_unis; ++u)
        free(names[u]);
    free(names);

    /* Map cells to universes */

    for (size_t c = 0; c < DATA.n_cells; c++)
    {
        Cell *cell = &DATA.cells[c];
        int uidx = cell->uni_idx;

        if (uidx < 0 || (size_t)uidx >= DATA.n_unis)
        {
            fprintf(stderr, "[ERROR] Cell \"%s\" references unknown universe \"%s\".\n",
                    cell->name, cell->uni_name);
            return EXIT_FAILURE;
        }

        Universe *uni = &DATA.unis[uidx];

        int *tmp = (int*)realloc(uni->cell_idxs, (uni->n_cells + 1) * sizeof(int));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed while mapping universe cells.\n");
            return EXIT_FAILURE;
        }
        uni->cell_idxs = tmp;
        uni->cell_idxs[uni->n_cells++] = (int)c;
    }

    if (VERBOSITY >= 1)
    {
        fprintf(stdout, "  Geometry consists %zu universe(s):\n", DATA.n_unis);
        for (size_t u = 0; u < DATA.n_unis; ++u)
            fprintf(stdout, "  %zu: %s%s | %zu cell(s)\n", u, DATA.unis[u].name,
                    (u == 0) ? " (root)" : "", DATA.unis[u].n_cells);
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;

fail:
    if (names)
    {
        for (size_t u = 0; u < n_unis; ++u)
            free(names[u]);
        free(names);
    }
    return EXIT_FAILURE;
}
