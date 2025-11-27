#include "header.h"

typedef bool (*MtGroupMatcher)(long mt);

typedef struct {
    long           response_id;
    MtGroupMatcher matches;
} MacroResponseGroup;

static bool matchCapture(long mt)     { return MT_IS_CAPTURE((int)mt); }
static bool matchElastic(long mt)     { return MT_IS_ELASTIC_SCATTER((int)mt); }
static bool matchInelastic(long mt)   { return MT_IS_INELASTIC_SCATTER((int)mt); }
static bool matchFission(long mt)     { return MT_IS_FISSION((int)mt); }

static bool macroscopicMtMatches(long response_id, long last_mt)
{
    if (response_id == DETECTOR_RESPONSE_MACRO_TOTAL)
        return true;

    if (last_mt <= 0)
        return false;

    /* Extend this table to add new grouped macroscopic responses. */
    static const MacroResponseGroup groups[] = {
        { DETECTOR_RESPONSE_MACRO_CAPTURE,   matchCapture },
        { DETECTOR_RESPONSE_MACRO_ELASTIC,   matchElastic },
        { DETECTOR_RESPONSE_MACRO_INELASTIC, matchInelastic },
        { DETECTOR_RESPONSE_MACRO_FISSION,   matchFission }
    };

    for (size_t i = 0; i < sizeof(groups) / sizeof(groups[0]); ++i)
    {
        if (response_id == groups[i].response_id)
            return groups[i].matches(last_mt);
    }

    long req_mt = labs(response_id);
    if (req_mt > 0 && req_mt != 1)
        return last_mt == req_mt;

    return true;
}

static bool macroscopicResponseMatches(const DetectorResponseBin *bin, const Neutron *n)
{
    if (!macroscopicMtMatches(bin->response_id, n->last_mt))
        return false;

    if (bin->rmat_id < 0 || n->mat_idx < 0)
        return false;

    return n->mat_idx == bin->rmat_id;
}

static bool microscopicResponseMatches(const DetectorResponseBin *bin, const Neutron *n, const Material *mat)
{
    long req_mt = bin->response_id;
    if (req_mt != 1 && n->last_mt != req_mt)
        return false;

    if (bin->rmat_id <= 0 || !mat)
        return false;

    if (n->last_nuc < 0 || (size_t)n->last_nuc >= mat->n_nucs)
        return false;

    long last_za = mat->nucs[n->last_nuc].nuc_data.ZA;
    return last_za == bin->rmat_id;
}

static bool responseMatchesBin(const DetectorResponseBin *bin, const Neutron *n, const Material *mat)
{
    if (bin->response_id < 0)
        return macroscopicResponseMatches(bin, n);

    if (bin->response_id > 0)
        return microscopicResponseMatches(bin, n, mat);

    return true;
}

static long detectorLocateBin(const double *edges, size_t n_bins, double value)
{
    if (!edges || n_bins == 0)
        return -1;

    if (!isfinite(value))
        return -1;

    if (value < edges[0] || value > edges[n_bins])
        return -1;

    if (value == edges[n_bins])
        return (long)n_bins - 1;

    size_t left = 0;
    size_t right = n_bins;

    while (left < right)
    {
        size_t mid = (left + right) / 2u;
        if (value < edges[mid])
        {
            right = mid;
            continue;
        }

        if (value >= edges[mid + 1u])
        {
            left = mid + 1u;
            continue;
        }

        return (long)mid;
    }

    return -1;
}

long computeDetectorBin(Neutron *n, size_t det_idx)
{
    if (!n || det_idx >= DATA.n_detectors)
        return -1;

    Detector *det = &DATA.detectors[det_idx];
    if (!det->responses.active || det->responses.n_bins == 0 || det->responses.bins == NULL)
        return -1;

    if (det->layout.total_bins == 0)
        return -1;

    if (det->has_material_filter)
    {
        if (n->mat_idx < 0 || (long)n->mat_idx != det->material_filter_index)
            return -1;
    }

    size_t bin_index = 0u;

    if (det->time.active)
    {
        long t_bin = detectorLocateBin(det->time.edges, det->time.n_bins, n->time);
        if (t_bin < 0)
            return -1;
        bin_index += (size_t)t_bin * det->layout.axes[DETECTOR_AXIS_TIME].stride;
    }

    if (det->energy.active)
    {
        long e_bin = detectorLocateBin(det->energy.edges, det->energy.n_bins, n->E);
        if (e_bin < 0)
            return -1;
        bin_index += (size_t)e_bin * det->layout.axes[DETECTOR_AXIS_ENERGY].stride;
    }

    const double coords[3] = { n->x, n->y, n->z };
    DetectorCartesianAxis *mesh_axes[3] = { &det->mesh_x, &det->mesh_y, &det->mesh_z };
    DetectorAxisKind mesh_layout[3] = {
        DETECTOR_AXIS_MESH_X,
        DETECTOR_AXIS_MESH_Y,
        DETECTOR_AXIS_MESH_Z
    };

    for (size_t a = 0; a < 3; a++)
    {
        DetectorCartesianAxis *axis = mesh_axes[a];
        if (!axis->active)
            continue;

        long m_bin = detectorLocateBin(axis->edges, axis->n_bins, coords[a]);
        if (m_bin < 0)
            return -1;
        bin_index += (size_t)m_bin * det->layout.axes[mesh_layout[a]].stride;
    }

    size_t response_stride = det->layout.axes[DETECTOR_AXIS_RESPONSE].stride;
    if (response_stride == 0)
        response_stride = 1u;

    Material *mat = NULL;
    if (n->mat_idx >= 0 && (size_t)n->mat_idx < DATA.n_mats)
        mat = &DATA.mats[n->mat_idx];

    for (size_t r = 0; r < det->responses.n_bins; r++)
    {
        DetectorResponseBin *bin = &det->responses.bins[r];
        if (!responseMatchesBin(bin, n, mat))
            continue;

        return (long)(bin_index + response_stride * r);
    }

    return -1;
}
