#ifndef DETECTORBUFFER_H
#define DETECTORBUFFER_H

#include <stdlib.h>
#include <stdio.h>

typedef struct {
    size_t *bins;
    double *values;
    size_t count;
    size_t capacity;
} DetectorHistoryBuffer;

static inline DetectorHistoryBuffer *createDetectorHistoryBuffers(size_t n)
{
    if (n == 0)
        return NULL;

    DetectorHistoryBuffer *bufs = (DetectorHistoryBuffer*)calloc(n, sizeof(DetectorHistoryBuffer));
    if (!bufs)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed for detector history buffers.\n");
        exit(EXIT_FAILURE);
    }
    return bufs;
}

static inline void destroyDetectorHistoryBuffers(DetectorHistoryBuffer *bufs, size_t n)
{
    if (!bufs)
        return;

    for (size_t i = 0; i < n; ++i)
    {
        free(bufs[i].bins);
        free(bufs[i].values);
    }
    free(bufs);
}

static inline void detectorHistoryBufferAccumulate(DetectorHistoryBuffer *buf, size_t bin, double value)
{
    for (size_t i = 0; i < buf->count; ++i)
    {
        if (buf->bins[i] == bin)
        {
            buf->values[i] += value;
            return;
        }
    }

    if (buf->count == buf->capacity)
    {
        size_t new_cap = buf->capacity ? buf->capacity * 2u : 8u;
        size_t *new_bins = (size_t*)realloc(buf->bins, new_cap * sizeof(size_t));
        double *new_vals = (double*)realloc(buf->values, new_cap * sizeof(double));
        if (!new_bins || !new_vals)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed for detector history buffer growth.\n");
            free(new_bins);
            free(new_vals);
            exit(EXIT_FAILURE);
        }
        buf->bins = new_bins;
        buf->values = new_vals;
        buf->capacity = new_cap;
    }

    buf->bins[buf->count] = bin;
    buf->values[buf->count] = value;
    buf->count++;
}

static inline void flushDetectorHistoryBuffers(DetectorHistoryBuffer *bufs, size_t n_detectors)
{
    if (!bufs)
        return;

    for (size_t d = 0; d < n_detectors; ++d)
    {
        DetectorHistoryBuffer *buf = &bufs[d];
        if (buf->count == 0)
            continue;

        Detector *det = &DATA.detectors[d];
        double *scores = det->scores;
        double *scores_sq = det->scores_sq;

        if (!scores || !scores_sq)
        {
            buf->count = 0;
            continue;
        }

        for (size_t i = 0; i < buf->count; ++i)
        {
            size_t bin = buf->bins[i];
            double value = buf->values[i];
            #pragma omp atomic
                scores[bin] += value;
            #pragma omp atomic
                scores_sq[bin] += value * value;
        }

        buf->count = 0;
    }
}

#endif
