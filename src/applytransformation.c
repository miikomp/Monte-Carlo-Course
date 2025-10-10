#include "header.h"

void applyTransformation(Transform* T, double *x, double *y, double *z, double *u, double *v, double *w) {
    if (!T || T->action == TRA_NONE)
        return;

    if ((T->action == TRA_TRANSLATION) || (T->action == TRA_BOTH))
    {
        if (x) *x -= T->translation[0];
        if (y) *y -= T->translation[1];
        if (z) *z -= T->translation[2];
    }

    if ((T->action == TRA_ROTATION) || (T->action == TRA_BOTH))
    {
        if (x && y && z)
        {
            double tx = *x, ty = *y, tz = *z;
            double nx = T->rotation[0][0]*tx + T->rotation[1][0]*ty + T->rotation[2][0]*tz;
            double ny = T->rotation[0][1]*tx + T->rotation[1][1]*ty + T->rotation[2][1]*tz;
            double nz = T->rotation[0][2]*tx + T->rotation[1][2]*ty + T->rotation[2][2]*tz;
            *x = nx; *y = ny; *z = nz;
        }

        if (u && v && w)
        {
            double tu = *u, tv = *v, tw = *w;
            double nu = T->rotation[0][0]*tu + T->rotation[1][0]*tv + T->rotation[2][0]*tw;
            double nv = T->rotation[0][1]*tu + T->rotation[1][1]*tv + T->rotation[2][1]*tw;
            double nw = T->rotation[0][2]*tu + T->rotation[1][2]*tv + T->rotation[2][2]*tw;
            *u = nu; *v = nv; *w = nw;
        }
    }
}