#pragma once

#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

struct euler2d_level1_init {
    long baseNx;
    long baseNy;
    long ref_factor;

    double coarse_x1;
    double coarse_y1;
    double coarse_x2;
    double coarse_y2;

    double refined_x1;
    double refined_y1;
    double refined_x2;
    double refined_y2;

    evalf_t eval;
    double gas_gamma;

    double cfl_frac;
    double t_end;
};

void euler2d_run_level1(int argc, char **argv, struct euler2d_level1_init* init);