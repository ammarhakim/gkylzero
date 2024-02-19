#pragma once

#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

void euler2d_run_level1(int argc, char **argv, evalf_t eval, double gas_gamma, int baseNx, int baseNy, int ref_factor, double refined_x1, double refined_y1,
    double refined_x2, double refined_y2, double coarse_x1, double coarse_y1, double coarse_x2, double coarse_y2, double cfl_frac, double t_end);