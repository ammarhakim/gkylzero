// Gkyl ------------------------------------------------------------------------
//
// Header file for FEM Poisson solver.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 
 
long fem_poisson_num_nodes_global_1x_ser_p1_periodicx(const int *numCells);
long fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(const int *numCells);

void fem_poisson_stiff_1x_ser_p1(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_1x_ser_p1(struct gkyl_mat *matout);
void fem_poisson_nodtomod_1x_ser_p1(struct gkyl_mat *matout);

void fem_poisson_local_to_global_1x_ser_p1_inx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_inx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);


long fem_poisson_num_nodes_global_2x_ser_p1_periodicx_periodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells);

void fem_poisson_stiff_2x_ser_p1(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_2x_ser_p1(struct gkyl_mat *matout);
void fem_poisson_nodtomod_2x_ser_p1(struct gkyl_mat *matout);

void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);


long fem_poisson_num_nodes_global_3x_ser_p1_periodicx_periodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_periodicx_periodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells);

void fem_poisson_stiff_3x_ser_p1(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_3x_ser_p1(struct gkyl_mat *matout);
void fem_poisson_nodtomod_3x_ser_p1(struct gkyl_mat *matout);

void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);



long fem_poisson_num_nodes_global_1x_ser_p2_periodicx(const int *numCells);
long fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(const int *numCells);

void fem_poisson_stiff_1x_ser_p2(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_1x_ser_p2(struct gkyl_mat *matout);
void fem_poisson_nodtomod_1x_ser_p2(struct gkyl_mat *matout);

void fem_poisson_local_to_global_1x_ser_p2_inx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_inx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);


long fem_poisson_num_nodes_global_2x_ser_p2_periodicx_periodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p2_periodicx_nonperiodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_periodicy(const int *numCells);
long fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_nonperiodicy(const int *numCells);

void fem_poisson_stiff_2x_ser_p2(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_2x_ser_p2(struct gkyl_mat *matout);
void fem_poisson_nodtomod_2x_ser_p2(struct gkyl_mat *matout);

void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);


long fem_poisson_num_nodes_global_3x_ser_p2_periodicx_periodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_periodicx_periodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_periodicx_nonperiodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_periodicx_nonperiodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_periodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_periodicy_nonperiodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_nonperiodicy_periodicz(const int *numCells);
long fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells);

void fem_poisson_stiff_3x_ser_p2(const double *dx, struct gkyl_mat *matout);
void fem_poisson_mass_times_modtonod_3x_ser_p2(struct gkyl_mat *matout);
void fem_poisson_nodtomod_3x_ser_p2(struct gkyl_mat *matout);

void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_periodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_periodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_nonperiodicy_inz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_nonperiodicy_inz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_periodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_periodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_periodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_periodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_nonperiodicy_upz_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);



EXTERN_C_END 
