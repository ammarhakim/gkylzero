#pragma once 
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
#include <gkyl_mat_triples.h>
EXTERN_C_BEG 

long fem_parproj_num_nodes_global_1x_ser_p1(const int numCellsPar);
void fem_parproj_mass_1x_ser_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_1x_ser_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_1x_ser_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_1x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_1x_ser_p2(const int numCellsPar);
void fem_parproj_mass_1x_ser_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_1x_ser_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_1x_ser_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_1x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_1x_ser_p3(const int numCellsPar);
void fem_parproj_mass_1x_ser_p3(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_1x_ser_p3(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_1x_ser_p3(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p3(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_1x_ser_p3(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_num_nodes_global_3x_ser_p1(const int numCellsPar);
void fem_parproj_mass_3x_ser_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_3x_ser_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_3x_ser_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_3x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_3x_ser_p2(const int numCellsPar);
void fem_parproj_mass_3x_ser_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_3x_ser_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_3x_ser_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_3x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_3x_ser_p3(const int numCellsPar);
void fem_parproj_mass_3x_ser_p3(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_3x_ser_p3(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_3x_ser_p3(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p3(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_3x_ser_p3(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_num_nodes_global_1x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_1x_tensor_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_1x_tensor_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_1x_tensor_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_1x_tensor_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_1x_tensor_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_1x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_1x_tensor_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_1x_tensor_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_1x_tensor_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_1x_tensor_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_1x_tensor_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_num_nodes_global_3x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_3x_tensor_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_3x_tensor_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_3x_tensor_p1(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_3x_tensor_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_3x_tensor_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_3x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_3x_tensor_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
void fem_parproj_lhs_stencil_3x_tensor_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_parproj_weighted_lhs_stencil_3x_tensor_p2(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);
GKYL_CU_DH void fem_parproj_src_stencil_3x_tensor_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_sol_stencil_3x_tensor_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


EXTERN_C_END 
