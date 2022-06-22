#pragma once 
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 

long fem_parproj_num_nodes_global_1x_ser_p1(const int numCellsPar);
void fem_parproj_mass_1x_ser_p1(struct gkyl_mat *matout);
void fem_parproj_nodtomod_1x_ser_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_1x_ser_p2(const int numCellsPar);
void fem_parproj_mass_1x_ser_p2(struct gkyl_mat *matout);
void fem_parproj_nodtomod_1x_ser_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_1x_ser_p3(const int numCellsPar);
void fem_parproj_mass_1x_ser_p3(struct gkyl_mat *matout);
void fem_parproj_nodtomod_1x_ser_p3(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p3(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);

long fem_parproj_num_nodes_global_3x_ser_p1(const int numCellsPar);
void fem_parproj_mass_3x_ser_p1(struct gkyl_mat *matout);
void fem_parproj_nodtomod_3x_ser_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_3x_ser_p2(const int numCellsPar);
void fem_parproj_mass_3x_ser_p2(struct gkyl_mat *matout);
void fem_parproj_nodtomod_3x_ser_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_3x_ser_p3(const int numCellsPar);
void fem_parproj_mass_3x_ser_p3(struct gkyl_mat *matout);
void fem_parproj_nodtomod_3x_ser_p3(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p3(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);

long fem_parproj_num_nodes_global_1x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_1x_tensor_p1(struct gkyl_mat *matout);
void fem_parproj_nodtomod_1x_tensor_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_1x_tensor_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_1x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_1x_tensor_p2(struct gkyl_mat *matout);
void fem_parproj_nodtomod_1x_tensor_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_1x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_1x_tensor_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);

long fem_parproj_num_nodes_global_3x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_3x_tensor_p1(struct gkyl_mat *matout);
void fem_parproj_nodtomod_3x_tensor_p1(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_3x_tensor_p1(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);
long fem_parproj_num_nodes_global_3x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_3x_tensor_p2(struct gkyl_mat *matout);
void fem_parproj_nodtomod_3x_tensor_p2(struct gkyl_mat *matout);
GKYL_CU_DH void fem_parproj_local_to_global_3x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_src_stencil_3x_tensor_p2(const double *rho, long nodeOff, const long *globalIdxs, double *bsrc);

EXTERN_C_END 
