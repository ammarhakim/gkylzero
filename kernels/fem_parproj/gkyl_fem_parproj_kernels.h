#pragma once 
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 

long fem_parproj_num_nodes_global_1x_ser_p1(const int numCellsPar);
void fem_parproj_mass_matrix_1x_ser_p1(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_1x_ser_p1(struct gkyl_mat *massout);
void fem_parproj_local_to_global_1x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_1x_ser_p2(const int numCellsPar);
void fem_parproj_mass_matrix_1x_ser_p2(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_1x_ser_p2(struct gkyl_mat *massout);
void fem_parproj_local_to_global_1x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_1x_ser_p3(const int numCellsPar);
void fem_parproj_mass_matrix_1x_ser_p3(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_1x_ser_p3(struct gkyl_mat *massout);
void fem_parproj_local_to_global_1x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);

long fem_parproj_num_nodes_global_3x_ser_p1(const int numCellsPar);
void fem_parproj_mass_matrix_3x_ser_p1(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_3x_ser_p1(struct gkyl_mat *massout);
void fem_parproj_local_to_global_3x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_3x_ser_p2(const int numCellsPar);
void fem_parproj_mass_matrix_3x_ser_p2(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_3x_ser_p2(struct gkyl_mat *massout);
void fem_parproj_local_to_global_3x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_3x_ser_p3(const int numCellsPar);
void fem_parproj_mass_matrix_3x_ser_p3(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_3x_ser_p3(struct gkyl_mat *massout);
void fem_parproj_local_to_global_3x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs);

long fem_parproj_num_nodes_global_1x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_matrix_1x_tensor_p1(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_1x_tensor_p1(struct gkyl_mat *massout);
void fem_parproj_local_to_global_1x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_1x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_matrix_1x_tensor_p2(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_1x_tensor_p2(struct gkyl_mat *massout);
void fem_parproj_local_to_global_1x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);

long fem_parproj_num_nodes_global_3x_tensor_p1(const int numCellsPar);
void fem_parproj_mass_matrix_3x_tensor_p1(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_3x_tensor_p1(struct gkyl_mat *massout);
void fem_parproj_local_to_global_3x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs);
long fem_parproj_num_nodes_global_3x_tensor_p2(const int numCellsPar);
void fem_parproj_mass_matrix_3x_tensor_p2(struct gkyl_mat *massM);
void fem_parproj_massmat_times_modtonod_3x_tensor_p2(struct gkyl_mat *massout);
void fem_parproj_local_to_global_3x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs);

EXTERN_C_END 
