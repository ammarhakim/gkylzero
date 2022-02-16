#pragma once 
#include <gkyl_util.h> 
EXTERN_C_BEG 

long fem_parproj_num_nodes_global_1x_ser_p1(const int numCellsPar);
long fem_parproj_num_nodes_global_1x_ser_p2(const int numCellsPar);
long fem_parproj_num_nodes_global_1x_ser_p3(const int numCellsPar);

long fem_parproj_num_nodes_global_3x_ser_p1(const int numCellsPar);
long fem_parproj_num_nodes_global_3x_ser_p2(const int numCellsPar);
long fem_parproj_num_nodes_global_3x_ser_p3(const int numCellsPar);

long fem_parproj_num_nodes_global_1x_tensor_p1(const int numCellsPar);
long fem_parproj_num_nodes_global_1x_tensor_p2(const int numCellsPar);

long fem_parproj_num_nodes_global_3x_tensor_p1(const int numCellsPar);
long fem_parproj_num_nodes_global_3x_tensor_p2(const int numCellsPar);

EXTERN_C_END 
