// Private header for fem_poisson updater.
#pragma once
#include <gkyl_fem_poisson_kernels.h>

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int *numCells, const int *idx,
  long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[3]; } local2global_kern_list;

//// Serendipity local-to-global kernels.
//GKYL_CU_D
//static const local2global_kern_list ser_loc2glob_list[] = {
//  { NULL, NULL, NULL }, // No 0D basis functions
//  { NULL, fem_poisson_local_to_global_1x_ser_p1, fem_poisson_local_to_global_1x_ser_p2 },
//  { NULL, NULL, NULL }, // No 2D basis functions
//  { NULL, fem_poisson_local_to_global_3x_ser_p1, fem_poisson_local_to_global_3x_ser_p2 }
//};
//
//// Tensor product local-to-global kernels.
//GKYL_CU_D
//static const local2global_kern_list ten_loc2glob_list[] = {
//  { NULL, NULL, NULL }, // No 0D basis functions
//  { NULL, fem_poisson_local_to_global_1x_tensor_p1, fem_poisson_local_to_global_1x_tensor_p2 },
//  { NULL, NULL, NULL }, // No 2D basis functions
//  { NULL, fem_poisson_local_to_global_3x_tensor_p1, fem_poisson_local_to_global_3x_tensor_p2 }
//};
//
//GKYL_CU_D
//static local2global_t
//choose_local2global_kern(const int dim, const int basis_type, const int poly_order)
//{
//  switch (basis_type) {
//    case GKYL_BASIS_MODAL_SERENDIPITY:
//      return ser_loc2glob_list[dim].kernels[poly_order];
//    case GKYL_BASIS_MODAL_TENSOR:
//      return ten_loc2glob_list[dim].kernels[poly_order];
//    default:
//      assert(false);
//      break;
//  }
//}
