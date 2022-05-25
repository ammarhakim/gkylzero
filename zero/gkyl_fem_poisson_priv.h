// Private header for fem_poisson updater.
#pragma once
#include <gkyl_fem_poisson_kernels.h>
#include <gkyl_basis.h>

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int *numCells, const int *idx,
  long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[2]; } local2global_kern_loc_list_1x;
typedef struct { local2global_kern_loc_list_1x list[3]; } local2global_kern_bcx_list_1x;

typedef struct { local2global_t kernels[4]; } local2global_kern_loc_list_2x;
typedef struct { local2global_kern_loc_list_2x list[3]; } local2global_kern_bcy_list_2x;
typedef struct { local2global_kern_bcy_list_2x list[2]; } local2global_kern_bcx_list_2x;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_bcx_list_1x ser_loc2glob_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_local_to_global_1x_ser_p1_inx_periodicx, fem_poisson_local_to_global_1x_ser_p1_upx_periodicx},
             {fem_poisson_local_to_global_1x_ser_p2_inx_periodicx, fem_poisson_local_to_global_1x_ser_p2_upx_periodicx}}, },
  // nonperiodicx
  { .list = {{NULL, NULL},
            {fem_poisson_local_to_global_1x_ser_p1_inx_nonperiodicx, fem_poisson_local_to_global_1x_ser_p1_upx_nonperiodicx},
            {fem_poisson_local_to_global_1x_ser_p2_inx_nonperiodicx, fem_poisson_local_to_global_1x_ser_p2_upx_nonperiodicx}}, }
};

GKYL_CU_D
static const local2global_kern_bcx_list_2x ser_loc2glob_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_periodicy}},
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_nonperiodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_nonperiodicy},}
    }}
  },
  // nonperiodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy}},
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy},}
    }}
  }
};

//// Tensor local-to-global kernels.
//GKYL_CU_D
//static const local2global_kern_bc_list_1x tensor_loc2glob_list_1x[] = {
//  // periodicx
//  { .list = {{NULL, NULL},
//             {fem_poisson_local_to_global_1x_tensor_p1_inx_periodicx, fem_poisson_local_to_global_1x_tensor_p1_upx_periodicx},
//             {fem_poisson_local_to_global_1x_tensor_p2_inx_periodicx, fem_poisson_local_to_global_1x_tensor_p2_upx_periodicx}}, },
//  // nonperiodicx
//  { .list = {{NULL, NULL},
//             {fem_poisson_local_to_global_1x_tensor_p1_inx_nonperiodicx, fem_poisson_local_to_global_1x_tensor_p1_upx_nonperiodicx},
//             {fem_poisson_local_to_global_1x_tensor_p2_inx_nonperiodicx, fem_poisson_local_to_global_1x_tensor_p2_upx_nonperiodicx}}, }
//};


// "Choose Kernel" based on cdim, vdim and polyorder
#define CK1(lst,poly_order,loc,bcx) lst[bcx].list[poly_order].kernels[loc]
#define CK2(lst,poly_order,loc,bcx,bcy) lst[bcx].list[bcy].list[poly_order].kernels[loc]

struct gkyl_fem_poisson {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int num_cells[GKYL_MAX_DIM];
  bool isdirperiodic[GKYL_MAX_DIM]; // =true if parallel direction is periodic.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range solve_range, solve_range_ext;
  struct gkyl_range_iter solve_iter;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_superlu_prob* prob;
  struct gkyl_mat *local_stiff; // local stiffness matrix.
  struct gkyl_mat *local_mass_modtonod; // local mass matrix times modal-to-nodal matrix.
  struct gkyl_mat *local_nodtomod; // local nodal-to-modal matrix.

  // Pointer to local-to-global kernels. 2^3, 2 (interior and upper) in each direction.
  local2global_t l2g[8];

  long *globalidx;
};


GKYL_CU_D
static void
choose_local2global_kernels(const struct gkyl_basis* basis, const bool *isdirperiodic, local2global_t *l2gout)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[GKYL_MAX_DIM] = {-1};
  for (int d=0; d<basis->ndim; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;

  int ki = 0;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
        for (int k=0; k<(int)(pow(dim,2)+0.5); k++) {
          if (dim == 1) {
            l2gout[k] = CK1(ser_loc2glob_list_1x, poly_order, k, bckey[0]);
          } else if (dim == 2) {
            l2gout[k] = CK2(ser_loc2glob_list_2x, poly_order, k, bckey[0], bckey[1]);
          }
        }
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}
