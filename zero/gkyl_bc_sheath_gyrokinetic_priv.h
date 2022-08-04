#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_kernels.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*sheath_reflectedf_t)(const double wv, const double dv,
  const double vlowerSq, const double vupperSq, const double q2Dm,
  const double *phi, const double *phiWall, const double *f, double *fRefl);

typedef struct { sheath_reflectedf_t kernels[3]; } sheath_reflectedf_kern_list;  // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const sheath_reflectedf_kern_list ser_sheath_reflect_list[] = {
//  { bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p2 },
  { NULL, NULL},
  { bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p2 },
  { NULL, NULL},
//  { bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p2 },
  { NULL, NULL},
};

// Serendipity  kernels.
GKYL_CU_D
static const sheath_reflectedf_kern_list tensor_sheath_reflect_list[] = {
//  { NULL, bc_sheath_gyrokinetic_reflectedf_upper_1x1v_tensor_p2 },
  { NULL, NULL},
  { NULL, bc_sheath_gyrokinetic_reflectedf_upper_1x2v_tensor_p2 },
  { NULL, NULL},
//  { NULL, bc_sheath_gyrokinetic_reflectedf_upper_3x2v_tensor_p2 },
  { NULL, NULL},
};

// Primary struct in this updater.
struct gkyl_bc_sheath_gyrokinetic {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  struct gkyl_range skin_r, ghost_r;
  const struct gkyl_basis *basis;
  bool use_gpu;
  double q2Dm; // charge-to-mass ratio times 2.
  sheath_reflectedf_t ker_reflectedf;  // reflectedf kernel.
  const struct gkyl_rect_grid *grid;
  struct gkyl_range *conf_r;
};

GKYL_CU_D
static sheath_reflectedf_t
bc_gksheath_choose_reflectedf_kernel(const int dim, const int basis_type, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sheath_reflect_list[dim-2].kernels[poly_order-1];
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sheath_reflect_list[dim-2].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.

 * @param dir Direction in which to apply BC .
 * @param cdim Number of configuration space dimensions.
 * @param bctype Type of BC .
 * @param basis Basis in which to expand coefficients in array we apply BC to.
 * @param num_comp Number of components (DOFs) within a cell.
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods.
 */
struct gkyl_array_copy_func* gkyl_bc_sheath_gyrokinetic_create_arr_copy_func_cu(int dir, int cdim,
  enum gkyl_bc_sheath_gyrokinetic_type bctype, const struct gkyl_basis *basis, int num_comp);

#endif

GKYL_CU_D
static void
bc_gksheath_reflect(int dir, const struct gkyl_basis *basis, int cdim, double *out, const double *inp)
{
  basis->flip_odd_sign(dir, inp, out);
  basis->flip_odd_sign(cdim+1, out, out); // cdim+1 is the vpar direction.
}

///* Modeled after gkyl_array_flip_copy_to_buffer_fn */
//void
//bc_gksheath_advance(const struct gkyl_array *phi, const struct gkyl_array *phi_wall,
//  const struct gkyl_array *arr, int dir, const struct gkyl_basis *basis, int cdim, 
//  struct gkyl_range skin_r, struct gkyl_range ghost_r)
//{
//#ifdef GKYL_HAVE_CUDA
//  if (gkyl_array_is_cu_dev(arr)) {
//    gkyl_array_flip_copy_to_buffer_fn_cu(buff, arr, dir, skin_r, cf);
//    return;
//  }
//#endif
//
//  struct gkyl_range_iter iter;
//  gkyl_range_iter_init(&iter, &skin_r);
//
//  int fidx[GKYL_MAX_DIM]; // Flipped index.
//
//  int vdir = cdim; 
//  int uplo = skin_r.upper[vdir]+skin_r.lower[vdir];
//
//  while (gkyl_range_iter_next(&iter)) {
//
//    gkyl_copy_int_arr(skin_r.ndim, iter.idx, fidx);
//    fidx[vdir] = uplo - iter.idx[vdir];
//    // Turn this skin fidx into a ghost fidx.
//    fidx[dir] = ghost_r.lower[dir];
//    
//
//    long skin_loc = gkyl_range_idx(&skin_r, iter.idx);
//    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);
//
//    const double *inp = gkyl_array_cfetch(arr, skin_loc);
//    double *out = gkyl_array_fetch(arr, ghost_loc);
//
//    // Calculate reflected distribution function fhat.
//    // note: reflected distribution can be
//    // 1) fhat=0 (no reflection, i.e. absorb),
//    // 2) fhat=f (full reflection)
//    // 3) fhat=c*f (partial reflection)
//    double fhat[basis->num_basis] = {0.};
//
//    // Reflect fhat into skin cells.
//    gksheath_bc_reflect(dir, basis, cdim, out, fhat)
//  }
//}

