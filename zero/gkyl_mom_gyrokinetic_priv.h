#pragma once

// Private header, not for direct use in user code

#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_gyrokinetic_kernels.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

typedef void (*gyrokinetic_momf_t)(const double *w, const double *dxv, 
  const int *idx, const double m_, const double *Bmag, 
  const double *f, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { gyrokinetic_momf_t kernels[3]; } gkyl_gyrokinetic_mom_kern_list;

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M0_1x1v_ser_p1, gyrokinetic_M0_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_M0_1x2v_ser_p1, gyrokinetic_M0_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M0_2x2v_ser_p1, gyrokinetic_M0_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M0_3x2v_ser_p1, NULL              }, // 3
};

// M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m1_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M1_1x1v_ser_p1, gyrokinetic_M1_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_M1_1x2v_ser_p1, gyrokinetic_M1_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M1_2x2v_ser_p1, gyrokinetic_M1_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M1_3x2v_ser_p1, NULL              }, // 3
};

// M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M2_1x1v_ser_p1, gyrokinetic_M2_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_M2_1x2v_ser_p1, gyrokinetic_M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_2x2v_ser_p1, gyrokinetic_M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_3x2v_ser_p1, NULL              }, // 3
};

// M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_par_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M2_par_1x1v_ser_p1, gyrokinetic_M2_par_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_M2_par_1x2v_ser_p1, gyrokinetic_M2_par_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_par_2x2v_ser_p1, gyrokinetic_M2_par_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_par_3x2v_ser_p1, NULL                  }, // 3
};

// M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, gyrokinetic_M2_perp_1x2v_ser_p1, gyrokinetic_M2_perp_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_perp_2x2v_ser_p1, gyrokinetic_M2_perp_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_perp_3x2v_ser_p1, NULL                   }, // 3
};

// M3 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3_par_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M3_par_1x1v_ser_p1, gyrokinetic_M3_par_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_M3_par_1x2v_ser_p1, gyrokinetic_M3_par_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M3_par_2x2v_ser_p1, gyrokinetic_M3_par_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M3_par_3x2v_ser_p1, NULL                  }, // 3
};

// M3 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, gyrokinetic_M3_perp_1x2v_ser_p1, gyrokinetic_M3_perp_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M3_perp_2x2v_ser_p1, gyrokinetic_M3_perp_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M3_perp_3x2v_ser_p1, NULL                   }, // 3
};

//
// Tensor-product basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m0_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M0_1x1v_ser_p1, gyrokinetic_M0_1x1v_tensor_p2 }, // 0
  { NULL, gyrokinetic_M0_1x2v_ser_p1, gyrokinetic_M0_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M0_2x2v_ser_p1, gyrokinetic_M0_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M0_3x2v_ser_p1, NULL                 }, // 3
};

// M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m1_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M1_1x1v_ser_p1, gyrokinetic_M1_1x1v_tensor_p2 }, // 0
  { NULL, gyrokinetic_M1_1x2v_ser_p1, gyrokinetic_M1_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M1_2x2v_ser_p1, gyrokinetic_M1_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M1_3x2v_ser_p1, NULL                 }, // 3
};

// M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m2_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M2_1x1v_ser_p1, gyrokinetic_M2_1x1v_tensor_p2 }, // 0
  { NULL, gyrokinetic_M2_1x2v_ser_p1, gyrokinetic_M2_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_2x2v_ser_p1, gyrokinetic_M2_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_3x2v_ser_p1, NULL                 }, // 3
};

// M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m2_par_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M2_par_1x1v_ser_p1, gyrokinetic_M2_par_1x1v_tensor_p2 }, // 0
  { NULL, gyrokinetic_M2_par_1x2v_ser_p1, gyrokinetic_M2_par_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_par_2x2v_ser_p1, gyrokinetic_M2_par_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_par_3x2v_ser_p1, NULL                     }, // 3
};

// M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m2_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, gyrokinetic_M2_perp_1x2v_ser_p1, gyrokinetic_M2_perp_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M2_perp_2x2v_ser_p1, gyrokinetic_M2_perp_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M2_perp_3x2v_ser_p1, NULL                      }, // 3
};

// M3 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m3_par_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_M3_par_1x1v_ser_p1, gyrokinetic_M3_par_1x1v_tensor_p2 }, // 0
  { NULL, gyrokinetic_M3_par_1x2v_ser_p1, gyrokinetic_M3_par_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M3_par_2x2v_ser_p1, gyrokinetic_M3_par_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M3_par_3x2v_ser_p1, NULL                     }, // 3
};

// M3 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ten_m3_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, gyrokinetic_M3_perp_1x2v_ser_p1, gyrokinetic_M3_perp_1x2v_tensor_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_M3_perp_2x2v_ser_p1, gyrokinetic_M3_perp_2x2v_tensor_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_M3_perp_3x2v_ser_p1, NULL                      }, // 3
};

struct mom_type_gyrokinetic {
  struct gkyl_mom_type momt;
  gyrokinetic_momf_t kernel; // moment calculation kernel
  double _m; // mass of species
  const double *Bmag; // pointer to magnitude of magnetic field
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_gk_mom_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gyrokinetic = container_of(momt, struct mom_type_gyrokinetic, momt);

  return mom_gyrokinetic->kernel(xc, dx, idx, mom_gyrokinetic->_m, mom_gyrokinetic->Bmag, f, out);
}
