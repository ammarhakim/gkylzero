#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_dg_prim_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

struct dg_prim_vars_type_gyrokinetic {
  struct gkyl_dg_prim_vars_type pvt;
};

// for use in kernel tables
typedef struct {
  pvf_t kernels[3];
} gkyl_dg_prim_vars_gyrokinetic_kern_list;

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_upar_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_upar_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_upar_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_upar_3x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_upar_3x2v_ser_p1(in, out);
}

// gyrokinetic u_parallel kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_gyrokinetic_kern_list ser_dg_prim_vars_gyrokinetic_upar_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_upar_1x1v_ser_p1, kernel_dg_prim_vars_gyrokinetic_upar_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_gyrokinetic_upar_1x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_upar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_upar_2x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_upar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_upar_3x2v_ser_p1, NULL }, // 4
};

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_vtSq_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_vtSq_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_vtSq_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_vtSq_3x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_vtSq_3x2v_ser_p1(in, out);
}

// gyrokinetic vth^2 kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_gyrokinetic_kern_list ser_dg_prim_vars_gyrokinetic_vtSq_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_vtSq_1x1v_ser_p1, kernel_dg_prim_vars_gyrokinetic_vtSq_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_gyrokinetic_vtSq_1x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_vtSq_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_vtSq_2x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_vtSq_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_vtSq_3x2v_ser_p1, NULL }, // 4
};

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gyrokinetic_3x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return gyrokinetic_prim_vars_3x2v_ser_p1(in, out);
}

// gyrokinetic primitive variable (u_parallel, vth^2) kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_gyrokinetic_kern_list ser_dg_prim_vars_gyrokinetic_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_1x1v_ser_p1, kernel_dg_prim_vars_gyrokinetic_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_gyrokinetic_1x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_2x2v_ser_p1, kernel_dg_prim_vars_gyrokinetic_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_dg_prim_vars_gyrokinetic_3x2v_ser_p1, NULL }, // 4
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_dg_prim_vars_gyrokinetic_free(const struct gkyl_ref_count *ref);
