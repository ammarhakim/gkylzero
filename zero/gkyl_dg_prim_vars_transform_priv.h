#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_dg_prim_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_array.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

struct dg_prim_vars_type_transform {
  struct gkyl_dg_prim_vars_type pvt;
  struct gkyl_range conf_range; // Configuration space range (for indexing fields)
  struct gkyl_dg_prim_vars_auxfields auxfields; // Auxiliary fields.
};

// for use in kernel tables
typedef struct {
  pvf_t kernels[3];
} gkyl_dg_prim_vars_transform_kern_list;

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_1x1v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_1x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_1x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_1x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_2x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_2x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_2x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_2x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_i_3x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_i_3x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}


// Transformation kernel for u_par -> u_par b_i for use of GK flow in Vlasov
GKYL_CU_D
static const gkyl_dg_prim_vars_transform_kern_list ser_dg_prim_vars_transform_u_par_i_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_i_1x1v_ser_p1, kernel_dg_prim_vars_transform_u_par_i_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_transform_u_par_i_1x2v_ser_p1, kernel_dg_prim_vars_transform_u_par_i_1x2v_ser_p2 }, // 1
  { NULL, kernel_dg_prim_vars_transform_u_par_i_1x3v_ser_p1, kernel_dg_prim_vars_transform_u_par_i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_i_2x2v_ser_p1, kernel_dg_prim_vars_transform_u_par_i_2x2v_ser_p2 }, // 3
  { NULL, kernel_dg_prim_vars_transform_u_par_i_2x3v_ser_p1, kernel_dg_prim_vars_transform_u_par_i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_i_3x3v_ser_p1, NULL }, // 5
};

// u_par
GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_1x1v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_1x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_2x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_transform_u_par_3x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_u_par_3x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

// Vlasov u_i kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_transform_kern_list ser_dg_prim_vars_transform_u_par_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_1x1v_ser_p1, kernel_dg_prim_vars_transform_u_par_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_transform_u_par_1x2v_ser_p1, kernel_dg_prim_vars_transform_u_par_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_2x2v_ser_p1, kernel_dg_prim_vars_transform_u_par_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_dg_prim_vars_transform_u_par_3x2v_ser_p1, NULL }, // 3
};


// prim vars GK

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_gk_1x1v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_gk_1x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_gk_2x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}


GKYL_CU_DH
static void
kernel_dg_prim_vars_gk_transform_3x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_gk_3x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

// Transformation kernel for u_i -> upar = u_i b_i and vth^2 -> vth_GK^2 = 1/vdim(M2/M0 - upar^2) 
GKYL_CU_D
static const gkyl_dg_prim_vars_transform_kern_list ser_dg_prim_vars_transform_gk_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_gk_transform_1x1v_ser_p1, kernel_dg_prim_vars_gk_transform_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_gk_transform_1x2v_ser_p1, kernel_dg_prim_vars_gk_transform_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_dg_prim_vars_gk_transform_2x2v_ser_p1, kernel_dg_prim_vars_gk_transform_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_dg_prim_vars_gk_transform_3x2v_ser_p1, NULL }, // 3
};

// prim vars Vlasov
GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_1x1v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_1x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_1x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_1x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_2x2v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_2x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_2x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_2x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_transform_3x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);

  long cidx = gkyl_range_idx(&prim->conf_range, idx);
  return transform_prim_vars_vlasov_3x3v_ser_p1((const double*) gkyl_array_cfetch(prim->auxfields.b_i, cidx), in, out);
}

// Transformation kernel for upar -> upar_i = upar*b_i and vth^2 
GKYL_CU_D
static const gkyl_dg_prim_vars_transform_kern_list ser_dg_prim_vars_transform_vlasov_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_vlasov_transform_1x1v_ser_p1, kernel_dg_prim_vars_vlasov_transform_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_vlasov_transform_1x2v_ser_p1, kernel_dg_prim_vars_vlasov_transform_1x2v_ser_p2 }, // 1
  { NULL, kernel_dg_prim_vars_vlasov_transform_1x3v_ser_p1, kernel_dg_prim_vars_vlasov_transform_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_dg_prim_vars_vlasov_transform_2x2v_ser_p1, kernel_dg_prim_vars_vlasov_transform_2x2v_ser_p2 }, // 3
  { NULL, kernel_dg_prim_vars_vlasov_transform_2x3v_ser_p1, kernel_dg_prim_vars_vlasov_transform_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_dg_prim_vars_vlasov_transform_3x3v_ser_p1, NULL }, // 5
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_dg_prim_vars_transform_free(const struct gkyl_ref_count *ref);
