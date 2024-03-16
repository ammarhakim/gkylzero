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
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

struct dg_prim_vars_type_vlasov {
  struct gkyl_dg_prim_vars_type pvt;
};

// for use in kernel tables
typedef struct {
  pvf_t kernels[3];
} gkyl_dg_prim_vars_vlasov_kern_list;

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_1x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_1x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_2x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_2x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_2x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_u_i_3x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_u_i_3x3v_ser_p1(in, out);
}

// Vlasov u_i kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_vlasov_kern_list ser_dg_prim_vars_vlasov_u_i_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_vlasov_u_i_1x1v_ser_p1, kernel_dg_prim_vars_vlasov_u_i_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_vlasov_u_i_1x2v_ser_p1, kernel_dg_prim_vars_vlasov_u_i_1x2v_ser_p2 }, // 1
  { NULL, kernel_dg_prim_vars_vlasov_u_i_1x3v_ser_p1, kernel_dg_prim_vars_vlasov_u_i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_dg_prim_vars_vlasov_u_i_2x2v_ser_p1, kernel_dg_prim_vars_vlasov_u_i_2x2v_ser_p2 }, // 3
  { NULL, kernel_dg_prim_vars_vlasov_u_i_2x3v_ser_p1, kernel_dg_prim_vars_vlasov_u_i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_dg_prim_vars_vlasov_u_i_3x3v_ser_p1, NULL }, // 5
};

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_1x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_1x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_2x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_2x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_2x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_vtSq_3x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_vtSq_3x3v_ser_p1(in, out);
}

// Vlasov vth^2 kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_vlasov_kern_list ser_dg_prim_vars_vlasov_vtSq_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_1x1v_ser_p1, kernel_dg_prim_vars_vlasov_vtSq_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_1x2v_ser_p1, kernel_dg_prim_vars_vlasov_vtSq_1x2v_ser_p2 }, // 1
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_1x3v_ser_p1, kernel_dg_prim_vars_vlasov_vtSq_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_2x2v_ser_p1, kernel_dg_prim_vars_vlasov_vtSq_2x2v_ser_p2 }, // 3
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_2x3v_ser_p1, kernel_dg_prim_vars_vlasov_vtSq_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_dg_prim_vars_vlasov_vtSq_3x3v_ser_p1, NULL }, // 5
};


GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x1v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_1x1v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x1v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_1x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_1x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_1x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_2x2v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_2x2v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_2x2v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_2x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_2x3v_ser_p1(in, out);
}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_2x3v_ser_p2(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{

}

GKYL_CU_DH
static void
kernel_dg_prim_vars_vlasov_3x3v_ser_p1(const struct gkyl_dg_prim_vars_type *pvt, 
  const int *idx, const double *in, double* out)
{
  return vlasov_prim_vars_3x3v_ser_p1(in, out);
}

// Vlasov primitive variable kernel list
GKYL_CU_D
static const gkyl_dg_prim_vars_vlasov_kern_list ser_dg_prim_vars_vlasov_kernels[] = {
  // 1x kernels
  { NULL, kernel_dg_prim_vars_vlasov_1x1v_ser_p1, kernel_dg_prim_vars_vlasov_1x1v_ser_p2 }, // 0
  { NULL, kernel_dg_prim_vars_vlasov_1x2v_ser_p1, kernel_dg_prim_vars_vlasov_1x2v_ser_p2 }, // 1
  { NULL, kernel_dg_prim_vars_vlasov_1x3v_ser_p1, kernel_dg_prim_vars_vlasov_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_dg_prim_vars_vlasov_2x2v_ser_p1, kernel_dg_prim_vars_vlasov_2x2v_ser_p2 }, // 3
  { NULL, kernel_dg_prim_vars_vlasov_2x3v_ser_p1, kernel_dg_prim_vars_vlasov_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_dg_prim_vars_vlasov_3x3v_ser_p1, NULL }, // 5
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_dg_prim_vars_vlasov_free(const struct gkyl_ref_count *ref);
