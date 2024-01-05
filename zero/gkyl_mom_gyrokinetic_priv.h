#pragma once

// Private header, not for direct use in user code

#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_gyrokinetic_kernels.h>


struct mom_type_gyrokinetic {
  struct gkyl_mom_type momt;
  double mass; // mass of species
  struct gkyl_range conf_range; // configuration space range
  const struct gkyl_array *bmag; // Pointer to magnitude of magnetic field
};

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_gyrokinetic_mom_kern_list;

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M0_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M0_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M1_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M1_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_par_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_par_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M2_perp_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M2_perp_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_par_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_par_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_M3_perp_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_M3_perp_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_1x1v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_1x1v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_1x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_1x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_2x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_2x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_3x2v_ser_p1(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);
}

GKYL_CU_DH
static void
kernel_gyrokinetic_three_moments_3x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  return gyrokinetic_three_moments_3x2v_ser_p2(xc, dx, idx, mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->bmag, cidx), f, out);
}

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_M0_1x1v_ser_p1, kernel_gyrokinetic_M0_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_M0_1x2v_ser_p1, kernel_gyrokinetic_M0_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M0_2x2v_ser_p1, kernel_gyrokinetic_M0_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M0_3x2v_ser_p1, kernel_gyrokinetic_M0_3x2v_ser_p2 }, // 3
};

// M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m1_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_M1_1x1v_ser_p1, kernel_gyrokinetic_M1_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_M1_1x2v_ser_p1, kernel_gyrokinetic_M1_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M1_2x2v_ser_p1, kernel_gyrokinetic_M1_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M1_3x2v_ser_p1, kernel_gyrokinetic_M1_3x2v_ser_p2 }, // 3
};

// M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_M2_1x1v_ser_p1, kernel_gyrokinetic_M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_M2_1x2v_ser_p1, kernel_gyrokinetic_M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M2_2x2v_ser_p1, kernel_gyrokinetic_M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M2_3x2v_ser_p1, kernel_gyrokinetic_M2_3x2v_ser_p2 }, // 3
};

// M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_par_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_M2_par_1x1v_ser_p1, kernel_gyrokinetic_M2_par_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_M2_par_1x2v_ser_p1, kernel_gyrokinetic_M2_par_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M2_par_2x2v_ser_p1, kernel_gyrokinetic_M2_par_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M2_par_3x2v_ser_p1, kernel_gyrokinetic_M2_par_3x2v_ser_p2 }, // 3
};

// M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_gyrokinetic_M2_perp_1x2v_ser_p1, kernel_gyrokinetic_M2_perp_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M2_perp_2x2v_ser_p1, kernel_gyrokinetic_M2_perp_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M2_perp_3x2v_ser_p1, kernel_gyrokinetic_M2_perp_3x2v_ser_p2 }, // 3
};

// M3 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3_par_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_M3_par_1x1v_ser_p1, kernel_gyrokinetic_M3_par_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_M3_par_1x2v_ser_p1, kernel_gyrokinetic_M3_par_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M3_par_2x2v_ser_p1, kernel_gyrokinetic_M3_par_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M3_par_3x2v_ser_p1, kernel_gyrokinetic_M3_par_3x2v_ser_p2 }, // 3
};

// M3 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3_perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_gyrokinetic_M3_perp_1x2v_ser_p1, kernel_gyrokinetic_M3_perp_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_M3_perp_2x2v_ser_p1, kernel_gyrokinetic_M3_perp_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_M3_perp_3x2v_ser_p1, kernel_gyrokinetic_M3_perp_3x2v_ser_p2 }, // 3
};

// Zeroth (density), First (parallel momentum), 
// and Second (total energy) computed together
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_three_moments_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_three_moments_1x1v_ser_p1, kernel_gyrokinetic_three_moments_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_three_moments_1x2v_ser_p1, kernel_gyrokinetic_three_moments_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_three_moments_2x2v_ser_p1, kernel_gyrokinetic_three_moments_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_three_moments_3x2v_ser_p1, kernel_gyrokinetic_three_moments_3x2v_ser_p2 }, // 3
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_gk_mom_free(const struct gkyl_ref_count *ref);
