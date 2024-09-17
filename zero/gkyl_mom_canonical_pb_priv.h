#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_canonical_pb {
  struct gkyl_mom_type momt;
  double mass; // mass of species
  struct gkyl_range phase_range; // phase space range
  struct gkyl_mom_canonical_pb_auxfields auxfields; // Auxiliary fields.
};

// The cv_index[cd].vdim[cd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0, -1, -1}, // 1x kernel indices
  {-1, -1,  1, -1}, // 2x kernel indices
  {-1, -1, -1,  2}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_canonical_pb_mom_kern_list;

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_3x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_mom_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_mom_1x1v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_mom_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_mom_1x1v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_mom_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_mom_2x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_mom_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_mom_2x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_mom_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_mom_3x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

//
// Serendipity basis kernels
//

// MEnergy kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list ser_menergy_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_MEnergy_1x1v_ser_p1, kernel_canonical_pb_MEnergy_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_MEnergy_2x2v_ser_p1, kernel_canonical_pb_MEnergy_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_canonical_pb_MEnergy_3x3v_ser_p1, NULL                            }, // 2
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list ser_int_mom_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_int_mom_1x1v_ser_p1, kernel_canonical_pb_int_mom_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_int_mom_2x2v_ser_p1, kernel_canonical_pb_int_mom_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_canonical_pb_int_mom_3x3v_ser_p1, NULL                                 }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_can_pb_free(const struct gkyl_ref_count *ref);
