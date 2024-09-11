#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_fpo_vlasov_kernels.h>

#include <assert.h>

struct mom_type_bcorr_fpo_vlasov {
  struct gkyl_mom_type momt;
  double vBoundary[6];
  struct gkyl_range phase_range; // velocity space range
  struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxfields; // Auxiliary fields.
};

// for use in kernel tables
typedef struct {
  momf_t kernels[4];
} gkyl_mom_bcorr_fpo_vlasov_kern_list;

// Fill pdim phase index given edge information and idx array of remaining directions
GKYL_CU_DH
static void
set_phase_idx(enum gkyl_vel_edge edge, const struct gkyl_range *phase_range, int cdim, int pdim, const int *idx, int fidx[GKYL_MAX_DIM])
{
  int vel_idx = 0;
  int cell = 1;

  switch (edge) {
    case GKYL_VX_LOWER:
      vel_idx = 1;
      break;
    case GKYL_VX_UPPER:
      vel_idx = 1;
      cell = phase_range->upper[vel_idx];
      break;
    case GKYL_VY_LOWER:
      vel_idx = 2;
      break;
    case GKYL_VY_UPPER:
      vel_idx = 2;
      cell = phase_range->upper[vel_idx];
      break;
    case GKYL_VZ_LOWER:
      vel_idx = 3;
      break;
    case GKYL_VZ_UPPER:
      vel_idx = 3;
      cell = phase_range->upper[vel_idx];
      break;
    default:
      assert(false);
      break;
  }
  
  for (int i=0; i<(cdim+vel_idx); ++i) fidx[i] = idx[i];
  fidx[vel_idx] = cell;
  for (int i=vel_idx+1; i<pdim; ++i) fidx[i] = idx[i-1];
}


GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  set_phase_idx(edge, &mom_fpo_vlasov->phase_range, 
    mom_fpo_vlasov->momt.cdim, mom_fpo_vlasov->momt.pdim, idx, fidx);

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  return mom_bcorr_fpo_vlasov_1x3v_ser_p1(xc, idx, edge, mom_fpo_vlasov->vBoundary, dx, 
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  set_phase_idx(edge, &mom_fpo_vlasov->phase_range, 
    mom_fpo_vlasov->momt.cdim, mom_fpo_vlasov->momt.pdim, idx, fidx);

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  return mom_bcorr_fpo_vlasov_1x3v_ser_p2(xc, idx, edge, mom_fpo_vlasov->vBoundary, dx, 
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  set_phase_idx(edge, &mom_fpo_vlasov->phase_range, 
    mom_fpo_vlasov->momt.cdim, mom_fpo_vlasov->momt.pdim, idx, fidx);

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  // return mom_bcorr_fpo_vlasov_2x3v_ser_p1(xc, idx, edge, mom_fpo_vlasov->vBoundary, dx, 
  //   (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
  //   f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  set_phase_idx(edge, &mom_fpo_vlasov->phase_range, 
    mom_fpo_vlasov->momt.cdim, mom_fpo_vlasov->momt.pdim, idx, fidx);

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  // return mom_bcorr_fpo_vlasov_2x3v_ser_p2(xc, idx, edge, mom_fpo_vlasov->vBoundary, dx, 
  //   (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
  //   f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  set_phase_idx(edge, &mom_fpo_vlasov->phase_range, 
    mom_fpo_vlasov->momt.cdim, mom_fpo_vlasov->momt.pdim, idx, fidx);

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  // return mom_bcorr_fpo_vlasov_3x3v_ser_p1(xc, idx, edge, mom_fpo_vlasov->vBoundary, dx, 
  //   (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
  //   f, out);
}

//
// Serendipity basis kernels
//

// FPO boundary corrections
GKYL_CU_D
static const gkyl_mom_bcorr_fpo_vlasov_kern_list ser_mom_bcorr_fpo_vlasov_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p1, kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p1, kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_mom_bcorr_fpo_vlasov_3x3v_ser_p1, NULL }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_bcorr_fpo_vlasov_free(const struct gkyl_ref_count *ref);
