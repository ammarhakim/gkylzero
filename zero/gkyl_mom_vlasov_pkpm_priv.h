#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_vlasov_pkpm {
  struct gkyl_mom_type momt;
  double mass; // mass of species
  struct gkyl_range conf_range; // velocity space range
  struct gkyl_mom_vlasov_pkpm_auxfields auxfields; // Auxiliary fields.
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_mom_vlasov_pkpm_kern_list;

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_1x1v_ser_p1(xc, dx, idx, mom_vlasov_pkpm->mass,
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_1x1v_ser_p2(xc, dx, idx, mom_vlasov_pkpm->mass, 
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_2x1v_ser_p1(xc, dx, idx, mom_vlasov_pkpm->mass, 
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_2x1v_ser_p2(xc, dx, idx, mom_vlasov_pkpm->mass, 
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_3x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_3x1v_ser_p1(xc, dx, idx, mom_vlasov_pkpm->mass, 
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_vlasov_pkpm_3x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = container_of(momt, struct mom_type_vlasov_pkpm, momt);
  long cidx = gkyl_range_idx(&mom_vlasov_pkpm->conf_range, idx);

  return mom_vlasov_pkpm_3x1v_ser_p2(xc, dx, idx, mom_vlasov_pkpm->mass, 
    (const double*) gkyl_array_cfetch(mom_vlasov_pkpm->auxfields.bvar, cidx), f, out);  
}

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_mom_vlasov_pkpm_kern_list ser_mom_vlasov_pkpm_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_mom_vlasov_pkpm_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, NULL, kernel_mom_vlasov_pkpm_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, kernel_mom_vlasov_pkpm_3x1v_ser_p2 }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_vlasov_pkpm_free(const struct gkyl_ref_count *ref);
