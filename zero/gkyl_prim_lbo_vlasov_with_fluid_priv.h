// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef void (*vlasov_with_fluid_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *fluid, const double *boundary_corrections);

typedef void (*vlasov_with_fluid_cross_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *greene, const double m_self, const double *moms_self, const double *u_self,
  const double *vtsq_self, const double m_other, const double *moms_other, const double *u_other,
  const double *vtsq_other, const double *fluid, const double *boundary_corrections);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vlasov_with_fluid_self_prim_t kernels[3]; } gkyl_prim_lbo_vlasov_with_fluid_self_kern_list;
typedef struct { vlasov_with_fluid_cross_prim_t kernels[3]; } gkyl_prim_lbo_vlasov_with_fluid_cross_kern_list;


//
// Serendipity basis kernels
//

// self-primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_vlasov_with_fluid_self_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, NULL, vlasov_with_fluid_self_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, NULL, vlasov_with_fluid_self_prim_moments_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, vlasov_with_fluid_self_prim_moments_2x2v_ser_p2 }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// cross-primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_vlasov_with_fluid_cross_kern_list ser_cross_prim_kernels[] = {
  // 1x kernels
  { NULL, NULL, vlasov_with_fluid_cross_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, NULL, vlasov_with_fluid_cross_prim_moments_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, vlasov_with_fluid_cross_prim_moments_2x2v_ser_p2 }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

struct prim_lbo_type_vlasov_with_fluid {
  struct gkyl_prim_lbo_type prim; // Base object
  vlasov_with_fluid_self_prim_t self_prim; // Self-primitive moments kernel
  vlasov_with_fluid_cross_prim_t cross_prim; // Cross-primitive moments kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_prim_lbo_vlasov_with_fluid_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_vlasov_with_fluid_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *moms, const double *boundary_corrections)
{
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid = container_of(prim, struct prim_lbo_type_vlasov_with_fluid, prim);

  long cidx = gkyl_range_idx(&prim_vlasov_with_fluid->conf_range, idx);
  return prim_vlasov_with_fluid->self_prim(A, rhs, 
    moms, (const double*) gkyl_array_cfetch(prim_vlasov_with_fluid->auxfields.fluid, cidx), boundary_corrections);
}

GKYL_CU_D
static void
cross_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *greene,
  const double m_self, const double *moms_self, const double *u_self, const double *vtsq_self,
  const double m_other, const double *moms_other, const double *u_other, const double *vtsq_other,
  const double *boundary_corrections)
{
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid = container_of(prim, struct prim_lbo_type_vlasov_with_fluid, prim);

  long cidx = gkyl_range_idx(&prim_vlasov_with_fluid->conf_range, idx);
  return prim_vlasov_with_fluid->cross_prim(A, rhs, 
    greene, m_self, moms_self, u_self, vtsq_self, m_other, moms_other, u_other, vtsq_other, 
    (const double*) gkyl_array_cfetch(prim_vlasov_with_fluid->auxfields.fluid, cidx), boundary_corrections);
}
