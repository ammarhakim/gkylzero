// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef void (*vlasov_pkpm_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *boundary_corrections);

// for use in kernel tables
typedef struct { vlasov_pkpm_self_prim_t kernels[3]; } gkyl_prim_lbo_vlasov_pkpm_self_kern_list;


//
// Serendipity basis kernels
//

// self-primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_vlasov_pkpm_self_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, vlasov_pkpm_self_prim_moments_1x1v_ser_p1, vlasov_pkpm_self_prim_moments_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, vlasov_pkpm_self_prim_moments_2x1v_ser_p1, vlasov_pkpm_self_prim_moments_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_self_prim_moments_3x1v_ser_p1, NULL }, // 2
};

struct prim_lbo_type_vlasov_pkpm {
  struct gkyl_prim_lbo_type prim; // Base object
  vlasov_pkpm_self_prim_t self_prim; // Self-primitive moments kernel
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_prim_lbo_vlasov_pkpm_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_vlasov_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *moms, const double *boundary_corrections)
{
  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm = container_of(prim, struct prim_lbo_type_vlasov_pkpm, prim);

  long cidx = gkyl_range_idx(&prim_vlasov_pkpm->conf_range, idx);
  return prim_vlasov_pkpm->self_prim(A, rhs, moms, boundary_corrections);
}
