#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_pkpm_auxfields { 
  const struct gkyl_array *bvar; // magnetic field unit vector and unit tensor
};

/**
 * Create new Vlasov parallel-kinetic-perpendicular-moment (pkpm) moment type object.
 * Computes mass density, parallel pressure, and parallel heat flux (q b_hat, 3 components)
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Range for indexing magnetic field unit vector and unit tensor
 * @param mass Mass of species (pkpm moments are mass weighted)
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_vlasov_pkpm_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  double mass, bool use_gpu);

/**
 * Create new Vlasov moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, double mass);

/**
 * Set the auxiliary fields bvar needed in computing moments. bvar contains magnetic field 
 * unit vector (first 3 components) and magnetic field unit tensor (last 6 components)
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_pkpm_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields bvar needed in computing moments. bvar contains magnetic field 
 * unit vector (first 3 components) and magnetic field unit tensor (last 6 components)
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_pkpm_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_pkpm_auxfields auxin);

#endif
