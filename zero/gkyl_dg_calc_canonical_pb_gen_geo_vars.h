#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_canonical_pb_gen_geo_vars gkyl_dg_calc_canonical_pb_gen_geo_vars;

/**
 * Create new updater to compute canonical_pb variables needed in 
 * updates for general geometry. Methods compute:
 * alpha_surf : Surface expansion of phase space flux alpha for streaming in general geometry (v^i = v . e^i)
 * cot_vec : Volume expansion of contangent vectors for volume term of streaming in general geometry
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_canonical_pb_gen_geo_vars* 
gkyl_dg_calc_canonical_pb_gen_geo_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,  bool use_gpu);

/**
 * Create new updater to compute canonical_pb general geometry variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_canonical_pb_gen_geo_vars* 
gkyl_dg_calc_canonical_pb_gen_geo_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis);

/**
 * Compute surface expansion of phase space flux alpha
 * Evaluates e^i = g^ij e_j at surfaces before computing v^i = v . e^i to insure continuity
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing general geometry canonical_pb variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param phase_ext_range Extended Phase space range (so we obtain alpha_surf at all the needed surfaces)
 * @param alpha_surf Output surface expansion in a cell on the *lower* edge in each direction 
 * @param sgn_alpha_surf Output sign(alpha) at quadrature points along a surface 
 * @param const_sgn_alpha Output boolean array for if sign(alpha) is a constant on the surface
 *                        If sign(alpha) is a constant, kernels are simpler and we exploit this fact.
 */
void gkyl_dg_calc_canonical_pb_gen_geo_vars_alpha_surf(struct gkyl_dg_calc_canonical_pb_gen_geo_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_ext_range, 
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha);

/**
 * Compute volume expansion of cotangent vectors e^i = g^ij e_j
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing general geometry canonical_pb variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param cot_vec Output volume expansion of cotangent vectors. 
 */
void gkyl_dg_calc_canonical_pb_gen_geo_vars_cot_vec(struct gkyl_dg_calc_canonical_pb_gen_geo_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* cot_vec);

/**
 * Delete pointer to updater to compute canonical_pb general geometry variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_canonical_pb_gen_geo_vars_release(struct gkyl_dg_calc_canonical_pb_gen_geo_vars *up);

/**
 * Host-side wrappers for canonical_pb general geometry variable operations on device
 */
void gkyl_dg_calc_canonical_pb_gen_geo_vars_alpha_surf_cu(struct gkyl_dg_calc_canonical_pb_gen_geo_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_ext_range, 
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha);

void gkyl_dg_calc_canonical_pb_gen_geo_vars_cot_vec_cu(struct gkyl_dg_calc_canonical_pb_gen_geo_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* cot_vec);