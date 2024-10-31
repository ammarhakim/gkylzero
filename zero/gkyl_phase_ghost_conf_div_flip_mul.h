#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_phase_ghost_conf_div_flip_mul gkyl_phase_ghost_conf_div_flip_mul;

/**
 * Create an updater which takes the phase-space quantity in the ghost cell,
 * divides it by the ghost cell of a conf-space quantity, and multiplies it
 * by the skin-cell value of the same ghost cell quantity but flipped in the
 * direction of the boundary (so that the conf-space quantity we multiply by
 * in the ghost cell has the same value at the boundary as it does in the skin
 * cell).
 *
 * @param conf_basis Configuration space basis object (on GPU if use_gpu=true).
 * @param phase_basis Phase space basis object.
 * @param use_gpu Whether to run on the GPU.
 * @return New phase_ghost_conf_div_flip_mul updater. 
 */
struct gkyl_phase_ghost_conf_div_flip_mul* gkyl_phase_ghost_conf_div_flip_mul_new(const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, bool use_gpu);

/**
 * Run the phase_ghost_conf_div_flip_mul updater.
 *
 * @param up Updater of type phase_ghost_conf_div_flip_mul. 
 * @param dir Direction perpendicular to the boundary.
 * @param edge Boundary edge (lower/upper).
 * @param conf_skin_r Configuration space skin range.
 * @param conf_ghost_r Configuration space ghost range.
 * @param phase_ghost_r Phase space ghost range.
 * @param jac Conf-space field (e.g. the jacobian).
 * @param jf Phase-space field (e.g. jacobian*distribution).
 */
void gkyl_phase_ghost_conf_div_flip_mul_advance(const struct gkyl_phase_ghost_conf_div_flip_mul *up,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf);

/**
 * Release the memory associated with phase_ghost_conf_div_flip_mul updater.
 *
 * @param up Updater of type phase_ghost_conf_div_flip_mul. 
 */
void gkyl_phase_ghost_conf_div_flip_mul_release(struct gkyl_phase_ghost_conf_div_flip_mul *up);
