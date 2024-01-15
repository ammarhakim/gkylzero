#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Struct definition. Used to store nodal values 
struct gkyl_nodal_ops {
  struct gkyl_array *nodes;
};

/**
 * Create new nodal_ops struct for transforming between nodal and modal representations.
 *
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param use_gpu Boolean for whether nodes are stored on device
 * Returns pointer to gkyl_nodal_ops struct.
 */
struct gkyl_nodal_ops* gkyl_nodal_ops_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, bool use_gpu)

/**
 * Transform nodal representation to modal representation
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld Input nodal representation
 * @param modal_fld Output modal representation
 */
void gkyl_nodal_ops_n2m(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld);

/**
 * Transform modal representation to nodal representation
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld Output nodal representation
 * @param modal_fld Input modal representation
 */
void gkyl_nodal_ops_m2n(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld);

/**
 * Delete pointer to gkyl_nodal_ops struct.
 *
 * @param up Struct to delete.
 */
void gkyl_nodal_ops_release(struct gkyl_nodal_ops *up);

/**
 * Host-side wrappers for nodal-modal transformation operations on device
 */

void gkyl_nodal_ops_n2m_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld);

void gkyl_nodal_ops_m2n_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld);