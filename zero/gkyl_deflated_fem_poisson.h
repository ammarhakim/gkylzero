#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>


// Object type
typedef struct gkyl_deflated_fem_poisson gkyl_deflated_fem_poisson;

/**
 * Create new updater to solve a poisson equation
 *   - nabla . (epsilon * nabla phi) = rho
 *   on lines or planes using FEM to ensure phi is continuous.
 * The inputs are a DG field rho and epsilon and the output is
 * a DG field phi.
 *
 * rho and epsilon are evaluated at surfaces that are constant
 * in the last coordinate and then the FEM poisson solver
 * is called to solve the poisson equation at each surface.
 * The surfaces are lines (1d) if the initial fields are 2d
 * and planes (2d) if the initial fields are 3d).
 * Finally, the surfaces are assembled and converted into
 * a solution with the initial dimenionality.
 * Free using gkyl_deflated_fem_poisson_release method.
 *
 *
 * @param grid Grid object
 * @param basis_on_dev Basis functions of the DG field stored on the gpu.
 * @param basis Basis functions of the DG field.
 * @param local range on which the poisson problem should be solver
 * @param global_sub_range local range as a sub-range of the global range
 * @param bcs Boundary conditions.
 * @param epsilon Permittivity tensor. Defined over the extended range.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_deflated_fem_poisson* gkyl_deflated_fem_poisson_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, struct gkyl_range local, struct gkyl_range global_sub_range, struct gkyl_array *epsilon, struct gkyl_poisson_bc poisson_bc, bool use_gpu);

/**
 * Solve the poisson equation for the given charge density
 *
 * @param up deflated FEM poisson updater to run.
 * @param field DG field to set as RHS source (charge density rho).
 * @param phi DG field solution to poison problem (phi).
 */
void gkyl_deflated_fem_poisson_advance(struct gkyl_deflated_fem_poisson* up, struct gkyl_array *field, struct gkyl_array* phi, double target_corner_bias);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson* up);
