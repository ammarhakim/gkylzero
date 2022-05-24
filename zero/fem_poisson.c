#include <gkyl_fem_poisson.h>
//#include <gkyl_fem_poisson_priv.h>

struct gkyl_fem_poisson {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int num_cells[GKYL_MAX_DIM];
  bool isdirperiodic[GKYL_MAX_DIM]; // =true if parallel direction is periodic.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range solve_range, solve_range_ext;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_superlu_prob* prob;
  struct gkyl_mat *local_mass; // local mass matrix.
  struct gkyl_mat *local_stiff; // local stiffness matrix.
  struct gkyl_mat *local_mass_modtonod; // local mass matrix times modal-to-nodal matrix.
  struct gkyl_mat *local_nodtomod; // local nodal-to-modal matrix.

//  local2global_t l2g;  // Pointer to local-to-global kernel:

  long *globalidx;
};


static long
global_num_nodes(const int dim, const int poly_order, const int basis_type, const int *num_cells)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(num_cells);
//    } else if (poly_order == 2) {
//      return fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(num_cells);
//    } else if (poly_order == 3) {
//      return fem_poisson_num_nodes_global_1x_ser_p3_nonperiodicx(num_cells);
    }
  }
  assert(false);  // Other dimensionalities not supported.
  return -1;
}

gkyl_fem_poisson*
gkyl_fem_poisson_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const bool *isdirperiodic, void *ctx)
{

  gkyl_fem_poisson *up = gkyl_malloc(sizeof(gkyl_fem_poisson));

  up->ctx = ctx;
  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  for (int d=0; d<up->ndim; d++) up->isdirperiodic[d] = isdirperiodic[d];

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  // Local and local-ext ranges for whole-grid arrays.
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<up->ndim; d++) ghost[d] = 1;
  gkyl_create_grid_ranges(grid, ghost, &up->local_range_ext, &up->local_range);
  // Range of cells we'll solve Poisson in, as
  // a sub-range of up->local_range_ext.
  int sublower[GKYL_MAX_DIM], subupper[GKYL_MAX_DIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->local_range.lower[d];
    subupper[d] = up->local_range.lower[d];
    if (up->isdirperiodic[d]) {
      // Include ghost cells if periodic in this direction. Assume that
      // periodic directions are sync-ed outside of Poisson solver.
      sublower[d] = up->local_range_ext.lower[d];
      subupper[d] = up->local_range_ext.upper[d];
    }
  }
  gkyl_sub_range_init(&up->solve_range, &up->local_range_ext, sublower, subupper);
  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range.upper[d]-up->solve_range.lower[d]+1;

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = global_num_nodes(up->ndim, up->poly_order, basis->b_type, &up->num_cells[0]);

  return up;
}

void gkyl_fem_poisson_release(gkyl_fem_poisson *up)
{
//  gkyl_mat_release(up->local_mass);
//  gkyl_mat_release(up->local_stiff);
//  gkyl_mat_release(up->local_mass_modtonod);
//  gkyl_mat_release(up->local_nodtomod);
//  gkyl_superlu_prob_release(up->prob);
  gkyl_free(up->globalidx);
  gkyl_free(up);
}
