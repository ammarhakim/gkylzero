#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_kernels.h>

struct gkyl_fem_parproj {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  int poly_order; // polynomial order of the basis.
  int pardir; // parallel (z) direction.
  int parnum_cells; // number of cells in parallel (z) direction.

  struct gkyl_range perp_range; // range of perpendicular cells.

  int numnodes_local, numnodes_global;
};

long
global_num_nodes(const int dim, const int poly_order, const int basis_type, const int parnum_cells)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_num_nodes_global_1x_ser_p1(parnum_cells);
    } else (poly_order == 2) {
      return fem_parproj_num_nodes_global_1x_ser_p2(parnum_cells);
    } else (poly_order == 3) {
      return fem_parproj_num_nodes_global_1x_ser_p3(parnum_cells);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_num_nodes_global_2x_ser_p1(parnum_cells);
      } else (poly_order == 2) {
        return fem_parproj_num_nodes_global_2x_ser_p2(parnum_cells);
      } else (poly_order == 3) {
        return fem_parproj_num_nodes_global_2x_ser_p3(parnum_cells);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_num_nodes_global_2x_tensor_p1(parnum_cells);
      } else (poly_order == 2) {
        return fem_parproj_num_nodes_global_2x_tensor_p2(parnum_cells);
      }
    }
  }
}


gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  void *ctx)
{
  gkyl_fem_parproj *up = gkyl_malloc(sizeof(gkyl_fem_parproj));

  up->ctx = ctx;
  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis = basis->num_basis;
  up->poly_order = basis->poly_order;
  up->pardir = grid->ndim; // Assume parallel direction is always the last.
  up->parnum_cells = grid->cells[up->pardir-1];

  // Create a range of perpendicular cells.
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, grid->ndim, grid->cells);
  gkyl_range_shorten(&up->perp_range, &range, up->pardir, 1);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = global_num_nodes(up->ndim, up->poly_order, basis->b_type, up->parnum_cells);

  return up;
}

//void gkyl_fem_parproj_begin_assembly(const gkyl_fem_proj *femproj,
//  double tm, const struct gkyl_array *src);

void gkyl_fem_parproj_release(gkyl_fem_parproj *up)
{
  gkyl_free(up);
}
