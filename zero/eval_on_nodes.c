#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>

struct gkyl_eval_on_nodes {
  struct gkyl_rect_grid grid;
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context
  
  int num_basis; // number of basis functions
  struct gkyl_array *nodes; // local nodal coordinates
};

gkyl_eval_on_nodes*
gkyl_eval_on_nodes_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_ret_vals, evalf_t eval, void *ctx)
{
  gkyl_eval_on_nodes *up = gkyl_malloc(sizeof(gkyl_eval_on_nodes));

  up->grid = *grid;
  up->num_ret_vals = num_ret_vals;
  up->eval = eval;
  up->ctx = ctx;
  up->num_basis = basis->num_basis;

  // initialize node local coordinates 
  up->nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, basis->num_basis);
  basis->node_list(gkyl_array_fetch(up->nodes, 0));

  return up;
}

static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

void
gkyl_eval_on_nodes_advance(const gkyl_eval_on_nodes *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int num_basis = up->num_basis;
  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, num_basis);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {
      comp_to_phys(up->grid.ndim, gkyl_array_cfetch(up->nodes, i),
        up->grid.dx, xc, xmu);
      up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }

    //long lidx = gkyl_range_idx(update_range, iter.idx);
    // TODO: transform to modal expansion

  }

  gkyl_array_release(fun_at_ords);
}

void
gkyl_eval_on_nodes_release(gkyl_eval_on_nodes* up)
{
  gkyl_array_release(up->nodes);
  gkyl_free(up);
}
