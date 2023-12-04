#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

static void
cubic_1d(void)
{
  double lower[] = { 0.0 }, upper[] = { 10.0 };
  int cells[] = { 16 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0] };
  int nc_cells[] = { cells[0] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, cells[0]+1);
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_1d(cells[0]);
  
  do {
    // initialize 1D nodal values
    for (int i=0; i<cells[0]+1; ++i) {
      double *pn = gkyl_array_fetch(psi_nodal, i);
      double xn = grid.lower[0] + i*grid.dx[0];
      pn[0] = -sq(xn) + 0.15*cub(xn);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_1d_a.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_1d_a.gkyl");
  } while (0);

  do {
    // initialize 1D nodal values
    for (int i=0; i<cells[0]+1; ++i) {
      double *pn = gkyl_array_fetch(psi_nodal, i);
      double xn = grid.lower[0] + i*grid.dx[0];
      pn[0] = sq(xn) + 10*sq(sin(xn));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_1d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_1d_b.gkyl");
  } while (0);  

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

static void
cubic_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 10.0, 10.0 };
  int cells[] = { 16, 16 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);
  
  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      double xn = grid.lower[0] + iter.idx[0]*grid.dx[0];
      double yn = grid.lower[1] + iter.idx[1]*grid.dx[1];
      pn[0] = (-sq(xn) + 0.15*cub(xn));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_2d_a.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_2d_a.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      double xn = grid.lower[0] + iter.idx[0]*grid.dx[0];
      double yn = grid.lower[1] + iter.idx[1]*grid.dx[1];
      pn[0] = (-sq(yn) + 0.15*cub(yn));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_2d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_2d_b.gkyl");
  } while (0);  

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

int
main(int argc, char **argv)
{
  cubic_1d();
  cubic_2d();
  
  return 0;
}
    
