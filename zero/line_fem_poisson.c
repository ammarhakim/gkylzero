#include <gkyl_line_fem_poisson.h>
#include <gkyl_line_fem_poisson_priv.h>


struct gkyl_line_fem_poisson* 
gkyl_line_fem_poisson_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, struct gkyl_range local, struct gkyl_range local_ext, 
  struct gkyl_array *epsilon, struct gkyl_poisson_bc poisson_bc, bool use_gpu)
{
  struct gkyl_line_fem_poisson *up = gkyl_malloc(sizeof(*up));
  up->grid = grid;
  up->basis = basis;
  up->basis_on_dev = basis_on_dev;
  up->local = local;
  up->local_ext = local_ext;
  up->poisson_bc = poisson_bc;
  up->num_solves_z = up->local.upper[1] - up->local.lower[1] + 2;
  up->d_fem_data = gkyl_malloc(sizeof(struct deflated_fem_data[up->num_solves_z]));

  int poly_order = up->basis.poly_order;

  // Create 2d nodal range nodal array to be populated
  int nodes[2];
  if (poly_order == 1){
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = up->grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = 2*(up->grid.cells[d]) + 1;
  }
  gkyl_range_init_from_shape(&up->nrange, up->grid.ndim, nodes);

  // Create deflated 1d grid, ranges, basis, and nodal range
  double deflated_lower[1] = { up->grid.lower[0]}, deflated_upper[1] = { up->grid.upper[0]};
  int deflated_cells[1] = { up->grid.cells[0] };
  gkyl_rect_grid_init(&up->deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&up->deflated_grid, deflated_nghost, &up->deflated_local_ext, &up->deflated_local);
  gkyl_cart_modal_serendip(&up->deflated_basis, 1, poly_order);

  if (up->use_gpu) {
    up->deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(up->deflated_basis_on_dev, up->deflated_grid.ndim, poly_order);
    up->nodal_fld = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->grid.ndim, up->nrange.volume);
  }
  else {
    up->deflated_basis_on_dev = &up->deflated_basis;
    up->nodal_fld = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, up->nrange.volume);
  }

  int deflated_nodes[1];
  if (poly_order == 1){
    for (int d=0; d<up->deflated_grid.ndim; ++d)
      deflated_nodes[d] = up->deflated_grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<up->deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(up->deflated_grid.cells[d]) + 1;
  }
  gkyl_range_init_from_shape(&up->deflated_nrange, up->deflated_grid.ndim, deflated_nodes);

  // Allocate a deflation object for the lower and upper edge
  up->deflator_lo = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 0, use_gpu);
  up->deflator_up = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 1, use_gpu);

  // Allocate the 1d and 2d nodal_ops updater to be used in the advance method
  up->n2m_2d = gkyl_nodal_ops_new(&up->basis, &up->grid, use_gpu);
  up->n2m_1d = gkyl_nodal_ops_new(&up->deflated_basis, &up->deflated_grid, use_gpu);

  // Allocate necessary fields and solvers for each z slice
  int ctr = 0;
  for (int zidx = up->local.lower[1]; zidx <= up->local.upper[1]+1; zidx++) {
    if (use_gpu) {
      up->d_fem_data[ctr].deflated_field = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_phi = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_epsilon = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_nodal_fld = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_grid.ndim, up->deflated_nrange.volume);
    }
    else {
      up->d_fem_data[ctr].deflated_field = gkyl_array_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_phi = gkyl_array_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_epsilon = gkyl_array_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, up->deflated_grid.ndim, up->deflated_nrange.volume);
    }
    if (zidx == up->local.upper[1] + 1 )
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx-1, &up->local, &up->deflated_local, epsilon, up->d_fem_data[ctr].deflated_epsilon, 1);
    else 
      gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, &up->local, &up->deflated_local, epsilon, up->d_fem_data[ctr].deflated_epsilon, 1);
    up->d_fem_data[ctr].fem_poisson = gkyl_fem_poisson_new(&up->deflated_local, &up->deflated_grid, up->deflated_basis, &up->poisson_bc, up->d_fem_data[ctr].deflated_epsilon, 0, false, use_gpu);
    ctr += 1;
  }

  return up;
}

void 
gkyl_line_fem_poisson_advance(struct gkyl_line_fem_poisson *up, struct gkyl_array *field, struct gkyl_array* phi)
{
  int ctr = 0;
  for(int zidx = up->local.lower[1]; zidx <= up->local.upper[1]; zidx++){
    // Deflate rho
    gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, &up->local, &up->deflated_local, field, up->d_fem_data[ctr].deflated_field, 1);
    // Do the poisson solve 
    gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_field);
    gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
    // Modal to Nodal in 1d -> Store the result in the 2d nodal field
    gkyl_nodal_ops_m2n_deflated(up->n2m_1d, up->deflated_basis_on_dev, &up->deflated_grid, &up->nrange, &up->deflated_local, 1, up->nodal_fld, up->d_fem_data[ctr].deflated_phi, zidx-1);
    ctr += 1;
    if (zidx == up->local.upper[1]) {
      // Deflate rho
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx, &up->local, &up->deflated_local, field, up->d_fem_data[ctr].deflated_field, 1);
      // Do the poisson solve 
      gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_field);
      gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
      // Modal to Nodal in 1d -> Store the result in the 2d nodal field
      gkyl_nodal_ops_m2n_deflated(up->n2m_1d, up->deflated_basis_on_dev, &up->deflated_grid, &up->nrange, &up->deflated_local, 1, up->nodal_fld, up->d_fem_data[ctr].deflated_phi, zidx);
    }
  }
  gkyl_nodal_ops_n2m(up->n2m_2d, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, phi);



}

void gkyl_line_fem_poisson_release(struct gkyl_line_fem_poisson* up){
  gkyl_array_release(up->nodal_fld);
  gkyl_nodal_ops_release(up->n2m_2d);
  gkyl_nodal_ops_release(up->n2m_1d);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  int ctr = 0;
  for (int zidx = up->local.lower[1]; zidx <= up->local.upper[1]; zidx++) {
    gkyl_array_release(up->d_fem_data[ctr].deflated_field);
    gkyl_array_release(up->d_fem_data[ctr].deflated_phi);
    gkyl_array_release(up->d_fem_data[ctr].deflated_epsilon);
    gkyl_array_release(up->d_fem_data[ctr].deflated_nodal_fld);
    gkyl_fem_poisson_release(up->d_fem_data[ctr].fem_poisson);
    ctr += 1;
  }
  if (up->use_gpu) {
    gkyl_free(up->deflated_basis_on_dev);
  }

  gkyl_free(up->d_fem_data);
  gkyl_free(up);
}


