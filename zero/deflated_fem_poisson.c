#include <assert.h>

#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_deflated_fem_poisson_priv.h>
#include <gkyl_array_average.h>

struct gkyl_deflated_fem_poisson* 
gkyl_deflated_fem_poisson_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, 
  struct gkyl_range local, struct gkyl_range global_sub_range, 
  struct gkyl_array *epsilon, struct gkyl_poisson_bc poisson_bc, bool use_gpu)
{
  struct gkyl_deflated_fem_poisson *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;
  up->grid = grid;
  up->basis = basis;
  up->basis_on_dev = basis_on_dev;

  // Assumes domain decomposition *only* in z.
  // local is range of potential; global_sub_range is range for indexing global charge density
  up->local = local;
  up->global_sub_range = global_sub_range;
  assert(up->local.volume == up->global_sub_range.volume);

  up->poisson_bc = poisson_bc;
  up->cdim = grid.ndim;
  up->num_solves_z = up->local.upper[up->cdim-1] - up->local.lower[up->cdim-1] + 2;
  up->d_fem_data = gkyl_malloc(sizeof(struct deflated_fem_data[up->num_solves_z]));

  int poly_order = up->basis.poly_order;

  // Create 2d/3d nodal range nodal array to be populated
  int nodes[GKYL_MAX_DIM];
  if (poly_order == 1) {
    for (int d=0; d<up->cdim; ++d)
      nodes[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->cdim; ++d)
      nodes[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }
  gkyl_range_init_from_shape(&up->nrange, up->cdim, nodes);

  // Create deflated 1d/2d grid, ranges, basis, and nodal range
  double deflated_lower[GKYL_MAX_DIM] = { 0.0 } ;
  double deflated_upper[GKYL_MAX_DIM] = { 0.0 };
  int deflated_cells[GKYL_MAX_DIM] = { 0 };
  for(int i = 0; i < up->cdim-1; i++){
    deflated_lower[i] = up->grid.lower[i];
    deflated_upper[i] = up->grid.upper[i];
    deflated_cells[i] = up->grid.cells[i];
  }

  gkyl_rect_grid_init(&up->deflated_grid, up->cdim-1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&up->deflated_grid, deflated_nghost, &up->deflated_local_ext, &up->deflated_local);
  gkyl_cart_modal_serendip(&up->deflated_basis, up->cdim-1, poly_order);

  if (up->use_gpu) {
    up->deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(up->deflated_basis_on_dev, up->cdim-1, poly_order);
    up->nodal_fld = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->cdim, up->nrange.volume);
  }
  else {
    up->deflated_basis_on_dev = &up->deflated_basis;
    up->nodal_fld = gkyl_array_new(GKYL_DOUBLE, up->cdim, up->nrange.volume);
  }

  int deflated_nodes[up->deflated_grid.ndim];
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
  up->n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, use_gpu);
  up->n2m_deflated = gkyl_nodal_ops_new(&up->deflated_basis, &up->deflated_grid, use_gpu);

  // Allocate necessary fields and solvers for each z slice
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1]+1; zidx++) {
    if (use_gpu) {
      up->d_fem_data[ctr].deflated_field = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_phi = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_epsilon = gkyl_array_cu_dev_new(GKYL_DOUBLE, (2*up->deflated_grid.ndim - 1)*up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_nodal_fld = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->deflated_grid.ndim, up->deflated_nrange.volume);
    }
    else {
      up->d_fem_data[ctr].deflated_field = gkyl_array_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_phi = gkyl_array_new(GKYL_DOUBLE, up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_epsilon = gkyl_array_new(GKYL_DOUBLE, (2*up->deflated_grid.ndim - 1)*up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_fem_data[ctr].deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, up->deflated_grid.ndim, up->deflated_nrange.volume);
    }
    if (zidx == up->local.upper[up->cdim-1] + 1 )
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx-1, &up->local, &up->deflated_local, epsilon, up->d_fem_data[ctr].deflated_epsilon,  2*up->deflated_grid.ndim-1);
    else 
      gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, &up->local, &up->deflated_local, epsilon, up->d_fem_data[ctr].deflated_epsilon, 2*up->deflated_grid.ndim-1);
    
    /// Here is the implementation of the target corner BC
    //---- check if we are at an extremal global index of z to apply target corner BC
    // get the global index of the z plane
    int global_zidx = zidx + up->global_sub_range.lower[up->cdim-1]; 
    // check if the global index is equal to the lower or upper limit of the z grid
    bool is_lower_z_edge = (zidx == up->local.lower[up->cdim-1])   && up->poisson_bc.contains_lower_z_edge;
    bool is_upper_z_edge = (zidx == up->local.upper[up->cdim-1]+1) && up->poisson_bc.contains_upper_z_edge;
    // store it in pisson_bc to pass it to gkyl_fem_poisson_new
    up->poisson_bc.is_z_edge = is_upper_z_edge || is_lower_z_edge;

    up->d_fem_data[ctr].fem_poisson = gkyl_fem_poisson_new(&up->deflated_local, &up->deflated_grid, up->deflated_basis, &up->poisson_bc, up->d_fem_data[ctr].deflated_epsilon, 0, false, use_gpu);
    ctr += 1;
  }

  return up;
}

void 
gkyl_deflated_fem_poisson_advance(struct gkyl_deflated_fem_poisson *up, struct gkyl_array *field, struct gkyl_array* phi, double target_corner_bias)
{
  int ctr = 0;
  int local_range_ctr = up->local.lower[up->cdim-1];
  for (int zidx = up->global_sub_range.lower[up->cdim-1]; zidx <= up->global_sub_range.upper[up->cdim-1]; zidx++) {
    // Deflate rho indexing global sub-range to fetch correct place in z
    gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, 
      &up->global_sub_range, &up->deflated_local, field, up->d_fem_data[ctr].deflated_field, 1);
    // Do the poisson solve 
    gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_field, target_corner_bias);
    gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
    // Modal to Nodal in 1d -> Store the result in the 2d nodal field
    gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
      &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
      up->nodal_fld, up->d_fem_data[ctr].deflated_phi, ctr);
    ctr += 1;
    local_range_ctr += 1;
    if (zidx == up->global_sub_range.upper[up->cdim-1]) {
      // Deflate rho indexing global sub-range to fetch correct place in z
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx, 
        &up->global_sub_range, &up->deflated_local, field, up->d_fem_data[ctr].deflated_field, 1);
      // Do the poisson solve 
      gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_field, target_corner_bias);
      gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
      // Modal to Nodal in 1d -> Store the result in the 2d nodal field
      gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
        &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
        up->nodal_fld, up->d_fem_data[ctr].deflated_phi, ctr);
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, phi);
}

void gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson* up){
  gkyl_array_release(up->nodal_fld);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_nodal_ops_release(up->n2m_deflated);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1] + 1; zidx++) {
    gkyl_array_release(up->d_fem_data[ctr].deflated_field);
    gkyl_array_release(up->d_fem_data[ctr].deflated_phi);
    gkyl_array_release(up->d_fem_data[ctr].deflated_epsilon);
    gkyl_array_release(up->d_fem_data[ctr].deflated_nodal_fld);
    gkyl_fem_poisson_release(up->d_fem_data[ctr].fem_poisson);
    ctr += 1;
  }
  if (up->use_gpu) {
    gkyl_cu_free(up->deflated_basis_on_dev);
  }

  gkyl_free(up->d_fem_data);
  gkyl_free(up);
}


