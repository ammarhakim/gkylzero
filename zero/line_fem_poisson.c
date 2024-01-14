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
#include <gkyl_line_fem_poisson.h>

void gkyl_line_fem_poisson_advance(struct gkyl_rect_grid grid, struct gkyl_basis basis, struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_array *epsilon, struct gkyl_array *field, struct gkyl_array* phi, struct gkyl_poisson_bc poisson_bc)
{

  int poly_order = basis.poly_order;
  // create deflated 1d grid, ranges, basis, and field
  // create xz grid
  double deflated_lower[1] = { grid.lower[0]}, deflated_upper[1] = { grid.upper[0]};
  int deflated_cells[1] = { grid.cells[0] };
  struct gkyl_rect_grid deflated_grid;
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);

  //ranges
  struct gkyl_range deflated_local, deflated_local_ext;
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);

  // basis function
  int deflated_poly_order = 1;
  struct gkyl_basis deflated_basis;
  gkyl_cart_modal_serendip(&deflated_basis, 1, deflated_poly_order);
  
  //field
  struct gkyl_array *deflated_field = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  struct gkyl_array *deflated_epsilon= gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  struct gkyl_array *deflated_phi = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);



  // create nrange and the 2d nodal array to be populated
  int nodes[2] = { 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = grid.cells[d] + 1;
  }
                   
  if (poly_order == 2){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = 2*(grid.cells[d]) + 1;
  }

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);


  //create the deflated nodal range and 1d array that will be used as an intermediate
  int deflated_nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = deflated_grid.cells[d] + 1;
  }
                   
  if (poly_order == 2){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(deflated_grid.cells[d]) + 1;
  }


  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);
  struct gkyl_array* deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, deflated_grid.ndim, deflated_nrange.volume);

  // Now for each z slice we want to deflate 
  // and then to a m2n to give us a nodal array at that z slice
  // Then fill the correct nodes in the 2d nodal array
  gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, &grid, &deflated_grid, 0, false);
  gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, &grid, &deflated_grid, 1, false);
  

  int nidx[2];
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    // first deflate epsilon
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, epsilon, deflated_epsilon, 1);
    // construct the fem poisson object
    struct gkyl_fem_poisson *fem_poisson = gkyl_fem_poisson_new(&deflated_local, &deflated_grid, deflated_basis, &poisson_bc, epsilon, 0, false, false);
    // then deflate deflate rho
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, field, deflated_field, 1);
    // do the poisson solve then free the solver
    gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field);
    gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
    gkyl_fem_poisson_release(fem_poisson);
    // then nodal to modal
    gkyl_nodal_ops_m2n(&deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi);
    // now loop through the 2d nrange and populate
    nidx[1] = zidx-1;
    for(int ix = 0; ix <=nrange.upper[0]; ix++){
      nidx[0] = ix;
      long lin_nidx = gkyl_range_idx(&nrange, nidx);
      long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
      const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
      double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
      output[0] = input[0];
    }
    if (zidx == local.upper[1]){
      // first deflate epsilon
      gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, epsilon, deflated_epsilon, 1);
      // construct the fem poisson object
      struct gkyl_fem_poisson *fem_poisson = gkyl_fem_poisson_new(&deflated_local, &deflated_grid, deflated_basis, &poisson_bc, epsilon, 0, false, false);
      // then deflate rho
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, field, deflated_field, 1);
      // do the poisson solve then free the solver
      gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field);
      gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
      gkyl_fem_poisson_release(fem_poisson);
      // then nodal to modal
      gkyl_nodal_ops_m2n(&deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi);
      // now loop through the 2d nrange and populate
      nidx[1] = zidx;
      for(int ix = 0; ix <= nrange.upper[0]; ix++){
        nidx[0] = ix;
        long lin_nidx = gkyl_range_idx(&nrange, nidx);
        long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
        const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
        double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
        output[0] = input[0];
      }
    }

  }
  gkyl_nodal_ops_n2m(&basis, &grid, &nrange, &local, 1, nodal_fld, phi);


  gkyl_array_release(deflated_field);
  gkyl_array_release(deflated_epsilon);
  gkyl_array_release(deflated_phi);
  gkyl_array_release(nodal_fld);
  gkyl_array_release(deflated_nodal_fld);
  gkyl_deflate_zsurf_release(deflator_lo);
  gkyl_deflate_zsurf_release(deflator_up);
}


