#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*cos(x);
}

void
proj_func2(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x);
}

void evalFunc1x_neumannx_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 5.0;
  double c0 = 0.;
  double c1 = a/12. - 1./2.;
  fout[0] = -(1.-a*pow(x,2));
}


void
test_deflate(){
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // project initial function
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, 0, field, "xzfunc.gkyl");

  // create deflated grid, ranges, basis, and field
  // create xz grid
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
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

  // now deflate
  int edge = 1; //lower = 0
  gkyl_deflate_zsurf *deflator = gkyl_deflate_zsurf_new(&basis, &deflated_basis, edge, false);

  int zidx = 1;
  gkyl_deflate_zsurf_advance(deflator, zidx, &local, &deflated_local, field, deflated_field, 1);

  gkyl_grid_sub_array_write(&deflated_grid, &deflated_local, 0,
    deflated_field, "xzfunc_deflated.gkyl");
}


void
test_reconstruct(){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // project initial function on 2d field
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, 0, field, "in_field.gkyl");




  // create deflated 1d grid, ranges, basis, and field
  // create xz grid
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
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

  printf("\n");
  for(int d=0; d<grid.ndim; d++){
    printf("\nd[%d] = %d\n", d, nodes[d]);
  }

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  printf("nrange lower = %d %d\n", nrange.lower[0], nrange.lower[1]);
  printf("nrange upper = %d %d\n", nrange.upper[0], nrange.upper[1]);
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

  for(int d=0; d<grid.ndim; d++){
    printf("deflated d[%d] = %d\n", d, nodes[d]);
  }

  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);
  printf("deflated_nrange lower = %d\n", deflated_nrange.lower[0]);
  printf("deflated_nrange upper = %d\n", deflated_nrange.upper[0]);
  struct gkyl_array* deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, deflated_grid.ndim, deflated_nrange.volume);

  printf("nrange.volume = %ld\n", nrange.volume);
  // Now for each z slice we want to deflate 
  // and then to a m2n to give us a nodal array at that z slice
  // Then fill the correct nodes in the 2d nodal array
  gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, false);
  gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, false);
  

  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&deflated_basis, &deflated_grid, false);
  int nidx[2];
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    // first deflate
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, field, deflated_field, 1);
    // then nodal to modal
    gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_field);
    // now loop through the 2d nrange and populate
    nidx[1] = zidx-1;
    printf("\niz = %d\n", zidx);
    for(int ix = 0; ix <=nrange.upper[0]; ix++){
      nidx[0] = ix;
      printf("nidx = %d, %d\n", nidx[0], nidx[1]);
      long lin_nidx = gkyl_range_idx(&nrange, nidx);
      long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
      printf("ix, lin_nidx_deflated, lin_nidx %d, %ld, %ld\n", ix, lin_nidx_deflated, lin_nidx);
      const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
      double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
      output[0] = input[0];
    }
    if (zidx == local.upper[1]){
      // first deflate
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, field, deflated_field, 1);
      // then nodal to modal
      gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_field);
      // now loop through the 2d nrange and populate
      nidx[1] = zidx;
      for(int ix = 0; ix <= nrange.upper[0]; ix++){
        nidx[0] = ix;
        printf("nidx = %d, %d\n", nidx[0], nidx[1]);
        long lin_nidx = gkyl_range_idx(&nrange, nidx);
        long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
        printf("doig upper ix, lin_nidx_deflated, lin_nidx %d, %ld, %ld\n", ix, lin_nidx_deflated, lin_nidx);
        const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
        double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
        output[0] = input[0];
      }
    }

  }


  struct gkyl_array *out_field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, out_field);
  gkyl_grid_sub_array_write(&grid, &local, 0, out_field, "out_field.gkyl");

}

void
test_poisson_slices(){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/2, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // project initial function on 2d field
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, field);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, 0, field, "in_field.gkyl");




  // create deflated 1d grid, ranges, basis, and field
  // create xz grid
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
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

  printf("\n");
  for(int d=0; d<grid.ndim; d++){
    printf("\nd[%d] = %d\n", d, nodes[d]);
  }

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  printf("nrange lower = %d %d\n", nrange.lower[0], nrange.lower[1]);
  printf("nrange upper = %d %d\n", nrange.upper[0], nrange.upper[1]);
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

  for(int d=0; d<grid.ndim; d++){
    printf("deflated d[%d] = %d\n", d, nodes[d]);
  }

  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);
  printf("deflated_nrange lower = %d\n", deflated_nrange.lower[0]);
  printf("deflated_nrange upper = %d\n", deflated_nrange.upper[0]);
  struct gkyl_array* deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, deflated_grid.ndim, deflated_nrange.volume);

  printf("nrange.volume = %ld\n", nrange.volume);
  // Now for each z slice we want to deflate 
  // and then to a m2n to give us a nodal array at that z slice
  // Then fill the correct nodes in the 2d nodal array
  gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, false);
  gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, false);
  
  struct gkyl_poisson_bc poisson_bc;
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(2.0), 0); // Sets weight=1.
  struct gkyl_fem_poisson *fem_poisson = gkyl_fem_poisson_new(&deflated_local, &deflated_grid,
    deflated_basis, &poisson_bc, NULL, epsilon, 0, false, false);


  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&basis, &grid, false);
  int nidx[2];
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    // first deflate
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, field, deflated_field, 1);
    // do the poisson solve
    gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field, NULL);
    gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
    // then nodal to modal
    gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi);
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
      // first deflate
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, field, deflated_field, 1);
      // do the poisson solve
      gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field, NULL);
      gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
      // then nodal to modal
      gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi);
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


  struct gkyl_array *out_field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, out_field);
  gkyl_nodal_ops_release(n2m);
  gkyl_grid_sub_array_write(&grid, &local, 0, out_field, "out_field.gkyl");

}



TEST_LIST = {
  //{ "test_deflate", test_deflate},
  //{ "test_reconstruct", test_reconstruct},
  { "test_poisson_slices", test_poisson_slices},
  { NULL, NULL },
};
