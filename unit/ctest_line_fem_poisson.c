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

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_line_fem_poisson.h>
#include <gkyl_fem_parproj.h>


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = cos(3*z)*cos(x);
}

void
proj_func2(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 4*z*cos(2*x);
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
test_1(){
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
  gkyl_grid_sub_array_write(&grid, &local, field, "in_field.gkyl");


  struct gkyl_poisson_bc poisson_bc;
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(2.0), 0); 
  struct gkyl_array *phi= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
                                            
  struct gkyl_line_fem_poisson* line_fem_poisson = gkyl_line_fem_poisson_new(grid, &basis, basis, local, local_ext, epsilon, poisson_bc, false);
  gkyl_line_fem_poisson_advance(line_fem_poisson, field, phi);
  gkyl_line_fem_poisson_release(line_fem_poisson);
  gkyl_grid_sub_array_write(&grid, &local, phi, "out_field.gkyl");

}

void
test_2(){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/4, 1.0 };
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
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func2, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, field);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, field, "in_field.gkyl");


  struct gkyl_poisson_bc poisson_bc;
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(2.0), 0); 
  struct gkyl_array *phi= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
                                            
  struct gkyl_line_fem_poisson* line_fem_poisson = gkyl_line_fem_poisson_new(grid, &basis, basis, local, local_ext, epsilon, poisson_bc, false);
  gkyl_line_fem_poisson_advance(line_fem_poisson, field, phi);
  gkyl_line_fem_poisson_release(line_fem_poisson);
  gkyl_grid_sub_array_write(&grid, &local, phi, "out_field.gkyl");

}





void
test_3(){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/4, 1.0 };
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
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func2, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, field_discont);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, field_discont, "in_field.gkyl");

  // smooth output
  struct gkyl_array *field= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // smooth
  struct gkyl_array *weight=0;
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &local_ext, &basis, GKYL_FEM_PARPROJ_DIRICHLET, weight, false);
  gkyl_fem_parproj_set_rhs(parproj, field_discont, field_discont);
  gkyl_fem_parproj_solve(parproj, field);
  gkyl_grid_sub_array_write(&grid, &local, field, "smooth_field.gkyl");


  struct gkyl_poisson_bc poisson_bc;
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(2.0), 0); 
  struct gkyl_array *phi= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
                                            
  struct gkyl_line_fem_poisson* line_fem_poisson = gkyl_line_fem_poisson_new(grid, &basis, basis, local, local_ext, epsilon, poisson_bc, false);
  gkyl_line_fem_poisson_advance(line_fem_poisson, field, phi);
  gkyl_line_fem_poisson_release(line_fem_poisson);
  gkyl_grid_sub_array_write(&grid, &local, phi, "out_field.gkyl");

}



TEST_LIST = {
  //{ "test_1", test_1},
  { "test_2", test_2},
  //{ "test_3", test_3},
  { NULL, NULL },
};
