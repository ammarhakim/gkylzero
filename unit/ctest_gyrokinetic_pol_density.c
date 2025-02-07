#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_dg_basis_ops.h>

#include <stdio.h>
#include <stdlib.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}


void evalFunc1x_1(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void evalFunc1x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = -pow(xn[0]-0.5, 2) + 1.0;
}

void evalFunc2x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0);
}

void evalFunc3x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0) * (-pow(xn[2]-0.5, 2) + 1.0);
}

void test_1x_flat( bool use_gpu ) 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 0, 0 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0]};
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0]};
  int nc_cells[] = { cells[0] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, ghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis_phi;
  gkyl_cart_modal_tensor(&basis_phi, 1, 3);

  struct gkyl_array *phi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1));
  struct gkyl_array *phi_cubic = gkyl_array_new(GKYL_DOUBLE, basis_phi.num_basis, localRange_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);
  double xn[1];
  // initialize 2D nodal values
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(phi_nodal, nidx);
    evalFunc1x_1(0.0, xn, pn, NULL);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
    phi_nodal, phi_cubic);

  // Print out the array to make sure it makes sense
  const char *fmt = "%s-%s.gkyl";
  char name[32] = "ctest_polarization";
  int sz = gkyl_calc_strlen(fmt, name, "phi_pol_projection");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "phi_pol_cubic");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_cubic, fileNm);
  sprintf(fileNm, fmt, name, "phi_pol_nm");
  gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, phi_nodal, fileNm);
    

  struct gkyl_array *phi_pol = use_gpu? mkarr(use_gpu, basis_phi.num_basis, localRange_ext.volume)
                                          : gkyl_array_acquire(phi_cubic);
  gkyl_array_set_range_to_range(phi_pol, 1.0, phi_cubic, &localRange, &localRange);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( npol_d[i] == 0.0 );
    }
  }
  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_array_release(phi_cubic);
  gkyl_array_release(phi_nodal);
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_1x_flat_cpu() 
{
  test_1x_flat(false);
}

void test_1x_flat_gpu() 
{
  test_1x_flat(true);
}

void test_1x_quad( bool use_gpu ) 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 0, 0 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);
  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0]};
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0]};
  int nc_cells[] = { cells[0] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, ghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis_phi;
  gkyl_cart_modal_tensor(&basis_phi, 1, 3);

  struct gkyl_array *phi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1));
  struct gkyl_array *phi_cubic = gkyl_array_new(GKYL_DOUBLE, basis_phi.num_basis, localRange_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);
  double xn[1];
  // initialize 2D nodal values
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(phi_nodal, nidx);
    evalFunc1x_quad(0.0, xn, pn, NULL);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
    phi_nodal, phi_cubic);

  // Print out the array to make sure it makes sense
  const char *fmt = "%s-%s.gkyl";
  char name[32] = "ctest_polarization";
  int sz = gkyl_calc_strlen(fmt, name, "phi_pol_projection");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "phi_pol_cubic");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_cubic, fileNm);
  sprintf(fileNm, fmt, name, "phi_pol_nm");
  gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, phi_nodal, fileNm);
    

  struct gkyl_array *phi_pol = use_gpu? mkarr(use_gpu, basis_phi.num_basis, localRange_ext.volume)
                                          : gkyl_array_acquire(phi_cubic);
  gkyl_array_set_range_to_range(phi_pol, 1.0, phi_cubic, &localRange, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    TEST_CHECK( gkyl_compare(npol_d[0], 2.8284271247462, 1e-14) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-12) );
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_array_release(phi_cubic);
  gkyl_array_release(phi_nodal);
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_1x_quad_cpu() 
{
  test_1x_quad(false);
}

void test_1x_quad_gpu() 
{
  test_1x_quad(true);
}

void
test_2x_quad1x( bool use_gpu ) 
{
  int cells[] = {7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 0, 0 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1};
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, dim, nc_lower, nc_upper, nc_cells);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, ghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis_phi;
  gkyl_cart_modal_tensor(&basis_phi, dim, 3);

  struct gkyl_array *phi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *phi_cubic = gkyl_array_new(GKYL_DOUBLE, basis_phi.num_basis, localRange_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);
  double xn[2];
  // initialize 2D nodal values
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(phi_nodal, nidx);
    evalFunc1x_quad(0.0, xn, pn, NULL);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
    phi_nodal, phi_cubic);

  // Print out the array to make sure it makes sense
  const char *fmt = "%s-%s.gkyl";
  char name[32] = "ctest_polarization";
  int sz = gkyl_calc_strlen(fmt, name, "phi_pol_projection");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "phi_pol_cubic");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_cubic, fileNm);
  sprintf(fileNm, fmt, name, "phi_pol_nm");
  gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, phi_nodal, fileNm);
    

  struct gkyl_array *phi_pol = use_gpu? mkarr(use_gpu, basis_phi.num_basis, localRange_ext.volume)
                                          : gkyl_array_acquire(phi_cubic);
  gkyl_array_set_range_to_range(phi_pol, 1.0, phi_cubic, &localRange, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  sprintf(fileNm, fmt, name, "npol");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, npol_ho, fileNm);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-12;
    TEST_CHECK( gkyl_compare(npol_d[0], 4.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_cubic);
  gkyl_array_release(phi_nodal);
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_2x_quad1x_cpu() 
{
  test_2x_quad1x(false);
}

void test_2x_quad1x_gpu() 
{
  test_2x_quad1x(true);
}

void
test_2x_quad( bool use_gpu ) 
{
  int cells[] = {7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int ghost[] = { 0, 0 };
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  struct gkyl_array *npol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = mkarr(false, basis.num_basis, localRange_ext.volume);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);

  struct gkyl_array *npol = use_gpu? mkarr(use_gpu, basis.num_basis, localRange_ext.volume)
                                      : gkyl_array_acquire(npol_ho);
  gkyl_array_copy(npol, npol_ho);
  struct gkyl_array *epsilon = use_gpu? mkarr(use_gpu, basis.num_basis, localRange_ext.volume)
                                      : gkyl_array_acquire(epsilon_ho);
  gkyl_array_copy(epsilon, epsilon_ho);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1};
  int nc_ghost[] = { 0, 0 };

  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, dim, nc_lower, nc_upper, nc_cells);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nc_ghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis_phi;
  int basis_phi_poly_order = 3;
  gkyl_cart_modal_tensor(&basis_phi, dim, basis_phi_poly_order);

  struct gkyl_array *phi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nc_cells[0]*nc_cells[1]);
  struct gkyl_array *phi_cubic = gkyl_array_new(GKYL_DOUBLE, basis_phi.num_basis, localRange_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);
  double xn[2];
  // initialize 2D nodal values
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(phi_nodal, nidx);
    evalFunc2x_quad(0.0, xn, pn, NULL);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
    phi_nodal, phi_cubic);

  // Print out the array to make sure it makes sense
  const char *fmt = "%s-%s.gkyl";
  char name[32] = "ctest_polarization";
  int sz = gkyl_calc_strlen(fmt, name, "phi_pol_projection");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "phi_pol_cubic");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_cubic, fileNm);
  sprintf(fileNm, fmt, name, "phi_pol_nm");
  gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, phi_nodal, fileNm);
    

  struct gkyl_array *phi_pol = use_gpu? mkarr(use_gpu, basis_phi.num_basis, localRange_ext.volume)
                                          : gkyl_array_acquire(phi_cubic);
  gkyl_array_set_range_to_range(phi_pol, 1.0, phi_cubic, &localRange, &localRange);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  sprintf(fileNm, fmt, name, "npol");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, npol_ho, fileNm);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-10;
    // Ignore the corners
    if (conf_iter.idx[0] == 1){
      if(conf_iter.idx[1] == 1){
        continue;
      } else if (conf_iter.idx[1] == 7){
        continue;
      }
    } else if (conf_iter.idx[0] == 7){
      if(conf_iter.idx[1] == 1){
        continue;
      } else if (conf_iter.idx[1] == 7){
        continue;
      }
    }
    if (conf_iter.idx[1] == 1) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.1413919026586975, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 2) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0942612684391317, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 3) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0471306342195649, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 4) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9931972789115648, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 5) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0471306342195649, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 6) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0942612684391317, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 7) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.1413919026586975, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    }
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_cubic);
  gkyl_array_release(phi_nodal);
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_2x_quad_cpu() 
{
  test_2x_quad(false);
}

void test_2x_quad_gpu() 
{
  test_2x_quad(true);
}

void
test_3x_flat( bool use_gpu ) 
{
  int cells[] = {7, 7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, dim, poly_order+1);

  struct gkyl_array *phi_pol = mkarr(use_gpu, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_pol_ho = use_gpu? mkarr(false, phi_pol->ncomp, phi_pol->size)
                                         : gkyl_array_acquire(phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, time, &localRange, phi_pol_ho);
  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_copy(phi_pol, phi_pol_ho);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-12;
    TEST_CHECK( gkyl_compare(npol_d[0], 5.6568542494927714, tol) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[4], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[5], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[6], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[7], 0.0, tol) );
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol_ho);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void
test_3x_flat_cpu() 
{
  test_3x_flat(false);
}

void
test_3x_flat_gpu() 
{
  test_3x_flat(true);
}

>>>>>>> gk-g0-app:unit/ctest_gyrokinetic_pol_density.c
TEST_LIST = {
  { "test_1x_flat_cpu", test_1x_flat_cpu },
  { "test_1x_quad_cpu", test_1x_quad_cpu },
  { "test_2x_quad1x_cpu", test_2x_quad1x_cpu },
  { "test_2x_quad_cpu", test_2x_quad_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_flat_gpu", test_1x_flat_gpu },
  { "test_1x_quad_gpu", test_1x_quad_gpu },
  { "test_2x_quad1x_gpu", test_2x_quad1x_gpu },
  { "test_2x_quad_gpu", test_2x_quad_gpu },
#endif
  { NULL, NULL },
};
