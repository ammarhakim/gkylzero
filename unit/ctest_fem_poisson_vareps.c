// Test the FEM Poisson solver with spatially varying permittivity.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson_vareps.h>

void evalFunc2x_periodicx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double gxx =  3.0;
  double gxy = 10.0;
  double gyy = -2.0;

  double x = xn[0], y = xn[1];
  double amn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  double bmn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  fout[0] = 0.;
  for (int m=0; m<3; m++) {
    for (int n=0; n<3; n++) {
      double a = amn[m*3+n];
      double b = bmn[m*3+n];
      double t1 = -(a*gxx*pow(m,2) - 2*b*gxy*m*n + a*gyy*pow(n,2))*cos(m*x)*cos(n*y);
      double t2 = -(b*gxx*pow(m,2) - 2*a*gxy*m*n + b*gyy*pow(n,2))*sin(m*x)*sin(n*y);
      fout[0] = fout[0] + (t1+t2);
    }
  }
  fout[0] = fout[0]/50.;
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// Apply periodic BCs in one direction.
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, sgr.lower_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(buff->data, fld, sgr.upper_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.lower_ghost[dir]);
}

void
test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double gxx =  3.0;
  double gxy = 10.0;
  double gyy = -2.0;

  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_periodicy, NULL);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // create DG field for permittivity tensor.
  int epsnum = dim+ceil((pow(3.,dim-1)-dim)/2);
  struct gkyl_array *eps = mkarr(epsnum*basis.num_basis, localRange_ext.volume);
  // device copies:
  struct gkyl_array *rho_cu, *phi_cu, *eps_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    eps_cu = mkarr_cu(epsnum*basis.num_basis, localRange_ext.volume);
  }

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  if (use_gpu) gkyl_array_copy(rho_cu, rho);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_vareps_2x_rho_1.gkyl");

  // project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, gxx*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps, gxy*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, gyy*dg0norm, 2*basis.num_basis);
  if (use_gpu) gkyl_array_copy(eps_cu, eps);
//  gkyl_grid_sub_array_write(&grid, &localRange, eps, "ctest_fem_poisson_vareps_2x_eps_1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson_vareps *poisson = gkyl_fem_poisson_vareps_new(&grid, basis, &bcs, eps, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_vareps_set_rhs(poisson, rho_cu);
  else
    gkyl_fem_poisson_vareps_set_rhs(poisson, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_vareps_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_vareps_solve(poisson, phi);
  }
//  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_vareps_2x_phi_1.gkyl");

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double));
    gkyl_array_clear(sol_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(basis, 0, sol_cellavg, 0, phi, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, localRange);
    gkyl_array_shiftc(phi, mavgfac*sol_avg[0], 0);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  gkyl_fem_poisson_vareps_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
}

void test_2x_p1_periodicx_periodicy() {
  int cells[] = {32,32};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, &cells[0], bc_tv, false);
}

TEST_LIST = {
  // 2x tests
  { "test_2x_p1_periodicx_periodicy", test_2x_p1_periodicx_periodicy },
//  { "test_2x_p2_periodicx_periodicy", test_2x_p2_periodicx_periodicy },
#ifdef GKYL_HAVE_CUDA
//  // 1x tests
//  { "gpu_test_1x_p1_periodicx", gpu_test_1x_p1_periodicx },
//  { "gpu_test_1x_p1_dirichletx", gpu_test_1x_p1_dirichletx },
//  { "gpu_test_1x_p1_neumannx_dirichletx", gpu_test_1x_p1_neumannx_dirichletx },
//  { "gpu_test_1x_p1_dirichletx_neumannx", gpu_test_1x_p1_dirichletx_neumannx },
//  { "gpu_test_1x_p2_periodicx", gpu_test_1x_p2_periodicx },
//  { "gpu_test_1x_p2_dirichletx", gpu_test_1x_p2_dirichletx },
//  { "gpu_test_1x_p2_neumannx_dirichletx", gpu_test_1x_p2_neumannx_dirichletx },
//  { "gpu_test_1x_p2_dirichletx_neumannx", gpu_test_1x_p2_dirichletx_neumannx },
//  // 2x tests
//  { "gpu_test_2x_p1_periodicx_periodicy", gpu_test_2x_p1_periodicx_periodicy },
//  { "gpu_test_2x_p1_dirichletx_dirichlety", gpu_test_2x_p1_dirichletx_dirichlety },
//  { "gpu_test_2x_p1_dirichletx_periodicy", gpu_test_2x_p1_dirichletx_periodicy },
//  { "gpu_test_2x_p1_periodicx_dirichlety", gpu_test_2x_p1_periodicx_dirichlety },
//  { "gpu_test_2x_p1_dirichletx_neumanny_dirichlety", gpu_test_2x_p1_dirichletx_neumanny_dirichlety },
//  { "gpu_test_2x_p1_dirichletx_dirichlety_neumanny", gpu_test_2x_p1_dirichletx_dirichlety_neumanny },
//  { "gpu_test_2x_p1_neumannx_dirichletx_dirichlety", gpu_test_2x_p1_neumannx_dirichletx_dirichlety },
//  { "gpu_test_2x_p1_dirichletx_neumannx_dirichlety", gpu_test_2x_p1_dirichletx_neumannx_dirichlety },
//  { "gpu_test_2x_p2_periodicx_periodicy", gpu_test_2x_p2_periodicx_periodicy },
//  { "gpu_test_2x_p2_dirichletx_dirichlety", gpu_test_2x_p2_dirichletx_dirichlety },
//  { "gpu_test_2x_p2_dirichletx_periodicy", gpu_test_2x_p2_dirichletx_periodicy },
//  { "gpu_test_2x_p2_dirichletx_neumanny_dirichlety", gpu_test_2x_p2_dirichletx_neumanny_dirichlety },
//  { "gpu_test_2x_p2_dirichletx_dirichlety_neumanny", gpu_test_2x_p2_dirichletx_dirichlety_neumanny },
//  { "gpu_test_2x_p2_neumannx_dirichletx_dirichlety", gpu_test_2x_p2_neumannx_dirichletx_dirichlety },
//  { "gpu_test_2x_p2_dirichletx_neumannx_dirichlety", gpu_test_2x_p2_dirichletx_neumannx_dirichlety },
#endif
  { NULL, NULL },
};
