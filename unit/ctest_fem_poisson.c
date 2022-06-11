// Test the FEM Poisson solver.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson.h>

void evalFunc1x_periodicx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x);
}
void evalFunc1x_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 2.0;
  double c0 = a/12. - 1./2.;
  double c1 = 0.;
  fout[0] = -(1.-a*pow(x,2));
}

void evalFunc2x_periodicx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x)*cos((2.*2.*M_PI/(2.*M_PI))*y);
}
void evalFunc2x_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 2.*M_PI;
  double kx = 2.*M_PI/Lx, ky = 2.*M_PI/Ly;
  double sig = 0.3*sqrt(Lx*Lx+Ly*Ly);
  fout[0] = exp(-(pow(kx*x,2)+pow(ky*y,2))/(2.*(sig*sig)));
}
void evalFunc2x_dirichletx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = (y/(2.*M_PI)+0.5)*sin(5.0*M_PI*(y/(2.*M_PI)+0.5))+exp(-(pow(x/(2.*M_PI)+0.5-0.5,2)+pow(y/(2.*M_PI)+0.5-0.5,2))/0.02);
}
void evalFunc2x_periodicx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = (x/(2.*M_PI)+0.5)*sin(5.0*M_PI*(x/(2.*M_PI)+0.5))+exp(-(pow(x/(2.*M_PI)+0.5-0.5,2)+pow(y/(2.*M_PI)+0.5-0.5,2))/0.02);
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate array (filled with zeros)
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
test_1x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs)
{
  double epsilon_0 = 1.0;
  double lower[] = {0.0}, upper[] = {1.0};
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    lower[0] = -M_PI;
    upper[0] =  M_PI;
  }
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
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_periodicx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx, NULL);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_1x_rho_1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&grid, basis, bcs, epsilon_0, NULL, false);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho);

  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_1x_phi_1.gkyl");

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double)); 
    gkyl_array_clear(sol_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(basis, 0, sol_cellavg, 0, phi, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, localRange);
    gkyl_array_shiftc0(phi, mavgfac*sol_avg[0]);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
	-2.7060255565796842e-01, -1.5623245835252897e-01,
	-7.7061088089325591e-01, -1.3244748281911814e-01,
	-1.1533006851175758e+00, -8.8498578665918093e-02,
	-1.3604109147294894e+00, -3.1076568152445378e-02,
	-1.3604109147294894e+00, 3.1076568152445333e-02,
	-1.1533006851175760e+00, 8.8498578665918051e-02,
	-7.7061088089325636e-01, 1.3244748281911811e-01,
	-2.7060255565796909e-01, 1.5623245835252891e-01,
	2.7060255565796831e-01, 1.5623245835252891e-01,
	7.7061088089325547e-01, 1.3244748281911814e-01,
	1.1533006851175758e+00, 8.8498578665918204e-02,
	1.3604109147294898e+00, 3.1076568152445395e-02,
	1.3604109147294898e+00, -3.1076568152445395e-02,
	1.1533006851175758e+00, -8.8498578665918023e-02,
	7.7061088089325658e-01, -1.3244748281911806e-01,
	2.7060255565796898e-01, -1.5623245835252900e-01,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
        -0.0133521216082515, -0.0077088510047765,
        -0.0373194046782545, -0.0061286663274337,
        -0.0558775050145713, -0.0045858578973869,
        -0.0691990561087026, -0.0031053432128332,
        -0.0775430081978993, -0.0017120397719699,
        -0.0812546282651626, -0.0004308650729942,
        -0.0807655000392442,  0.0007132633858965,
        -0.0765935239946458,  0.0016954281065051,
        -0.0693429173516197,  0.0024907115906342,
        -0.0597042140761683,  0.0030741963400866,
        -0.0484542648800444,  0.0034209648566651,
        -0.0364562372207511,  0.0035060996421724,
        -0.024659615301542 ,  0.0033046831984112,
        -0.0141002000714206,  0.0027917980271844,
        -0.0059001092251411,  0.0019425266302945,
        -0.0012677772032077,  0.0007319515095444,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
      }
    } else {
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
	-5.2739372946099083e-01, -2.8867496150229432e-01, 1.2250982323291082e-02, 
	-1.2732410943752506e+00, -1.1957308417178132e-01, 2.9576487677282437e-02, 
	-1.2732410943752503e+00, 1.1957308417178146e-01, 2.9576487677282326e-02, 
	-5.2739372946099017e-01, 2.8867496150229438e-01, 1.2250982323291082e-02, 
	5.2739372946099117e-01, 2.8867496150229432e-01, -1.2250982323291082e-02, 
	1.2732410943752508e+00, 1.1957308417178118e-01, -2.9576487677282381e-02, 
	1.2732410943752503e+00, -1.1957308417178146e-01, -2.9576487677282326e-02, 
	5.2739372946099039e-01, -2.8867496150229427e-01, -1.2250982323291082e-02,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
        -2.5791443630192343e-02, -1.3837517332209948e-02,  8.1578940289174257e-04,
        -6.2965188799992161e-02, -7.6912011102198800e-03,  7.6432003993326847e-04,
        -7.9768181972718816e-02, -2.1429048449640056e-03,  6.6138131401632689e-04,
        -7.8962559012382361e-02,  2.4086914924015437e-03,  5.0697322514091727e-04,
        -6.4691523714997728e-02,  5.5649079307207141e-03,  3.0109577330704072e-04,
        -4.2479347808584755e-02,  6.9270644988373735e-03,  4.3748958514695155e-05,
        -1.9231370953168247e-02,  6.0964812255954216e-03, -2.6506721923612130e-04,
        -3.2340007407779619e-03,  2.6744781398387418e-03, -6.2535275994540806e-04,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
      }
    } else {
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(perbuff);
}

void
gpu_test_1x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs)
{
  double epsilon_0 = 1.0;
  double lower[] = {0.0}, upper[] = {1.0};
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    lower[0] = -M_PI;
    upper[0] =  M_PI;
  }
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
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_periodicx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx, NULL);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_1x_rho_1.gkyl");
  gkyl_array_copy(rho_cu, rho);

  // FEM poisson solver.
  gkyl_fem_poisson *poisson_cu = gkyl_fem_poisson_new(&grid, basis, bcs, epsilon_0, NULL, true);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson_cu, rho_cu);

  // Solve the problem.
  gkyl_fem_poisson_solve(poisson_cu, phi_cu);
  gkyl_array_copy(phi, phi_cu);

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif

  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_1x_phi_1.gkyl");

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double)); 
    gkyl_array_clear(sol_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(basis, 0, sol_cellavg, 0, phi, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, localRange);
    gkyl_array_shiftc0(phi, mavgfac*sol_avg[0]);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
	-2.7060255565796842e-01, -1.5623245835252897e-01,
	-7.7061088089325591e-01, -1.3244748281911814e-01,
	-1.1533006851175758e+00, -8.8498578665918093e-02,
	-1.3604109147294894e+00, -3.1076568152445378e-02,
	-1.3604109147294894e+00, 3.1076568152445333e-02,
	-1.1533006851175760e+00, 8.8498578665918051e-02,
	-7.7061088089325636e-01, 1.3244748281911811e-01,
	-2.7060255565796909e-01, 1.5623245835252891e-01,
	2.7060255565796831e-01, 1.5623245835252891e-01,
	7.7061088089325547e-01, 1.3244748281911814e-01,
	1.1533006851175758e+00, 8.8498578665918204e-02,
	1.3604109147294898e+00, 3.1076568152445395e-02,
	1.3604109147294898e+00, -3.1076568152445395e-02,
	1.1533006851175758e+00, -8.8498578665918023e-02,
	7.7061088089325658e-01, -1.3244748281911806e-01,
	2.7060255565796898e-01, -1.5623245835252900e-01,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
        -0.0133521216082515, -0.0077088510047765,
        -0.0373194046782545, -0.0061286663274337,
        -0.0558775050145713, -0.0045858578973869,
        -0.0691990561087026, -0.0031053432128332,
        -0.0775430081978993, -0.0017120397719699,
        -0.0812546282651626, -0.0004308650729942,
        -0.0807655000392442,  0.0007132633858965,
        -0.0765935239946458,  0.0016954281065051,
        -0.0693429173516197,  0.0024907115906342,
        -0.0597042140761683,  0.0030741963400866,
        -0.0484542648800444,  0.0034209648566651,
        -0.0364562372207511,  0.0035060996421724,
        -0.024659615301542 ,  0.0033046831984112,
        -0.0141002000714206,  0.0027917980271844,
        -0.0059001092251411,  0.0019425266302945,
        -0.0012677772032077,  0.0007319515095444,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else {
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
	-5.2739372946099083e-01, -2.8867496150229432e-01, 1.2250982323291082e-02, 
	-1.2732410943752506e+00, -1.1957308417178132e-01, 2.9576487677282437e-02, 
	-1.2732410943752503e+00, 1.1957308417178146e-01, 2.9576487677282326e-02, 
	-5.2739372946099017e-01, 2.8867496150229438e-01, 1.2250982323291082e-02, 
	5.2739372946099117e-01, 2.8867496150229432e-01, -1.2250982323291082e-02, 
	1.2732410943752508e+00, 1.1957308417178118e-01, -2.9576487677282381e-02, 
	1.2732410943752503e+00, -1.1957308417178146e-01, -2.9576487677282326e-02, 
	5.2739372946099039e-01, -2.8867496150229427e-01, -1.2250982323291082e-02,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
        -2.5791443630192343e-02, -1.3837517332209948e-02,  8.1578940289174257e-04,
        -6.2965188799992161e-02, -7.6912011102198800e-03,  7.6432003993326847e-04,
        -7.9768181972718816e-02, -2.1429048449640056e-03,  6.6138131401632689e-04,
        -7.8962559012382361e-02,  2.4086914924015437e-03,  5.0697322514091727e-04,
        -6.4691523714997728e-02,  5.5649079307207141e-03,  3.0109577330704072e-04,
        -4.2479347808584755e-02,  6.9270644988373735e-03,  4.3748958514695155e-05,
        -1.9231370953168247e-02,  6.0964812255954216e-03, -2.6506721923612130e-04,
        -3.2340007407779619e-03,  2.6744781398387418e-03, -6.2535275994540806e-04,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
      }
    } else {
    }
  }

  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_cu);
  gkyl_array_release(phi_cu);
  gkyl_array_release(perbuff);
  gkyl_fem_poisson_release(poisson_cu);
}

void
test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs)
{
  double epsilon_0 = 1.0;
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
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_dirichlety, NULL);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_2x_rho_1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&grid, basis, bcs, epsilon_0, NULL, false);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho);

  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_2x_phi_1.gkyl");

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
    gkyl_array_shiftc0(phi, mavgfac*sol_avg[0]);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
	-7.6861958028898086e-02, -4.4376272158426031e-02, 4.4376272158426031e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426100e-02, 4.4376272158426051e-02, 2.5620652676299389e-02,
	7.6861958028898059e-02, 4.4376272158426072e-02, -4.4376272158426079e-02, -2.5620652676299396e-02,
	-7.6861958028898225e-02, -4.4376272158426079e-02, -4.4376272158426065e-02, -2.5620652676299392e-02,
	-7.6861958028898170e-02, -4.4376272158426079e-02, 4.4376272158426100e-02, 2.5620652676299396e-02,
	7.6861958028898197e-02, 4.4376272158426065e-02, 4.4376272158426086e-02, 2.5620652676299392e-02,
	7.6861958028898253e-02, 4.4376272158426079e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898031e-02, -4.4376272158426031e-02, -4.4376272158426072e-02, -2.5620652676299389e-02,
	-1.8556118150391737e-01, -1.8381253775579646e-02, 1.0713379809243175e-01, 1.0612421815373736e-02,
	1.8556118150391748e-01, 1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373724e-02,
	1.8556118150391740e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373727e-02,
	-1.8556118150391759e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, -1.0612421815373736e-02,
	-1.8556118150391751e-01, -1.8381253775579632e-02, 1.0713379809243183e-01, 1.0612421815373741e-02,
	1.8556118150391762e-01, 1.8381253775579681e-02, 1.0713379809243183e-01, 1.0612421815373745e-02,
	1.8556118150391770e-01, 1.8381253775579684e-02, -1.0713379809243177e-01, -1.0612421815373741e-02,
	-1.8556118150391729e-01, -1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373745e-02,
	-1.8556118150391737e-01, 1.8381253775579646e-02, 1.0713379809243175e-01, -1.0612421815373736e-02,
	1.8556118150391746e-01, -1.8381253775579646e-02, 1.0713379809243173e-01, -1.0612421815373734e-02,
	1.8556118150391737e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, 1.0612421815373736e-02,
	-1.8556118150391757e-01, 1.8381253775579663e-02, -1.0713379809243176e-01, 1.0612421815373745e-02,
	-1.8556118150391748e-01, 1.8381253775579656e-02, 1.0713379809243181e-01, -1.0612421815373750e-02,
	1.8556118150391762e-01, -1.8381253775579656e-02, 1.0713379809243183e-01, -1.0612421815373741e-02,
	1.8556118150391770e-01, -1.8381253775579670e-02, -1.0713379809243177e-01, 1.0612421815373731e-02,
	-1.8556118150391729e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, 1.0612421815373745e-02,
	-7.6861958028898114e-02, 4.4376272158426017e-02, 4.4376272158426044e-02, -2.5620652676299378e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426031e-02, -2.5620652676299378e-02,
	7.6861958028898045e-02, -4.4376272158426058e-02, -4.4376272158426065e-02, 2.5620652676299389e-02,
	-7.6861958028898197e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, 2.5620652676299378e-02,
	-7.6861958028898142e-02, 4.4376272158426051e-02, 4.4376272158426079e-02, -2.5620652676299389e-02,
	7.6861958028898197e-02, -4.4376272158426107e-02, 4.4376272158426100e-02, -2.5620652676299396e-02,
	7.6861958028898253e-02, -4.4376272158426114e-02, -4.4376272158426065e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426017e-02, -4.4376272158426072e-02, 2.5620652676299378e-02,
	7.6861958028898031e-02, 4.4376272158426031e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426051e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426044e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426065e-02, 4.4376272158426065e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898114e-02, -4.4376272158426058e-02, -4.4376272158426051e-02, -2.5620652676299392e-02,
	-7.6861958028898114e-02, -4.4376272158426079e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426017e-02, 4.4376272158426031e-02, 2.5620652676299378e-02,
	1.8556118150391737e-01, 1.8381253775579677e-02, -1.0713379809243177e-01, -1.0612421815373757e-02,
	-1.8556118150391754e-01, -1.8381253775579670e-02, -1.0713379809243175e-01, -1.0612421815373741e-02,
	-1.8556118150391748e-01, -1.8381253775579663e-02, 1.0713379809243177e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, 1.8381253775579656e-02, 1.0713379809243177e-01, 1.0612421815373731e-02,
	1.8556118150391748e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, -1.0612421815373731e-02,
	-1.8556118150391740e-01, -1.8381253775579625e-02, -1.0713379809243176e-01, -1.0612421815373736e-02,
	-1.8556118150391743e-01, -1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373731e-02,
	1.8556118150391743e-01, 1.8381253775579687e-02, 1.0713379809243175e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, -1.8381253775579642e-02, -1.0713379809243181e-01, 1.0612421815373738e-02,
	-1.8556118150391765e-01, 1.8381253775579615e-02, -1.0713379809243181e-01, 1.0612421815373708e-02,
	-1.8556118150391762e-01, 1.8381253775579601e-02, 1.0713379809243183e-01, -1.0612421815373717e-02,
	1.8556118150391746e-01, -1.8381253775579642e-02, 1.0713379809243180e-01, -1.0612421815373722e-02,
	1.8556118150391751e-01, -1.8381253775579632e-02, -1.0713379809243176e-01, 1.0612421815373727e-02,
	-1.8556118150391737e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, 1.0612421815373750e-02,
	-1.8556118150391740e-01, 1.8381253775579663e-02, 1.0713379809243173e-01, -1.0612421815373745e-02,
	1.8556118150391748e-01, -1.8381253775579642e-02, 1.0713379809243177e-01, -1.0612421815373736e-02,
	7.6861958028898114e-02, -4.4376272158426051e-02, -4.4376272158426100e-02, 2.5620652676299389e-02,
	-7.6861958028898281e-02, 4.4376272158426114e-02, -4.4376272158426114e-02, 2.5620652676299406e-02,
	-7.6861958028898308e-02, 4.4376272158426107e-02, 4.4376272158426107e-02, -2.5620652676299410e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426086e-02, -2.5620652676299396e-02,
	7.6861958028898156e-02, -4.4376272158426079e-02, -4.4376272158426051e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426044e-02, -4.4376272158426031e-02, 2.5620652676299378e-02,
	-7.6861958028898059e-02, 4.4376272158426058e-02, 4.4376272158426044e-02, -2.5620652676299368e-02,
	7.6861958028898170e-02, -4.4376272158426051e-02, 4.4376272158426058e-02, -2.5620652676299389e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        2.3681559834188498e-01,  1.3672554945099072e-01,  1.3672554945098905e-01,  7.8938532780630122e-02,
        6.3388180537571182e-01,  3.6597183096807612e-01,  9.2520732066094527e-02,  5.3416869563979862e-02,
        8.9067264750613684e-01,  5.1423009279750498e-01,  5.5737529763337020e-02,  3.2180077812826891e-02,
        1.0196938687973565e+00,  5.8872052964116539e-01,  1.8752907080322980e-02,  1.0826995950913283e-02,
        1.0196938687973551e+00,  5.8872052964116617e-01, -1.8752907080323782e-02, -1.0826995950912827e-02,
        8.9067264750613595e-01,  5.1423009279750564e-01, -5.5737529763335895e-02, -3.2180077812827390e-02,
        6.3388180537571470e-01,  3.6597183096807495e-01, -9.2520732066093542e-02, -5.3416869563980424e-02,
        2.3681559834188717e-01,  1.3672554945099089e-01, -1.3672554945099044e-01, -7.8938532780628831e-02,
                                                        
        6.3388180537571459e-01,  9.2520732066094485e-02,  3.6597183096807473e-01,  5.3416869563980014e-02,
        1.7090814385197843e+00,  2.5479496666024115e-01,  2.5479496666024093e-01,  4.0272203461491046e-02,
        2.4167368175137374e+00,  3.6684346655702615e-01,  1.5377005689540441e-01,  2.4419028116211553e-02,
        2.7726579759463506e+00,  4.2335376950105846e-01,  5.1721119735948180e-02,  8.2072105005127904e-03,
        2.7726579759463510e+00,  4.2335376950105857e-01, -5.1721119735948097e-02, -8.2072105005127401e-03,
        2.4167368175137383e+00,  3.6684346655702638e-01, -1.5377005689540432e-01, -2.4419028116211581e-02,
        1.7090814385197852e+00,  2.5479496666024121e-01, -2.5479496666024098e-01, -4.0272203461491102e-02,
        6.3388180537571515e-01,  9.2520732066093431e-02, -3.6597183096807484e-01, -5.3416869563980611e-02,
                                                        
        8.9067264750613706e-01,  5.5737529763335694e-02,  5.1423009279750498e-01,  3.2180077812827605e-02,
        2.4167368175137374e+00,  1.5377005689540435e-01,  3.6684346655702599e-01,  2.4419028116211518e-02,
        3.4373208199425331e+00,  2.2239098197586474e-01,  2.2239098197586490e-01,  1.5199281451033411e-02,
        3.9521264825956046e+00,  2.5761269031357442e-01,  7.4832205937229848e-02,  5.1359813390617843e-03,
        3.9521264825956055e+00,  2.5761269031357470e-01, -7.4832205937229418e-02, -5.1359813390615995e-03,
        3.4373208199425345e+00,  2.2239098197586496e-01, -2.2239098197586502e-01, -1.5199281451033622e-02,
        2.4167368175137383e+00,  1.5377005689540424e-01, -3.6684346655702615e-01, -2.4419028116211505e-02,
        8.9067264750613750e-01,  5.5737529763336638e-02, -5.1423009279750520e-01, -3.2180077812827015e-02,
                                                        
        1.0196938687973565e+00,  1.8752907080324045e-02,  5.8872052964116561e-01,  1.0826995950912709e-02,
        2.7726579759463501e+00,  5.1721119735948284e-02,  4.2335376950105824e-01,  8.2072105005128095e-03,
        3.9521264825956042e+00,  7.4832205937229584e-02,  2.5761269031357442e-01,  5.1359813390616359e-03,
        4.5485900655688685e+00,  8.6755719877848164e-02,  8.6755719877848719e-02,  1.7480626442406752e-03,
        4.5485900655688685e+00,  8.6755719877847776e-02, -8.6755719877848275e-02, -1.7480626442409205e-03,
        3.9521264825956046e+00,  7.4832205937228988e-02, -2.5761269031357459e-01, -5.1359813390615587e-03,
        2.7726579759463506e+00,  5.1721119735947924e-02, -4.2335376950105819e-01, -8.2072105005127505e-03,
        1.0196938687973565e+00,  1.8752907080322814e-02, -5.8872052964116572e-01, -1.0826995950913227e-02,
                                                        
        1.0196938687973578e+00, -1.8752907080323244e-02,  5.8872052964116495e-01, -1.0826995950913076e-02,
        2.7726579759463510e+00, -5.1721119735948159e-02,  4.2335376950105824e-01, -8.2072105005127211e-03,
        3.9521264825956046e+00, -7.4832205937229349e-02,  2.5761269031357431e-01, -5.1359813390616368e-03,
        4.5485900655688685e+00, -8.6755719877848525e-02,  8.6755719877847970e-02, -1.7480626442410686e-03,
        4.5485900655688667e+00, -8.6755719877848803e-02, -8.6755719877848247e-02,  1.7480626442408095e-03,
        3.9521264825956033e+00, -7.4832205937229904e-02, -2.5761269031357409e-01,  5.1359813390617765e-03,
        2.7726579759463501e+00, -5.1721119735948243e-02, -4.2335376950105824e-01,  8.2072105005128078e-03,
        1.0196938687973569e+00, -1.8752907080322453e-02, -5.8872052964116517e-01,  1.0826995950913517e-02,
                                                        
        8.9067264750613750e-01, -5.5737529763336964e-02,  5.1423009279750509e-01, -3.2180077812826779e-02,
        2.4167368175137387e+00, -1.5377005689540407e-01,  3.6684346655702621e-01, -2.4419028116211449e-02,
        3.4373208199425340e+00, -2.2239098197586465e-01,  2.2239098197586460e-01, -1.5199281451033603e-02,
        3.9521264825956037e+00, -2.5761269031357431e-01,  7.4832205937229085e-02, -5.1359813390615770e-03,
        3.9521264825956024e+00, -2.5761269031357420e-01, -7.4832205937229654e-02,  5.1359813390616099e-03,
        3.4373208199425314e+00, -2.2239098197586477e-01, -2.2239098197586432e-01,  1.5199281451033496e-02,
        2.4167368175137374e+00, -1.5377005689540416e-01, -3.6684346655702588e-01,  2.4419028116211574e-02,
        8.9067264750613806e-01, -5.5737529763336978e-02, -5.1423009279750453e-01,  3.2180077812826766e-02,
                                                        
        6.3388180537571459e-01, -9.2520732066093431e-02,  3.6597183096807528e-01, -5.3416869563980569e-02,
        1.7090814385197859e+00, -2.5479496666024121e-01,  2.5479496666024121e-01, -4.0272203461491109e-02,
        2.4167368175137387e+00, -3.6684346655702615e-01,  1.5377005689540399e-01, -2.4419028116211421e-02,
        2.7726579759463506e+00, -4.2335376950105813e-01,  5.1721119735947743e-02, -8.2072105005127610e-03,
        2.7726579759463492e+00, -4.2335376950105796e-01, -5.1721119735948368e-02,  8.2072105005128824e-03,
        2.4167368175137369e+00, -3.6684346655702577e-01, -1.5377005689540402e-01,  2.4419028116211414e-02,
        1.7090814385197848e+00, -2.5479496666024104e-01, -2.5479496666024082e-01,  4.0272203461491039e-02,
        6.3388180537571470e-01, -9.2520732066093667e-02, -3.6597183096807501e-01,  5.3416869563980458e-02,
                                                        
        2.3681559834188717e-01, -1.3672554945099052e-01,  1.3672554945099086e-01, -7.8938532780628887e-02,
        6.3388180537571537e-01, -3.6597183096807490e-01,  9.2520732066093569e-02, -5.3416869563980528e-02,
        8.9067264750613795e-01, -5.1423009279750498e-01,  5.5737529763336596e-02, -3.2180077812827008e-02,
        1.0196938687973571e+00, -5.8872052964116517e-01,  1.8752907080323126e-02, -1.0826995950912990e-02,
        1.0196938687973565e+00, -5.8872052964116517e-01, -1.8752907080323424e-02,  1.0826995950912984e-02,
        8.9067264750613695e-01, -5.1423009279750498e-01, -5.5737529763336513e-02,  3.2180077812827092e-02,
        6.3388180537571459e-01, -3.6597183096807484e-01, -9.2520732066093569e-02,  5.3416869563980382e-02,
        2.3681559834188676e-01, -1.3672554945099064e-01, -1.3672554945099064e-01,  7.8938532780628942e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        3.3943303570130012e-01,  1.9597175453399679e-01, -2.5382106732311723e-02, -1.4654366154499068e-02,
        2.5852815051768180e-01,  1.4926129729448101e-01, -2.1328350507206096e-02, -1.2313928906705882e-02,
        2.3704347887786692e-01,  1.3685711633978268e-01,  8.9241695525082501e-03,  5.1523716934339632e-03,
        2.9838358920906116e-01,  1.7227184555161776e-01,  2.6490559659328117e-02,  1.5294331750296710e-02,
        2.7855144864018777e-01,  1.6082175385557235e-01, -3.7940651355373674e-02, -2.1905045273254424e-02,
        2.0381615307077769e-01,  1.1767331084060774e-01, -5.2077916595914715e-03, -3.0067199165486271e-03,
        3.1270423370037509e-01,  1.8053987350364709e-01,  6.8074354322631203e-02,  3.9302746793080885e-02,
        4.0700430412220923e-01,  2.3498404454629240e-01, -1.3630183279984597e-02, -7.8693899858035734e-03,
                                                        
        9.2471721182926658e-01,  1.4194222210591301e-01, -5.7372732217365302e-02, -3.8154300808411140e-03,
        7.5028336345656188e-01,  1.3465370727118317e-01, -4.3336697096398154e-02, -3.9259592098273228e-04,
        7.2217959325319470e-01,  1.4323634988842188e-01,  2.7110977804240665e-02,  5.3477869464038522e-03,
        8.6863423587078081e-01,  1.5696251880408427e-01,  5.7444649535093165e-02,  2.5770203719961103e-03,
        8.2578551259476063e-01,  1.5512398027833324e-01, -8.2183371452935167e-02, -3.6385010847539575e-03,
        6.5796549632452850e-01,  1.4452993473923709e-01, -1.4707560116090753e-02, -2.4779739590503548e-03,
        8.7197741194623291e-01,  1.4235664650713659e-01,  1.3826739721006934e-01,  1.2232254132204976e-03,
        1.0677766339715837e+00,  1.4651303796530019e-01, -2.5222663666613812e-02,  1.1764683140077245e-03,
                                                        
        1.3236042788217255e+00,  8.8355333398443453e-02, -6.3319731576856511e-02,  3.8192839976832317e-04,
        1.1500636120564152e+00,  9.6159526874637705e-02, -3.6874019095442509e-02,  4.1238248045203435e-03,
        1.1823773590563609e+00,  1.2245895407850058e-01,  5.5530369624386594e-02,  1.1060156571162853e-02,
        1.4109181535507069e+00,  1.5612524704562614e-01,  7.6417719597735217e-02,  8.3770867360238831e-03,
        1.3654745745531423e+00,  1.5646564492203485e-01, -1.0265458216491942e-01, -8.1805579304477273e-03,
        1.1227306086903228e+00,  1.2380232799510095e-01, -3.7493711870136587e-02, -1.0677616889943637e-02,
        1.2888335450678960e+00,  9.8315354164002000e-02,  1.3339328690091568e-01, -4.0372943122702881e-03,
        1.4765773832016793e+00,  8.9508184647616132e-02, -2.4999331415682426e-02, -1.0475273788137696e-03,
                                                        
        1.5295214682804199e+00,  3.0531011366304865e-02, -6.1912475511358124e-02,  4.3055126846589144e-04,
        1.3790145293808918e+00,  3.6025346873860498e-02, -2.4982746177190984e-02,  2.7416048158394697e-03,
        1.4901027388499875e+00,  5.5206376781829561e-02,  8.9119553817967004e-02,  8.3325679648607077e-03,
        1.8232429688045224e+00,  8.1930596034726758e-02,  1.0321904762417754e-01,  7.0966672146820734e-03,
        1.7786475086544757e+00,  8.2079859803237842e-02, -1.2896624854644198e-01, -7.0104897377519498e-03,
        1.4335217115489933e+00,  5.5632998902093111e-02, -7.0292223316692631e-02, -8.2586125227445510e-03,
        1.5224137527023140e+00,  3.6542241591097972e-02,  1.2161406720537776e-01, -2.7634413497923073e-03,
        1.6849062587527901e+00,  3.0770547731789839e-02, -2.7798975095838528e-02, -5.6884765355937322e-04,
                                                        
        1.5295214682804197e+00, -3.0531011366304792e-02, -6.1912475511358096e-02, -4.3055126846587653e-04,
        1.3790145293808918e+00, -3.6025346873860394e-02, -2.4982746177190963e-02, -2.7416048158394866e-03,
        1.4901027388499877e+00, -5.5206376781829401e-02,  8.9119553817967018e-02, -8.3325679648606973e-03,
        1.8232429688045226e+00, -8.1930596034726688e-02,  1.0321904762417747e-01, -7.0966672146820916e-03,
        1.7786475086544760e+00, -8.2079859803237690e-02, -1.2896624854644190e-01,  7.0104897377519706e-03,
        1.4335217115489931e+00, -5.5632998902093125e-02, -7.0292223316692812e-02,  8.2586125227444539e-03,
        1.5224137527023138e+00, -3.6542241591098125e-02,  1.2161406720537783e-01,  2.7634413497923224e-03,
        1.6849062587527901e+00, -3.0770547731789867e-02, -2.7798975095838487e-02,  5.6884765355936563e-04,
                                                        
        1.3236042788217257e+00, -8.8355333398443481e-02, -6.3319731576856414e-02, -3.8192839976835028e-04,
        1.1500636120564156e+00, -9.6159526874637732e-02, -3.6874019095442509e-02, -4.1238248045203391e-03,
        1.1823773590563613e+00, -1.2245895407850059e-01,  5.5530369624386615e-02, -1.1060156571162844e-02,
        1.4109181535507074e+00, -1.5612524704562611e-01,  7.6417719597735093e-02, -8.3770867360238987e-03,
        1.3654745745531423e+00, -1.5646564492203499e-01, -1.0265458216491942e-01,  8.1805579304476995e-03,
        1.1227306086903228e+00, -1.2380232799510095e-01, -3.7493711870136739e-02,  1.0677616889943734e-02,
        1.2888335450678954e+00, -9.8315354164002056e-02,  1.3339328690091565e-01,  4.0372943122701779e-03,
        1.4765773832016789e+00, -8.9508184647616243e-02, -2.4999331415682322e-02,  1.0475273788138081e-03,
                                                        
        9.2471721182926658e-01, -1.4194222210591306e-01, -5.7372732217365198e-02,  3.8154300808411362e-03,
        7.5028336345656221e-01, -1.3465370727118314e-01, -4.3336697096398057e-02,  3.9259592098277505e-04,
        7.2217959325319514e-01, -1.4323634988842177e-01,  2.7110977804240696e-02, -5.3477869464038444e-03,
        8.6863423587078115e-01, -1.5696251880408424e-01,  5.7444649535093005e-02, -2.5770203719961463e-03,
        8.2578551259476063e-01, -1.5512398027833330e-01, -8.2183371452935264e-02,  3.6385010847539753e-03,
        6.5796549632452839e-01, -1.4452993473923711e-01, -1.4707560116090757e-02,  2.4779739590503531e-03,
        8.7197741194623257e-01, -1.4235664650713659e-01,  1.3826739721006917e-01, -1.2232254132204556e-03,
        1.0677766339715833e+00, -1.4651303796530019e-01, -2.5222663666613656e-02, -1.1764683140077696e-03,
                                                        
        3.3943303570129962e-01, -1.9597175453399704e-01, -2.5382106732311140e-02,  1.4654366154499322e-02,
        2.5852815051768296e-01, -1.4926129729448057e-01, -2.1328350507205690e-02,  1.2313928906706016e-02,
        2.3704347887786834e-01, -1.3685711633978229e-01,  8.9241695525079847e-03, -5.1523716934341332e-03,
        2.9838358920906127e-01, -1.7227184555161790e-01,  2.6490559659327607e-02, -1.5294331750296876e-02,
        2.7855144864018744e-01, -1.6082175385557249e-01, -3.7940651355373403e-02,  2.1905045273254604e-02,
        2.0381615307077783e-01, -1.1767331084060755e-01, -5.2077916595914967e-03,  3.0067199165486245e-03,
        3.1270423370037476e-01, -1.8053987350364703e-01,  6.8074354322630953e-02, -3.9302746793080955e-02,
        4.0700430412220806e-01, -2.3498404454629279e-01, -1.3630183279984807e-02,  7.8693899858033965e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against dirichletx-periodicy,
      // rotated by pgkyl):
      const double sol[256] = {
         3.3943303570129929e-01, -2.5382106732310627e-02,  1.9597175453399684e-01, -1.4654366154499540e-02,
         9.2471721182926570e-01, -5.7372732217365149e-02,  1.4194222210591295e-01, -3.8154300808411951e-03,
         1.3236042788217244e+00, -6.3319731576856414e-02,  8.8355333398443286e-02,  3.8192839976837950e-04,
         1.5295214682804186e+00, -6.1912475511358069e-02,  3.0531011366304959e-02,  4.3055126846582503e-04,
         1.5295214682804188e+00, -6.1912475511358138e-02, -3.0531011366304928e-02, -4.3055126846586200e-04,
         1.3236042788217248e+00, -6.3319731576856483e-02, -8.8355333398443120e-02, -3.8192839976834253e-04,
         9.2471721182926625e-01, -5.7372732217365198e-02, -1.4194222210591306e-01,  3.8154300808411765e-03,
         3.3943303570129940e-01, -2.5382106732310755e-02, -1.9597175453399698e-01,  1.4654366154499521e-02,
                                                         
         2.5852815051768291e-01, -2.1328350507206047e-02,  1.4926129729448026e-01, -1.2313928906705863e-02,
         7.5028336345656133e-01, -4.3336697096398091e-02,  1.3465370727118300e-01, -3.9259592098274870e-04,
         1.1500636120564143e+00, -3.6874019095442460e-02,  9.6159526874637649e-02,  4.1238248045203573e-03,
         1.3790145293808904e+00, -2.4982746177191012e-02,  3.6025346873860394e-02,  2.7416048158394237e-03,
         1.3790145293808904e+00, -2.4982746177190981e-02, -3.6025346873860421e-02, -2.7416048158394055e-03,
         1.1500636120564147e+00, -3.6874019095442495e-02, -9.6159526874637552e-02, -4.1238248045204129e-03,
         7.5028336345656177e-01, -4.3336697096398126e-02, -1.3465370727118303e-01,  3.9259592098280421e-04,
         2.5852815051768296e-01, -2.1328350507205975e-02, -1.4926129729448043e-01,  1.2313928906705860e-02,
                                                         
         2.3704347887786717e-01,  8.9241695525077245e-03,  1.3685711633978231e-01,  5.1523716934341627e-03,
         7.2217959325319403e-01,  2.7110977804240554e-02,  1.4323634988842177e-01,  5.3477869464038947e-03,
         1.1823773590563600e+00,  5.5530369624386580e-02,  1.2245895407850053e-01,  1.1060156571162867e-02,
         1.4901027388499859e+00,  8.9119553817966893e-02,  5.5206376781829269e-02,  8.3325679648606505e-03,
         1.4901027388499863e+00,  8.9119553817966962e-02, -5.5206376781829172e-02, -8.3325679648606332e-03,
         1.1823773590563604e+00,  5.5530369624386740e-02, -1.2245895407850046e-01, -1.1060156571162811e-02,
         7.2217959325319470e-01,  2.7110977804240727e-02, -1.4323634988842171e-01, -5.3477869464039407e-03,
         2.3704347887786764e-01,  8.9241695525078789e-03, -1.3685711633978243e-01, -5.1523716934341384e-03,
                                                         
         2.9838358920905961e-01,  2.6490559659327586e-02,  1.7227184555161815e-01,  1.5294331750296947e-02,
         8.6863423587077970e-01,  5.7444649535092998e-02,  1.5696251880408424e-01,  2.5770203719960869e-03,
         1.4109181535507060e+00,  7.6417719597735009e-02,  1.5612524704562616e-01,  8.3770867360239074e-03,
         1.8232429688045209e+00,  1.0321904762417757e-01,  8.1930596034726591e-02,  7.0966672146821827e-03,
         1.8232429688045211e+00,  1.0321904762417759e-01, -8.1930596034726424e-02, -7.0966672146822009e-03,
         1.4109181535507065e+00,  7.6417719597734982e-02, -1.5612524704562597e-01, -8.3770867360238883e-03,
         8.6863423587078070e-01,  5.7444649535092970e-02, -1.5696251880408427e-01, -2.5770203719961051e-03,
         2.9838358920906094e-01,  2.6490559659327881e-02, -1.7227184555161787e-01, -1.5294331750296758e-02,
                                                         
         2.7855144864018666e-01, -3.7940651355372897e-02,  1.6082175385557240e-01, -2.1905045273254868e-02,
         8.2578551259475919e-01, -8.2183371452935167e-02,  1.5512398027833305e-01, -3.6385010847539805e-03,
         1.3654745745531405e+00, -1.0265458216491943e-01,  1.5646564492203499e-01, -8.1805579304476978e-03,
         1.7786475086544749e+00, -1.2896624854644190e-01,  8.2079859803237828e-02, -7.0104897377519255e-03,
         1.7786475086544753e+00, -1.2896624854644170e-01, -8.2079859803237634e-02,  7.0104897377520183e-03,
         1.3654745745531418e+00, -1.0265458216491913e-01, -1.5646564492203482e-01,  8.1805579304476423e-03,
         8.2578551259476041e-01, -8.2183371452935056e-02, -1.5512398027833324e-01,  3.6385010847539249e-03,
         2.7855144864018744e-01, -3.7940651355373500e-02, -1.6082175385557240e-01,  2.1905045273254500e-02,
                                                         
         2.0381615307077769e-01, -5.2077916595916198e-03,  1.1767331084060717e-01, -3.0067199165485217e-03,
         6.5796549632452717e-01, -1.4707560116090692e-02,  1.4452993473923695e-01, -2.4779739590503344e-03,
         1.1227306086903213e+00, -3.7493711870136663e-02,  1.2380232799510092e-01, -1.0677616889943762e-02,
         1.4335217115489918e+00, -7.0292223316692895e-02,  5.5632998902093299e-02, -8.2586125227444972e-03,
         1.4335217115489924e+00, -7.0292223316692923e-02, -5.5632998902092813e-02,  8.2586125227444782e-03,
         1.1227306086903228e+00, -3.7493711870136760e-02, -1.2380232799510088e-01,  1.0677616889943752e-02,
         6.5796549632452839e-01, -1.4707560116090772e-02, -1.4452993473923714e-01,  2.4779739590503622e-03,
         2.0381615307077774e-01, -5.2077916595914247e-03, -1.1767331084060761e-01,  3.0067199165486531e-03,
                                                         
         3.1270423370037409e-01,  6.8074354322630801e-02,  1.8053987350364703e-01,  3.9302746793081066e-02,
         8.7197741194623157e-01,  1.3826739721006912e-01,  1.4235664650713634e-01,  1.2232254132204063e-03,
         1.2888335450678943e+00,  1.3339328690091576e-01,  9.8315354164002167e-02, -4.0372943122700409e-03,
         1.5224137527023127e+00,  1.2161406720537798e-01,  3.6542241591098042e-02, -2.7634413497924326e-03,
         1.5224137527023132e+00,  1.2161406720537785e-01, -3.6542241591097750e-02,  2.7634413497923402e-03,
         1.2888335450678956e+00,  1.3339328690091576e-01, -9.8315354164001792e-02,  4.0372943122702351e-03,
         8.7197741194623291e-01,  1.3826739721006936e-01, -1.4235664650713664e-01, -1.2232254132204989e-03,
         3.1270423370037487e-01,  6.8074354322631009e-02, -1.8053987350364709e-01, -3.9302746793080982e-02,
                                                         
         4.0700430412220695e-01, -1.3630183279984911e-02,  2.3498404454629299e-01, -7.8693899858033878e-03,
         1.0677766339715822e+00, -2.5222663666613573e-02,  1.4651303796530005e-01,  1.1764683140078709e-03,
         1.4765773832016775e+00, -2.4999331415682360e-02,  8.9508184647616285e-02, -1.0475273788140100e-03,
         1.6849062587527890e+00, -2.7798975095838598e-02,  3.0770547731789843e-02, -5.6884765355922739e-04,
         1.6849062587527892e+00, -2.7798975095838566e-02, -3.0770547731789812e-02,  5.6884765355924582e-04,
         1.4765773832016786e+00, -2.4999331415682616e-02, -8.9508184647615868e-02,  1.0475273788138251e-03,
         1.0677766339715835e+00, -2.5222663666613909e-02, -1.4651303796530027e-01, -1.1764683140077321e-03,
         4.0700430412220773e-01, -1.3630183279985131e-02, -2.3498404454629301e-01,  7.8693899858033340e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else {
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
	-9.5101835527468453e-02, -5.2237861001368492e-02, 4.3177375058357795e-02, 2.3387400381385207e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, -1.1937067064604241e-03, 5.2456792020110785e-03, 
	9.5101835527468453e-02, 5.2237861001368610e-02, 4.3177375058357802e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905626e-03, -1.1937067064604172e-03, -5.2456792020110915e-03, 
	9.5101835527468384e-02, 5.2237861001368582e-02, -4.3177375058357843e-02, -2.3387400381385235e-02, -2.0675606649251488e-03, -9.0857828980905643e-03, 1.1937067064604181e-03, -5.2456792020110924e-03, 
	-9.5101835527468620e-02, -5.2237861001368589e-02, -4.3177375058357788e-02, -2.3387400381385214e-02, 2.0675606649251388e-03, 9.0857828980905712e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, 2.3387400381385225e-02, 2.0675606649251596e-03, 9.0857828980905470e-03, -1.1937067064604137e-03, 5.2456792020110837e-03, 
	9.5101835527468453e-02, 5.2237861001368548e-02, 4.3177375058357816e-02, 2.3387400381385211e-02, -2.0675606649251488e-03, -9.0857828980905591e-03, -1.1937067064604259e-03, -5.2456792020110759e-03, 
	9.5101835527468537e-02, 5.2237861001368562e-02, -4.3177375058357795e-02, -2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905782e-03, 1.1937067064604241e-03, -5.2456792020110915e-03, 
	-9.5101835527468398e-02, -5.2237861001368506e-02, -4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251596e-03, 9.0857828980905504e-03, 1.1937067064604207e-03, 5.2456792020110776e-03, 
	-2.2959614113698959e-01, -2.1637630496127446e-02, 1.0423940445355719e-01, 9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347764e-02, -2.8818629202324655e-03, 2.1728314693314739e-03, 
	2.2959614113698990e-01, 2.1637630496127457e-02, 1.0423940445355717e-01, 9.6873784266194498e-03, -4.9915329982914251e-03, -2.1935020297347764e-02, -2.8818629202324387e-03, -2.1728314693314635e-03, 
	2.2959614113698973e-01, 2.1637630496127453e-02, -1.0423940445355724e-01, -9.6873784266194377e-03, -4.9915329982914043e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, -2.1728314693314518e-03, 
	-2.2959614113698987e-01, -2.1637630496127373e-02, -1.0423940445355719e-01, -9.6873784266194654e-03, 4.9915329982914147e-03, 2.1935020297347754e-02, 2.8818629202324413e-03, 2.1728314693314496e-03, 
	-2.2959614113698981e-01, -2.1637630496127425e-02, 1.0423940445355721e-01, 9.6873784266194429e-03, 4.9915329982914147e-03, 2.1935020297347736e-02, -2.8818629202324378e-03, 2.1728314693314679e-03, 
	2.2959614113698976e-01, 2.1637630496127457e-02, 1.0423940445355720e-01, 9.6873784266194567e-03, -4.9915329982914113e-03, -2.1935020297347757e-02, -2.8818629202324551e-03, -2.1728314693314783e-03, 
	2.2959614113698987e-01, 2.1637630496127443e-02, -1.0423940445355716e-01, -9.6873784266194429e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324551e-03, -2.1728314693314600e-03, 
	-2.2959614113698959e-01, -2.1637630496127443e-02, -1.0423940445355719e-01, -9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324560e-03, 2.1728314693314809e-03, 
	-2.2959614113698965e-01, 2.1637630496127432e-02, 1.0423940445355719e-01, -9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347771e-02, -2.8818629202324655e-03, -2.1728314693314657e-03, 
	2.2959614113698987e-01, -2.1637630496127495e-02, 1.0423940445355716e-01, -9.6873784266194689e-03, -4.9915329982914321e-03, -2.1935020297347757e-02, -2.8818629202324499e-03, 2.1728314693314679e-03, 
	2.2959614113698973e-01, -2.1637630496127470e-02, -1.0423940445355720e-01, 9.6873784266194654e-03, -4.9915329982914217e-03, -2.1935020297347750e-02, 2.8818629202324586e-03, 2.1728314693314540e-03, 
	-2.2959614113698976e-01, 2.1637630496127436e-02, -1.0423940445355717e-01, 9.6873784266194793e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324551e-03, -2.1728314693314583e-03, 
	-2.2959614113698976e-01, 2.1637630496127460e-02, 1.0423940445355719e-01, -9.6873784266194758e-03, 4.9915329982914286e-03, 2.1935020297347740e-02, -2.8818629202324586e-03, -2.1728314693314653e-03, 
	2.2959614113698981e-01, -2.1637630496127446e-02, 1.0423940445355721e-01, -9.6873784266194377e-03, -4.9915329982914251e-03, -2.1935020297347761e-02, -2.8818629202324577e-03, 2.1728314693314722e-03, 
	2.2959614113698987e-01, -2.1637630496127453e-02, -1.0423940445355717e-01, 9.6873784266194238e-03, -4.9915329982914321e-03, -2.1935020297347747e-02, 2.8818629202324586e-03, 2.1728314693314683e-03, 
	-2.2959614113698959e-01, 2.1637630496127436e-02, -1.0423940445355719e-01, 9.6873784266194619e-03, 4.9915329982914147e-03, 2.1935020297347750e-02, 2.8818629202324603e-03, -2.1728314693314774e-03, 
	-9.5101835527468509e-02, 5.2237861001368499e-02, 4.3177375058357795e-02, -2.3387400381385211e-02, 2.0675606649251735e-03, 9.0857828980905678e-03, -1.1937067064604276e-03, -5.2456792020110863e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, 4.3177375058357781e-02, -2.3387400381385200e-02, -2.0675606649251449e-03, -9.0857828980905608e-03, -1.1937067064604207e-03, 5.2456792020110811e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, -4.3177375058357788e-02, 2.3387400381385211e-02, -2.0675606649251509e-03, -9.0857828980905660e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, 5.2237861001368534e-02, -4.3177375058357767e-02, 2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604215e-03, -5.2456792020110932e-03, 
	-9.5101835527468426e-02, 5.2237861001368534e-02, 4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251457e-03, 9.0857828980905556e-03, -1.1937067064604207e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, -5.2237861001368527e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251423e-03, -9.0857828980905574e-03, -1.1937067064604128e-03, 5.2456792020110846e-03, 
	9.5101835527468495e-02, -5.2237861001368555e-02, -4.3177375058357816e-02, 2.3387400381385221e-02, -2.0675606649251466e-03, -9.0857828980905504e-03, 1.1937067064604102e-03, 5.2456792020110854e-03, 
	-9.5101835527468453e-02, 5.2237861001368499e-02, -4.3177375058357823e-02, 2.3387400381385200e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, 1.1937067064604337e-03, -5.2456792020110776e-03, 
	9.5101835527468231e-02, 5.2237861001368499e-02, -4.3177375058357809e-02, -2.3387400381385207e-02, -2.0675606649251488e-03, -9.0857828980905678e-03, 1.1937067064604233e-03, -5.2456792020110932e-03, 
	-9.5101835527468676e-02, -5.2237861001368610e-02, -4.3177375058357767e-02, -2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905591e-03, 1.1937067064604137e-03, 5.2456792020110846e-03, 
	-9.5101835527468565e-02, -5.2237861001368562e-02, 4.3177375058357836e-02, 2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905678e-03, -1.1937067064604172e-03, 5.2456792020110837e-03, 
	9.5101835527468426e-02, 5.2237861001368582e-02, 4.3177375058357823e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604154e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, 5.2237861001368589e-02, -4.3177375058357753e-02, -2.3387400381385214e-02, -2.0675606649251423e-03, -9.0857828980905435e-03, 1.1937067064604198e-03, -5.2456792020110820e-03, 
	-9.5101835527468231e-02, -5.2237861001368478e-02, -4.3177375058357809e-02, -2.3387400381385228e-02, 2.0675606649251596e-03, 9.0857828980905383e-03, 1.1937067064604172e-03, 5.2456792020110707e-03, 
	-9.5101835527468370e-02, -5.2237861001368541e-02, 4.3177375058357739e-02, 2.3387400381385197e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, -1.1937067064604207e-03, 5.2456792020110724e-03, 
	9.5101835527468315e-02, 5.2237861001368499e-02, 4.3177375058357753e-02, 2.3387400381385200e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604215e-03, -5.2456792020110898e-03, 
	2.2959614113698942e-01, 2.1637630496127418e-02, -1.0423940445355720e-01, -9.6873784266194619e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324586e-03, -2.1728314693314574e-03, 
	-2.2959614113699003e-01, -2.1637630496127470e-02, -1.0423940445355714e-01, -9.6873784266194619e-03, 4.9915329982914008e-03, 2.1935020297347740e-02, 2.8818629202324387e-03, 2.1728314693314609e-03, 
	-2.2959614113698998e-01, -2.1637630496127495e-02, 1.0423940445355724e-01, 9.6873784266194567e-03, 4.9915329982914008e-03, 2.1935020297347771e-02, -2.8818629202324447e-03, 2.1728314693314791e-03, 
	2.2959614113698981e-01, 2.1637630496127450e-02, 1.0423940445355726e-01, 9.6873784266194706e-03, -4.9915329982914286e-03, -2.1935020297347750e-02, -2.8818629202324638e-03, -2.1728314693314713e-03, 
	2.2959614113699001e-01, 2.1637630496127470e-02, -1.0423940445355717e-01, -9.6873784266194706e-03, -4.9915329982914355e-03, -2.1935020297347764e-02, 2.8818629202324621e-03, -2.1728314693314817e-03, 
	-2.2959614113698937e-01, -2.1637630496127439e-02, -1.0423940445355716e-01, -9.6873784266194151e-03, 4.9915329982914147e-03, 2.1935020297347722e-02, 2.8818629202324586e-03, 2.1728314693314783e-03, 
	-2.2959614113698959e-01, -2.1637630496127418e-02, 1.0423940445355706e-01, 9.6873784266194290e-03, 4.9915329982914286e-03, 2.1935020297347764e-02, -2.8818629202324586e-03, 2.1728314693314869e-03, 
	2.2959614113698953e-01, 2.1637630496127457e-02, 1.0423940445355712e-01, 9.6873784266194498e-03, -4.9915329982914182e-03, -2.1935020297347757e-02, -2.8818629202324560e-03, -2.1728314693314627e-03, 
	2.2959614113698945e-01, -2.1637630496127394e-02, -1.0423940445355723e-01, 9.6873784266194342e-03, -4.9915329982914147e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, 2.1728314693314613e-03, 
	-2.2959614113699015e-01, 2.1637630496127450e-02, -1.0423940445355720e-01, 9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324577e-03, -2.1728314693314653e-03, 
	-2.2959614113699009e-01, 2.1637630496127460e-02, 1.0423940445355727e-01, -9.6873784266194654e-03, 4.9915329982914425e-03, 2.1935020297347778e-02, -2.8818629202324517e-03, -2.1728314693314800e-03, 
	2.2959614113698970e-01, -2.1637630496127484e-02, 1.0423940445355724e-01, -9.6873784266194619e-03, -4.9915329982914008e-03, -2.1935020297347764e-02, -2.8818629202324517e-03, 2.1728314693314661e-03, 
	2.2959614113698992e-01, -2.1637630496127498e-02, -1.0423940445355714e-01, 9.6873784266194654e-03, -4.9915329982914182e-03, -2.1935020297347764e-02, 2.8818629202324447e-03, 2.1728314693314761e-03, 
	-2.2959614113698948e-01, 2.1637630496127401e-02, -1.0423940445355717e-01, 9.6873784266194342e-03, 4.9915329982914564e-03, 2.1935020297347750e-02, 2.8818629202324829e-03, -2.1728314693314635e-03, 
	-2.2959614113698953e-01, 2.1637630496127484e-02, 1.0423940445355713e-01, -9.6873784266194203e-03, 4.9915329982914425e-03, 2.1935020297347750e-02, -2.8818629202324794e-03, -2.1728314693314913e-03, 
	2.2959614113698965e-01, -2.1637630496127401e-02, 1.0423940445355709e-01, -9.6873784266194446e-03, -4.9915329982914251e-03, -2.1935020297347743e-02, -2.8818629202324517e-03, 2.1728314693314670e-03, 
	9.5101835527468259e-02, -5.2237861001368513e-02, -4.3177375058357836e-02, 2.3387400381385239e-02, -2.0675606649251405e-03, -9.0857828980905626e-03, 1.1937067064604181e-03, 5.2456792020110898e-03, 
	-9.5101835527468676e-02, 5.2237861001368631e-02, -4.3177375058357767e-02, 2.3387400381385214e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604068e-03, -5.2456792020110776e-03, 
	-9.5101835527468620e-02, 5.2237861001368603e-02, 4.3177375058357795e-02, -2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905608e-03, -1.1937067064604102e-03, -5.2456792020110828e-03, 
	9.5101835527468356e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251492e-03, -9.0857828980905608e-03, -1.1937067064604163e-03, 5.2456792020110880e-03, 
	9.5101835527468453e-02, -5.2237861001368582e-02, -4.3177375058357760e-02, 2.3387400381385211e-02, -2.0675606649251466e-03, -9.0857828980905591e-03, 1.1937067064604172e-03, 5.2456792020110811e-03, 
	-9.5101835527468342e-02, 5.2237861001368485e-02, -4.3177375058357788e-02, 2.3387400381385207e-02, 2.0675606649251596e-03, 9.0857828980905574e-03, 1.1937067064604267e-03, -5.2456792020110863e-03, 
	-9.5101835527468287e-02, 5.2237861001368478e-02, 4.3177375058357788e-02, -2.3387400381385197e-02, 2.0675606649251735e-03, 9.0857828980905383e-03, -1.1937067064604241e-03, -5.2456792020110768e-03, 
	9.5101835527468398e-02, -5.2237861001368541e-02, 4.3177375058357746e-02, -2.3387400381385211e-02, -2.0675606649251379e-03, -9.0857828980905539e-03, -1.1937067064604189e-03, 5.2456792020110811e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.6174232732491770e-01,  1.4239797833436429e-01,  1.4239797833437121e-01,  7.7179579805169304e-02,
        -6.7537275729811939e-03, -6.7537275729829148e-03, -3.8992664322967511e-03, -3.8992664322970148e-03,
         6.6919163252610048e-01,  3.6607316345625468e-01,  9.5063998582039566e-02,  5.3241649763932503e-02,
        -1.5712542948565247e-02, -5.0333399450138384e-03, -1.2731080364168593e-03, -2.9060001721712362e-03,
         9.3261285108184055e-01,  5.1325052577226338e-01,  5.7325647937302174e-02,  3.1906370575560429e-02,
        -1.9514998592151618e-02, -4.7983776055374530e-03, -9.2224075298933392e-04, -2.7703446022362605e-03,
         1.0655299621632952e+00,  5.8724826482558512e-01,  1.9276481795971818e-02,  1.0736800514483670e-02,
        -2.1638935540250832e-02, -4.9049622075441517e-03, -3.0401481573749179e-04, -2.8318812508971261e-03,
         1.0655299621633028e+00,  5.8724826482558323e-01, -1.9276481795972602e-02, -1.0736800514480182e-02,
        -2.1638935540251626e-02, -4.9049622075482769e-03,  3.0401481573705659e-04, -2.8318812508936436e-03,
         9.3261285108186132e-01,  5.1325052577225383e-01, -5.7325647937298004e-02, -3.1906370575565557e-02,
        -1.9514998592151112e-02, -4.7983776055450077e-03,  9.2224075299059171e-04, -2.7703446022306617e-03,
         6.6919163252609015e-01,  3.6607316345626151e-01, -9.5063998582045506e-02, -5.3241649763929401e-02,
        -1.5712542948564903e-02, -5.0333399450088745e-03,  1.2731080364155342e-03, -2.9060001721745977e-03,
         2.6174232732491648e-01,  1.4239797833437010e-01, -1.4239797833436521e-01, -7.7179579805168486e-02,
        -6.7537275729838489e-03, -6.7537275729819659e-03,  3.8992664322963539e-03, -3.8992664322969480e-03,
                                                         
         6.6919163252608815e-01,  9.5063998582044021e-02,  3.6607316345626351e-01,  5.3241649763928853e-02,
        -5.0333399450085128e-03, -1.5712542948565542e-02, -2.9060001721746753e-03, -1.2731080364157741e-03,
         1.7423637960024951e+00,  2.5663155896154988e-01,  2.5663155896154877e-01,  3.9269506879060931e-02,
        -1.3304720861538551e-02, -1.3304720861538329e-02, -1.8694838265539103e-03, -1.8694838265541740e-03,
         2.4546088336295506e+00,  3.6679544630615912e-01,  1.5482377206844999e-01,  2.4122215546182844e-02,
        -1.8491615532474651e-02, -1.3118456058247588e-02, -1.1251712079693174e-03, -2.0332549321147595e-03,
         2.8136363432486524e+00,  4.2275290286990563e-01,  5.2063204566974094e-02,  8.1141391202813138e-03,
        -2.1072335114916858e-02, -1.3426476507161618e-02, -3.6480793765665813e-04, -2.0880173238903453e-03,
         2.8136363432486524e+00,  4.2275290286990524e-01, -5.2063204566972256e-02, -8.1141391202823199e-03,
        -2.1072335114916522e-02, -1.3426476507160893e-02,  3.6480793765687319e-04, -2.0880173238910131e-03,
         2.4546088336295533e+00,  3.6679544630615807e-01, -1.5482377206845077e-01, -2.4122215546182161e-02,
        -1.8491615532474186e-02, -1.3118456058246705e-02,  1.1251712079691587e-03, -2.0332549321155132e-03,
         1.7423637960024967e+00,  2.5663155896154877e-01, -2.5663155896155010e-01, -3.9269506879060175e-02,
        -1.3304720861539127e-02, -1.3304720861538711e-02,  1.8694838265534816e-03, -1.8694838265539018e-03,
         6.6919163252608682e-01,  9.5063998582043396e-02, -3.6607316345626262e-01, -5.3241649763929991e-02,
        -5.0333399450070591e-03, -1.5712542948565018e-02,  2.9060001721762860e-03, -1.2731080364161030e-03,
                                                         
         9.3261285108185710e-01,  5.7325647937300842e-02,  5.1325052577225627e-01,  3.1906370575564183e-02,
        -4.7983776055426555e-03, -1.9514998592151105e-02, -2.7703446022324602e-03, -9.2224075298994704e-04,
         2.4546088336295533e+00,  1.5482377206845188e-01,  3.6679544630615785e-01,  2.4122215546181956e-02,
        -1.3118456058246891e-02, -1.8491615532474508e-02, -2.0332549321152001e-03, -1.1251712079690640e-03,
         3.4763948230693575e+00,  2.2273928406436580e-01,  2.2273928406436491e-01,  1.4934123050414328e-02,
        -1.8796563154148561e-02, -1.8796563154148505e-02, -1.2450017281912052e-03, -1.2450017281912329e-03,
         3.9932287149419077e+00,  2.5750472735140589e-01,  7.4989196262071997e-02,  5.0536521254172630e-03,
        -2.1676786288648296e-02, -1.9312217015521736e-02, -4.1789587383843661e-04, -1.3101165429916409e-03,
         3.9932287149419086e+00,  2.5750472735140550e-01, -7.4989196262071983e-02, -5.0536521254171546e-03,
        -2.1676786288648310e-02, -1.9312217015521806e-02,  4.1789587383839465e-04, -1.3101165429915658e-03,
         3.4763948230693589e+00,  2.2273928406436544e-01, -2.2273928406436508e-01, -1.4934123050414191e-02,
        -1.8796563154148464e-02, -1.8796563154148630e-02,  1.2450017281912470e-03, -1.2450017281910922e-03,
         2.4546088336295528e+00,  1.5482377206845233e-01, -3.6679544630615840e-01, -2.4122215546181762e-02,
        -1.3118456058247220e-02, -1.8491615532474616e-02,  2.0332549321149230e-03, -1.1251712079691938e-03,
         9.3261285108185132e-01,  5.7325647937300522e-02, -5.1325052577225905e-01, -3.1906370575564613e-02,
        -4.7983776055399753e-03, -1.9514998592151157e-02,  2.7703446022344595e-03, -9.2224075299001112e-04,
                                                         
         1.0655299621633085e+00,  1.9276481795971488e-02,  5.8724826482558201e-01,  1.0736800514482708e-02,
        -4.9049622075492960e-03, -2.1638935540251314e-02, -2.8318812508928769e-03, -3.0401481573748951e-04,
         2.8136363432486569e+00,  5.2063204566973623e-02,  4.2275290286990502e-01,  8.1141391202819799e-03,
        -1.3426476507160820e-02, -2.1072335114916560e-02, -2.0880173238911450e-03, -3.6480793765673593e-04,
         3.9932287149419117e+00,  7.4989196262073329e-02,  2.5750472735140501e-01,  5.0536521254175961e-03,
        -1.9312217015521948e-02, -2.1676786288648508e-02, -1.3101165429914947e-03, -4.1789587383842935e-04,
         4.5910338085415718e+00,  8.6771187773204306e-02,  8.6771187773202557e-02,  1.7165201446689979e-03,
        -2.2348374449956559e-02, -2.2348374449956337e-02, -4.4280976908144193e-04, -4.4280976908130315e-04,
         4.5910338085415718e+00,  8.6771187773204195e-02, -8.6771187773203362e-02, -1.7165201446691534e-03,
        -2.2348374449956399e-02, -2.2348374449956399e-02,  4.4280976908152368e-04, -4.4280976908143911e-04,
         3.9932287149419095e+00,  7.4989196262072677e-02, -2.5750472735140451e-01, -5.0536521254171502e-03,
        -1.9312217015521778e-02, -2.1676786288648116e-02,  1.3101165429914377e-03, -4.1789587383835627e-04,
         2.8136363432486573e+00,  5.2063204566974289e-02, -4.2275290286990530e-01, -8.1141391202814977e-03,
        -1.3426476507160697e-02, -2.1072335114916685e-02,  2.0880173238911949e-03, -3.6480793765670535e-04,
         1.0655299621633088e+00,  1.9276481795971901e-02, -5.8724826482558179e-01, -1.0736800514483127e-02,
        -4.9049622075491390e-03, -2.1638935540251369e-02,  2.8318812508928517e-03, -3.0401481573744788e-04,
                                                         
         1.0655299621633059e+00, -1.9276481795971124e-02,  5.8724826482558301e-01, -1.0736800514483365e-02,
        -4.9049622075478649e-03, -2.1638935540251245e-02, -2.8318812508938908e-03,  3.0401481573757262e-04,
         2.8136363432486569e+00, -5.2063204566973290e-02,  4.2275290286990536e-01, -8.1141391202814873e-03,
        -1.3426476507160872e-02, -2.1072335114916612e-02, -2.0880173238910075e-03,  3.6480793765668437e-04,
         3.9932287149419130e+00, -7.4989196262072108e-02,  2.5750472735140517e-01, -5.0536521254173029e-03,
        -1.9312217015521708e-02, -2.1676786288648268e-02, -1.3101165429914407e-03,  4.1789587383857528e-04,
         4.5910338085415718e+00, -8.6771187773202793e-02,  8.6771187773202904e-02, -1.7165201446689685e-03,
        -2.2348374449956222e-02, -2.2348374449956222e-02, -4.4280976908148416e-04,  4.4280976908131210e-04,
         4.5910338085415727e+00, -8.6771187773202474e-02, -8.6771187773202779e-02,  1.7165201446694099e-03,
        -2.2348374449956278e-02, -2.2348374449956500e-02,  4.4280976908150861e-04,  4.4280976908139759e-04,
         3.9932287149419117e+00, -7.4989196262071275e-02, -2.5750472735140512e-01,  5.0536521254170332e-03,
        -1.9312217015521788e-02, -2.1676786288648237e-02,  1.3101165429913700e-03,  4.1789587383833242e-04,
         2.8136363432486591e+00, -5.2063204566973408e-02, -4.2275290286990502e-01,  8.1141391202817024e-03,
        -1.3426476507160754e-02, -2.1072335114916716e-02,  2.0880173238911632e-03,  3.6480793765672633e-04,
         1.0655299621633081e+00, -1.9276481795971284e-02, -5.8724826482558301e-01,  1.0736800514483185e-02,
        -4.9049622075484347e-03, -2.1638935540251134e-02,  2.8318812508933353e-03,  3.0401481573757126e-04,
                                                         
         9.3261285108185266e-01, -5.7325647937301265e-02,  5.1325052577225927e-01, -3.1906370575562747e-02,
        -4.7983776055406909e-03, -1.9514998592151191e-02, -2.7703446022337964e-03,  9.2224075298978842e-04,
         2.4546088336295555e+00, -1.5482377206845083e-01,  3.6679544630615885e-01, -2.4122215546182133e-02,
        -1.3118456058247074e-02, -1.8491615532474526e-02, -2.0332549321151190e-03,  1.1251712079692335e-03,
         3.4763948230693624e+00, -2.2273928406436483e-01,  2.2273928406436511e-01, -1.4934123050414449e-02,
        -1.8796563154148675e-02, -1.8796563154148564e-02, -1.2450017281912206e-03,  1.2450017281911187e-03,
         3.9932287149419121e+00, -2.5750472735140495e-01,  7.4989196262071955e-02, -5.0536521254169655e-03,
        -2.1676786288648216e-02, -1.9312217015521767e-02, -4.1789587383831805e-04,  1.3101165429915314e-03,
         3.9932287149419126e+00, -2.5750472735140467e-01, -7.4989196262072025e-02,  5.0536521254170704e-03,
        -2.1676786288648092e-02, -1.9312217015522087e-02,  4.1789587383837931e-04,  1.3101165429914169e-03,
         3.4763948230693611e+00, -2.2273928406436524e-01, -2.2273928406436566e-01,  1.4934123050414329e-02,
        -1.8796563154148464e-02, -1.8796563154148713e-02,  1.2450017281911525e-03,  1.2450017281911802e-03,
         2.4546088336295551e+00, -1.5482377206845158e-01, -3.6679544630615818e-01,  2.4122215546182102e-02,
        -1.3118456058246773e-02, -1.8491615532474696e-02,  2.0332549321152413e-03,  1.1251712079691884e-03,
         9.3261285108185532e-01, -5.7325647937301619e-02, -5.1325052577225805e-01,  3.1906370575563073e-02,
        -4.7983776055419217e-03, -1.9514998592151160e-02,  2.7703446022327573e-03,  9.2224075298971849e-04,
                                                         
         6.6919163252609370e-01, -9.5063998582043147e-02,  3.6607316345626068e-01, -5.3241649763929991e-02,
        -5.0333399450107766e-03, -1.5712542948565014e-02, -2.9060001721731752e-03,  1.2731080364163324e-03,
         1.7423637960024991e+00, -2.5663155896154932e-01,  2.5663155896154954e-01, -3.9269506879060966e-02,
        -1.3304720861538491e-02, -1.3304720861538602e-02, -1.8694838265540792e-03,  1.8694838265539704e-03,
         2.4546088336295564e+00, -3.6679544630615801e-01,  1.5482377206845152e-01, -2.4122215546181613e-02,
        -1.8491615532474564e-02, -1.3118456058247055e-02, -1.1251712079691414e-03,  2.0332549321152114e-03,
         2.8136363432486586e+00, -4.2275290286990402e-01,  5.2063204566973983e-02, -8.1141391202815896e-03,
        -2.1072335114916560e-02, -1.3426476507160459e-02, -3.6480793765670513e-04,  2.0880173238911949e-03,
         2.8136363432486595e+00, -4.2275290286990463e-01, -5.2063204566974060e-02,  8.1141391202815445e-03,
        -2.1072335114916383e-02, -1.3426476507160671e-02,  3.6480793765677414e-04,  2.0880173238913636e-03,
         2.4546088336295546e+00, -3.6679544630615818e-01, -1.5482377206845166e-01,  2.4122215546182019e-02,
        -1.8491615532474439e-02, -1.3118456058246849e-02,  1.1251712079690349e-03,  2.0332549321153233e-03,
         1.7423637960024978e+00, -2.5663155896154966e-01, -2.5663155896154949e-01,  3.9269506879060584e-02,
        -1.3304720861538579e-02, -1.3304720861538607e-02,  1.8694838265540319e-03,  1.8694838265540597e-03,
         6.6919163252609004e-01, -9.5063998582043216e-02, -3.6607316345626228e-01,  5.3241649763930428e-02,
        -5.0333399450094695e-03, -1.5712542948565080e-02,  2.9060001721740326e-03,  1.2731080364163240e-03,
                                                         
         2.6174232732492186e-01, -1.4239797833436652e-01,  1.4239797833436699e-01, -7.7179579805169304e-02,
        -6.7537275729839755e-03, -6.7537275729838663e-03, -3.8992664322955498e-03,  3.8992664322956201e-03,
         6.6919163252609259e-01, -3.6607316345626134e-01,  9.5063998582043535e-02, -5.3241649763929630e-02,
        -1.5712542948564886e-02, -5.0333399450095380e-03, -1.2731080364162786e-03,  2.9060001721740491e-03,
         9.3261285108185321e-01, -5.1325052577226005e-01,  5.7325647937300730e-02, -3.1906370575564294e-02,
        -1.9514998592151306e-02, -4.7983776055403804e-03, -9.2224075298998217e-04,  2.7703446022338254e-03,
         1.0655299621633147e+00, -5.8724826482558012e-01,  1.9276481795972494e-02, -1.0736800514482444e-02,
        -2.1638935540251224e-02, -4.9049622075506205e-03, -3.0401481573732217e-04,  2.8318812508918382e-03,
         1.0655299621633119e+00, -5.8724826482558135e-01, -1.9276481795972668e-02,  1.0736800514482532e-02,
        -2.1638935540251137e-02, -4.9049622075495979e-03,  3.0401481573736000e-04,  2.8318812508924367e-03,
         9.3261285108185554e-01, -5.1325052577225794e-01, -5.7325647937301508e-02,  3.1906370575563600e-02,
        -1.9514998592151313e-02, -4.7983776055423831e-03,  9.2224075298983992e-04,  2.7703446022324342e-03,
         6.6919163252608904e-01, -3.6607316345626240e-01, -9.5063998582043160e-02,  5.3241649763930317e-02,
        -1.5712542948565038e-02, -5.0333399450093533e-03,  1.2731080364162960e-03,  2.9060001721740885e-03,
         2.6174232732492109e-01, -1.4239797833436624e-01, -1.4239797833436613e-01,  7.7179579805169332e-02,
        -6.7537275729838767e-03, -6.7537275729839148e-03,  3.8992664322956569e-03,  3.8992664322956318e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174043031e-01,  1.9147096450810874e-01, -2.8048110764785784e-02, -1.3769366894035592e-02,
        -7.1159030934106020e-03,  2.2827046419852215e-03,  1.8777907284159622e-03,  1.3179201395296886e-03,
         2.6018094650847862e-01,  1.5014541606507584e-01, -2.4623528318145276e-02, -1.1377222883164622e-02,
        -5.4317382757281126e-05,  5.5564633406491573e-04,  2.1992176825353575e-03,  3.2080256054490733e-04,
         2.2061875008365231e-01,  1.3141643256470523e-01,  9.6631517600851008e-03,  5.1140605928803942e-03,
         3.1310264886016711e-03,  6.6601405906152981e-03, -3.6015854094473347e-04,  3.8452339628319973e-03,
         3.1728477488866691e-01,  1.7738358224001630e-01,  3.1338192874284236e-02,  1.2875214739491091e-02,
        -4.4933329941175878e-03, -4.8107149369187568e-03, -4.0417674588017032e-03, -2.7774675638256754e-03,
         2.9856868259017194e-01,  1.6904725884004743e-01, -4.4358254629947795e-02, -1.8966632025180653e-02,
        -2.5805304029333908e-03, -6.5259204756552055e-03,  5.1461245497285120e-03, -3.7677419433303117e-03,
         1.6472836762602577e-01,  1.0615968567242218e-01, -6.4622802488493898e-03, -2.0693364430828479e-03,
         8.5621733627070755e-03,  1.3963870956120502e-02,  1.2871184688642739e-03,  8.0620446554442865e-03,
         3.2567896066607988e-01,  1.8138329229246208e-01,  7.9694807939026877e-02,  3.4130327869920304e-02,
        -5.1491649123717653e-03, -1.2897504262264271e-03, -9.2033633129310552e-03, -7.4463775577016827e-04,
         4.5758071423995317e-01,  2.4387817433520742e-01, -1.7203978611667874e-02, -5.9370449568280142e-03,
        -1.5729094905946293e-02, -1.1872819639637478e-02,  3.0950378831334125e-03, -6.8547756149879950e-03,
                                                         
         9.2555564511589195e-01,  1.4283166346900319e-01, -5.7087203270645710e-02, -4.2107167636395056e-03,
        -6.6585943733794726e-03,  4.3810045656668519e-03,  9.3715478409904616e-04, -1.0646611375156228e-04,
         7.5643742119286139e-01,  1.3194920980432359e-01, -4.3151229134044317e-02, -9.9766889501246210e-04,
        -3.4774021071348250e-03,  2.3686176801918867e-03,  8.9950742716126829e-04,  7.2591693417456003e-04,
         7.0097359003902937e-01,  1.4059912442543343e-01,  2.5969251836348434e-02,  5.4595811792827933e-03,
        -9.8784417722106539e-04,  1.3863694644921702e-02,  5.3783951383761853e-04,  3.1373990954381765e-04,
         8.9395264296783827e-01,  1.6114696623034061e-01,  5.3336355793392402e-02,  4.9642255681518815e-03,
        -1.6240840546830664e-04, -1.1009301949378427e-02, -6.1273948817407056e-05, -8.0128831641327464e-04,
         8.5858667148480283e-01,  1.5815426669652310e-01, -7.7423345159838022e-02, -6.2531566299462403e-03,
         4.2102309020289271e-04, -1.3850860419405031e-02,  3.9811827989687796e-04, -4.6131410499152816e-04,
         6.0862841609240770e-01,  1.3886189996819126e-01, -1.1804099479237351e-02, -3.6257448687969356e-03,
        -1.6305922506315459e-04,  2.8818588277168135e-02, -7.3533836184461665e-04,  5.1433038859834352e-04,
         8.8573073934651869e-01,  1.4220407629839985e-01,  1.3180886462596331e-01,  5.2114498588328763e-03,
        -4.9622366030997230e-03, -2.1500466509340923e-03, -2.0354679892535497e-03,  2.4794549885235803e-04,
         1.1259773284051051e+00,  1.5150225754918686e-01, -2.1648595211938813e-02, -5.4796944887231718e-04,
        -8.3847823258310743e-03, -2.4261809779577561e-02,  5.9460294920793655e-05, -2.9801117729407584e-04,
                                                         
         1.3294757533414434e+00,  8.9788893137442261e-02, -6.2911497442604103e-02, -3.9614170768240779e-05,
        -7.1100499232079391e-03,  4.6679346543523435e-03,  2.4956613702374012e-04,  2.7212527769270804e-04,
         1.1558500674443988e+00,  9.5187518915131650e-02, -3.6801896226935131e-02,  3.8875289024224544e-03,
        -6.1608141969413492e-03,  5.0780393003260101e-03,  2.9847536506069259e-04,  8.3836836755803482e-04,
         1.1608685714961731e+00,  1.2037284313135879e-01,  5.3537968971594352e-02,  1.0607159984611206e-02,
        -4.5111155521692318e-03,  1.5797402263497122e-02,  6.5397859158025575e-04,  8.0268670457471416e-04,
         1.4515046707409818e+00,  1.5802986276729328e-01,  7.3764923340376418e-02,  7.7325221486257015e-03,
        -2.2739467835626078e-03, -1.5570616292196532e-02,  6.3765139919740651e-04, -1.8321877472712517e-03,
         1.4115464704823018e+00,  1.5776040893168652e-01, -9.9620570722380883e-02, -7.8173884689785901e-03,
        -2.1636612802878562e-03, -1.7728464227987777e-02, -5.7397803419402017e-04, -1.7774214977044778e-03,
         1.0684562036485115e+00,  1.2101078886542831e-01, -3.5201132471140130e-02, -9.9394198301684059e-03,
        -4.5077782114690034e-03,  3.1273725348626311e-02, -7.7939850703535667e-04,  9.0314366050543921e-04,
         1.3060178631336954e+00,  9.8310629669391275e-02,  1.3131146816260594e-01, -3.4378558002396307e-03,
        -6.6184622280532198e-03, -5.2033263638921626e-04, -4.3920547811377179e-04,  6.9297032614724626e-04,
         1.5457214754324096e+00,  9.2030059599753342e-02, -2.4079263611516431e-02, -9.9293276550431963e-04,
        -7.4607497917453056e-03, -2.4657338466129757e-02, -4.7089473518925407e-05,  6.9652583640923910e-05,
                                                         
         1.5410607675732244e+00,  3.1259713850908775e-02, -6.1609756908925520e-02,  4.0278881947996382e-04,
        -7.9698971022144766e-03,  5.4767589025950247e-03, -5.1276054205655116e-05,  1.9484961975728570e-04,
         1.3895291108288024e+00,  3.6049179432957129e-02, -2.3292866819651249e-02,  2.8177193660516302e-03,
        -9.0097499591924313e-03,  7.4784863973146746e-03, -5.4908327268816275e-04,  5.4753041006383430e-04,
         1.4814229261146523e+00,  5.3722083116000623e-02,  9.0436045201508081e-02,  7.5759200922297047e-03,
        -1.3014031031221703e-02,  1.8391205195993725e-02, -1.7627894821588457e-03,  6.9484611672705331e-04,
         1.9036696745231794e+00,  8.2014576705607065e-02,  1.0223309369326436e-01,  6.0293706608141534e-03,
        -1.8550660598603910e-02, -2.1202408094434324e-02, -1.4337850889725194e-03, -1.4193287657706831e-03,
         1.8635410509569001e+00,  8.2038435164303550e-02, -1.2804640133826681e-01, -6.0050462696767748e-03,
        -1.8554396836116659e-02, -2.3251244190818931e-02,  1.4316279712387684e-03, -1.4111570005111908e-03,
         1.3908819207316627e+00,  5.4053805794842834e-02, -7.0885965719392860e-02, -7.4371696780804368e-03,
        -1.3096490503921039e-02,  3.4035995955905070e-02,  1.7194957188661293e-03,  6.9165401818151314e-04,
         1.5458079707862993e+00,  3.6785251742612025e-02,  1.1844442680216077e-01, -2.7477367760387942e-03,
        -9.2110386667161315e-03,  1.5895856233272834e-03,  5.2377094526747942e-04,  5.2519154906815752e-04,
         1.7614161838416025e+00,  3.1685800667273945e-02, -2.7278574910696574e-02, -6.3584621477922470e-04,
        -8.0924625745515209e-03, -2.4216573332765313e-02,  1.2203926265282081e-04,  1.8482328475642772e-04,
                                                         
         1.5410607675732240e+00, -3.1259713850908713e-02, -6.1609756908925596e-02, -4.0278881947999954e-04,
        -7.9698971022143933e-03,  5.4767589025950534e-03, -5.1276054205712274e-05, -1.9484961975725637e-04,
         1.3895291108288028e+00, -3.6049179432956852e-02, -2.3292866819651072e-02, -2.8177193660514754e-03,
        -9.0097499591924035e-03,  7.4784863973146321e-03, -5.4908327268812123e-04, -5.4753041006383484e-04,
         1.4814229261146523e+00, -5.3722083116000553e-02,  9.0436045201507997e-02, -7.5759200922298131e-03,
        -1.3014031031221738e-02,  1.8391205195993847e-02, -1.7627894821589177e-03, -6.9484611672699509e-04,
         1.9036696745231787e+00, -8.2014576705607023e-02,  1.0223309369326418e-01, -6.0293706608139782e-03,
        -1.8550660598603833e-02, -2.1202408094434275e-02, -1.4337850889724136e-03,  1.4193287657707004e-03,
         1.8635410509568999e+00, -8.2038435164303550e-02, -1.2804640133826664e-01,  6.0050462696767818e-03,
        -1.8554396836116590e-02, -2.3251244190819000e-02,  1.4316279712386628e-03,  1.4111570005111837e-03,
         1.3908819207316627e+00, -5.4053805794842855e-02, -7.0885965719392902e-02,  7.4371696780805912e-03,
        -1.3096490503921063e-02,  3.4035995955905098e-02,  1.7194957188662130e-03, -6.9165401818148885e-04,
         1.5458079707862988e+00, -3.6785251742611901e-02,  1.1844442680216073e-01,  2.7477367760388098e-03,
        -9.2110386667160222e-03,  1.5895856233272539e-03,  5.2377094526745579e-04, -5.2519154906813107e-04,
         1.7614161838416023e+00, -3.1685800667273778e-02, -2.7278574910696542e-02,  6.3584621477932846e-04,
        -8.0924625745514220e-03, -2.4216573332765355e-02,  1.2203926265283721e-04, -1.8482328475641070e-04,
                                                         
         1.3294757533414430e+00, -8.9788893137442191e-02, -6.2911497442604464e-02,  3.9614170768292408e-05,
        -7.1100499232079313e-03,  4.6679346543524623e-03,  2.4956613702379482e-04, -2.7212527769268630e-04,
         1.1558500674443986e+00, -9.5187518915131830e-02, -3.6801896226934798e-02, -3.8875289024225706e-03,
        -6.1608141969413084e-03,  5.0780393003259120e-03,  2.9847536506066310e-04, -8.3836836755803027e-04,
         1.1608685714961731e+00, -1.2037284313135867e-01,  5.3537968971594171e-02, -1.0607159984610976e-02,
        -4.5111155521692517e-03,  1.5797402263497157e-02,  6.5397859158026453e-04, -8.0268670457479710e-04,
         1.4515046707409809e+00, -1.5802986276729319e-01,  7.3764923340376501e-02, -7.7325221486257874e-03,
        -2.2739467835625805e-03, -1.5570616292196462e-02,  6.3765139919740119e-04,  1.8321877472712586e-03,
         1.4115464704823009e+00, -1.5776040893168661e-01, -9.9620570722381063e-02,  7.8173884689786872e-03,
        -2.1636612802878735e-03, -1.7728464227987822e-02, -5.7397803419405682e-04,  1.7774214977044585e-03,
         1.0684562036485103e+00, -1.2101078886542833e-01, -3.5201132471139970e-02,  9.9394198301683799e-03,
        -4.5077782114689643e-03,  3.1273725348626304e-02, -7.7939850703529119e-04, -9.0314366050547369e-04,
         1.3060178631336949e+00, -9.8310629669391220e-02,  1.3131146816260625e-01,  3.4378558002397842e-03,
        -6.6184622280531843e-03, -5.2033263638922255e-04, -4.3920547811384877e-04, -6.9297032614728052e-04,
         1.5457214754324098e+00, -9.2030059599753106e-02, -2.4079263611516445e-02,  9.9293276550434868e-04,
        -7.4607497917453776e-03, -2.4657338466129785e-02, -4.7089473518925374e-05, -6.9652583640937327e-05,
                                                         
         9.2555564511589228e-01, -1.4283166346900242e-01, -5.7087203270645939e-02,  4.2107167636395108e-03,
        -6.6585943733794293e-03,  4.3810045656671034e-03,  9.3715478409907424e-04,  1.0646611375162348e-04,
         7.5643742119286161e-01, -1.3194920980432301e-01, -4.3151229134044268e-02,  9.9766889501250069e-04,
        -3.4774021071347677e-03,  2.3686176801920480e-03,  8.9950742716123533e-04, -7.2591693417441388e-04,
         7.0097359003902948e-01, -1.4059912442543290e-01,  2.5969251836348382e-02, -5.4595811792828132e-03,
        -9.8784417722099275e-04,  1.3863694644921844e-02,  5.3783951383765062e-04, -3.1373990954366011e-04,
         8.9395264296783794e-01, -1.6114696623034006e-01,  5.3336355793392118e-02, -4.9642255681519613e-03,
        -1.6240840546821416e-04, -1.1009301949378315e-02, -6.1273948817443689e-05,  8.0128831641332842e-04,
         8.5858667148480194e-01, -1.5815426669652288e-01, -7.7423345159838203e-02,  6.2531566299462507e-03,
         4.2102309020295776e-04, -1.3850860419404911e-02,  3.9811827989688620e-04,  4.6131410499162032e-04,
         6.0862841609240648e-01, -1.3886189996819115e-01, -1.1804099479237122e-02,  3.6257448687969642e-03,
        -1.6305922506312725e-04,  2.8818588277168312e-02, -7.3533836184465481e-04, -5.1433038859821081e-04,
         8.8573073934651858e-01, -1.4220407629839935e-01,  1.3180886462596375e-01, -5.2114498588326360e-03,
        -4.9622366030995808e-03, -2.1500466509339223e-03, -2.0354679892534604e-03, -2.4794549885219248e-04,
         1.1259773284051056e+00, -1.5150225754918639e-01, -2.1648595211938699e-02,  5.4796944887228096e-04,
        -8.3847823258309216e-03, -2.4261809779577131e-02,  5.9460294920692505e-05,  2.9801117729436375e-04,
                                                         
         3.4754908174043120e-01, -1.9147096450810924e-01, -2.8048110764785136e-02,  1.3769366894036281e-02,
        -7.1159030934106497e-03,  2.2827046419849114e-03,  1.8777907284160782e-03, -1.3179201395300791e-03,
         2.6018094650848222e-01, -1.5014541606507417e-01, -2.4623528318145321e-02,  1.1377222883164596e-02,
        -5.4317382757140756e-05,  5.5564633406385939e-04,  2.1992176825353384e-03, -3.2080256054575182e-04,
         2.2061875008365400e-01, -1.3141643256470487e-01,  9.6631517600846636e-03, -5.1140605928807351e-03,
         3.1310264886016273e-03,  6.6601405906146901e-03, -3.6015854094482351e-04, -3.8452339628325693e-03,
         3.1728477488866830e-01, -1.7738358224001571e-01,  3.1338192874284528e-02, -1.2875214739490477e-02,
        -4.4933329941175626e-03, -4.8107149369193128e-03, -4.0417674588015843e-03,  2.7774675638252209e-03,
         2.9856868259017150e-01, -1.6904725884004720e-01, -4.4358254629948746e-02,  1.8966632025180070e-02,
        -2.5805304029333561e-03, -6.5259204756554336e-03,  5.1461245497283941e-03,  3.7677419433300198e-03,
         1.6472836762602511e-01, -1.0615968567242212e-01, -6.4622802488490489e-03,  2.0693364430829689e-03,
         8.5621733627069246e-03,  1.3963870956119893e-02,  1.2871184688642873e-03, -8.0620446554448624e-03,
         3.2567896066608065e-01, -1.8138329229246203e-01,  7.9694807939027446e-02, -3.4130327869920409e-02,
        -5.1491649123716690e-03, -1.2897504262269553e-03, -9.2033633129309355e-03,  7.4463775576960654e-04,
         4.5758071423995794e-01, -2.4387817433520531e-01, -1.7203978611668343e-02,  5.9370449568276907e-03,
        -1.5729094905946251e-02, -1.1872819639639669e-02,  3.0950378831332307e-03,  6.8547756149861657e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174042964e-01, -2.8048110764784386e-02,  1.9147096450811057e-01, -1.3769366894036857e-02,  2.2827046419850081e-03, -7.1159030934107581e-03,  1.3179201395300958e-03,  1.8777907284162560e-03,
         9.2555564511589383e-01, -5.7087203270645925e-02,  1.4283166346900303e-01, -4.2107167636395723e-03,  4.3810045656670115e-03, -6.6585943733795568e-03, -1.0646611375170161e-04,  9.3715478409909473e-04,
         1.3294757533414452e+00, -6.2911497442604436e-02,  8.9788893137442344e-02, -3.9614170768272533e-05,  4.6679346543522706e-03, -7.1100499232081170e-03,  2.7212527769271232e-04,  2.4956613702375823e-04,
         1.5410607675732264e+00, -6.1609756908925783e-02,  3.1259713850908692e-02,  4.0278881947991574e-04,  5.4767589025949666e-03, -7.9698971022144141e-03,  1.9484961975729643e-04, -5.1276054205724837e-05,
         1.5410607675732255e+00, -6.1609756908925964e-02, -3.1259713850908942e-02, -4.0278881947993422e-04,  5.4767589025949874e-03, -7.9698971022144575e-03, -1.9484961975728526e-04, -5.1276054205658064e-05,
         1.3294757533414445e+00, -6.2911497442604505e-02, -8.9788893137442580e-02,  3.9614170768346550e-05,  4.6679346543523348e-03, -7.1100499232081196e-03, -2.7212527769269334e-04,  2.4956613702381553e-04,
         9.2555564511589206e-01, -5.7087203270645800e-02, -1.4283166346900306e-01,  4.2107167636395082e-03,  4.3810045656669195e-03, -6.6585943733794492e-03,  1.0646611375155487e-04,  9.3715478409906871e-04,
         3.4754908174042604e-01, -2.8048110764785483e-02, -1.9147096450811166e-01,  1.3769366894035920e-02,  2.2827046419861504e-03, -7.1159030934108405e-03, -1.3179201395291872e-03,  1.8777907284159925e-03,
                                                         
         2.6018094650848411e-01, -2.4623528318145346e-02,  1.5014541606507359e-01, -1.1377222883164428e-02,  5.5564633406300114e-04, -5.4317382757119831e-05,  3.2080256054629734e-04,  2.1992176825352508e-03,
         7.5643742119286284e-01, -4.3151229134044358e-02,  1.3194920980432326e-01, -9.9766889501257138e-04,  2.3686176801920706e-03, -3.4774021071348232e-03,  7.2591693417436444e-04,  8.9950742716123674e-04,
         1.1558500674444006e+00, -3.6801896226935180e-02,  9.5187518915131983e-02,  3.8875289024226109e-03,  5.0780393003259051e-03, -6.1608141969413475e-03,  8.3836836755805694e-04,  2.9847536506072300e-04,
         1.3895291108288046e+00, -2.3292866819651030e-02,  3.6049179432956914e-02,  2.8177193660515751e-03,  7.4784863973144733e-03, -9.0097499591925388e-03,  5.4753041006382476e-04, -5.4908327268816091e-04,
         1.3895291108288044e+00, -2.3292866819650922e-02, -3.6049179432957101e-02, -2.8177193660515378e-03,  7.4784863973144950e-03, -9.0097499591925163e-03, -5.4753041006382487e-04, -5.4908327268817989e-04,
         1.1558500674444003e+00, -3.6801896226934715e-02, -9.5187518915131886e-02, -3.8875289024224817e-03,  5.0780393003259268e-03, -6.1608141969413258e-03, -8.3836836755802658e-04,  2.9847536506067546e-04,
         7.5643742119286161e-01, -4.3151229134043845e-02, -1.3194920980432370e-01,  9.9766889501253452e-04,  2.3686176801920966e-03, -3.4774021071348306e-03, -7.2591693417442061e-04,  8.9950742716124173e-04,
         2.6018094650848284e-01, -2.4623528318144457e-02, -1.5014541606507342e-01,  1.1377222883164968e-02,  5.5564633406298997e-04, -5.4317382757213696e-05, -3.2080256054627907e-04,  2.1992176825354850e-03,
                                                         
         2.2061875008365345e-01,  9.6631517600846272e-03,  1.3141643256470570e-01,  5.1140605928805252e-03,  6.6601405906149182e-03,  3.1310264886015996e-03,  3.8452339628324505e-03, -3.6015854094478215e-04,
         7.0097359003903059e-01,  2.5969251836348413e-02,  1.4059912442543326e-01,  5.4595811792830344e-03,  1.3863694644921848e-02, -9.8784417722103503e-04,  3.1373990954367497e-04,  5.3783951383765680e-04,
         1.1608685714961748e+00,  5.3537968971594505e-02,  1.2037284313135892e-01,  1.0607159984610891e-02,  1.5797402263497132e-02, -4.5111155521692725e-03,  8.0268670457474148e-04,  6.5397859158021802e-04,
         1.4814229261146545e+00,  9.0436045201507859e-02,  5.3722083116000588e-02,  7.5759200922298322e-03,  1.8391205195993732e-02, -1.3014031031221736e-02,  6.9484611672702631e-04, -1.7627894821588743e-03,
         1.4814229261146548e+00,  9.0436045201508108e-02, -5.3722083116000581e-02, -7.5759200922297211e-03,  1.8391205195993788e-02, -1.3014031031221748e-02, -6.9484611672701677e-04, -1.7627894821588838e-03,
         1.1608685714961748e+00,  5.3537968971594407e-02, -1.2037284313135886e-01, -1.0607159984611132e-02,  1.5797402263497066e-02, -4.5111155521693393e-03, -8.0268670457477877e-04,  6.5397859158024664e-04,
         7.0097359003903081e-01,  2.5969251836348378e-02, -1.4059912442543340e-01, -5.4595811792828314e-03,  1.3863694644921699e-02, -9.8784417722108534e-04, -3.1373990954372105e-04,  5.3783951383762341e-04,
         2.2061875008365178e-01,  9.6631517600844936e-03, -1.3141643256470645e-01, -5.1140605928809068e-03,  6.6601405906156555e-03,  3.1310264886016746e-03, -3.8452339628318906e-03, -3.6015854094490277e-04,
                                                         
         3.1728477488866852e-01,  3.1338192874284548e-02,  1.7738358224001591e-01,  1.2875214739490581e-02, -4.8107149369194429e-03, -4.4933329941175887e-03, -2.7774675638251212e-03, -4.0417674588016208e-03,
         8.9395264296783927e-01,  5.3336355793392055e-02,  1.6114696623034069e-01,  4.9642255681518581e-03, -1.1009301949378283e-02, -1.6240840546819605e-04, -8.0128831641333893e-04, -6.1273948817444122e-05,
         1.4515046707409838e+00,  7.3764923340376473e-02,  1.5802986276729358e-01,  7.7325221486260034e-03, -1.5570616292196622e-02, -2.2739467835626980e-03, -1.8321877472712773e-03,  6.3765139919743480e-04,
         1.9036696745231823e+00,  1.0223309369326505e-01,  8.2014576705607287e-02,  6.0293706608141795e-03, -2.1202408094434320e-02, -1.8550660598603982e-02, -1.4193287657706390e-03, -1.4337850889724951e-03,
         1.9036696745231827e+00,  1.0223309369326498e-01, -8.2014576705607231e-02, -6.0293706608142350e-03, -2.1202408094434320e-02, -1.8550660598604118e-02,  1.4193287657706867e-03, -1.4337850889725031e-03,
         1.4515046707409842e+00,  7.3764923340376751e-02, -1.5802986276729358e-01, -7.7325221486257995e-03, -1.5570616292196521e-02, -2.2739467835627301e-03,  1.8321877472712252e-03,  6.3765139919739501e-04,
         8.9395264296783961e-01,  5.3336355793392465e-02, -1.6114696623034064e-01, -4.9642255681518858e-03, -1.1009301949378327e-02, -1.6240840546830605e-04,  8.0128831641336322e-04, -6.1273948817420066e-05,
         3.1728477488866941e-01,  3.1338192874284715e-02, -1.7738358224001569e-01, -1.2875214739490671e-02, -4.8107149369198323e-03, -4.4933329941176468e-03,  2.7774675638248801e-03, -4.0417674588015757e-03,
                                                         
         2.9856868259017222e-01, -4.4358254629948406e-02,  1.6904725884004743e-01, -1.8966632025180216e-02, -6.5259204756554969e-03, -2.5805304029333379e-03, -3.7677419433299738e-03,  5.1461245497284617e-03,
         8.5858667148480361e-01, -7.7423345159837897e-02,  1.5815426669652335e-01, -6.2531566299461995e-03, -1.3850860419404974e-02,  4.2102309020287818e-04, -4.6131410499162384e-04,  3.9811827989688772e-04,
         1.4115464704823038e+00, -9.9620570722380813e-02,  1.5776040893168708e-01, -7.8173884689787462e-03, -1.7728464227987965e-02, -2.1636612802880054e-03, -1.7774214977044575e-03, -5.7397803419405899e-04,
         1.8635410509569048e+00, -1.2804640133826678e-01,  8.2038435164304063e-02, -6.0050462696768147e-03, -2.3251244190819049e-02, -1.8554396836116767e-02, -1.4111570005111832e-03,  1.4316279712387229e-03,
         1.8635410509569053e+00, -1.2804640133826684e-01, -8.2038435164303827e-02,  6.0050462696768702e-03, -2.3251244190819070e-02, -1.8554396836116722e-02,  1.4111570005111642e-03,  1.4316279712387932e-03,
         1.4115464704823049e+00, -9.9620570722380869e-02, -1.5776040893168702e-01,  7.8173884689786716e-03, -1.7728464227987951e-02, -2.1636612802879928e-03,  1.7774214977044958e-03, -5.7397803419400651e-04,
         8.5858667148480450e-01, -7.7423345159837939e-02, -1.5815426669652330e-01,  6.2531566299462082e-03, -1.3850860419404939e-02,  4.2102309020288029e-04,  4.6131410499163864e-04,  3.9811827989687552e-04,
         2.9856868259017277e-01, -4.4358254629948260e-02, -1.6904725884004765e-01,  1.8966632025180313e-02, -6.5259204756555290e-03, -2.5805304029333370e-03,  3.7677419433299339e-03,  5.1461245497284478e-03,
                                                         
         1.6472836762602650e-01, -6.4622802488492528e-03,  1.0615968567242214e-01, -2.0693364430828132e-03,  1.3963870956119633e-02,  8.5621733627069419e-03,  8.0620446554450047e-03,  1.2871184688642244e-03,
         6.0862841609240859e-01, -1.1804099479237192e-02,  1.3886189996819157e-01, -3.6257448687969846e-03,  2.8818588277168294e-02, -1.6305922506323095e-04,  5.1433038859822794e-04, -7.3533836184467519e-04,
         1.0684562036485135e+00, -3.5201132471140248e-02,  1.2101078886542867e-01, -9.9394198301683955e-03,  3.1273725348626248e-02, -4.5077782114690086e-03,  9.0314366050544658e-04, -7.7939850703530583e-04,
         1.3908819207316665e+00, -7.0885965719393249e-02,  5.4053805794843257e-02, -7.4371696780805782e-03,  3.4035995955904931e-02, -1.3096490503921115e-02,  6.9165401818148148e-04,  1.7194957188661963e-03,
         1.3908819207316667e+00, -7.0885965719393138e-02, -5.4053805794842896e-02,  7.4371696780806146e-03,  3.4035995955904966e-02, -1.3096490503921015e-02, -6.9165401818150035e-04,  1.7194957188661772e-03,
         1.0684562036485146e+00, -3.5201132471140234e-02, -1.2101078886542882e-01,  9.9394198301683487e-03,  3.1273725348626262e-02, -4.5077782114689956e-03, -9.0314366050543834e-04, -7.7939850703529628e-04,
         6.0862841609240925e-01, -1.1804099479237469e-02, -1.3886189996819157e-01,  3.6257448687969478e-03,  2.8818588277168287e-02, -1.6305922506316907e-04, -5.1433038859822100e-04, -7.3533836184460841e-04,
         1.6472836762602747e-01, -6.4622802488492458e-03, -1.0615968567242191e-01,  2.0693364430829832e-03,  1.3963870956119570e-02,  8.5621733627069974e-03, -8.0620446554450394e-03,  1.2871184688642598e-03,
                                                         
         3.2567896066608137e-01,  7.9694807939027321e-02,  1.8138329229246233e-01,  3.4130327869920436e-02, -1.2897504262272103e-03, -5.1491649123716135e-03, -7.4463775576943090e-04, -9.2033633129308540e-03,
         8.8573073934652091e-01,  1.3180886462596381e-01,  1.4220407629839993e-01,  5.2114498588327497e-03, -2.1500466509339136e-03, -4.9622366030995175e-03,  2.4794549885218521e-04, -2.0354679892533663e-03,
         1.3060178631336981e+00,  1.3131146816260650e-01,  9.8310629669391358e-02, -3.4378558002398189e-03, -5.2033263638920748e-04, -6.6184622280531470e-03,  6.9297032614725450e-04, -4.3920547811381597e-04,
         1.5458079707863022e+00,  1.1844442680216062e-01,  3.6785251742611789e-02, -2.7477367760390123e-03,  1.5895856233272446e-03, -9.2110386667161575e-03,  5.2519154906816750e-04,  5.2377094526740299e-04,
         1.5458079707863026e+00,  1.1844442680216057e-01, -3.6785251742611727e-02,  2.7477367760389013e-03,  1.5895856233273001e-03, -9.2110386667162338e-03, -5.2519154906811025e-04,  5.2377094526734563e-04,
         1.3060178631336978e+00,  1.3131146816260572e-01, -9.8310629669391830e-02,  3.4378558002395136e-03, -5.2033263638935265e-04, -6.6184622280532918e-03, -6.9297032614738189e-04, -4.3920547811387062e-04,
         8.8573073934651991e-01,  1.3180886462596297e-01, -1.4220407629839971e-01, -5.2114498588324946e-03, -2.1500466509341257e-03, -4.9622366030995643e-03, -2.4794549885213566e-04, -2.0354679892534569e-03,
         3.2567896066608099e-01,  7.9694807939026988e-02, -1.8138329229246228e-01, -3.4130327869920402e-02, -1.2897504262269091e-03, -5.1491649123716274e-03,  7.4463775576965577e-04, -9.2033633129309320e-03,
                                                         
         4.5758071423995789e-01, -1.7203978611669109e-02,  2.4387817433520595e-01, -5.9370449568272301e-03, -1.1872819639640254e-02, -1.5729094905946386e-02, -6.8547756149857598e-03,  3.0950378831330659e-03,
         1.1259773284051073e+00, -2.1648595211938900e-02,  1.5150225754918717e-01, -5.4796944887231913e-04, -2.4261809779577124e-02, -8.3847823258308921e-03, -2.9801117729441520e-04,  5.9460294920611305e-05,
         1.5457214754324125e+00, -2.4079263611516789e-02,  9.2030059599753189e-02, -9.9293276550425479e-04, -2.4657338466129764e-02, -7.4607497917453663e-03,  6.9652583640970205e-05, -4.7089473518950642e-05,
         1.7614161838416049e+00, -2.7278574910696709e-02,  3.1685800667273785e-02, -6.3584621477909753e-04, -2.4216573332765473e-02, -8.0924625745515469e-03,  1.8482328475636194e-04,  1.2203926265291815e-04,
         1.7614161838416049e+00, -2.7278574910696782e-02, -3.1685800667273903e-02,  6.3584621477904201e-04, -2.4216573332765518e-02, -8.0924625745517915e-03, -1.8482328475635242e-04,  1.2203926265290857e-04,
         1.5457214754324113e+00, -2.4079263611516525e-02, -9.2030059599753952e-02,  9.9293276550453235e-04, -2.4657338466129847e-02, -7.4607497917455797e-03, -6.9652583640970558e-05, -4.7089473518941840e-05,
         1.1259773284051051e+00, -2.1648595211938775e-02, -1.5150225754918667e-01,  5.4796944887201382e-04, -2.4261809779577238e-02, -8.3847823258309407e-03,  2.9801117729432342e-04,  5.9460294920663889e-05,
         4.5758071423995744e-01, -1.7203978611668745e-02, -2.4387817433520526e-01,  5.9370449568277878e-03, -1.1872819639639924e-02, -1.5729094905946272e-02,  6.8547756149861354e-03,  3.0950378831332199e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(perbuff);
}

void
gpu_test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs)
{
  double epsilon_0 = 1.0;
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
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_dirichlety, NULL);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest-gpu_fem_poisson_2x_rho_1.gkyl");
  gkyl_array_copy(rho_cu, rho);

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif

  // FEM poisson solver.
  gkyl_fem_poisson *poisson_cu = gkyl_fem_poisson_new(&grid, basis, bcs, epsilon_0, NULL, true);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson_cu, rho_cu);

  // Solve the problem.
  gkyl_fem_poisson_solve(poisson_cu, phi_cu);
  gkyl_array_copy(phi, phi_cu);

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif

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
    gkyl_array_shiftc0(phi, mavgfac*sol_avg[0]);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest-gpu_fem_poisson_2x_phi_1.gkyl");

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
    gkyl_array_shiftc0(phi, mavgfac*sol_avg[0]);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
	-7.6861958028898086e-02, -4.4376272158426031e-02, 4.4376272158426031e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426100e-02, 4.4376272158426051e-02, 2.5620652676299389e-02,
	7.6861958028898059e-02, 4.4376272158426072e-02, -4.4376272158426079e-02, -2.5620652676299396e-02,
	-7.6861958028898225e-02, -4.4376272158426079e-02, -4.4376272158426065e-02, -2.5620652676299392e-02,
	-7.6861958028898170e-02, -4.4376272158426079e-02, 4.4376272158426100e-02, 2.5620652676299396e-02,
	7.6861958028898197e-02, 4.4376272158426065e-02, 4.4376272158426086e-02, 2.5620652676299392e-02,
	7.6861958028898253e-02, 4.4376272158426079e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898031e-02, -4.4376272158426031e-02, -4.4376272158426072e-02, -2.5620652676299389e-02,
	-1.8556118150391737e-01, -1.8381253775579646e-02, 1.0713379809243175e-01, 1.0612421815373736e-02,
	1.8556118150391748e-01, 1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373724e-02,
	1.8556118150391740e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373727e-02,
	-1.8556118150391759e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, -1.0612421815373736e-02,
	-1.8556118150391751e-01, -1.8381253775579632e-02, 1.0713379809243183e-01, 1.0612421815373741e-02,
	1.8556118150391762e-01, 1.8381253775579681e-02, 1.0713379809243183e-01, 1.0612421815373745e-02,
	1.8556118150391770e-01, 1.8381253775579684e-02, -1.0713379809243177e-01, -1.0612421815373741e-02,
	-1.8556118150391729e-01, -1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373745e-02,
	-1.8556118150391737e-01, 1.8381253775579646e-02, 1.0713379809243175e-01, -1.0612421815373736e-02,
	1.8556118150391746e-01, -1.8381253775579646e-02, 1.0713379809243173e-01, -1.0612421815373734e-02,
	1.8556118150391737e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, 1.0612421815373736e-02,
	-1.8556118150391757e-01, 1.8381253775579663e-02, -1.0713379809243176e-01, 1.0612421815373745e-02,
	-1.8556118150391748e-01, 1.8381253775579656e-02, 1.0713379809243181e-01, -1.0612421815373750e-02,
	1.8556118150391762e-01, -1.8381253775579656e-02, 1.0713379809243183e-01, -1.0612421815373741e-02,
	1.8556118150391770e-01, -1.8381253775579670e-02, -1.0713379809243177e-01, 1.0612421815373731e-02,
	-1.8556118150391729e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, 1.0612421815373745e-02,
	-7.6861958028898114e-02, 4.4376272158426017e-02, 4.4376272158426044e-02, -2.5620652676299378e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426031e-02, -2.5620652676299378e-02,
	7.6861958028898045e-02, -4.4376272158426058e-02, -4.4376272158426065e-02, 2.5620652676299389e-02,
	-7.6861958028898197e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, 2.5620652676299378e-02,
	-7.6861958028898142e-02, 4.4376272158426051e-02, 4.4376272158426079e-02, -2.5620652676299389e-02,
	7.6861958028898197e-02, -4.4376272158426107e-02, 4.4376272158426100e-02, -2.5620652676299396e-02,
	7.6861958028898253e-02, -4.4376272158426114e-02, -4.4376272158426065e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426017e-02, -4.4376272158426072e-02, 2.5620652676299378e-02,
	7.6861958028898031e-02, 4.4376272158426031e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426051e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426044e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426065e-02, 4.4376272158426065e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898114e-02, -4.4376272158426058e-02, -4.4376272158426051e-02, -2.5620652676299392e-02,
	-7.6861958028898114e-02, -4.4376272158426079e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426017e-02, 4.4376272158426031e-02, 2.5620652676299378e-02,
	1.8556118150391737e-01, 1.8381253775579677e-02, -1.0713379809243177e-01, -1.0612421815373757e-02,
	-1.8556118150391754e-01, -1.8381253775579670e-02, -1.0713379809243175e-01, -1.0612421815373741e-02,
	-1.8556118150391748e-01, -1.8381253775579663e-02, 1.0713379809243177e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, 1.8381253775579656e-02, 1.0713379809243177e-01, 1.0612421815373731e-02,
	1.8556118150391748e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, -1.0612421815373731e-02,
	-1.8556118150391740e-01, -1.8381253775579625e-02, -1.0713379809243176e-01, -1.0612421815373736e-02,
	-1.8556118150391743e-01, -1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373731e-02,
	1.8556118150391743e-01, 1.8381253775579687e-02, 1.0713379809243175e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, -1.8381253775579642e-02, -1.0713379809243181e-01, 1.0612421815373738e-02,
	-1.8556118150391765e-01, 1.8381253775579615e-02, -1.0713379809243181e-01, 1.0612421815373708e-02,
	-1.8556118150391762e-01, 1.8381253775579601e-02, 1.0713379809243183e-01, -1.0612421815373717e-02,
	1.8556118150391746e-01, -1.8381253775579642e-02, 1.0713379809243180e-01, -1.0612421815373722e-02,
	1.8556118150391751e-01, -1.8381253775579632e-02, -1.0713379809243176e-01, 1.0612421815373727e-02,
	-1.8556118150391737e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, 1.0612421815373750e-02,
	-1.8556118150391740e-01, 1.8381253775579663e-02, 1.0713379809243173e-01, -1.0612421815373745e-02,
	1.8556118150391748e-01, -1.8381253775579642e-02, 1.0713379809243177e-01, -1.0612421815373736e-02,
	7.6861958028898114e-02, -4.4376272158426051e-02, -4.4376272158426100e-02, 2.5620652676299389e-02,
	-7.6861958028898281e-02, 4.4376272158426114e-02, -4.4376272158426114e-02, 2.5620652676299406e-02,
	-7.6861958028898308e-02, 4.4376272158426107e-02, 4.4376272158426107e-02, -2.5620652676299410e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426086e-02, -2.5620652676299396e-02,
	7.6861958028898156e-02, -4.4376272158426079e-02, -4.4376272158426051e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426044e-02, -4.4376272158426031e-02, 2.5620652676299378e-02,
	-7.6861958028898059e-02, 4.4376272158426058e-02, 4.4376272158426044e-02, -2.5620652676299368e-02,
	7.6861958028898170e-02, -4.4376272158426051e-02, 4.4376272158426058e-02, -2.5620652676299389e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        2.3681559834188498e-01,  1.3672554945099072e-01,  1.3672554945098905e-01,  7.8938532780630122e-02,
        6.3388180537571182e-01,  3.6597183096807612e-01,  9.2520732066094527e-02,  5.3416869563979862e-02,
        8.9067264750613684e-01,  5.1423009279750498e-01,  5.5737529763337020e-02,  3.2180077812826891e-02,
        1.0196938687973565e+00,  5.8872052964116539e-01,  1.8752907080322980e-02,  1.0826995950913283e-02,
        1.0196938687973551e+00,  5.8872052964116617e-01, -1.8752907080323782e-02, -1.0826995950912827e-02,
        8.9067264750613595e-01,  5.1423009279750564e-01, -5.5737529763335895e-02, -3.2180077812827390e-02,
        6.3388180537571470e-01,  3.6597183096807495e-01, -9.2520732066093542e-02, -5.3416869563980424e-02,
        2.3681559834188717e-01,  1.3672554945099089e-01, -1.3672554945099044e-01, -7.8938532780628831e-02,
                                                        
        6.3388180537571459e-01,  9.2520732066094485e-02,  3.6597183096807473e-01,  5.3416869563980014e-02,
        1.7090814385197843e+00,  2.5479496666024115e-01,  2.5479496666024093e-01,  4.0272203461491046e-02,
        2.4167368175137374e+00,  3.6684346655702615e-01,  1.5377005689540441e-01,  2.4419028116211553e-02,
        2.7726579759463506e+00,  4.2335376950105846e-01,  5.1721119735948180e-02,  8.2072105005127904e-03,
        2.7726579759463510e+00,  4.2335376950105857e-01, -5.1721119735948097e-02, -8.2072105005127401e-03,
        2.4167368175137383e+00,  3.6684346655702638e-01, -1.5377005689540432e-01, -2.4419028116211581e-02,
        1.7090814385197852e+00,  2.5479496666024121e-01, -2.5479496666024098e-01, -4.0272203461491102e-02,
        6.3388180537571515e-01,  9.2520732066093431e-02, -3.6597183096807484e-01, -5.3416869563980611e-02,
                                                        
        8.9067264750613706e-01,  5.5737529763335694e-02,  5.1423009279750498e-01,  3.2180077812827605e-02,
        2.4167368175137374e+00,  1.5377005689540435e-01,  3.6684346655702599e-01,  2.4419028116211518e-02,
        3.4373208199425331e+00,  2.2239098197586474e-01,  2.2239098197586490e-01,  1.5199281451033411e-02,
        3.9521264825956046e+00,  2.5761269031357442e-01,  7.4832205937229848e-02,  5.1359813390617843e-03,
        3.9521264825956055e+00,  2.5761269031357470e-01, -7.4832205937229418e-02, -5.1359813390615995e-03,
        3.4373208199425345e+00,  2.2239098197586496e-01, -2.2239098197586502e-01, -1.5199281451033622e-02,
        2.4167368175137383e+00,  1.5377005689540424e-01, -3.6684346655702615e-01, -2.4419028116211505e-02,
        8.9067264750613750e-01,  5.5737529763336638e-02, -5.1423009279750520e-01, -3.2180077812827015e-02,
                                                        
        1.0196938687973565e+00,  1.8752907080324045e-02,  5.8872052964116561e-01,  1.0826995950912709e-02,
        2.7726579759463501e+00,  5.1721119735948284e-02,  4.2335376950105824e-01,  8.2072105005128095e-03,
        3.9521264825956042e+00,  7.4832205937229584e-02,  2.5761269031357442e-01,  5.1359813390616359e-03,
        4.5485900655688685e+00,  8.6755719877848164e-02,  8.6755719877848719e-02,  1.7480626442406752e-03,
        4.5485900655688685e+00,  8.6755719877847776e-02, -8.6755719877848275e-02, -1.7480626442409205e-03,
        3.9521264825956046e+00,  7.4832205937228988e-02, -2.5761269031357459e-01, -5.1359813390615587e-03,
        2.7726579759463506e+00,  5.1721119735947924e-02, -4.2335376950105819e-01, -8.2072105005127505e-03,
        1.0196938687973565e+00,  1.8752907080322814e-02, -5.8872052964116572e-01, -1.0826995950913227e-02,
                                                        
        1.0196938687973578e+00, -1.8752907080323244e-02,  5.8872052964116495e-01, -1.0826995950913076e-02,
        2.7726579759463510e+00, -5.1721119735948159e-02,  4.2335376950105824e-01, -8.2072105005127211e-03,
        3.9521264825956046e+00, -7.4832205937229349e-02,  2.5761269031357431e-01, -5.1359813390616368e-03,
        4.5485900655688685e+00, -8.6755719877848525e-02,  8.6755719877847970e-02, -1.7480626442410686e-03,
        4.5485900655688667e+00, -8.6755719877848803e-02, -8.6755719877848247e-02,  1.7480626442408095e-03,
        3.9521264825956033e+00, -7.4832205937229904e-02, -2.5761269031357409e-01,  5.1359813390617765e-03,
        2.7726579759463501e+00, -5.1721119735948243e-02, -4.2335376950105824e-01,  8.2072105005128078e-03,
        1.0196938687973569e+00, -1.8752907080322453e-02, -5.8872052964116517e-01,  1.0826995950913517e-02,
                                                        
        8.9067264750613750e-01, -5.5737529763336964e-02,  5.1423009279750509e-01, -3.2180077812826779e-02,
        2.4167368175137387e+00, -1.5377005689540407e-01,  3.6684346655702621e-01, -2.4419028116211449e-02,
        3.4373208199425340e+00, -2.2239098197586465e-01,  2.2239098197586460e-01, -1.5199281451033603e-02,
        3.9521264825956037e+00, -2.5761269031357431e-01,  7.4832205937229085e-02, -5.1359813390615770e-03,
        3.9521264825956024e+00, -2.5761269031357420e-01, -7.4832205937229654e-02,  5.1359813390616099e-03,
        3.4373208199425314e+00, -2.2239098197586477e-01, -2.2239098197586432e-01,  1.5199281451033496e-02,
        2.4167368175137374e+00, -1.5377005689540416e-01, -3.6684346655702588e-01,  2.4419028116211574e-02,
        8.9067264750613806e-01, -5.5737529763336978e-02, -5.1423009279750453e-01,  3.2180077812826766e-02,
                                                        
        6.3388180537571459e-01, -9.2520732066093431e-02,  3.6597183096807528e-01, -5.3416869563980569e-02,
        1.7090814385197859e+00, -2.5479496666024121e-01,  2.5479496666024121e-01, -4.0272203461491109e-02,
        2.4167368175137387e+00, -3.6684346655702615e-01,  1.5377005689540399e-01, -2.4419028116211421e-02,
        2.7726579759463506e+00, -4.2335376950105813e-01,  5.1721119735947743e-02, -8.2072105005127610e-03,
        2.7726579759463492e+00, -4.2335376950105796e-01, -5.1721119735948368e-02,  8.2072105005128824e-03,
        2.4167368175137369e+00, -3.6684346655702577e-01, -1.5377005689540402e-01,  2.4419028116211414e-02,
        1.7090814385197848e+00, -2.5479496666024104e-01, -2.5479496666024082e-01,  4.0272203461491039e-02,
        6.3388180537571470e-01, -9.2520732066093667e-02, -3.6597183096807501e-01,  5.3416869563980458e-02,
                                                        
        2.3681559834188717e-01, -1.3672554945099052e-01,  1.3672554945099086e-01, -7.8938532780628887e-02,
        6.3388180537571537e-01, -3.6597183096807490e-01,  9.2520732066093569e-02, -5.3416869563980528e-02,
        8.9067264750613795e-01, -5.1423009279750498e-01,  5.5737529763336596e-02, -3.2180077812827008e-02,
        1.0196938687973571e+00, -5.8872052964116517e-01,  1.8752907080323126e-02, -1.0826995950912990e-02,
        1.0196938687973565e+00, -5.8872052964116517e-01, -1.8752907080323424e-02,  1.0826995950912984e-02,
        8.9067264750613695e-01, -5.1423009279750498e-01, -5.5737529763336513e-02,  3.2180077812827092e-02,
        6.3388180537571459e-01, -3.6597183096807484e-01, -9.2520732066093569e-02,  5.3416869563980382e-02,
        2.3681559834188676e-01, -1.3672554945099064e-01, -1.3672554945099064e-01,  7.8938532780628942e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        3.3943303570130012e-01,  1.9597175453399679e-01, -2.5382106732311723e-02, -1.4654366154499068e-02,
        2.5852815051768180e-01,  1.4926129729448101e-01, -2.1328350507206096e-02, -1.2313928906705882e-02,
        2.3704347887786692e-01,  1.3685711633978268e-01,  8.9241695525082501e-03,  5.1523716934339632e-03,
        2.9838358920906116e-01,  1.7227184555161776e-01,  2.6490559659328117e-02,  1.5294331750296710e-02,
        2.7855144864018777e-01,  1.6082175385557235e-01, -3.7940651355373674e-02, -2.1905045273254424e-02,
        2.0381615307077769e-01,  1.1767331084060774e-01, -5.2077916595914715e-03, -3.0067199165486271e-03,
        3.1270423370037509e-01,  1.8053987350364709e-01,  6.8074354322631203e-02,  3.9302746793080885e-02,
        4.0700430412220923e-01,  2.3498404454629240e-01, -1.3630183279984597e-02, -7.8693899858035734e-03,
                                                        
        9.2471721182926658e-01,  1.4194222210591301e-01, -5.7372732217365302e-02, -3.8154300808411140e-03,
        7.5028336345656188e-01,  1.3465370727118317e-01, -4.3336697096398154e-02, -3.9259592098273228e-04,
        7.2217959325319470e-01,  1.4323634988842188e-01,  2.7110977804240665e-02,  5.3477869464038522e-03,
        8.6863423587078081e-01,  1.5696251880408427e-01,  5.7444649535093165e-02,  2.5770203719961103e-03,
        8.2578551259476063e-01,  1.5512398027833324e-01, -8.2183371452935167e-02, -3.6385010847539575e-03,
        6.5796549632452850e-01,  1.4452993473923709e-01, -1.4707560116090753e-02, -2.4779739590503548e-03,
        8.7197741194623291e-01,  1.4235664650713659e-01,  1.3826739721006934e-01,  1.2232254132204976e-03,
        1.0677766339715837e+00,  1.4651303796530019e-01, -2.5222663666613812e-02,  1.1764683140077245e-03,
                                                        
        1.3236042788217255e+00,  8.8355333398443453e-02, -6.3319731576856511e-02,  3.8192839976832317e-04,
        1.1500636120564152e+00,  9.6159526874637705e-02, -3.6874019095442509e-02,  4.1238248045203435e-03,
        1.1823773590563609e+00,  1.2245895407850058e-01,  5.5530369624386594e-02,  1.1060156571162853e-02,
        1.4109181535507069e+00,  1.5612524704562614e-01,  7.6417719597735217e-02,  8.3770867360238831e-03,
        1.3654745745531423e+00,  1.5646564492203485e-01, -1.0265458216491942e-01, -8.1805579304477273e-03,
        1.1227306086903228e+00,  1.2380232799510095e-01, -3.7493711870136587e-02, -1.0677616889943637e-02,
        1.2888335450678960e+00,  9.8315354164002000e-02,  1.3339328690091568e-01, -4.0372943122702881e-03,
        1.4765773832016793e+00,  8.9508184647616132e-02, -2.4999331415682426e-02, -1.0475273788137696e-03,
                                                        
        1.5295214682804199e+00,  3.0531011366304865e-02, -6.1912475511358124e-02,  4.3055126846589144e-04,
        1.3790145293808918e+00,  3.6025346873860498e-02, -2.4982746177190984e-02,  2.7416048158394697e-03,
        1.4901027388499875e+00,  5.5206376781829561e-02,  8.9119553817967004e-02,  8.3325679648607077e-03,
        1.8232429688045224e+00,  8.1930596034726758e-02,  1.0321904762417754e-01,  7.0966672146820734e-03,
        1.7786475086544757e+00,  8.2079859803237842e-02, -1.2896624854644198e-01, -7.0104897377519498e-03,
        1.4335217115489933e+00,  5.5632998902093111e-02, -7.0292223316692631e-02, -8.2586125227445510e-03,
        1.5224137527023140e+00,  3.6542241591097972e-02,  1.2161406720537776e-01, -2.7634413497923073e-03,
        1.6849062587527901e+00,  3.0770547731789839e-02, -2.7798975095838528e-02, -5.6884765355937322e-04,
                                                        
        1.5295214682804197e+00, -3.0531011366304792e-02, -6.1912475511358096e-02, -4.3055126846587653e-04,
        1.3790145293808918e+00, -3.6025346873860394e-02, -2.4982746177190963e-02, -2.7416048158394866e-03,
        1.4901027388499877e+00, -5.5206376781829401e-02,  8.9119553817967018e-02, -8.3325679648606973e-03,
        1.8232429688045226e+00, -8.1930596034726688e-02,  1.0321904762417747e-01, -7.0966672146820916e-03,
        1.7786475086544760e+00, -8.2079859803237690e-02, -1.2896624854644190e-01,  7.0104897377519706e-03,
        1.4335217115489931e+00, -5.5632998902093125e-02, -7.0292223316692812e-02,  8.2586125227444539e-03,
        1.5224137527023138e+00, -3.6542241591098125e-02,  1.2161406720537783e-01,  2.7634413497923224e-03,
        1.6849062587527901e+00, -3.0770547731789867e-02, -2.7798975095838487e-02,  5.6884765355936563e-04,
                                                        
        1.3236042788217257e+00, -8.8355333398443481e-02, -6.3319731576856414e-02, -3.8192839976835028e-04,
        1.1500636120564156e+00, -9.6159526874637732e-02, -3.6874019095442509e-02, -4.1238248045203391e-03,
        1.1823773590563613e+00, -1.2245895407850059e-01,  5.5530369624386615e-02, -1.1060156571162844e-02,
        1.4109181535507074e+00, -1.5612524704562611e-01,  7.6417719597735093e-02, -8.3770867360238987e-03,
        1.3654745745531423e+00, -1.5646564492203499e-01, -1.0265458216491942e-01,  8.1805579304476995e-03,
        1.1227306086903228e+00, -1.2380232799510095e-01, -3.7493711870136739e-02,  1.0677616889943734e-02,
        1.2888335450678954e+00, -9.8315354164002056e-02,  1.3339328690091565e-01,  4.0372943122701779e-03,
        1.4765773832016789e+00, -8.9508184647616243e-02, -2.4999331415682322e-02,  1.0475273788138081e-03,
                                                        
        9.2471721182926658e-01, -1.4194222210591306e-01, -5.7372732217365198e-02,  3.8154300808411362e-03,
        7.5028336345656221e-01, -1.3465370727118314e-01, -4.3336697096398057e-02,  3.9259592098277505e-04,
        7.2217959325319514e-01, -1.4323634988842177e-01,  2.7110977804240696e-02, -5.3477869464038444e-03,
        8.6863423587078115e-01, -1.5696251880408424e-01,  5.7444649535093005e-02, -2.5770203719961463e-03,
        8.2578551259476063e-01, -1.5512398027833330e-01, -8.2183371452935264e-02,  3.6385010847539753e-03,
        6.5796549632452839e-01, -1.4452993473923711e-01, -1.4707560116090757e-02,  2.4779739590503531e-03,
        8.7197741194623257e-01, -1.4235664650713659e-01,  1.3826739721006917e-01, -1.2232254132204556e-03,
        1.0677766339715833e+00, -1.4651303796530019e-01, -2.5222663666613656e-02, -1.1764683140077696e-03,
                                                        
        3.3943303570129962e-01, -1.9597175453399704e-01, -2.5382106732311140e-02,  1.4654366154499322e-02,
        2.5852815051768296e-01, -1.4926129729448057e-01, -2.1328350507205690e-02,  1.2313928906706016e-02,
        2.3704347887786834e-01, -1.3685711633978229e-01,  8.9241695525079847e-03, -5.1523716934341332e-03,
        2.9838358920906127e-01, -1.7227184555161790e-01,  2.6490559659327607e-02, -1.5294331750296876e-02,
        2.7855144864018744e-01, -1.6082175385557249e-01, -3.7940651355373403e-02,  2.1905045273254604e-02,
        2.0381615307077783e-01, -1.1767331084060755e-01, -5.2077916595914967e-03,  3.0067199165486245e-03,
        3.1270423370037476e-01, -1.8053987350364703e-01,  6.8074354322630953e-02, -3.9302746793080955e-02,
        4.0700430412220806e-01, -2.3498404454629279e-01, -1.3630183279984807e-02,  7.8693899858033965e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against dirichletx-periodicy,
      // rotated by pgkyl):
      const double sol[256] = {
         3.3943303570129929e-01, -2.5382106732310627e-02,  1.9597175453399684e-01, -1.4654366154499540e-02,
         9.2471721182926570e-01, -5.7372732217365149e-02,  1.4194222210591295e-01, -3.8154300808411951e-03,
         1.3236042788217244e+00, -6.3319731576856414e-02,  8.8355333398443286e-02,  3.8192839976837950e-04,
         1.5295214682804186e+00, -6.1912475511358069e-02,  3.0531011366304959e-02,  4.3055126846582503e-04,
         1.5295214682804188e+00, -6.1912475511358138e-02, -3.0531011366304928e-02, -4.3055126846586200e-04,
         1.3236042788217248e+00, -6.3319731576856483e-02, -8.8355333398443120e-02, -3.8192839976834253e-04,
         9.2471721182926625e-01, -5.7372732217365198e-02, -1.4194222210591306e-01,  3.8154300808411765e-03,
         3.3943303570129940e-01, -2.5382106732310755e-02, -1.9597175453399698e-01,  1.4654366154499521e-02,
                                                         
         2.5852815051768291e-01, -2.1328350507206047e-02,  1.4926129729448026e-01, -1.2313928906705863e-02,
         7.5028336345656133e-01, -4.3336697096398091e-02,  1.3465370727118300e-01, -3.9259592098274870e-04,
         1.1500636120564143e+00, -3.6874019095442460e-02,  9.6159526874637649e-02,  4.1238248045203573e-03,
         1.3790145293808904e+00, -2.4982746177191012e-02,  3.6025346873860394e-02,  2.7416048158394237e-03,
         1.3790145293808904e+00, -2.4982746177190981e-02, -3.6025346873860421e-02, -2.7416048158394055e-03,
         1.1500636120564147e+00, -3.6874019095442495e-02, -9.6159526874637552e-02, -4.1238248045204129e-03,
         7.5028336345656177e-01, -4.3336697096398126e-02, -1.3465370727118303e-01,  3.9259592098280421e-04,
         2.5852815051768296e-01, -2.1328350507205975e-02, -1.4926129729448043e-01,  1.2313928906705860e-02,
                                                         
         2.3704347887786717e-01,  8.9241695525077245e-03,  1.3685711633978231e-01,  5.1523716934341627e-03,
         7.2217959325319403e-01,  2.7110977804240554e-02,  1.4323634988842177e-01,  5.3477869464038947e-03,
         1.1823773590563600e+00,  5.5530369624386580e-02,  1.2245895407850053e-01,  1.1060156571162867e-02,
         1.4901027388499859e+00,  8.9119553817966893e-02,  5.5206376781829269e-02,  8.3325679648606505e-03,
         1.4901027388499863e+00,  8.9119553817966962e-02, -5.5206376781829172e-02, -8.3325679648606332e-03,
         1.1823773590563604e+00,  5.5530369624386740e-02, -1.2245895407850046e-01, -1.1060156571162811e-02,
         7.2217959325319470e-01,  2.7110977804240727e-02, -1.4323634988842171e-01, -5.3477869464039407e-03,
         2.3704347887786764e-01,  8.9241695525078789e-03, -1.3685711633978243e-01, -5.1523716934341384e-03,
                                                         
         2.9838358920905961e-01,  2.6490559659327586e-02,  1.7227184555161815e-01,  1.5294331750296947e-02,
         8.6863423587077970e-01,  5.7444649535092998e-02,  1.5696251880408424e-01,  2.5770203719960869e-03,
         1.4109181535507060e+00,  7.6417719597735009e-02,  1.5612524704562616e-01,  8.3770867360239074e-03,
         1.8232429688045209e+00,  1.0321904762417757e-01,  8.1930596034726591e-02,  7.0966672146821827e-03,
         1.8232429688045211e+00,  1.0321904762417759e-01, -8.1930596034726424e-02, -7.0966672146822009e-03,
         1.4109181535507065e+00,  7.6417719597734982e-02, -1.5612524704562597e-01, -8.3770867360238883e-03,
         8.6863423587078070e-01,  5.7444649535092970e-02, -1.5696251880408427e-01, -2.5770203719961051e-03,
         2.9838358920906094e-01,  2.6490559659327881e-02, -1.7227184555161787e-01, -1.5294331750296758e-02,
                                                         
         2.7855144864018666e-01, -3.7940651355372897e-02,  1.6082175385557240e-01, -2.1905045273254868e-02,
         8.2578551259475919e-01, -8.2183371452935167e-02,  1.5512398027833305e-01, -3.6385010847539805e-03,
         1.3654745745531405e+00, -1.0265458216491943e-01,  1.5646564492203499e-01, -8.1805579304476978e-03,
         1.7786475086544749e+00, -1.2896624854644190e-01,  8.2079859803237828e-02, -7.0104897377519255e-03,
         1.7786475086544753e+00, -1.2896624854644170e-01, -8.2079859803237634e-02,  7.0104897377520183e-03,
         1.3654745745531418e+00, -1.0265458216491913e-01, -1.5646564492203482e-01,  8.1805579304476423e-03,
         8.2578551259476041e-01, -8.2183371452935056e-02, -1.5512398027833324e-01,  3.6385010847539249e-03,
         2.7855144864018744e-01, -3.7940651355373500e-02, -1.6082175385557240e-01,  2.1905045273254500e-02,
                                                         
         2.0381615307077769e-01, -5.2077916595916198e-03,  1.1767331084060717e-01, -3.0067199165485217e-03,
         6.5796549632452717e-01, -1.4707560116090692e-02,  1.4452993473923695e-01, -2.4779739590503344e-03,
         1.1227306086903213e+00, -3.7493711870136663e-02,  1.2380232799510092e-01, -1.0677616889943762e-02,
         1.4335217115489918e+00, -7.0292223316692895e-02,  5.5632998902093299e-02, -8.2586125227444972e-03,
         1.4335217115489924e+00, -7.0292223316692923e-02, -5.5632998902092813e-02,  8.2586125227444782e-03,
         1.1227306086903228e+00, -3.7493711870136760e-02, -1.2380232799510088e-01,  1.0677616889943752e-02,
         6.5796549632452839e-01, -1.4707560116090772e-02, -1.4452993473923714e-01,  2.4779739590503622e-03,
         2.0381615307077774e-01, -5.2077916595914247e-03, -1.1767331084060761e-01,  3.0067199165486531e-03,
                                                         
         3.1270423370037409e-01,  6.8074354322630801e-02,  1.8053987350364703e-01,  3.9302746793081066e-02,
         8.7197741194623157e-01,  1.3826739721006912e-01,  1.4235664650713634e-01,  1.2232254132204063e-03,
         1.2888335450678943e+00,  1.3339328690091576e-01,  9.8315354164002167e-02, -4.0372943122700409e-03,
         1.5224137527023127e+00,  1.2161406720537798e-01,  3.6542241591098042e-02, -2.7634413497924326e-03,
         1.5224137527023132e+00,  1.2161406720537785e-01, -3.6542241591097750e-02,  2.7634413497923402e-03,
         1.2888335450678956e+00,  1.3339328690091576e-01, -9.8315354164001792e-02,  4.0372943122702351e-03,
         8.7197741194623291e-01,  1.3826739721006936e-01, -1.4235664650713664e-01, -1.2232254132204989e-03,
         3.1270423370037487e-01,  6.8074354322631009e-02, -1.8053987350364709e-01, -3.9302746793080982e-02,
                                                         
         4.0700430412220695e-01, -1.3630183279984911e-02,  2.3498404454629299e-01, -7.8693899858033878e-03,
         1.0677766339715822e+00, -2.5222663666613573e-02,  1.4651303796530005e-01,  1.1764683140078709e-03,
         1.4765773832016775e+00, -2.4999331415682360e-02,  8.9508184647616285e-02, -1.0475273788140100e-03,
         1.6849062587527890e+00, -2.7798975095838598e-02,  3.0770547731789843e-02, -5.6884765355922739e-04,
         1.6849062587527892e+00, -2.7798975095838566e-02, -3.0770547731789812e-02,  5.6884765355924582e-04,
         1.4765773832016786e+00, -2.4999331415682616e-02, -8.9508184647615868e-02,  1.0475273788138251e-03,
         1.0677766339715835e+00, -2.5222663666613909e-02, -1.4651303796530027e-01, -1.1764683140077321e-03,
         4.0700430412220773e-01, -1.3630183279985131e-02, -2.3498404454629301e-01,  7.8693899858033340e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else {
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
	-9.5101835527468453e-02, -5.2237861001368492e-02, 4.3177375058357795e-02, 2.3387400381385207e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, -1.1937067064604241e-03, 5.2456792020110785e-03, 
	9.5101835527468453e-02, 5.2237861001368610e-02, 4.3177375058357802e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905626e-03, -1.1937067064604172e-03, -5.2456792020110915e-03, 
	9.5101835527468384e-02, 5.2237861001368582e-02, -4.3177375058357843e-02, -2.3387400381385235e-02, -2.0675606649251488e-03, -9.0857828980905643e-03, 1.1937067064604181e-03, -5.2456792020110924e-03, 
	-9.5101835527468620e-02, -5.2237861001368589e-02, -4.3177375058357788e-02, -2.3387400381385214e-02, 2.0675606649251388e-03, 9.0857828980905712e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, 2.3387400381385225e-02, 2.0675606649251596e-03, 9.0857828980905470e-03, -1.1937067064604137e-03, 5.2456792020110837e-03, 
	9.5101835527468453e-02, 5.2237861001368548e-02, 4.3177375058357816e-02, 2.3387400381385211e-02, -2.0675606649251488e-03, -9.0857828980905591e-03, -1.1937067064604259e-03, -5.2456792020110759e-03, 
	9.5101835527468537e-02, 5.2237861001368562e-02, -4.3177375058357795e-02, -2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905782e-03, 1.1937067064604241e-03, -5.2456792020110915e-03, 
	-9.5101835527468398e-02, -5.2237861001368506e-02, -4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251596e-03, 9.0857828980905504e-03, 1.1937067064604207e-03, 5.2456792020110776e-03, 
	-2.2959614113698959e-01, -2.1637630496127446e-02, 1.0423940445355719e-01, 9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347764e-02, -2.8818629202324655e-03, 2.1728314693314739e-03, 
	2.2959614113698990e-01, 2.1637630496127457e-02, 1.0423940445355717e-01, 9.6873784266194498e-03, -4.9915329982914251e-03, -2.1935020297347764e-02, -2.8818629202324387e-03, -2.1728314693314635e-03, 
	2.2959614113698973e-01, 2.1637630496127453e-02, -1.0423940445355724e-01, -9.6873784266194377e-03, -4.9915329982914043e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, -2.1728314693314518e-03, 
	-2.2959614113698987e-01, -2.1637630496127373e-02, -1.0423940445355719e-01, -9.6873784266194654e-03, 4.9915329982914147e-03, 2.1935020297347754e-02, 2.8818629202324413e-03, 2.1728314693314496e-03, 
	-2.2959614113698981e-01, -2.1637630496127425e-02, 1.0423940445355721e-01, 9.6873784266194429e-03, 4.9915329982914147e-03, 2.1935020297347736e-02, -2.8818629202324378e-03, 2.1728314693314679e-03, 
	2.2959614113698976e-01, 2.1637630496127457e-02, 1.0423940445355720e-01, 9.6873784266194567e-03, -4.9915329982914113e-03, -2.1935020297347757e-02, -2.8818629202324551e-03, -2.1728314693314783e-03, 
	2.2959614113698987e-01, 2.1637630496127443e-02, -1.0423940445355716e-01, -9.6873784266194429e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324551e-03, -2.1728314693314600e-03, 
	-2.2959614113698959e-01, -2.1637630496127443e-02, -1.0423940445355719e-01, -9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324560e-03, 2.1728314693314809e-03, 
	-2.2959614113698965e-01, 2.1637630496127432e-02, 1.0423940445355719e-01, -9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347771e-02, -2.8818629202324655e-03, -2.1728314693314657e-03, 
	2.2959614113698987e-01, -2.1637630496127495e-02, 1.0423940445355716e-01, -9.6873784266194689e-03, -4.9915329982914321e-03, -2.1935020297347757e-02, -2.8818629202324499e-03, 2.1728314693314679e-03, 
	2.2959614113698973e-01, -2.1637630496127470e-02, -1.0423940445355720e-01, 9.6873784266194654e-03, -4.9915329982914217e-03, -2.1935020297347750e-02, 2.8818629202324586e-03, 2.1728314693314540e-03, 
	-2.2959614113698976e-01, 2.1637630496127436e-02, -1.0423940445355717e-01, 9.6873784266194793e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324551e-03, -2.1728314693314583e-03, 
	-2.2959614113698976e-01, 2.1637630496127460e-02, 1.0423940445355719e-01, -9.6873784266194758e-03, 4.9915329982914286e-03, 2.1935020297347740e-02, -2.8818629202324586e-03, -2.1728314693314653e-03, 
	2.2959614113698981e-01, -2.1637630496127446e-02, 1.0423940445355721e-01, -9.6873784266194377e-03, -4.9915329982914251e-03, -2.1935020297347761e-02, -2.8818629202324577e-03, 2.1728314693314722e-03, 
	2.2959614113698987e-01, -2.1637630496127453e-02, -1.0423940445355717e-01, 9.6873784266194238e-03, -4.9915329982914321e-03, -2.1935020297347747e-02, 2.8818629202324586e-03, 2.1728314693314683e-03, 
	-2.2959614113698959e-01, 2.1637630496127436e-02, -1.0423940445355719e-01, 9.6873784266194619e-03, 4.9915329982914147e-03, 2.1935020297347750e-02, 2.8818629202324603e-03, -2.1728314693314774e-03, 
	-9.5101835527468509e-02, 5.2237861001368499e-02, 4.3177375058357795e-02, -2.3387400381385211e-02, 2.0675606649251735e-03, 9.0857828980905678e-03, -1.1937067064604276e-03, -5.2456792020110863e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, 4.3177375058357781e-02, -2.3387400381385200e-02, -2.0675606649251449e-03, -9.0857828980905608e-03, -1.1937067064604207e-03, 5.2456792020110811e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, -4.3177375058357788e-02, 2.3387400381385211e-02, -2.0675606649251509e-03, -9.0857828980905660e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, 5.2237861001368534e-02, -4.3177375058357767e-02, 2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604215e-03, -5.2456792020110932e-03, 
	-9.5101835527468426e-02, 5.2237861001368534e-02, 4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251457e-03, 9.0857828980905556e-03, -1.1937067064604207e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, -5.2237861001368527e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251423e-03, -9.0857828980905574e-03, -1.1937067064604128e-03, 5.2456792020110846e-03, 
	9.5101835527468495e-02, -5.2237861001368555e-02, -4.3177375058357816e-02, 2.3387400381385221e-02, -2.0675606649251466e-03, -9.0857828980905504e-03, 1.1937067064604102e-03, 5.2456792020110854e-03, 
	-9.5101835527468453e-02, 5.2237861001368499e-02, -4.3177375058357823e-02, 2.3387400381385200e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, 1.1937067064604337e-03, -5.2456792020110776e-03, 
	9.5101835527468231e-02, 5.2237861001368499e-02, -4.3177375058357809e-02, -2.3387400381385207e-02, -2.0675606649251488e-03, -9.0857828980905678e-03, 1.1937067064604233e-03, -5.2456792020110932e-03, 
	-9.5101835527468676e-02, -5.2237861001368610e-02, -4.3177375058357767e-02, -2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905591e-03, 1.1937067064604137e-03, 5.2456792020110846e-03, 
	-9.5101835527468565e-02, -5.2237861001368562e-02, 4.3177375058357836e-02, 2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905678e-03, -1.1937067064604172e-03, 5.2456792020110837e-03, 
	9.5101835527468426e-02, 5.2237861001368582e-02, 4.3177375058357823e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604154e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, 5.2237861001368589e-02, -4.3177375058357753e-02, -2.3387400381385214e-02, -2.0675606649251423e-03, -9.0857828980905435e-03, 1.1937067064604198e-03, -5.2456792020110820e-03, 
	-9.5101835527468231e-02, -5.2237861001368478e-02, -4.3177375058357809e-02, -2.3387400381385228e-02, 2.0675606649251596e-03, 9.0857828980905383e-03, 1.1937067064604172e-03, 5.2456792020110707e-03, 
	-9.5101835527468370e-02, -5.2237861001368541e-02, 4.3177375058357739e-02, 2.3387400381385197e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, -1.1937067064604207e-03, 5.2456792020110724e-03, 
	9.5101835527468315e-02, 5.2237861001368499e-02, 4.3177375058357753e-02, 2.3387400381385200e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604215e-03, -5.2456792020110898e-03, 
	2.2959614113698942e-01, 2.1637630496127418e-02, -1.0423940445355720e-01, -9.6873784266194619e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324586e-03, -2.1728314693314574e-03, 
	-2.2959614113699003e-01, -2.1637630496127470e-02, -1.0423940445355714e-01, -9.6873784266194619e-03, 4.9915329982914008e-03, 2.1935020297347740e-02, 2.8818629202324387e-03, 2.1728314693314609e-03, 
	-2.2959614113698998e-01, -2.1637630496127495e-02, 1.0423940445355724e-01, 9.6873784266194567e-03, 4.9915329982914008e-03, 2.1935020297347771e-02, -2.8818629202324447e-03, 2.1728314693314791e-03, 
	2.2959614113698981e-01, 2.1637630496127450e-02, 1.0423940445355726e-01, 9.6873784266194706e-03, -4.9915329982914286e-03, -2.1935020297347750e-02, -2.8818629202324638e-03, -2.1728314693314713e-03, 
	2.2959614113699001e-01, 2.1637630496127470e-02, -1.0423940445355717e-01, -9.6873784266194706e-03, -4.9915329982914355e-03, -2.1935020297347764e-02, 2.8818629202324621e-03, -2.1728314693314817e-03, 
	-2.2959614113698937e-01, -2.1637630496127439e-02, -1.0423940445355716e-01, -9.6873784266194151e-03, 4.9915329982914147e-03, 2.1935020297347722e-02, 2.8818629202324586e-03, 2.1728314693314783e-03, 
	-2.2959614113698959e-01, -2.1637630496127418e-02, 1.0423940445355706e-01, 9.6873784266194290e-03, 4.9915329982914286e-03, 2.1935020297347764e-02, -2.8818629202324586e-03, 2.1728314693314869e-03, 
	2.2959614113698953e-01, 2.1637630496127457e-02, 1.0423940445355712e-01, 9.6873784266194498e-03, -4.9915329982914182e-03, -2.1935020297347757e-02, -2.8818629202324560e-03, -2.1728314693314627e-03, 
	2.2959614113698945e-01, -2.1637630496127394e-02, -1.0423940445355723e-01, 9.6873784266194342e-03, -4.9915329982914147e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, 2.1728314693314613e-03, 
	-2.2959614113699015e-01, 2.1637630496127450e-02, -1.0423940445355720e-01, 9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324577e-03, -2.1728314693314653e-03, 
	-2.2959614113699009e-01, 2.1637630496127460e-02, 1.0423940445355727e-01, -9.6873784266194654e-03, 4.9915329982914425e-03, 2.1935020297347778e-02, -2.8818629202324517e-03, -2.1728314693314800e-03, 
	2.2959614113698970e-01, -2.1637630496127484e-02, 1.0423940445355724e-01, -9.6873784266194619e-03, -4.9915329982914008e-03, -2.1935020297347764e-02, -2.8818629202324517e-03, 2.1728314693314661e-03, 
	2.2959614113698992e-01, -2.1637630496127498e-02, -1.0423940445355714e-01, 9.6873784266194654e-03, -4.9915329982914182e-03, -2.1935020297347764e-02, 2.8818629202324447e-03, 2.1728314693314761e-03, 
	-2.2959614113698948e-01, 2.1637630496127401e-02, -1.0423940445355717e-01, 9.6873784266194342e-03, 4.9915329982914564e-03, 2.1935020297347750e-02, 2.8818629202324829e-03, -2.1728314693314635e-03, 
	-2.2959614113698953e-01, 2.1637630496127484e-02, 1.0423940445355713e-01, -9.6873784266194203e-03, 4.9915329982914425e-03, 2.1935020297347750e-02, -2.8818629202324794e-03, -2.1728314693314913e-03, 
	2.2959614113698965e-01, -2.1637630496127401e-02, 1.0423940445355709e-01, -9.6873784266194446e-03, -4.9915329982914251e-03, -2.1935020297347743e-02, -2.8818629202324517e-03, 2.1728314693314670e-03, 
	9.5101835527468259e-02, -5.2237861001368513e-02, -4.3177375058357836e-02, 2.3387400381385239e-02, -2.0675606649251405e-03, -9.0857828980905626e-03, 1.1937067064604181e-03, 5.2456792020110898e-03, 
	-9.5101835527468676e-02, 5.2237861001368631e-02, -4.3177375058357767e-02, 2.3387400381385214e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604068e-03, -5.2456792020110776e-03, 
	-9.5101835527468620e-02, 5.2237861001368603e-02, 4.3177375058357795e-02, -2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905608e-03, -1.1937067064604102e-03, -5.2456792020110828e-03, 
	9.5101835527468356e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251492e-03, -9.0857828980905608e-03, -1.1937067064604163e-03, 5.2456792020110880e-03, 
	9.5101835527468453e-02, -5.2237861001368582e-02, -4.3177375058357760e-02, 2.3387400381385211e-02, -2.0675606649251466e-03, -9.0857828980905591e-03, 1.1937067064604172e-03, 5.2456792020110811e-03, 
	-9.5101835527468342e-02, 5.2237861001368485e-02, -4.3177375058357788e-02, 2.3387400381385207e-02, 2.0675606649251596e-03, 9.0857828980905574e-03, 1.1937067064604267e-03, -5.2456792020110863e-03, 
	-9.5101835527468287e-02, 5.2237861001368478e-02, 4.3177375058357788e-02, -2.3387400381385197e-02, 2.0675606649251735e-03, 9.0857828980905383e-03, -1.1937067064604241e-03, -5.2456792020110768e-03, 
	9.5101835527468398e-02, -5.2237861001368541e-02, 4.3177375058357746e-02, -2.3387400381385211e-02, -2.0675606649251379e-03, -9.0857828980905539e-03, -1.1937067064604189e-03, 5.2456792020110811e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.6174232732491770e-01,  1.4239797833436429e-01,  1.4239797833437121e-01,  7.7179579805169304e-02,
        -6.7537275729811939e-03, -6.7537275729829148e-03, -3.8992664322967511e-03, -3.8992664322970148e-03,
         6.6919163252610048e-01,  3.6607316345625468e-01,  9.5063998582039566e-02,  5.3241649763932503e-02,
        -1.5712542948565247e-02, -5.0333399450138384e-03, -1.2731080364168593e-03, -2.9060001721712362e-03,
         9.3261285108184055e-01,  5.1325052577226338e-01,  5.7325647937302174e-02,  3.1906370575560429e-02,
        -1.9514998592151618e-02, -4.7983776055374530e-03, -9.2224075298933392e-04, -2.7703446022362605e-03,
         1.0655299621632952e+00,  5.8724826482558512e-01,  1.9276481795971818e-02,  1.0736800514483670e-02,
        -2.1638935540250832e-02, -4.9049622075441517e-03, -3.0401481573749179e-04, -2.8318812508971261e-03,
         1.0655299621633028e+00,  5.8724826482558323e-01, -1.9276481795972602e-02, -1.0736800514480182e-02,
        -2.1638935540251626e-02, -4.9049622075482769e-03,  3.0401481573705659e-04, -2.8318812508936436e-03,
         9.3261285108186132e-01,  5.1325052577225383e-01, -5.7325647937298004e-02, -3.1906370575565557e-02,
        -1.9514998592151112e-02, -4.7983776055450077e-03,  9.2224075299059171e-04, -2.7703446022306617e-03,
         6.6919163252609015e-01,  3.6607316345626151e-01, -9.5063998582045506e-02, -5.3241649763929401e-02,
        -1.5712542948564903e-02, -5.0333399450088745e-03,  1.2731080364155342e-03, -2.9060001721745977e-03,
         2.6174232732491648e-01,  1.4239797833437010e-01, -1.4239797833436521e-01, -7.7179579805168486e-02,
        -6.7537275729838489e-03, -6.7537275729819659e-03,  3.8992664322963539e-03, -3.8992664322969480e-03,
                                                         
         6.6919163252608815e-01,  9.5063998582044021e-02,  3.6607316345626351e-01,  5.3241649763928853e-02,
        -5.0333399450085128e-03, -1.5712542948565542e-02, -2.9060001721746753e-03, -1.2731080364157741e-03,
         1.7423637960024951e+00,  2.5663155896154988e-01,  2.5663155896154877e-01,  3.9269506879060931e-02,
        -1.3304720861538551e-02, -1.3304720861538329e-02, -1.8694838265539103e-03, -1.8694838265541740e-03,
         2.4546088336295506e+00,  3.6679544630615912e-01,  1.5482377206844999e-01,  2.4122215546182844e-02,
        -1.8491615532474651e-02, -1.3118456058247588e-02, -1.1251712079693174e-03, -2.0332549321147595e-03,
         2.8136363432486524e+00,  4.2275290286990563e-01,  5.2063204566974094e-02,  8.1141391202813138e-03,
        -2.1072335114916858e-02, -1.3426476507161618e-02, -3.6480793765665813e-04, -2.0880173238903453e-03,
         2.8136363432486524e+00,  4.2275290286990524e-01, -5.2063204566972256e-02, -8.1141391202823199e-03,
        -2.1072335114916522e-02, -1.3426476507160893e-02,  3.6480793765687319e-04, -2.0880173238910131e-03,
         2.4546088336295533e+00,  3.6679544630615807e-01, -1.5482377206845077e-01, -2.4122215546182161e-02,
        -1.8491615532474186e-02, -1.3118456058246705e-02,  1.1251712079691587e-03, -2.0332549321155132e-03,
         1.7423637960024967e+00,  2.5663155896154877e-01, -2.5663155896155010e-01, -3.9269506879060175e-02,
        -1.3304720861539127e-02, -1.3304720861538711e-02,  1.8694838265534816e-03, -1.8694838265539018e-03,
         6.6919163252608682e-01,  9.5063998582043396e-02, -3.6607316345626262e-01, -5.3241649763929991e-02,
        -5.0333399450070591e-03, -1.5712542948565018e-02,  2.9060001721762860e-03, -1.2731080364161030e-03,
                                                         
         9.3261285108185710e-01,  5.7325647937300842e-02,  5.1325052577225627e-01,  3.1906370575564183e-02,
        -4.7983776055426555e-03, -1.9514998592151105e-02, -2.7703446022324602e-03, -9.2224075298994704e-04,
         2.4546088336295533e+00,  1.5482377206845188e-01,  3.6679544630615785e-01,  2.4122215546181956e-02,
        -1.3118456058246891e-02, -1.8491615532474508e-02, -2.0332549321152001e-03, -1.1251712079690640e-03,
         3.4763948230693575e+00,  2.2273928406436580e-01,  2.2273928406436491e-01,  1.4934123050414328e-02,
        -1.8796563154148561e-02, -1.8796563154148505e-02, -1.2450017281912052e-03, -1.2450017281912329e-03,
         3.9932287149419077e+00,  2.5750472735140589e-01,  7.4989196262071997e-02,  5.0536521254172630e-03,
        -2.1676786288648296e-02, -1.9312217015521736e-02, -4.1789587383843661e-04, -1.3101165429916409e-03,
         3.9932287149419086e+00,  2.5750472735140550e-01, -7.4989196262071983e-02, -5.0536521254171546e-03,
        -2.1676786288648310e-02, -1.9312217015521806e-02,  4.1789587383839465e-04, -1.3101165429915658e-03,
         3.4763948230693589e+00,  2.2273928406436544e-01, -2.2273928406436508e-01, -1.4934123050414191e-02,
        -1.8796563154148464e-02, -1.8796563154148630e-02,  1.2450017281912470e-03, -1.2450017281910922e-03,
         2.4546088336295528e+00,  1.5482377206845233e-01, -3.6679544630615840e-01, -2.4122215546181762e-02,
        -1.3118456058247220e-02, -1.8491615532474616e-02,  2.0332549321149230e-03, -1.1251712079691938e-03,
         9.3261285108185132e-01,  5.7325647937300522e-02, -5.1325052577225905e-01, -3.1906370575564613e-02,
        -4.7983776055399753e-03, -1.9514998592151157e-02,  2.7703446022344595e-03, -9.2224075299001112e-04,
                                                         
         1.0655299621633085e+00,  1.9276481795971488e-02,  5.8724826482558201e-01,  1.0736800514482708e-02,
        -4.9049622075492960e-03, -2.1638935540251314e-02, -2.8318812508928769e-03, -3.0401481573748951e-04,
         2.8136363432486569e+00,  5.2063204566973623e-02,  4.2275290286990502e-01,  8.1141391202819799e-03,
        -1.3426476507160820e-02, -2.1072335114916560e-02, -2.0880173238911450e-03, -3.6480793765673593e-04,
         3.9932287149419117e+00,  7.4989196262073329e-02,  2.5750472735140501e-01,  5.0536521254175961e-03,
        -1.9312217015521948e-02, -2.1676786288648508e-02, -1.3101165429914947e-03, -4.1789587383842935e-04,
         4.5910338085415718e+00,  8.6771187773204306e-02,  8.6771187773202557e-02,  1.7165201446689979e-03,
        -2.2348374449956559e-02, -2.2348374449956337e-02, -4.4280976908144193e-04, -4.4280976908130315e-04,
         4.5910338085415718e+00,  8.6771187773204195e-02, -8.6771187773203362e-02, -1.7165201446691534e-03,
        -2.2348374449956399e-02, -2.2348374449956399e-02,  4.4280976908152368e-04, -4.4280976908143911e-04,
         3.9932287149419095e+00,  7.4989196262072677e-02, -2.5750472735140451e-01, -5.0536521254171502e-03,
        -1.9312217015521778e-02, -2.1676786288648116e-02,  1.3101165429914377e-03, -4.1789587383835627e-04,
         2.8136363432486573e+00,  5.2063204566974289e-02, -4.2275290286990530e-01, -8.1141391202814977e-03,
        -1.3426476507160697e-02, -2.1072335114916685e-02,  2.0880173238911949e-03, -3.6480793765670535e-04,
         1.0655299621633088e+00,  1.9276481795971901e-02, -5.8724826482558179e-01, -1.0736800514483127e-02,
        -4.9049622075491390e-03, -2.1638935540251369e-02,  2.8318812508928517e-03, -3.0401481573744788e-04,
                                                         
         1.0655299621633059e+00, -1.9276481795971124e-02,  5.8724826482558301e-01, -1.0736800514483365e-02,
        -4.9049622075478649e-03, -2.1638935540251245e-02, -2.8318812508938908e-03,  3.0401481573757262e-04,
         2.8136363432486569e+00, -5.2063204566973290e-02,  4.2275290286990536e-01, -8.1141391202814873e-03,
        -1.3426476507160872e-02, -2.1072335114916612e-02, -2.0880173238910075e-03,  3.6480793765668437e-04,
         3.9932287149419130e+00, -7.4989196262072108e-02,  2.5750472735140517e-01, -5.0536521254173029e-03,
        -1.9312217015521708e-02, -2.1676786288648268e-02, -1.3101165429914407e-03,  4.1789587383857528e-04,
         4.5910338085415718e+00, -8.6771187773202793e-02,  8.6771187773202904e-02, -1.7165201446689685e-03,
        -2.2348374449956222e-02, -2.2348374449956222e-02, -4.4280976908148416e-04,  4.4280976908131210e-04,
         4.5910338085415727e+00, -8.6771187773202474e-02, -8.6771187773202779e-02,  1.7165201446694099e-03,
        -2.2348374449956278e-02, -2.2348374449956500e-02,  4.4280976908150861e-04,  4.4280976908139759e-04,
         3.9932287149419117e+00, -7.4989196262071275e-02, -2.5750472735140512e-01,  5.0536521254170332e-03,
        -1.9312217015521788e-02, -2.1676786288648237e-02,  1.3101165429913700e-03,  4.1789587383833242e-04,
         2.8136363432486591e+00, -5.2063204566973408e-02, -4.2275290286990502e-01,  8.1141391202817024e-03,
        -1.3426476507160754e-02, -2.1072335114916716e-02,  2.0880173238911632e-03,  3.6480793765672633e-04,
         1.0655299621633081e+00, -1.9276481795971284e-02, -5.8724826482558301e-01,  1.0736800514483185e-02,
        -4.9049622075484347e-03, -2.1638935540251134e-02,  2.8318812508933353e-03,  3.0401481573757126e-04,
                                                         
         9.3261285108185266e-01, -5.7325647937301265e-02,  5.1325052577225927e-01, -3.1906370575562747e-02,
        -4.7983776055406909e-03, -1.9514998592151191e-02, -2.7703446022337964e-03,  9.2224075298978842e-04,
         2.4546088336295555e+00, -1.5482377206845083e-01,  3.6679544630615885e-01, -2.4122215546182133e-02,
        -1.3118456058247074e-02, -1.8491615532474526e-02, -2.0332549321151190e-03,  1.1251712079692335e-03,
         3.4763948230693624e+00, -2.2273928406436483e-01,  2.2273928406436511e-01, -1.4934123050414449e-02,
        -1.8796563154148675e-02, -1.8796563154148564e-02, -1.2450017281912206e-03,  1.2450017281911187e-03,
         3.9932287149419121e+00, -2.5750472735140495e-01,  7.4989196262071955e-02, -5.0536521254169655e-03,
        -2.1676786288648216e-02, -1.9312217015521767e-02, -4.1789587383831805e-04,  1.3101165429915314e-03,
         3.9932287149419126e+00, -2.5750472735140467e-01, -7.4989196262072025e-02,  5.0536521254170704e-03,
        -2.1676786288648092e-02, -1.9312217015522087e-02,  4.1789587383837931e-04,  1.3101165429914169e-03,
         3.4763948230693611e+00, -2.2273928406436524e-01, -2.2273928406436566e-01,  1.4934123050414329e-02,
        -1.8796563154148464e-02, -1.8796563154148713e-02,  1.2450017281911525e-03,  1.2450017281911802e-03,
         2.4546088336295551e+00, -1.5482377206845158e-01, -3.6679544630615818e-01,  2.4122215546182102e-02,
        -1.3118456058246773e-02, -1.8491615532474696e-02,  2.0332549321152413e-03,  1.1251712079691884e-03,
         9.3261285108185532e-01, -5.7325647937301619e-02, -5.1325052577225805e-01,  3.1906370575563073e-02,
        -4.7983776055419217e-03, -1.9514998592151160e-02,  2.7703446022327573e-03,  9.2224075298971849e-04,
                                                         
         6.6919163252609370e-01, -9.5063998582043147e-02,  3.6607316345626068e-01, -5.3241649763929991e-02,
        -5.0333399450107766e-03, -1.5712542948565014e-02, -2.9060001721731752e-03,  1.2731080364163324e-03,
         1.7423637960024991e+00, -2.5663155896154932e-01,  2.5663155896154954e-01, -3.9269506879060966e-02,
        -1.3304720861538491e-02, -1.3304720861538602e-02, -1.8694838265540792e-03,  1.8694838265539704e-03,
         2.4546088336295564e+00, -3.6679544630615801e-01,  1.5482377206845152e-01, -2.4122215546181613e-02,
        -1.8491615532474564e-02, -1.3118456058247055e-02, -1.1251712079691414e-03,  2.0332549321152114e-03,
         2.8136363432486586e+00, -4.2275290286990402e-01,  5.2063204566973983e-02, -8.1141391202815896e-03,
        -2.1072335114916560e-02, -1.3426476507160459e-02, -3.6480793765670513e-04,  2.0880173238911949e-03,
         2.8136363432486595e+00, -4.2275290286990463e-01, -5.2063204566974060e-02,  8.1141391202815445e-03,
        -2.1072335114916383e-02, -1.3426476507160671e-02,  3.6480793765677414e-04,  2.0880173238913636e-03,
         2.4546088336295546e+00, -3.6679544630615818e-01, -1.5482377206845166e-01,  2.4122215546182019e-02,
        -1.8491615532474439e-02, -1.3118456058246849e-02,  1.1251712079690349e-03,  2.0332549321153233e-03,
         1.7423637960024978e+00, -2.5663155896154966e-01, -2.5663155896154949e-01,  3.9269506879060584e-02,
        -1.3304720861538579e-02, -1.3304720861538607e-02,  1.8694838265540319e-03,  1.8694838265540597e-03,
         6.6919163252609004e-01, -9.5063998582043216e-02, -3.6607316345626228e-01,  5.3241649763930428e-02,
        -5.0333399450094695e-03, -1.5712542948565080e-02,  2.9060001721740326e-03,  1.2731080364163240e-03,
                                                         
         2.6174232732492186e-01, -1.4239797833436652e-01,  1.4239797833436699e-01, -7.7179579805169304e-02,
        -6.7537275729839755e-03, -6.7537275729838663e-03, -3.8992664322955498e-03,  3.8992664322956201e-03,
         6.6919163252609259e-01, -3.6607316345626134e-01,  9.5063998582043535e-02, -5.3241649763929630e-02,
        -1.5712542948564886e-02, -5.0333399450095380e-03, -1.2731080364162786e-03,  2.9060001721740491e-03,
         9.3261285108185321e-01, -5.1325052577226005e-01,  5.7325647937300730e-02, -3.1906370575564294e-02,
        -1.9514998592151306e-02, -4.7983776055403804e-03, -9.2224075298998217e-04,  2.7703446022338254e-03,
         1.0655299621633147e+00, -5.8724826482558012e-01,  1.9276481795972494e-02, -1.0736800514482444e-02,
        -2.1638935540251224e-02, -4.9049622075506205e-03, -3.0401481573732217e-04,  2.8318812508918382e-03,
         1.0655299621633119e+00, -5.8724826482558135e-01, -1.9276481795972668e-02,  1.0736800514482532e-02,
        -2.1638935540251137e-02, -4.9049622075495979e-03,  3.0401481573736000e-04,  2.8318812508924367e-03,
         9.3261285108185554e-01, -5.1325052577225794e-01, -5.7325647937301508e-02,  3.1906370575563600e-02,
        -1.9514998592151313e-02, -4.7983776055423831e-03,  9.2224075298983992e-04,  2.7703446022324342e-03,
         6.6919163252608904e-01, -3.6607316345626240e-01, -9.5063998582043160e-02,  5.3241649763930317e-02,
        -1.5712542948565038e-02, -5.0333399450093533e-03,  1.2731080364162960e-03,  2.9060001721740885e-03,
         2.6174232732492109e-01, -1.4239797833436624e-01, -1.4239797833436613e-01,  7.7179579805169332e-02,
        -6.7537275729838767e-03, -6.7537275729839148e-03,  3.8992664322956569e-03,  3.8992664322956318e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++)
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174043031e-01,  1.9147096450810874e-01, -2.8048110764785784e-02, -1.3769366894035592e-02,
        -7.1159030934106020e-03,  2.2827046419852215e-03,  1.8777907284159622e-03,  1.3179201395296886e-03,
         2.6018094650847862e-01,  1.5014541606507584e-01, -2.4623528318145276e-02, -1.1377222883164622e-02,
        -5.4317382757281126e-05,  5.5564633406491573e-04,  2.1992176825353575e-03,  3.2080256054490733e-04,
         2.2061875008365231e-01,  1.3141643256470523e-01,  9.6631517600851008e-03,  5.1140605928803942e-03,
         3.1310264886016711e-03,  6.6601405906152981e-03, -3.6015854094473347e-04,  3.8452339628319973e-03,
         3.1728477488866691e-01,  1.7738358224001630e-01,  3.1338192874284236e-02,  1.2875214739491091e-02,
        -4.4933329941175878e-03, -4.8107149369187568e-03, -4.0417674588017032e-03, -2.7774675638256754e-03,
         2.9856868259017194e-01,  1.6904725884004743e-01, -4.4358254629947795e-02, -1.8966632025180653e-02,
        -2.5805304029333908e-03, -6.5259204756552055e-03,  5.1461245497285120e-03, -3.7677419433303117e-03,
         1.6472836762602577e-01,  1.0615968567242218e-01, -6.4622802488493898e-03, -2.0693364430828479e-03,
         8.5621733627070755e-03,  1.3963870956120502e-02,  1.2871184688642739e-03,  8.0620446554442865e-03,
         3.2567896066607988e-01,  1.8138329229246208e-01,  7.9694807939026877e-02,  3.4130327869920304e-02,
        -5.1491649123717653e-03, -1.2897504262264271e-03, -9.2033633129310552e-03, -7.4463775577016827e-04,
         4.5758071423995317e-01,  2.4387817433520742e-01, -1.7203978611667874e-02, -5.9370449568280142e-03,
        -1.5729094905946293e-02, -1.1872819639637478e-02,  3.0950378831334125e-03, -6.8547756149879950e-03,
                                                         
         9.2555564511589195e-01,  1.4283166346900319e-01, -5.7087203270645710e-02, -4.2107167636395056e-03,
        -6.6585943733794726e-03,  4.3810045656668519e-03,  9.3715478409904616e-04, -1.0646611375156228e-04,
         7.5643742119286139e-01,  1.3194920980432359e-01, -4.3151229134044317e-02, -9.9766889501246210e-04,
        -3.4774021071348250e-03,  2.3686176801918867e-03,  8.9950742716126829e-04,  7.2591693417456003e-04,
         7.0097359003902937e-01,  1.4059912442543343e-01,  2.5969251836348434e-02,  5.4595811792827933e-03,
        -9.8784417722106539e-04,  1.3863694644921702e-02,  5.3783951383761853e-04,  3.1373990954381765e-04,
         8.9395264296783827e-01,  1.6114696623034061e-01,  5.3336355793392402e-02,  4.9642255681518815e-03,
        -1.6240840546830664e-04, -1.1009301949378427e-02, -6.1273948817407056e-05, -8.0128831641327464e-04,
         8.5858667148480283e-01,  1.5815426669652310e-01, -7.7423345159838022e-02, -6.2531566299462403e-03,
         4.2102309020289271e-04, -1.3850860419405031e-02,  3.9811827989687796e-04, -4.6131410499152816e-04,
         6.0862841609240770e-01,  1.3886189996819126e-01, -1.1804099479237351e-02, -3.6257448687969356e-03,
        -1.6305922506315459e-04,  2.8818588277168135e-02, -7.3533836184461665e-04,  5.1433038859834352e-04,
         8.8573073934651869e-01,  1.4220407629839985e-01,  1.3180886462596331e-01,  5.2114498588328763e-03,
        -4.9622366030997230e-03, -2.1500466509340923e-03, -2.0354679892535497e-03,  2.4794549885235803e-04,
         1.1259773284051051e+00,  1.5150225754918686e-01, -2.1648595211938813e-02, -5.4796944887231718e-04,
        -8.3847823258310743e-03, -2.4261809779577561e-02,  5.9460294920793655e-05, -2.9801117729407584e-04,
                                                         
         1.3294757533414434e+00,  8.9788893137442261e-02, -6.2911497442604103e-02, -3.9614170768240779e-05,
        -7.1100499232079391e-03,  4.6679346543523435e-03,  2.4956613702374012e-04,  2.7212527769270804e-04,
         1.1558500674443988e+00,  9.5187518915131650e-02, -3.6801896226935131e-02,  3.8875289024224544e-03,
        -6.1608141969413492e-03,  5.0780393003260101e-03,  2.9847536506069259e-04,  8.3836836755803482e-04,
         1.1608685714961731e+00,  1.2037284313135879e-01,  5.3537968971594352e-02,  1.0607159984611206e-02,
        -4.5111155521692318e-03,  1.5797402263497122e-02,  6.5397859158025575e-04,  8.0268670457471416e-04,
         1.4515046707409818e+00,  1.5802986276729328e-01,  7.3764923340376418e-02,  7.7325221486257015e-03,
        -2.2739467835626078e-03, -1.5570616292196532e-02,  6.3765139919740651e-04, -1.8321877472712517e-03,
         1.4115464704823018e+00,  1.5776040893168652e-01, -9.9620570722380883e-02, -7.8173884689785901e-03,
        -2.1636612802878562e-03, -1.7728464227987777e-02, -5.7397803419402017e-04, -1.7774214977044778e-03,
         1.0684562036485115e+00,  1.2101078886542831e-01, -3.5201132471140130e-02, -9.9394198301684059e-03,
        -4.5077782114690034e-03,  3.1273725348626311e-02, -7.7939850703535667e-04,  9.0314366050543921e-04,
         1.3060178631336954e+00,  9.8310629669391275e-02,  1.3131146816260594e-01, -3.4378558002396307e-03,
        -6.6184622280532198e-03, -5.2033263638921626e-04, -4.3920547811377179e-04,  6.9297032614724626e-04,
         1.5457214754324096e+00,  9.2030059599753342e-02, -2.4079263611516431e-02, -9.9293276550431963e-04,
        -7.4607497917453056e-03, -2.4657338466129757e-02, -4.7089473518925407e-05,  6.9652583640923910e-05,
                                                         
         1.5410607675732244e+00,  3.1259713850908775e-02, -6.1609756908925520e-02,  4.0278881947996382e-04,
        -7.9698971022144766e-03,  5.4767589025950247e-03, -5.1276054205655116e-05,  1.9484961975728570e-04,
         1.3895291108288024e+00,  3.6049179432957129e-02, -2.3292866819651249e-02,  2.8177193660516302e-03,
        -9.0097499591924313e-03,  7.4784863973146746e-03, -5.4908327268816275e-04,  5.4753041006383430e-04,
         1.4814229261146523e+00,  5.3722083116000623e-02,  9.0436045201508081e-02,  7.5759200922297047e-03,
        -1.3014031031221703e-02,  1.8391205195993725e-02, -1.7627894821588457e-03,  6.9484611672705331e-04,
         1.9036696745231794e+00,  8.2014576705607065e-02,  1.0223309369326436e-01,  6.0293706608141534e-03,
        -1.8550660598603910e-02, -2.1202408094434324e-02, -1.4337850889725194e-03, -1.4193287657706831e-03,
         1.8635410509569001e+00,  8.2038435164303550e-02, -1.2804640133826681e-01, -6.0050462696767748e-03,
        -1.8554396836116659e-02, -2.3251244190818931e-02,  1.4316279712387684e-03, -1.4111570005111908e-03,
         1.3908819207316627e+00,  5.4053805794842834e-02, -7.0885965719392860e-02, -7.4371696780804368e-03,
        -1.3096490503921039e-02,  3.4035995955905070e-02,  1.7194957188661293e-03,  6.9165401818151314e-04,
         1.5458079707862993e+00,  3.6785251742612025e-02,  1.1844442680216077e-01, -2.7477367760387942e-03,
        -9.2110386667161315e-03,  1.5895856233272834e-03,  5.2377094526747942e-04,  5.2519154906815752e-04,
         1.7614161838416025e+00,  3.1685800667273945e-02, -2.7278574910696574e-02, -6.3584621477922470e-04,
        -8.0924625745515209e-03, -2.4216573332765313e-02,  1.2203926265282081e-04,  1.8482328475642772e-04,
                                                         
         1.5410607675732240e+00, -3.1259713850908713e-02, -6.1609756908925596e-02, -4.0278881947999954e-04,
        -7.9698971022143933e-03,  5.4767589025950534e-03, -5.1276054205712274e-05, -1.9484961975725637e-04,
         1.3895291108288028e+00, -3.6049179432956852e-02, -2.3292866819651072e-02, -2.8177193660514754e-03,
        -9.0097499591924035e-03,  7.4784863973146321e-03, -5.4908327268812123e-04, -5.4753041006383484e-04,
         1.4814229261146523e+00, -5.3722083116000553e-02,  9.0436045201507997e-02, -7.5759200922298131e-03,
        -1.3014031031221738e-02,  1.8391205195993847e-02, -1.7627894821589177e-03, -6.9484611672699509e-04,
         1.9036696745231787e+00, -8.2014576705607023e-02,  1.0223309369326418e-01, -6.0293706608139782e-03,
        -1.8550660598603833e-02, -2.1202408094434275e-02, -1.4337850889724136e-03,  1.4193287657707004e-03,
         1.8635410509568999e+00, -8.2038435164303550e-02, -1.2804640133826664e-01,  6.0050462696767818e-03,
        -1.8554396836116590e-02, -2.3251244190819000e-02,  1.4316279712386628e-03,  1.4111570005111837e-03,
         1.3908819207316627e+00, -5.4053805794842855e-02, -7.0885965719392902e-02,  7.4371696780805912e-03,
        -1.3096490503921063e-02,  3.4035995955905098e-02,  1.7194957188662130e-03, -6.9165401818148885e-04,
         1.5458079707862988e+00, -3.6785251742611901e-02,  1.1844442680216073e-01,  2.7477367760388098e-03,
        -9.2110386667160222e-03,  1.5895856233272539e-03,  5.2377094526745579e-04, -5.2519154906813107e-04,
         1.7614161838416023e+00, -3.1685800667273778e-02, -2.7278574910696542e-02,  6.3584621477932846e-04,
        -8.0924625745514220e-03, -2.4216573332765355e-02,  1.2203926265283721e-04, -1.8482328475641070e-04,
                                                         
         1.3294757533414430e+00, -8.9788893137442191e-02, -6.2911497442604464e-02,  3.9614170768292408e-05,
        -7.1100499232079313e-03,  4.6679346543524623e-03,  2.4956613702379482e-04, -2.7212527769268630e-04,
         1.1558500674443986e+00, -9.5187518915131830e-02, -3.6801896226934798e-02, -3.8875289024225706e-03,
        -6.1608141969413084e-03,  5.0780393003259120e-03,  2.9847536506066310e-04, -8.3836836755803027e-04,
         1.1608685714961731e+00, -1.2037284313135867e-01,  5.3537968971594171e-02, -1.0607159984610976e-02,
        -4.5111155521692517e-03,  1.5797402263497157e-02,  6.5397859158026453e-04, -8.0268670457479710e-04,
         1.4515046707409809e+00, -1.5802986276729319e-01,  7.3764923340376501e-02, -7.7325221486257874e-03,
        -2.2739467835625805e-03, -1.5570616292196462e-02,  6.3765139919740119e-04,  1.8321877472712586e-03,
         1.4115464704823009e+00, -1.5776040893168661e-01, -9.9620570722381063e-02,  7.8173884689786872e-03,
        -2.1636612802878735e-03, -1.7728464227987822e-02, -5.7397803419405682e-04,  1.7774214977044585e-03,
         1.0684562036485103e+00, -1.2101078886542833e-01, -3.5201132471139970e-02,  9.9394198301683799e-03,
        -4.5077782114689643e-03,  3.1273725348626304e-02, -7.7939850703529119e-04, -9.0314366050547369e-04,
         1.3060178631336949e+00, -9.8310629669391220e-02,  1.3131146816260625e-01,  3.4378558002397842e-03,
        -6.6184622280531843e-03, -5.2033263638922255e-04, -4.3920547811384877e-04, -6.9297032614728052e-04,
         1.5457214754324098e+00, -9.2030059599753106e-02, -2.4079263611516445e-02,  9.9293276550434868e-04,
        -7.4607497917453776e-03, -2.4657338466129785e-02, -4.7089473518925374e-05, -6.9652583640937327e-05,
                                                         
         9.2555564511589228e-01, -1.4283166346900242e-01, -5.7087203270645939e-02,  4.2107167636395108e-03,
        -6.6585943733794293e-03,  4.3810045656671034e-03,  9.3715478409907424e-04,  1.0646611375162348e-04,
         7.5643742119286161e-01, -1.3194920980432301e-01, -4.3151229134044268e-02,  9.9766889501250069e-04,
        -3.4774021071347677e-03,  2.3686176801920480e-03,  8.9950742716123533e-04, -7.2591693417441388e-04,
         7.0097359003902948e-01, -1.4059912442543290e-01,  2.5969251836348382e-02, -5.4595811792828132e-03,
        -9.8784417722099275e-04,  1.3863694644921844e-02,  5.3783951383765062e-04, -3.1373990954366011e-04,
         8.9395264296783794e-01, -1.6114696623034006e-01,  5.3336355793392118e-02, -4.9642255681519613e-03,
        -1.6240840546821416e-04, -1.1009301949378315e-02, -6.1273948817443689e-05,  8.0128831641332842e-04,
         8.5858667148480194e-01, -1.5815426669652288e-01, -7.7423345159838203e-02,  6.2531566299462507e-03,
         4.2102309020295776e-04, -1.3850860419404911e-02,  3.9811827989688620e-04,  4.6131410499162032e-04,
         6.0862841609240648e-01, -1.3886189996819115e-01, -1.1804099479237122e-02,  3.6257448687969642e-03,
        -1.6305922506312725e-04,  2.8818588277168312e-02, -7.3533836184465481e-04, -5.1433038859821081e-04,
         8.8573073934651858e-01, -1.4220407629839935e-01,  1.3180886462596375e-01, -5.2114498588326360e-03,
        -4.9622366030995808e-03, -2.1500466509339223e-03, -2.0354679892534604e-03, -2.4794549885219248e-04,
         1.1259773284051056e+00, -1.5150225754918639e-01, -2.1648595211938699e-02,  5.4796944887228096e-04,
        -8.3847823258309216e-03, -2.4261809779577131e-02,  5.9460294920692505e-05,  2.9801117729436375e-04,
                                                         
         3.4754908174043120e-01, -1.9147096450810924e-01, -2.8048110764785136e-02,  1.3769366894036281e-02,
        -7.1159030934106497e-03,  2.2827046419849114e-03,  1.8777907284160782e-03, -1.3179201395300791e-03,
         2.6018094650848222e-01, -1.5014541606507417e-01, -2.4623528318145321e-02,  1.1377222883164596e-02,
        -5.4317382757140756e-05,  5.5564633406385939e-04,  2.1992176825353384e-03, -3.2080256054575182e-04,
         2.2061875008365400e-01, -1.3141643256470487e-01,  9.6631517600846636e-03, -5.1140605928807351e-03,
         3.1310264886016273e-03,  6.6601405906146901e-03, -3.6015854094482351e-04, -3.8452339628325693e-03,
         3.1728477488866830e-01, -1.7738358224001571e-01,  3.1338192874284528e-02, -1.2875214739490477e-02,
        -4.4933329941175626e-03, -4.8107149369193128e-03, -4.0417674588015843e-03,  2.7774675638252209e-03,
         2.9856868259017150e-01, -1.6904725884004720e-01, -4.4358254629948746e-02,  1.8966632025180070e-02,
        -2.5805304029333561e-03, -6.5259204756554336e-03,  5.1461245497283941e-03,  3.7677419433300198e-03,
         1.6472836762602511e-01, -1.0615968567242212e-01, -6.4622802488490489e-03,  2.0693364430829689e-03,
         8.5621733627069246e-03,  1.3963870956119893e-02,  1.2871184688642873e-03, -8.0620446554448624e-03,
         3.2567896066608065e-01, -1.8138329229246203e-01,  7.9694807939027446e-02, -3.4130327869920409e-02,
        -5.1491649123716690e-03, -1.2897504262269553e-03, -9.2033633129309355e-03,  7.4463775576960654e-04,
         4.5758071423995794e-01, -2.4387817433520531e-01, -1.7203978611668343e-02,  5.9370449568276907e-03,
        -1.5729094905946251e-02, -1.1872819639639669e-02,  3.0950378831332307e-03,  6.8547756149861657e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174042964e-01, -2.8048110764784386e-02,  1.9147096450811057e-01, -1.3769366894036857e-02,  2.2827046419850081e-03, -7.1159030934107581e-03,  1.3179201395300958e-03,  1.8777907284162560e-03,
         9.2555564511589383e-01, -5.7087203270645925e-02,  1.4283166346900303e-01, -4.2107167636395723e-03,  4.3810045656670115e-03, -6.6585943733795568e-03, -1.0646611375170161e-04,  9.3715478409909473e-04,
         1.3294757533414452e+00, -6.2911497442604436e-02,  8.9788893137442344e-02, -3.9614170768272533e-05,  4.6679346543522706e-03, -7.1100499232081170e-03,  2.7212527769271232e-04,  2.4956613702375823e-04,
         1.5410607675732264e+00, -6.1609756908925783e-02,  3.1259713850908692e-02,  4.0278881947991574e-04,  5.4767589025949666e-03, -7.9698971022144141e-03,  1.9484961975729643e-04, -5.1276054205724837e-05,
         1.5410607675732255e+00, -6.1609756908925964e-02, -3.1259713850908942e-02, -4.0278881947993422e-04,  5.4767589025949874e-03, -7.9698971022144575e-03, -1.9484961975728526e-04, -5.1276054205658064e-05,
         1.3294757533414445e+00, -6.2911497442604505e-02, -8.9788893137442580e-02,  3.9614170768346550e-05,  4.6679346543523348e-03, -7.1100499232081196e-03, -2.7212527769269334e-04,  2.4956613702381553e-04,
         9.2555564511589206e-01, -5.7087203270645800e-02, -1.4283166346900306e-01,  4.2107167636395082e-03,  4.3810045656669195e-03, -6.6585943733794492e-03,  1.0646611375155487e-04,  9.3715478409906871e-04,
         3.4754908174042604e-01, -2.8048110764785483e-02, -1.9147096450811166e-01,  1.3769366894035920e-02,  2.2827046419861504e-03, -7.1159030934108405e-03, -1.3179201395291872e-03,  1.8777907284159925e-03,
                                                         
         2.6018094650848411e-01, -2.4623528318145346e-02,  1.5014541606507359e-01, -1.1377222883164428e-02,  5.5564633406300114e-04, -5.4317382757119831e-05,  3.2080256054629734e-04,  2.1992176825352508e-03,
         7.5643742119286284e-01, -4.3151229134044358e-02,  1.3194920980432326e-01, -9.9766889501257138e-04,  2.3686176801920706e-03, -3.4774021071348232e-03,  7.2591693417436444e-04,  8.9950742716123674e-04,
         1.1558500674444006e+00, -3.6801896226935180e-02,  9.5187518915131983e-02,  3.8875289024226109e-03,  5.0780393003259051e-03, -6.1608141969413475e-03,  8.3836836755805694e-04,  2.9847536506072300e-04,
         1.3895291108288046e+00, -2.3292866819651030e-02,  3.6049179432956914e-02,  2.8177193660515751e-03,  7.4784863973144733e-03, -9.0097499591925388e-03,  5.4753041006382476e-04, -5.4908327268816091e-04,
         1.3895291108288044e+00, -2.3292866819650922e-02, -3.6049179432957101e-02, -2.8177193660515378e-03,  7.4784863973144950e-03, -9.0097499591925163e-03, -5.4753041006382487e-04, -5.4908327268817989e-04,
         1.1558500674444003e+00, -3.6801896226934715e-02, -9.5187518915131886e-02, -3.8875289024224817e-03,  5.0780393003259268e-03, -6.1608141969413258e-03, -8.3836836755802658e-04,  2.9847536506067546e-04,
         7.5643742119286161e-01, -4.3151229134043845e-02, -1.3194920980432370e-01,  9.9766889501253452e-04,  2.3686176801920966e-03, -3.4774021071348306e-03, -7.2591693417442061e-04,  8.9950742716124173e-04,
         2.6018094650848284e-01, -2.4623528318144457e-02, -1.5014541606507342e-01,  1.1377222883164968e-02,  5.5564633406298997e-04, -5.4317382757213696e-05, -3.2080256054627907e-04,  2.1992176825354850e-03,
                                                         
         2.2061875008365345e-01,  9.6631517600846272e-03,  1.3141643256470570e-01,  5.1140605928805252e-03,  6.6601405906149182e-03,  3.1310264886015996e-03,  3.8452339628324505e-03, -3.6015854094478215e-04,
         7.0097359003903059e-01,  2.5969251836348413e-02,  1.4059912442543326e-01,  5.4595811792830344e-03,  1.3863694644921848e-02, -9.8784417722103503e-04,  3.1373990954367497e-04,  5.3783951383765680e-04,
         1.1608685714961748e+00,  5.3537968971594505e-02,  1.2037284313135892e-01,  1.0607159984610891e-02,  1.5797402263497132e-02, -4.5111155521692725e-03,  8.0268670457474148e-04,  6.5397859158021802e-04,
         1.4814229261146545e+00,  9.0436045201507859e-02,  5.3722083116000588e-02,  7.5759200922298322e-03,  1.8391205195993732e-02, -1.3014031031221736e-02,  6.9484611672702631e-04, -1.7627894821588743e-03,
         1.4814229261146548e+00,  9.0436045201508108e-02, -5.3722083116000581e-02, -7.5759200922297211e-03,  1.8391205195993788e-02, -1.3014031031221748e-02, -6.9484611672701677e-04, -1.7627894821588838e-03,
         1.1608685714961748e+00,  5.3537968971594407e-02, -1.2037284313135886e-01, -1.0607159984611132e-02,  1.5797402263497066e-02, -4.5111155521693393e-03, -8.0268670457477877e-04,  6.5397859158024664e-04,
         7.0097359003903081e-01,  2.5969251836348378e-02, -1.4059912442543340e-01, -5.4595811792828314e-03,  1.3863694644921699e-02, -9.8784417722108534e-04, -3.1373990954372105e-04,  5.3783951383762341e-04,
         2.2061875008365178e-01,  9.6631517600844936e-03, -1.3141643256470645e-01, -5.1140605928809068e-03,  6.6601405906156555e-03,  3.1310264886016746e-03, -3.8452339628318906e-03, -3.6015854094490277e-04,
                                                         
         3.1728477488866852e-01,  3.1338192874284548e-02,  1.7738358224001591e-01,  1.2875214739490581e-02, -4.8107149369194429e-03, -4.4933329941175887e-03, -2.7774675638251212e-03, -4.0417674588016208e-03,
         8.9395264296783927e-01,  5.3336355793392055e-02,  1.6114696623034069e-01,  4.9642255681518581e-03, -1.1009301949378283e-02, -1.6240840546819605e-04, -8.0128831641333893e-04, -6.1273948817444122e-05,
         1.4515046707409838e+00,  7.3764923340376473e-02,  1.5802986276729358e-01,  7.7325221486260034e-03, -1.5570616292196622e-02, -2.2739467835626980e-03, -1.8321877472712773e-03,  6.3765139919743480e-04,
         1.9036696745231823e+00,  1.0223309369326505e-01,  8.2014576705607287e-02,  6.0293706608141795e-03, -2.1202408094434320e-02, -1.8550660598603982e-02, -1.4193287657706390e-03, -1.4337850889724951e-03,
         1.9036696745231827e+00,  1.0223309369326498e-01, -8.2014576705607231e-02, -6.0293706608142350e-03, -2.1202408094434320e-02, -1.8550660598604118e-02,  1.4193287657706867e-03, -1.4337850889725031e-03,
         1.4515046707409842e+00,  7.3764923340376751e-02, -1.5802986276729358e-01, -7.7325221486257995e-03, -1.5570616292196521e-02, -2.2739467835627301e-03,  1.8321877472712252e-03,  6.3765139919739501e-04,
         8.9395264296783961e-01,  5.3336355793392465e-02, -1.6114696623034064e-01, -4.9642255681518858e-03, -1.1009301949378327e-02, -1.6240840546830605e-04,  8.0128831641336322e-04, -6.1273948817420066e-05,
         3.1728477488866941e-01,  3.1338192874284715e-02, -1.7738358224001569e-01, -1.2875214739490671e-02, -4.8107149369198323e-03, -4.4933329941176468e-03,  2.7774675638248801e-03, -4.0417674588015757e-03,
                                                         
         2.9856868259017222e-01, -4.4358254629948406e-02,  1.6904725884004743e-01, -1.8966632025180216e-02, -6.5259204756554969e-03, -2.5805304029333379e-03, -3.7677419433299738e-03,  5.1461245497284617e-03,
         8.5858667148480361e-01, -7.7423345159837897e-02,  1.5815426669652335e-01, -6.2531566299461995e-03, -1.3850860419404974e-02,  4.2102309020287818e-04, -4.6131410499162384e-04,  3.9811827989688772e-04,
         1.4115464704823038e+00, -9.9620570722380813e-02,  1.5776040893168708e-01, -7.8173884689787462e-03, -1.7728464227987965e-02, -2.1636612802880054e-03, -1.7774214977044575e-03, -5.7397803419405899e-04,
         1.8635410509569048e+00, -1.2804640133826678e-01,  8.2038435164304063e-02, -6.0050462696768147e-03, -2.3251244190819049e-02, -1.8554396836116767e-02, -1.4111570005111832e-03,  1.4316279712387229e-03,
         1.8635410509569053e+00, -1.2804640133826684e-01, -8.2038435164303827e-02,  6.0050462696768702e-03, -2.3251244190819070e-02, -1.8554396836116722e-02,  1.4111570005111642e-03,  1.4316279712387932e-03,
         1.4115464704823049e+00, -9.9620570722380869e-02, -1.5776040893168702e-01,  7.8173884689786716e-03, -1.7728464227987951e-02, -2.1636612802879928e-03,  1.7774214977044958e-03, -5.7397803419400651e-04,
         8.5858667148480450e-01, -7.7423345159837939e-02, -1.5815426669652330e-01,  6.2531566299462082e-03, -1.3850860419404939e-02,  4.2102309020288029e-04,  4.6131410499163864e-04,  3.9811827989687552e-04,
         2.9856868259017277e-01, -4.4358254629948260e-02, -1.6904725884004765e-01,  1.8966632025180313e-02, -6.5259204756555290e-03, -2.5805304029333370e-03,  3.7677419433299339e-03,  5.1461245497284478e-03,
                                                         
         1.6472836762602650e-01, -6.4622802488492528e-03,  1.0615968567242214e-01, -2.0693364430828132e-03,  1.3963870956119633e-02,  8.5621733627069419e-03,  8.0620446554450047e-03,  1.2871184688642244e-03,
         6.0862841609240859e-01, -1.1804099479237192e-02,  1.3886189996819157e-01, -3.6257448687969846e-03,  2.8818588277168294e-02, -1.6305922506323095e-04,  5.1433038859822794e-04, -7.3533836184467519e-04,
         1.0684562036485135e+00, -3.5201132471140248e-02,  1.2101078886542867e-01, -9.9394198301683955e-03,  3.1273725348626248e-02, -4.5077782114690086e-03,  9.0314366050544658e-04, -7.7939850703530583e-04,
         1.3908819207316665e+00, -7.0885965719393249e-02,  5.4053805794843257e-02, -7.4371696780805782e-03,  3.4035995955904931e-02, -1.3096490503921115e-02,  6.9165401818148148e-04,  1.7194957188661963e-03,
         1.3908819207316667e+00, -7.0885965719393138e-02, -5.4053805794842896e-02,  7.4371696780806146e-03,  3.4035995955904966e-02, -1.3096490503921015e-02, -6.9165401818150035e-04,  1.7194957188661772e-03,
         1.0684562036485146e+00, -3.5201132471140234e-02, -1.2101078886542882e-01,  9.9394198301683487e-03,  3.1273725348626262e-02, -4.5077782114689956e-03, -9.0314366050543834e-04, -7.7939850703529628e-04,
         6.0862841609240925e-01, -1.1804099479237469e-02, -1.3886189996819157e-01,  3.6257448687969478e-03,  2.8818588277168287e-02, -1.6305922506316907e-04, -5.1433038859822100e-04, -7.3533836184460841e-04,
         1.6472836762602747e-01, -6.4622802488492458e-03, -1.0615968567242191e-01,  2.0693364430829832e-03,  1.3963870956119570e-02,  8.5621733627069974e-03, -8.0620446554450394e-03,  1.2871184688642598e-03,
                                                         
         3.2567896066608137e-01,  7.9694807939027321e-02,  1.8138329229246233e-01,  3.4130327869920436e-02, -1.2897504262272103e-03, -5.1491649123716135e-03, -7.4463775576943090e-04, -9.2033633129308540e-03,
         8.8573073934652091e-01,  1.3180886462596381e-01,  1.4220407629839993e-01,  5.2114498588327497e-03, -2.1500466509339136e-03, -4.9622366030995175e-03,  2.4794549885218521e-04, -2.0354679892533663e-03,
         1.3060178631336981e+00,  1.3131146816260650e-01,  9.8310629669391358e-02, -3.4378558002398189e-03, -5.2033263638920748e-04, -6.6184622280531470e-03,  6.9297032614725450e-04, -4.3920547811381597e-04,
         1.5458079707863022e+00,  1.1844442680216062e-01,  3.6785251742611789e-02, -2.7477367760390123e-03,  1.5895856233272446e-03, -9.2110386667161575e-03,  5.2519154906816750e-04,  5.2377094526740299e-04,
         1.5458079707863026e+00,  1.1844442680216057e-01, -3.6785251742611727e-02,  2.7477367760389013e-03,  1.5895856233273001e-03, -9.2110386667162338e-03, -5.2519154906811025e-04,  5.2377094526734563e-04,
         1.3060178631336978e+00,  1.3131146816260572e-01, -9.8310629669391830e-02,  3.4378558002395136e-03, -5.2033263638935265e-04, -6.6184622280532918e-03, -6.9297032614738189e-04, -4.3920547811387062e-04,
         8.8573073934651991e-01,  1.3180886462596297e-01, -1.4220407629839971e-01, -5.2114498588324946e-03, -2.1500466509341257e-03, -4.9622366030995643e-03, -2.4794549885213566e-04, -2.0354679892534569e-03,
         3.2567896066608099e-01,  7.9694807939026988e-02, -1.8138329229246228e-01, -3.4130327869920402e-02, -1.2897504262269091e-03, -5.1491649123716274e-03,  7.4463775576965577e-04, -9.2033633129309320e-03,
                                                         
         4.5758071423995789e-01, -1.7203978611669109e-02,  2.4387817433520595e-01, -5.9370449568272301e-03, -1.1872819639640254e-02, -1.5729094905946386e-02, -6.8547756149857598e-03,  3.0950378831330659e-03,
         1.1259773284051073e+00, -2.1648595211938900e-02,  1.5150225754918717e-01, -5.4796944887231913e-04, -2.4261809779577124e-02, -8.3847823258308921e-03, -2.9801117729441520e-04,  5.9460294920611305e-05,
         1.5457214754324125e+00, -2.4079263611516789e-02,  9.2030059599753189e-02, -9.9293276550425479e-04, -2.4657338466129764e-02, -7.4607497917453663e-03,  6.9652583640970205e-05, -4.7089473518950642e-05,
         1.7614161838416049e+00, -2.7278574910696709e-02,  3.1685800667273785e-02, -6.3584621477909753e-04, -2.4216573332765473e-02, -8.0924625745515469e-03,  1.8482328475636194e-04,  1.2203926265291815e-04,
         1.7614161838416049e+00, -2.7278574910696782e-02, -3.1685800667273903e-02,  6.3584621477904201e-04, -2.4216573332765518e-02, -8.0924625745517915e-03, -1.8482328475635242e-04,  1.2203926265290857e-04,
         1.5457214754324113e+00, -2.4079263611516525e-02, -9.2030059599753952e-02,  9.9293276550453235e-04, -2.4657338466129847e-02, -7.4607497917455797e-03, -6.9652583640970558e-05, -4.7089473518941840e-05,
         1.1259773284051051e+00, -2.1648595211938775e-02, -1.5150225754918667e-01,  5.4796944887201382e-04, -2.4261809779577238e-02, -8.3847823258309407e-03,  2.9801117729432342e-04,  5.9460294920663889e-05,
         4.5758071423995744e-01, -1.7203978611668745e-02, -2.4387817433520526e-01,  5.9370449568277878e-03, -1.1872819639639924e-02, -1.5729094905946272e-02,  6.8547756149861354e-03,  3.0950378831332199e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
    }
  }

  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_cu);
  gkyl_array_release(phi_cu);
  gkyl_array_release(perbuff);
  gkyl_fem_poisson_release(poisson_cu);
}

void test_1x_p1_periodicx() {
  int cells[] = {16}; // MF 2022/05/31: N=8 is broken.
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(1, &cells[0], bc_tv);
}
void test_1x_p1_dirichletx() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv);
}

void test_1x_p2_periodicx() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(2, &cells[0], bc_tv);
}
void test_1x_p2_dirichletx() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv);
}

void test_2x_p1_periodicx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, &cells[0], bc_tv);
}
void test_2x_p1_dirichletx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv);
}
void test_2x_p1_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv);
}
void test_2x_p1_periodicx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv);
}
void test_2x_p2_periodicx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(2, &cells[0], bc_tv);
}
void test_2x_p2_dirichletx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv);
}
void test_2x_p2_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv);
}
void test_2x_p2_periodicx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv);
}

void gpu_test_1x_p1_periodicx() {
  int cells[] = {16}; // MF 2022/05/31: N=8 is broken.
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  gpu_test_1x(1, &cells[0], bc_tv);
}
void gpu_test_1x_p1_dirichletx() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  gpu_test_1x(1, &cells[0], bc_tv);
}

void gpu_test_1x_p2_periodicx() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  gpu_test_1x(2, &cells[0], bc_tv);
}
void gpu_test_1x_p2_dirichletx() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  gpu_test_1x(2, &cells[0], bc_tv);
}

void gpu_test_2x_p1_periodicx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  gpu_test_2x(1, &cells[0], bc_tv);
}
void gpu_test_2x_p1_dirichletx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  gpu_test_2x(1, &cells[0], bc_tv);
}
void gpu_test_2x_p1_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  gpu_test_2x(1, &cells[0], bc_tv);
}
void gpu_test_2x_p1_periodicx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  gpu_test_2x(1, &cells[0], bc_tv);
}
void gpu_test_2x_p2_periodicx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  gpu_test_2x(2, &cells[0], bc_tv);
}
void gpu_test_2x_p2_dirichletx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  gpu_test_2x(2, &cells[0], bc_tv);
}
void gpu_test_2x_p2_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  gpu_test_2x(2, &cells[0], bc_tv);
}
void gpu_test_2x_p2_periodicx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  gpu_test_2x(2, &cells[0], bc_tv);
}


TEST_LIST = {
  { "test_1x_p1_periodicx", test_1x_p1_periodicx },
  { "test_1x_p1_dirichletx", test_1x_p1_dirichletx },
  { "test_1x_p2_periodicx", test_1x_p2_periodicx },
  { "test_1x_p2_dirichletx", test_1x_p2_dirichletx },
  { "test_2x_p1_periodicx_periodicy", test_2x_p1_periodicx_periodicy },
  { "test_2x_p1_dirichletx_dirichlety", test_2x_p1_dirichletx_dirichlety },
  { "test_2x_p1_dirichletx_periodicy", test_2x_p1_dirichletx_periodicy },
  { "test_2x_p1_periodicx_dirichlety", test_2x_p1_periodicx_dirichlety },
  { "test_2x_p2_periodicx_periodicy", test_2x_p2_periodicx_periodicy },
  { "test_2x_p2_dirichletx_dirichlety", test_2x_p2_dirichletx_dirichlety },
  { "test_2x_p2_dirichletx_periodicy", test_2x_p2_dirichletx_periodicy },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_1x_p1_periodicx", gpu_test_1x_p1_periodicx },
  { "gpu_test_1x_p1_dirichletx", gpu_test_1x_p1_dirichletx },
  { "gpu_test_1x_p2_periodicx", gpu_test_1x_p2_periodicx },
  { "gpu_test_1x_p2_dirichletx", gpu_test_1x_p2_dirichletx },
  { "gpu_test_2x_p1_periodicx_periodicy", gpu_test_2x_p1_periodicx_periodicy },
  { "gpu_test_2x_p1_dirichletx_dirichlety", gpu_test_2x_p1_dirichletx_dirichlety },
  { "gpu_test_2x_p1_dirichletx_periodicy", gpu_test_2x_p1_dirichletx_periodicy },
  { "gpu_test_2x_p1_periodicx_dirichlety", gpu_test_2x_p1_periodicx_dirichlety },
  { "gpu_test_2x_p2_periodicx_periodicy", gpu_test_2x_p2_periodicx_periodicy },
  { "gpu_test_2x_p2_dirichletx_dirichlety", gpu_test_2x_p2_dirichletx_dirichlety },
  { "gpu_test_2x_p2_dirichletx_periodicy", gpu_test_2x_p2_dirichletx_periodicy },
#endif
  { NULL, NULL },
};
