// Test the perpendicular FEM Helmholtz/Poisson solver.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson_perp.h>

#define PERP_DIM 2

void evalFunc_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 2.*M_PI;
  double kx = 2.*M_PI/Lx, ky = 2.*M_PI/Ly;
  double sig = 0.3*sqrt(Lx*Lx+Ly*Ly);
  fout[0] = exp(-(pow(kx*x,2)+pow(ky*y,2))/(2.*(sig*sig)));
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
test_fem_poisson_perp(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI,-M_PI}, upper[] = {M_PI,M_PI,M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
      (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_dirichletx_dirichlety, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = PERP_DIM+ceil((pow(3.,PERP_DIM-1)-PERP_DIM)/2);
  struct gkyl_array *eps = use_gpu? mkarr_cu(epsnum*basis.num_basis, localRange_ext.volume)
                                  : mkarr(   epsnum*basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  if (use_gpu) gkyl_array_copy(rho_cu, rho);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_perp_2x_rho_1.gkyl");

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps,        0.*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 2*basis.num_basis);

  // FEM poisson solver.
  gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&grid, basis, &bcs, eps, NULL, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_perp_set_rhs(poisson, rho_cu);
  else
    gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_perp_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_perp_solve(poisson, phi);
  }
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_perp_phi_1.gkyl");

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(perbuff);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }

}

void test_p1_dirichletx_dirichlety() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp(1, cells, bc_tv, false);
}

TEST_LIST = {
  { "test_p1_dirichletx_dirichlety", test_p1_dirichletx_dirichlety },
  { NULL, NULL },
};

