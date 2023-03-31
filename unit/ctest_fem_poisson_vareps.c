// Test the FEM Poisson solver with spatially varying permittivity.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson.h>

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
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&grid, basis, &bcs, 0, eps, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_set_rhs(poisson, rho_cu);
  else
    gkyl_fem_poisson_set_rhs(poisson, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_solve(poisson, phi);
  }

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
//  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_vareps_2x_phi_1.gkyl");

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
         4.8206219913678056e-01,  8.2198705544354114e-02, -3.3829630915641395e-02, -1.7439753345153505e-15,
         2.8200758611238580e-01,  8.2198705544354322e-02, -8.1671953766619371e-02,  1.8549976369778664e-15,
        -9.1236084201473666e-04,  8.2198705544354267e-02, -8.1671953766618413e-02, -1.8792837656415417e-15,
        -2.0096697386640955e-01,  8.2198705544354267e-02, -3.3829630915642360e-02,  1.8850661772281309e-15,
        -2.0096697386640949e-01,  8.2198705544354322e-02,  3.3829630915642388e-02, -1.8596235662471376e-15,
        -9.1236084201488932e-04,  8.2198705544354156e-02,  8.1671953766618288e-02,  1.7740438747656152e-15,
         2.8200758611238541e-01,  8.2198705544354073e-02,  8.1671953766619398e-02, -1.8179902028236942e-15,
         4.8206219913678033e-01,  8.2198705544354031e-02,  3.3829630915641465e-02,  1.7902346272080653e-15,
         6.8394153384520473e-01,  3.4356382693377228e-02, -3.3829630915643179e-02,  7.2164496600635195e-16,
         4.8388692082081014e-01,  3.4356382693377083e-02, -8.1671953766617539e-02, -8.0491169285323869e-16,
         2.0096697386640971e-01,  3.4356382693377173e-02, -8.1671953766620189e-02,  8.5579691481522510e-16,
         9.1236084201498646e-04,  3.4356382693377222e-02, -3.3829630915640535e-02, -8.3382375078618551e-16,
         9.1236084201504197e-04,  3.4356382693377187e-02,  3.3829630915640577e-02,  8.1647651602641744e-16,
         2.0096697386640971e-01,  3.4356382693377353e-02,  8.1671953766620092e-02, -7.1933200137171624e-16,
         4.8388692082080992e-01,  3.4356382693377381e-02,  8.1671953766617511e-02,  7.4014868308343792e-16,
         6.8394153384520462e-01,  3.4356382693377326e-02,  3.3829630915643213e-02, -7.7715611723760978e-16,
         4.0102158689080447e-01, -1.9770029022661495e-01, -3.3829630915640618e-02,  7.4940054162198086e-16,
         2.0096697386640955e-01, -1.9770029022661495e-01, -8.1671953766620231e-02, -7.4940054162198086e-16,
        -8.1952973087990932e-02, -1.9770029022661506e-01, -8.1671953766617539e-02,  6.6613381477509412e-16,
        -2.8200758611238574e-01, -1.9770029022661514e-01, -3.3829630915643193e-02, -6.9620235502535874e-16,
        -2.8200758611238558e-01, -1.9770029022661503e-01,  3.3829630915643276e-02,  7.4130516540075573e-16,
        -8.1952973087990738e-02, -1.9770029022661514e-01,  8.1671953766617511e-02, -7.8640797577615281e-16,
         2.0096697386640960e-01, -1.9770029022661506e-01,  8.1671953766620176e-02,  8.0028576358396727e-16,
         4.0102158689080453e-01, -1.9770029022661495e-01,  3.3829630915640618e-02, -7.3089682454489489e-16,
        -2.0096697386640988e-01, -1.4985796737563822e-01, -3.3829630915642478e-02, -1.8341809552661447e-15,
        -4.0102158689080469e-01, -1.4985796737563817e-01, -8.1671953766618302e-02,  1.8815967302761773e-15,
        -6.8394153384520517e-01, -1.4985796737563800e-01, -8.1671953766619468e-02, -1.7902346272080653e-15,
        -8.8399614686960004e-01, -1.4985796737563797e-01, -3.3829630915641298e-02,  1.7994864857466085e-15,
        -8.8399614686960004e-01, -1.4985796737563817e-01,  3.3829630915641257e-02, -1.8966310004013097e-15,
        -6.8394153384520540e-01, -1.4985796737563817e-01,  8.1671953766619398e-02,  1.8827532125934953e-15,
        -4.0102158689080503e-01, -1.4985796737563822e-01,  8.1671953766618288e-02, -1.9058828589398526e-15,
        -2.0096697386641005e-01, -1.4985796737563828e-01,  3.3829630915642582e-02,  1.8723448717376344e-15,
        -2.0096697386640988e-01,  1.4985796737563822e-01, -3.3829630915642443e-02,  1.8549976369778664e-15,
        -4.0102158689080469e-01,  1.4985796737563817e-01, -8.1671953766618330e-02, -1.8873791418627665e-15,
        -6.8394153384520517e-01,  1.4985796737563803e-01, -8.1671953766619440e-02,  1.8133642735544230e-15,
        -8.8399614686960015e-01,  1.4985796737563789e-01, -3.3829630915641430e-02, -1.8873791418627665e-15,
        -8.8399614686960004e-01,  1.4985796737563817e-01,  3.3829630915641465e-02,  2.0354088784794543e-15,
        -6.8394153384520517e-01,  1.4985796737563836e-01,  8.1671953766619385e-02, -1.9058828589398526e-15,
        -4.0102158689080475e-01,  1.4985796737563839e-01,  8.1671953766618288e-02,  1.9151347174783954e-15,
        -2.0096697386640994e-01,  1.4985796737563836e-01,  3.3829630915642485e-02, -1.9243865760169387e-15,
         4.0102158689080475e-01,  1.9770029022661512e-01, -3.3829630915640639e-02, -8.0953762212251021e-16,
         2.0096697386640977e-01,  1.9770029022661501e-01, -8.1671953766620259e-02,  7.7021722333370255e-16,
        -8.1952973087990835e-02,  1.9770029022661506e-01, -8.1671953766617594e-02, -7.4477461235270934e-16,
        -2.8200758611238586e-01,  1.9770029022661520e-01, -3.3829630915643297e-02,  8.0028576358396727e-16,
        -2.8200758611238563e-01,  1.9770029022661501e-01,  3.3829630915643435e-02, -8.9743027823866847e-16,
        -8.1952973087990530e-02,  1.9770029022661489e-01,  8.1671953766617497e-02,  8.1416355139178173e-16,
         2.0096697386640983e-01,  1.9770029022661492e-01,  8.1671953766620217e-02, -7.9565983431469575e-16,
         4.0102158689080480e-01,  1.9770029022661501e-01,  3.3829630915640632e-02,  8.5232746786327148e-16,
         6.8394153384520517e-01, -3.4356382693377353e-02, -3.3829630915643276e-02, -7.2164496600635195e-16,
         4.8388692082081025e-01, -3.4356382693377215e-02, -8.1671953766617511e-02,  8.0491169285323869e-16,
         2.0096697386640983e-01, -3.4356382693377201e-02, -8.1671953766620273e-02, -8.0028576358396727e-16,
         9.1236084201498646e-04, -3.4356382693377215e-02, -3.3829630915640535e-02,  7.9565983431469575e-16,
         9.1236084201505585e-04, -3.4356382693377138e-02,  3.3829630915640570e-02, -7.5171350625661657e-16,
         2.0096697386640983e-01, -3.4356382693377187e-02,  8.1671953766620176e-02,  7.2395793064098766e-16,
         4.8388692082081014e-01, -3.4356382693377222e-02,  8.1671953766617580e-02, -7.4477461235270934e-16,
         6.8394153384520517e-01, -3.4356382693377305e-02,  3.3829630915643324e-02,  6.9388939039072303e-16,
         4.8206219913678056e-01, -8.2198705544354142e-02, -3.3829630915641465e-02,  1.7763568394002509e-15,
         2.8200758611238569e-01, -8.2198705544354295e-02, -8.1671953766619343e-02, -1.8596235662471376e-15,
        -9.1236084201470891e-04, -8.2198705544354281e-02, -8.1671953766618413e-02,  1.8688754247856808e-15,
        -2.0096697386640955e-01, -8.2198705544354281e-02, -3.3829630915642395e-02, -1.8654059778337272e-15,
        -2.0096697386640949e-01, -8.2198705544354350e-02,  3.3829630915642436e-02,  1.8336027141074855e-15,
        -9.1236084201469503e-04, -8.2198705544354267e-02,  8.1671953766618316e-02, -1.7948605564773369e-15,
         2.8200758611238569e-01, -8.2198705544354239e-02,  8.1671953766619398e-02,  1.8133642735544230e-15,
         4.8206219913678056e-01, -8.2198705544354170e-02,  3.3829630915641451e-02, -1.7763568394002509e-15,
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
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
}

void test_2x_p1_periodicx_periodicy() {
  int cells[] = {8,8};
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
