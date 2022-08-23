// Test the projection onto an FEM basis that is continuous in the parallel
// direction.
//
#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_parproj.h>

void evalFunc1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin(2.*M_PI*x);
}

void evalFunc3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double mu[2] = {.2, 0.2};
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu[0],2)+pow(y-mu[1],2))/(2.0*sig*sig))*sin(2.*M_PI*z);
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

// Apply periodic BCs along parallel direction
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, sgr.lower_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(buff->data, fld, sgr.upper_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.lower_ghost[dir]);
}

void
test_1x(int poly_order, const bool isperiodic, bool use_gpu)
{
  double lower[] = {-0.5}, upper[] = {0.5};
  int cells[] = {4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc1x, NULL);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // create an array to use as the multiplicative weight.
  struct gkyl_array *wgt = mkarr(basis.num_basis, localRange_ext.volume);
  gkyl_array_shiftc0(wgt, sqrt(2.)); // Sets wgt=1.
  // device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  }

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *parbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  if (isperiodic) apply_periodic_bc(parbuff, rho, dim-1, skin_ghost);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_parproj_1x_p2_rho_1.gkyl");
  if (use_gpu) gkyl_array_copy(rho_cu, rho);

  // parallel FEM projection method.
  gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&grid, basis, isperiodic,
    false, wgt, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_parproj_set_rhs(parproj, rho_cu);
  else
    gkyl_fem_parproj_set_rhs(parproj, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_parproj_solve(parproj, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_parproj_solve(parproj, phi);
  }

  if (isperiodic)
    apply_periodic_bc(parbuff, phi, dim-1, skin_ghost);
//  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_parproj_1x_p2_phi_1.gkyl");

  if (poly_order == 1) {
    if (!isperiodic) {
      // Solution (checked visually, also checked that phi is actually continuous,
      // and checked that visually looks like results in g2):
      const double sol[8] = {-0.9089542445638024, -0.4554124667453318,
                             -0.8488758876834943,  0.4900987222626481,
                              0.8488758876834943,  0.490098722262648 ,
                              0.9089542445638024, -0.4554124667453318};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 2);
      TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 3);
      TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 4);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
    } else {
      // Solution (checked visually against g2):
      const double sol[8] = {-0.8638954769035714, -0.498770286141977,
                             -0.8638954769035713,  0.498770286141977,
                              0.8638954769035713,  0.498770286141977,
                              0.8638954769035713, -0.498770286141977};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi, 0);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 2);
      TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 3);
      TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 4);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 5);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
    }
  } if (poly_order == 2) {
    if (!isperiodic) {
      // Solution (checked visually against g2):
      const double sol[12] = {-0.9010465429057769, -0.4272439810948228,  0.0875367707148495,
                              -0.9039382020247494,  0.4172269800703625,  0.08107082435707  ,
                               0.9039382020247495,  0.4172269800703625, -0.0810708243570699,
                               0.9010465429057768, -0.4272439810948229, -0.0875367707148495};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 2);
      TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 3);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 4);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
    } else {
      // Solution (checked visually against g2):
      const double sol[12] = {-0.9044201452112453, -0.418896480241106,   0.0799931666307734,
                              -0.9044201452112451,  0.418896480241106,   0.0799931666307734,
                               0.904420145211245 ,  0.418896480241106,  -0.0799931666307734,
                               0.9044201452112451, -0.418896480241106,  -0.0799931666307734};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi, 0);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 2);
      TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 3);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 4);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi, 5);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
    }
  }

  gkyl_fem_parproj_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(wgt);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
  gkyl_array_release(parbuff);

}

void
test_3x(const int poly_order, const bool isperiodic, bool use_gpu)
{
  double lower[] = {-2., -2., -0.5}, upper[] = {2., 2., 0.5};
  int cells[] = {3, 3, 4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc3x, NULL);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // create an array to use as the multiplicative weight.
  struct gkyl_array *wgt = mkarr(basis.num_basis, localRange_ext.volume);
  gkyl_array_shiftc0(wgt, 2.*sqrt(2.)); // Sets wgt=1.
  // device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  }

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *parbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  if (isperiodic) apply_periodic_bc(parbuff, rho, dim-1, skin_ghost);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_parproj_3x_p2_rho_1.gkyl");
  if (use_gpu) gkyl_array_copy(rho_cu, rho);

  // parallel FEM projection method.
  gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&grid, basis, isperiodic,
    false, wgt, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_parproj_set_rhs(parproj, rho_cu);
  else
    gkyl_fem_parproj_set_rhs(parproj, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_parproj_solve(parproj, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_parproj_solve(parproj, phi);
  }

  if (isperiodic)
    apply_periodic_bc(parbuff, phi, dim-1, skin_ghost);
//  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_parproj_3x_p1_phi_1.gkyl");

  if (poly_order == 1) {
    if (!isperiodic) {
      // Solution (checked visually, also checked that phi is actually continuous,
      // and checked that visually looks like results in g2):
      const double sol[96] = {
        // idx = [0,1,:]
        -2.9175130738000619e-04, -2.9175013189013625e-04,
        -2.0243038272814358e-04, -1.4617587558971766e-04,
        -2.0242956711958172e-04, -1.4617528663516419e-04,
        -1.0142349903065383e-04, -1.0142309038708445e-04,
        -2.7246767537113221e-04, -2.7246657757659220e-04,
        -1.8905051806532438e-04,  1.5730928572099475e-04,
        -1.8904975636533129e-04,  1.5730865190902163e-04,
         1.0914836749545332e-04,  1.0914792772775724e-04,
         2.7246767537113210e-04,  2.7246657757659220e-04,
         1.8905051806532432e-04,  1.5730928572099478e-04,
         1.8904975636533134e-04,  1.5730865190902163e-04,
         1.0914836749545326e-04,  1.0914792772775718e-04,
         2.9175130738000619e-04,  2.9175013189013631e-04,
         2.0243038272814358e-04, -1.4617587558971761e-04,
         2.0242956711958167e-04, -1.4617528663516419e-04,
        -1.0142349903065386e-04, -1.0142309038708456e-04,
        // idx = [1,0,:]
        -2.9175130738000619e-04, -2.0243038272814342e-04,
        -2.9175013189013631e-04, -1.4617587558971766e-04,
        -2.0242956711958167e-04, -1.0142349903065388e-04,
        -1.4617528663516419e-04, -1.0142309038708458e-04,
        -2.7246767537113221e-04, -1.8905051806532432e-04,
        -2.7246657757659220e-04,  1.5730928572099475e-04,
        -1.8904975636533140e-04,  1.0914836749545325e-04,
         1.5730865190902168e-04,  1.0914792772775730e-04,
         2.7246767537113210e-04,  1.8905051806532430e-04,
         2.7246657757659220e-04,  1.5730928572099478e-04,
         1.8904975636533134e-04,  1.0914836749545328e-04,
         1.5730865190902160e-04,  1.0914792772775720e-04,
         2.9175130738000625e-04,  2.0243038272814347e-04,
         2.9175013189013631e-04, -1.4617587558971761e-04,
         2.0242956711958167e-04, -1.0142349903065384e-04,
        -1.4617528663516416e-04, -1.0142309038708461e-04,
        // idx = [1,2,:]
        -1.9756145041229078e-02, -1.3707715786581542e-02,
         1.9753708637470873e-02, -9.8984022577751872e-03,
         1.3706025298371970e-02, -6.8679635934884792e-03,
         9.8971815487547514e-03,  6.8671166098145447e-03,
        -1.8450340332725093e-02, -1.2801688837518459e-02,
         1.8448064965823664e-02,  1.0652308957739044e-02,
         1.2800110084035459e-02,  7.3910584964231056e-03,
        -1.0650995274046001e-02, -7.3901470026747126e-03,
         1.8450340332725097e-02,  1.2801688837518463e-02,
        -1.8448064965823671e-02,  1.0652308957739044e-02,
        -1.2800110084035459e-02,  7.3910584964231134e-03,
        -1.0650995274045996e-02, -7.3901470026747221e-03,
         1.9756145041229078e-02,  1.3707715786581542e-02,
        -1.9753708637470870e-02, -9.8984022577751907e-03,
        -1.3706025298371980e-02, -6.8679635934884914e-03,
         9.8971815487547583e-03,  6.8671166098145560e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx1[] = {2,1,k+1};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[32+k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[32+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx2[] = {2,3,k+1};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[64+k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[64+k*basis.num_basis+m], idx2[0], idx0[2], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    } else {
      // Solution (checked visually against g2):
      const double sol[96] = {
        // idx = [0,0,:]
        -1.8667872535731117e-07, -1.8667797321346457e-07,
        -1.8667797321346459e-07, -1.0777901233701990e-07,
        -1.8667722107264865e-07, -1.0777857808656757e-07,
        -1.0777857808656758e-07, -1.0777814383786501e-07,
        -1.8667872535731120e-07, -1.8667797321346457e-07,
        -1.8667797321346459e-07,  1.0777901233701988e-07,
        -1.8667722107264867e-07,  1.0777857808656758e-07,
         1.0777857808656763e-07,  1.0777814383786508e-07,
         1.8667872535731120e-07,  1.8667797321346457e-07,
         1.8667797321346459e-07,  1.0777901233701986e-07,
         1.8667722107264859e-07,  1.0777857808656754e-07,
         1.0777857808656751e-07,  1.0777814383786489e-07,
         1.8667872535731122e-07,  1.8667797321346457e-07,
         1.8667797321346459e-07, -1.0777901233701986e-07,
         1.8667722107264859e-07, -1.0777857808656754e-07,
        -1.0777857808656757e-07, -1.0777814383786493e-07,
        // ...
        // idx = [1,1,:]
        -4.1187852725065882e-01, -2.8578013465508395e-01,
        -2.8578013465508390e-01, -2.3779817858159463e-01,
        -1.9828731035977376e-01, -1.6499523767216026e-01,
        -1.6499523767216037e-01, -1.1448123201310233e-01,
        -4.1187852725065888e-01, -2.8578013465508395e-01,
        -2.8578013465508401e-01,  2.3779817858159458e-01,
        -1.9828731035977384e-01,  1.6499523767216023e-01,
         1.6499523767216037e-01,  1.1448123201310237e-01,
         4.1187852725065888e-01,  2.8578013465508384e-01,
         2.8578013465508390e-01,  2.3779817858159458e-01,
         1.9828731035977371e-01,  1.6499523767216021e-01,
         1.6499523767216018e-01,  1.1448123201310215e-01,
         4.1187852725065888e-01,  2.8578013465508390e-01,
         2.8578013465508401e-01, -2.3779817858159455e-01,
         1.9828731035977379e-01, -1.6499523767216021e-01,
        -1.6499523767216021e-01, -1.1448123201310219e-01,
        // ...
        // idx = [2,1,:]
        -1.8776791509851089e-02,  1.8774475883735472e-02,
        -1.3028195574784233e-02, -1.0840785632730009e-02,
         1.3026588887619593e-02,  1.0839448705368813e-02,
        -7.5218322221567678e-03,  7.5209046008897530e-03,
        -1.8776791509851089e-02,  1.8774475883735469e-02,
        -1.3028195574784237e-02,  1.0840785632730009e-02,
         1.3026588887619593e-02, -1.0839448705368815e-02,
         7.5218322221567661e-03, -7.5209046008897521e-03,
         1.8776791509851089e-02, -1.8774475883735479e-02,
         1.3028195574784240e-02,  1.0840785632730013e-02,
        -1.3026588887619596e-02, -1.0839448705368820e-02,
         7.5218322221567817e-03, -7.5209046008897652e-03,
         1.8776791509851089e-02, -1.8774475883735483e-02,
         1.3028195574784242e-02, -1.0840785632730015e-02,
        -1.3026588887619596e-02,  1.0839448705368822e-02,
        -7.5218322221567791e-03,  7.5209046008897643e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,1,k+1};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-14) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx1[] = {2,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[32+k*basis.num_basis+m], phi_p[m], 1e-14) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[32+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx2[] = {3,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[64+k*basis.num_basis+m], phi_p[m], 1e-14) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[64+k*basis.num_basis+m], idx2[0], idx2[1], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  } if (poly_order == 2) {
    if (!isperiodic) {
      // Solution (checked visually against g2):
      const double sol[240] = {
        // idx = [0,1,:]
        -8.5122112791349098e-04, -1.1408199279301246e-03, -3.0793150106132720e-04, -4.0361855482954289e-04, -4.1269463518755719e-04, -5.4093592784822922e-04,
        -1.4601008289058379e-04, -7.6870970263647495e-04, 3.9096141311460275e-04, 8.2696226169962582e-05, -1.9568500684258495e-04, -2.7808277408889083e-04,
        5.2397262768387270e-04, -3.7603430197459076e-04, 1.9124892215530208e-04, 1.1083078143344461e-04, 2.9915579185682800e-05, -1.3603140624741703e-04,
        2.5631480990695021e-04, 4.0093329184926548e-05, -8.5395288617411710e-04, -1.1444810732657105e-03, -3.0891972185865852e-04, 3.9415546662673637e-04,
        -4.1401906422464335e-04, 5.2825334837797448e-04, 1.4258678563043659e-04, -7.1001010952258183e-04, 3.6110713159018460e-04, 7.6587829001101642e-05,
        1.9109705007607910e-04, -2.5684804056722943e-04, 4.8396145058753455e-04, 4.0992452786021841e-04, -2.0848529962988769e-04, 1.0264421158757013e-04,
        2.7705850306071661e-05, 1.4829128536231821e-04, -2.7941527377411485e-04, 3.7131815826626999e-05, 8.5395288617411884e-04, 1.1444810732657098e-03,
        3.0891972185865814e-04, 3.9415546662673609e-04, 4.1401906422464346e-04, 5.2825334837797405e-04, 1.4258678563043765e-04, 7.1001010952258291e-04,
        -3.6110713159018493e-04, -7.6587829001100328e-05, 1.9109705007607869e-04, 2.5684804056722867e-04, -4.8396145058753423e-04, 4.0992452786021858e-04,
        -2.0848529962988585e-04, -1.0264421158756809e-04, -2.7705850306071424e-05, 1.4829128536231699e-04, -2.7941527377411458e-04, -3.7131815826627277e-05,
        8.5122112791348968e-04, 1.1408199279301255e-03, 3.0793150106132769e-04, -4.0361855482954419e-04, 4.1269463518755649e-04, -5.4093592784823020e-04,
        -1.4601008289058344e-04, 7.6870970263647463e-04, -3.9096141311460134e-04, -8.2696226169961281e-05, -1.9568500684258492e-04, 2.7808277408888969e-04,
        -5.2397262768387270e-04, -3.7603430197459179e-04, 1.9124892215530127e-04, -1.1083078143344442e-04, -2.9915579185681851e-05, -1.3603140624741602e-04,
        2.5631480990694972e-04, -4.0093329184926751e-05,
        // idx = [1,0,:]
        -8.5122112791349153e-04, -3.0793150106132720e-04, -1.1408199279301255e-03, -4.0361855482954289e-04, -4.1269463518755741e-04, -1.4601008289058415e-04,
        -5.4093592784822966e-04, 3.9096141311460253e-04, -7.6870970263647528e-04, 8.2696226169961824e-05, -1.9568500684258636e-04, 5.2397262768387226e-04,
        -2.7808277408888986e-04, 1.9124892215530219e-04, -3.7603430197459168e-04, 2.9915579185682231e-05, 1.1083078143344385e-04, 2.5631480990695086e-04,
        -1.3603140624741619e-04, 4.0093329184925192e-05, -8.5395288617411776e-04, -3.0891972185865824e-04, -1.1444810732657109e-03, 3.9415546662673685e-04,
        -4.1401906422464319e-04, 1.4258678563043716e-04, 5.2825334837797470e-04, 3.6110713159018574e-04, -7.1001010952258345e-04, 7.6587829001101182e-05,
        1.9109705007608045e-04, 4.8396145058753482e-04, -2.5684804056722797e-04, -2.0848529962988704e-04, 4.0992452786021874e-04, 2.7705850306071156e-05,
        1.0264421158756911e-04, -2.7941527377411507e-04, 1.4829128536231762e-04, 3.7131815826625535e-05, 8.5395288617411971e-04, 3.0891972185865786e-04,
        1.1444810732657100e-03, 3.9415546662673626e-04, 4.1401906422464341e-04, 1.4258678563043749e-04, 5.2825334837797459e-04, -3.6110713159018422e-04,
        7.1001010952258291e-04, -7.6587829001100924e-05, 1.9109705007608010e-04, -4.8396145058753303e-04, 2.5684804056722867e-04, -2.0848529962988666e-04,
        4.0992452786021901e-04, -2.7705850306071342e-05, -1.0264421158756869e-04, -2.7941527377411393e-04, 1.4829128536231675e-04, -3.7131815826626545e-05,
        8.5122112791348979e-04, 3.0793150106132747e-04, 1.1408199279301255e-03, -4.0361855482954446e-04, 4.1269463518755681e-04, -1.4601008289058388e-04,
        -5.4093592784823074e-04, -3.9096141311460242e-04, 7.6870970263647539e-04, -8.2696226169961674e-05, -1.9568500684258601e-04, -5.2397262768387172e-04,
        2.7808277408888893e-04, 1.9124892215530097e-04, -3.7603430197459184e-04, -2.9915579185682275e-05, -1.1083078143344484e-04, 2.5631480990694907e-04,
        -1.3603140624741627e-04, -4.0093329184925863e-05,
        // idx = [1,2,:]
        -3.2433083445698618e-02, -1.1732753971887125e-02, 4.3058932372746953e-02, -1.5378605910671510e-02, 1.5576683008490425e-02, -5.5632514830916188e-03,
        2.0417002688116338e-02, 1.4896345637796328e-02, -2.8669545501892997e-02, 3.1508776228205857e-03, 7.3859048827129778e-03, -1.9776742488664067e-02,
        -1.0371284137159195e-02, 7.2869340853499104e-03, -1.4024452265605638e-02, 1.1398383384044293e-03, -4.1831800144129828e-03, -9.6743068697453770e-03,
        -5.0733828098187637e-03, -1.5132764669568534e-03, -3.2537168436913032e-02, -1.1770406993565723e-02, 4.3197118080720860e-02, 1.5018044924493976e-02,
        1.5626672054915903e-02, 5.4328175898799803e-03, -1.9938313353934981e-02, 1.3758842852513904e-02, -2.6480304687123233e-02, 2.9181365070276312e-03,
        -7.2127377462512807e-03, -1.8266566757537843e-02, -9.5793204650059353e-03, -7.9436716246366681e-03, 1.5288411039333878e-02, 1.0556436223728178e-03,
        -3.8741873778643924e-03, 1.0546207234634739e-02, 5.5306232491248520e-03, -1.4014975610190871e-03, 3.2537168436913053e-02, 1.1770406993565723e-02,
        -4.3197118080720881e-02, 1.5018044924493958e-02, -1.5626672054915910e-02, 5.4328175898799742e-03, -1.9938313353934967e-02, -1.3758842852513912e-02,
        2.6480304687123216e-02, -2.9181365070276317e-03, -7.2127377462513102e-03, 1.8266566757537846e-02, 9.5793204650059301e-03, -7.9436716246366473e-03,
        1.5288411039333913e-02, -1.0556436223727863e-03, 3.8741873778644024e-03, 1.0546207234634748e-02, 5.5306232491248806e-03, 1.4014975610190752e-03,
        3.2433083445698618e-02, 1.1732753971887144e-02, -4.3058932372746966e-02, -1.5378605910671526e-02, -1.5576683008490423e-02, -5.5632514830916543e-03,
        2.0417002688116349e-02, -1.4896345637796302e-02, 2.8669545501893032e-02, -3.1508776228206014e-03, 7.3859048827130082e-03, 1.9776742488664040e-02,
        1.0371284137159235e-02, 7.2869340853499148e-03, -1.4024452265605647e-02, -1.1398383384044384e-03, 4.1831800144129993e-03, -9.6743068697453978e-03,
        -5.0733828098187663e-03, 1.5132764669568367e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx1[] = {2,1,k+1};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[80+k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[80+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx2[] = {2,3,k+1};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[160+k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[160+k*basis.num_basis+m], idx2[0], idx2[1], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else {
      // Solution (checked visually against g2):
      const double sol[240] = {
        // idx = [0,0,:]
        -1.4301948349911598e-06, -1.9167695855717431e-06, -1.9167695855717423e-06, -6.6241733513894548e-07, -2.5688847101699803e-06, -8.8778211883114730e-07,
        -8.8778211883114889e-07, -1.2130510689887676e-06, -1.2130510689887695e-06, 1.2649631297541619e-07, -1.1898195423146482e-06, -1.6257500991445802e-06,
        -1.6257500991445825e-06, -7.0035536122142517e-07, -7.0035536122142835e-07, 1.6953234584975017e-07, 1.6953234584974911e-07, -9.3862725737618055e-07,
        -9.3862725737618489e-07, 2.2720991318463379e-07, -1.4301948349911524e-06, -1.9167695855717435e-06, -1.9167695855717427e-06, 6.6241733513894823e-07,
        -2.5688847101699791e-06, 8.8778211883114804e-07, 8.8778211883114931e-07, -1.2130510689887666e-06, -1.2130510689887695e-06, 1.2649631297541566e-07,
        1.1898195423146473e-06, -1.6257500991445804e-06, -1.6257500991445833e-06, 7.0035536122142591e-07, 7.0035536122142835e-07, 1.6953234584975067e-07,
        1.6953234584974940e-07, 9.3862725737618055e-07, 9.3862725737618447e-07, 2.2720991318463252e-07, 1.4301948349911598e-06, 1.9167695855717423e-06,
        1.9167695855717431e-06, 6.6241733513894505e-07, 2.5688847101699820e-06, 8.8778211883114942e-07, 8.8778211883114804e-07, 1.2130510689887683e-06,
        1.2130510689887704e-06, -1.2649631297541619e-07, 1.1898195423146507e-06, 1.6257500991445850e-06, 1.6257500991445838e-06, 7.0035536122142941e-07,
        7.0035536122142994e-07, -1.6953234584974699e-07, -1.6953234584974848e-07, 9.3862725737619072e-07, 9.3862725737618786e-07, -2.2720991318463040e-07,
        1.4301948349911547e-06, 1.9167695855717457e-06, 1.9167695855717448e-06, -6.6241733513894759e-07, 2.5688847101699812e-06, -8.8778211883115005e-07,
        -8.8778211883114836e-07, 1.2130510689887672e-06, 1.2130510689887702e-06, -1.2649631297541658e-07, -1.1898195423146498e-06, 1.6257500991445850e-06,
        1.6257500991445846e-06, -7.0035536122143015e-07, -7.0035536122142983e-07, -1.6953234584974898e-07, -1.6953234584974961e-07, -9.3862725737619093e-07,
        -9.3862725737618743e-07, -2.2720991318462924e-07,
        // idx = [1,1,:]
        -5.1042929176735119e-01, -1.8464915031522100e-01, -1.8464915031522100e-01, -2.3641339134851103e-01, -6.6797319946272804e-02, -8.5523171455331198e-02,
        -8.5523171455331504e-02, 2.2018665896396955e-01, 2.2018665896396972e-01, 4.5145893317134277e-02, -3.0938234141719069e-02, 7.9653107970486733e-02,
        7.9653107970487164e-02, 1.2712482682481197e-01, 1.2712482682481213e-01, 1.6331646666214235e-02, 1.6331646666214214e-02, 4.5987743328550920e-02,
        4.5987743328551274e-02, 5.9080165045449981e-03, -5.1042929176735075e-01, -1.8464915031522106e-01, -1.8464915031522114e-01, 2.3641339134851089e-01,
        -6.6797319946272957e-02, 8.5523171455331171e-02, 8.5523171455331337e-02, 2.2018665896396969e-01, 2.2018665896396958e-01, 4.5145893317134027e-02,
        3.0938234141719135e-02, 7.9653107970487066e-02, 7.9653107970487289e-02, -1.2712482682481188e-01, -1.2712482682481221e-01, 1.6331646666214231e-02,
        1.6331646666214113e-02, -4.5987743328550733e-02, -4.5987743328551191e-02, 5.9080165045451117e-03, 5.1042929176735097e-01, 1.8464915031522072e-01,
        1.8464915031522092e-01, 2.3641339134851092e-01, 6.6797319946272721e-02, 8.5523171455331018e-02, 8.5523171455331462e-02, -2.2018665896396972e-01,
        -2.2018665896396961e-01, -4.5145893317134471e-02, 3.0938234141718597e-02, -7.9653107970486553e-02, -7.9653107970487039e-02, -1.2712482682481249e-01,
        -1.2712482682481208e-01, -1.6331646666214422e-02, -1.6331646666214193e-02, -4.5987743328551031e-02, -4.5987743328550976e-02, -5.9080165045452219e-03,
        5.1042929176735052e-01, 1.8464915031522081e-01, 1.8464915031522106e-01, -2.3641339134851078e-01, 6.6797319946272776e-02, -8.5523171455331004e-02,
        -8.5523171455331282e-02, -2.2018665896396991e-01, -2.2018665896396952e-01, -4.5145893317134263e-02, -3.0938234141718670e-02, -7.9653107970486886e-02,
        -7.9653107970487164e-02, 1.2712482682481238e-01, 1.2712482682481210e-01, -1.6331646666214450e-02, -1.6331646666214103e-02, 4.5987743328550851e-02,
        4.5987743328550892e-02, -5.9080165045452965e-03,
        // idx = [2,1,:]
        -3.2554515935448780e-02, 4.3220149032049821e-02, -1.1776682497178807e-02, -1.5078138422190220e-02, 1.5635003562653463e-02, 2.0018094909631873e-02,
        -5.4545565720819306e-03, -2.7027614890815656e-02, 1.4043218548834502e-02, 2.8793463210621123e-03, 7.2415989356616299e-03, -9.7773113830442026e-03,
        -1.8644110690319519e-02, -1.5604400732765941e-02, 8.1078560094583526e-03, -3.8226886051062874e-03, 1.0416111697008567e-03, -5.6449333589513725e-03,
        -1.0764182325857134e-02, -1.3828677433627428e-03, -3.2554515935448794e-02, 4.3220149032049848e-02, -1.1776682497178824e-02, 1.5078138422190237e-02,
        1.5635003562653481e-02, -2.0018094909631842e-02, 5.4545565720819245e-03, -2.7027614890815649e-02, 1.4043218548834510e-02, 2.8793463210621357e-03,
        -7.2415989356616030e-03, -9.7773113830442199e-03, -1.8644110690319463e-02, 1.5604400732765944e-02, -8.1078560094583457e-03, -3.8226886051062748e-03,
        1.0416111697008567e-03, 5.6449333589513638e-03, 1.0764182325857158e-02, -1.3828677433627308e-03, 3.2554515935448752e-02, -4.3220149032049814e-02,
        1.1776682497178821e-02, 1.5078138422190192e-02, -1.5635003562653467e-02, -2.0018094909631825e-02, 5.4545565720819332e-03, 2.7027614890815625e-02,
        -1.4043218548834505e-02, -2.8793463210621470e-03, -7.2415989356616230e-03, 9.7773113830442251e-03, 1.8644110690319488e-02, 1.5604400732765936e-02,
        -8.1078560094583474e-03, 3.8226886051063516e-03, -1.0416111697008641e-03, 5.6449333589513847e-03, 1.0764182325857099e-02, 1.3828677433627831e-03,
        3.2554515935448794e-02, -4.3220149032049821e-02, 1.1776682497178829e-02, -1.5078138422190215e-02, -1.5635003562653474e-02, 2.0018094909631794e-02,
        -5.4545565720819254e-03, 2.7027614890815618e-02, -1.4043218548834509e-02, -2.8793463210621708e-03, 7.2415989356615969e-03, 9.7773113830442425e-03,
        1.8644110690319435e-02, -1.5604400732765939e-02, 8.1078560094583405e-03, 3.8226886051063334e-03, -1.0416111697008599e-03, -5.6449333589513760e-03,
        -1.0764182325857125e-02, 1.3828677433627668e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,1,k+1};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx1[] = {2,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[80+k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[80+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}

        int idx2[] = {3,2,k+1};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[160+k*basis.num_basis+m], phi_p[m], 1e-12) );
	  TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[160+k*basis.num_basis+m], idx2[0], idx2[1], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  }

  gkyl_fem_parproj_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
  gkyl_array_release(parbuff);

}

void test_1x_p1_nonperiodic() {test_1x(1, false, false);}
void test_1x_p1_periodic() {test_1x(1, true, false);}

void test_1x_p2_nonperiodic() {test_1x(2, false, false);}
void test_1x_p2_periodic() {test_1x(2, true, false);}

void test_3x_p1_nonperiodic() {test_3x(1, false, false);}
void test_3x_p1_periodic() {test_3x(1, true, false);}

void test_3x_p2_nonperiodic() {test_3x(2, false, false);}
void test_3x_p2_periodic() {test_3x(2, true, false);}

#ifdef GKYL_HAVE_CUDA
// ......... GPU tests ............ //
void gpu_test_1x_p1_nonperiodic() {test_1x(1, false, true);}
void gpu_test_1x_p1_periodic() {test_1x(1, true, true);}

void gpu_test_1x_p2_nonperiodic() {test_1x(2, false, true);}
void gpu_test_1x_p2_periodic() {test_1x(2, true, true);}

void gpu_test_3x_p1_nonperiodic() {test_3x(1, false, true);}
void gpu_test_3x_p1_periodic() {test_3x(1, true, true);}

void gpu_test_3x_p2_nonperiodic() {test_3x(2, false, true);}
void gpu_test_3x_p2_periodic() {test_3x(2, true, true);}
#endif


TEST_LIST = {
  { "test_1x_p1_nonperiodic", test_1x_p1_nonperiodic },
  { "test_1x_p1_periodic", test_1x_p1_periodic },
  { "test_1x_p2_nonperiodic", test_1x_p2_nonperiodic },
  { "test_1x_p2_periodic", test_1x_p2_periodic },
  { "test_3x_p1_nonperiodic", test_3x_p1_nonperiodic },
  { "test_3x_p1_periodic", test_3x_p1_periodic },
  { "test_3x_p2_nonperiodic", test_3x_p2_nonperiodic },
  { "test_3x_p2_periodic", test_3x_p2_periodic },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_1x_p1_nonperiodic", gpu_test_1x_p1_nonperiodic },
  { "gpu_test_1x_p1_periodic", gpu_test_1x_p1_periodic },
  { "gpu_test_1x_p2_nonperiodic", gpu_test_1x_p2_nonperiodic },
  { "gpu_test_1x_p2_periodic", gpu_test_1x_p2_periodic },
  { "gpu_test_3x_p1_nonperiodic", gpu_test_3x_p1_nonperiodic },
  { "gpu_test_3x_p1_periodic", gpu_test_3x_p1_periodic },
  { "gpu_test_3x_p2_nonperiodic", gpu_test_3x_p2_nonperiodic },
  { "gpu_test_3x_p2_periodic", gpu_test_3x_p2_periodic },
#endif
  { NULL, NULL },
};

