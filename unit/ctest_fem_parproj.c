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
test_1x(int poly_order, const bool isperiodic)
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

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *parbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  if (isperiodic) apply_periodic_bc(parbuff, rho, dim-1, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_parproj_1x_p2_rho_1.gkyl");

  // parallel FEM projection method.
  gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&grid, &basis, isperiodic, NULL);

  // Set the RHS source.
  gkyl_fem_parproj_set_rhs(parproj, rho);

  // Solve the problem.
  gkyl_fem_parproj_solve(parproj, phi);
  if (isperiodic) {
    apply_periodic_bc(parbuff, phi, dim-1, skin_ghost);
  }
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_parproj_1x_p2_phi_1.gkyl");

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
  gkyl_array_release(parbuff);

}

void
test_3x(const int poly_order, const bool isperiodic)
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

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *parbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  if (isperiodic) apply_periodic_bc(parbuff, rho, dim-1, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_parproj_3x_p1_rho_1.gkyl");

  // parallel FEM projection method.
  gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&grid, &basis, isperiodic, NULL);

  // Set the RHS source.
  gkyl_fem_parproj_set_rhs(parproj, rho);

  // Solve the problem.
  gkyl_fem_parproj_solve(parproj, phi);
  if (isperiodic) {
    apply_periodic_bc(parbuff, phi, dim-1, skin_ghost);
  }
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_parproj_3x_p1_phi_1.gkyl");

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
      const double *phi_p;
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        int idx0[] = {0,1,k};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-14) );

        int idx1[] = {1,0,k};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[32+k*basis.num_basis+m], phi_p[m], 1e-14) );

        int idx2[] = {1,2,k};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[64+k*basis.num_basis+m], phi_p[m], 1e-14) );
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
      const double *phi_p;
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        int idx0[] = {0,0,k};
        linidx= gkyl_range_idx(&localRange, idx0); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-14) );

        int idx1[] = {1,1,k};
        linidx= gkyl_range_idx(&localRange, idx1); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[32+k*basis.num_basis+m], phi_p[m], 1e-14) );

        int idx2[] = {2,1,k};
        linidx= gkyl_range_idx(&localRange, idx2); 
        phi_p = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++)
          TEST_CHECK( gkyl_compare(sol[64+k*basis.num_basis+m], phi_p[m], 1e-14) );
      }
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
  gkyl_array_release(parbuff);

}

void test_1x_p1_nonperiodic() {test_1x(1, false);}
void test_1x_p1_periodic() {test_1x(1, true);}

void test_1x_p2_nonperiodic() {test_1x(2, false);}
void test_1x_p2_periodic() {test_1x(2, true);}

void test_3x_p1_nonperiodic() {test_3x(1, false);}
void test_3x_p1_periodic() {test_3x(1, true);}


TEST_LIST = {
  { "test_1x_p1_nonperiodic", test_1x_p1_nonperiodic },
  { "test_1x_p1_periodic", test_1x_p1_periodic },
  { "test_1x_p2_nonperiodic", test_1x_p2_nonperiodic },
  { "test_1x_p2_periodic", test_1x_p2_periodic },
  { "test_3x_p1_nonperiodic", test_3x_p1_nonperiodic },
  { "test_3x_p1_periodic", test_3x_p1_periodic },
  { NULL, NULL },
};

