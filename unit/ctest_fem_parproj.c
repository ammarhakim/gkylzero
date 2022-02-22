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

// Apply periodic BCs along parallerl direction
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
  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_parproj_1x1v_p2_rho_1.gkyl");

  // parallel FEM projection method.
  gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&grid, &basis, isperiodic, NULL);

  // Set the RHS source.
  gkyl_fem_parproj_set_rhs(parproj, rho);

  // Solve the problem.
  gkyl_fem_parproj_solve(parproj, phi);
  if (isperiodic) {
    apply_periodic_bc(parbuff, phi, dim-1, skin_ghost);
  }
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_parproj_1x1v_p2_phi_1.gkyl");

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

void test_1x_p1_nonperiodic() {
  test_1x(1, false);
}

void test_1x_p1_periodic() {
  test_1x(1, true);
}

void test_1x_p2_nonperiodic() {
  test_1x(2, false);
}

void test_1x_p2_periodic() {
  test_1x(2, true);
}

TEST_LIST = {
  { "test_1x_p1_nonperiodic", test_1x_p1_nonperiodic },
  { "test_1x_p1_periodic", test_1x_p1_periodic },
  { "test_1x_p2_nonperiodic", test_1x_p2_nonperiodic },
  { "test_1x_p2_periodic", test_1x_p2_periodic },
  { NULL, NULL },
};

