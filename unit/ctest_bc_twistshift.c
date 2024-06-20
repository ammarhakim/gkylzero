// Test creation and deallocation of updater that applies the
// twist shift BCs.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>
#include <gkyl_bc_twistshift.h>

void shift_2x_linear(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.6*x+1.8;
}

void
test_bc_twistshift_2x(bool use_gpu)
{
  assert(!use_gpu); // 2x test only available on CPUs.

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.5}, upper[] = {2.0, 1.5};
  const int cells[] = {12, 10};
  const int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Parameters expected by the updater, but not used in this 2x test.
  int bc_dir = 0;
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE;
  int cdim = ndim;

  struct gkyl_bc_twistshift_inp tsinp = {
    .dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .local_r = local,
    .local_ext_r = local_ext,
    .basis = basis,
    .grid = grid,
    .shift_func = shift_2x_linear,
    .shift_func_ctx = NULL,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  gkyl_bc_twistshift_release(tsup);

}

void test_bc_twistshift_2x_ho(){ test_bc_twistshift_2x(false); }

TEST_LIST = {
  { "test_bc_twistshift_2x_ho", test_bc_twistshift_2x_ho },
  { NULL, NULL },
};
