// Test creation and deallocation of updater that applies the
// conducting sheath BC for gyrokinetics.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>

void
test_bc_sheath_gyrokinetic(bool use_gpu)
{

  double charge = -1., mass = 1.;

  int poly_order = 1;
  double lower[] = {-2.0, -2.0, 0.}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {6, 8, 4};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = 1;
  int dir = 0;
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  int vdim = ndim-cdim;
  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  // Basis functions.
  struct gkyl_basis basis;
  if (poly_order == 1)  // Force gkhybrid for p=1.
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double q2Dm = 2.*charge/mass;

  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    &local_ext, ghost, &basis, &grid, cdim, q2Dm, use_gpu);

  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
}

void test_bc_gksheath(){ test_bc_sheath_gyrokinetic(false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_gksheath_cu(){ test_bc_sheath_gyrokinetic(true); }
#endif

TEST_LIST = {
  { "test_bc_gksheath", test_bc_gksheath },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_gksheath_cu", test_bc_gksheath_cu },
#endif
  { NULL, NULL },
};
