#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_gyrokinetic.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_3x2v_p1(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 3, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0};
  double lower[] = {0., 0., 0., -1., 0.};
  double upper[] = {1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_dg_updater_gyrokinetic* up;
  up = gkyl_dg_updater_gyrokinetic_new(&phaseGrid, &confBasis, &basis, &confRange, &phaseRange, 0, 1.0, 1.0, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *bmag, *jacobtot_inv, *cmag, *b_i, *phi, *apar, *apardot;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  // Initialize gyrokinetic variables
  bmag = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  jacobtot_inv = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  cmag = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  b_i = mkarr1(use_gpu, 3*confBasis.num_basis, confRange_ext.volume);
  phi = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  apar = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  apardot = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);

  struct timespec tm = gkyl_wall_clock();
  // run hyper_dg_advance
  int nrep = 1;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_updater_gyrokinetic_advance(up, &phaseRange, bmag, jacobtot_inv, cmag, b_i, phi, apar, apardot, fin, cflrate, rhs);
  }
  double gk_tm = gkyl_time_diff_now_sec(tm);

  printf("\ngyrokinetic update on (%d, %d, %d, %d, %d) took %g sec\n", cells[0], cells[1], cells[2], cells[3], cells[4], gk_tm); 

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(bmag);
  gkyl_array_release(jacobtot_inv);
  gkyl_array_release(cmag);
  gkyl_array_release(b_i);
  gkyl_array_release(phi);
  gkyl_array_release(apar);
  gkyl_array_release(apardot);

  gkyl_dg_updater_gyrokinetic_release(up);
}

void
test_gyrokinetic_3x2v_p1()
{
  test_3x2v_p1(false);
}

void
test_gyrokinetic_3x2v_p1_cu()
{
  test_3x2v_p1(true);
}

TEST_LIST = {
  { "test_gyrokinetic_3x2v_p1", test_gyrokinetic_3x2v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "test_gyrokinetic_3x2v_p1_cu", test_gyrokinetic_3x2v_p1_cu },
#endif
  { NULL, NULL },
};

