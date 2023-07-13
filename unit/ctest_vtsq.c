#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_vtsq.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void moms_vlasov_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double m0 = 1.5, m1 = 2., m2 = 9.;
  fout[0] = m0;
  fout[1] = m1;
  fout[2] = m2;
}
void moms_vlasov_1x1v_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double m0 = 1.5, m1 = 2., m2 = 9.;
  fout[0] = (m0*m2 - m1*m1)/(1.*m0*m0);
}

void
test_vtsq_vlasov_1x1v(int poly_order, bool use_gpu)
{
  int m1comps = 1, vdim_phys = 1;
  double Lx = 1.;
  double lower[] = {-Lx/2.}, upper[] = {Lx/2.};
  int cells[] = {32};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create moments and v_t^2 arrays.
  struct gkyl_array *moms_ho, *vtsq_ho,  *moms, *vtsq;
  moms_ho = mkarr((m1comps+2)*basis.num_basis, local_ext.volume);
  vtsq_ho = mkarr(basis.num_basis, local_ext.volume);
  if (use_gpu) { // Create device copies
    moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, moms_ho->ncomp, moms_ho->size);
    vtsq = gkyl_array_cu_dev_new(GKYL_DOUBLE, vtsq_ho->ncomp, vtsq_ho->size);
  } else {
    moms = moms_ho;
    vtsq = vtsq_ho;
  }

  gkyl_proj_on_basis *proj_moms = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, m1comps+2, moms_vlasov_1x1v, NULL);

  gkyl_proj_on_basis_advance(proj_moms, 0.0, &local, moms);

  // Copy host array to device.
  gkyl_array_copy(moms, moms_ho);

  // Create and call updater that computes v_t^2.
  struct gkyl_vtsq *vtsq_up = gkyl_vtsq_new(&basis, &local, m1comps, vdim_phys, use_gpu); 

  gkyl_vtsq_advance(vtsq_up, basis, moms, &local, vtsq);
  gkyl_array_copy(vtsq_ho, vtsq);

  // Project expected vtsq and compare.
  struct gkyl_array *vtsqA;
  vtsqA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_sol = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, moms_vlasov_1x1v_sol, NULL);
  gkyl_proj_on_basis_advance(proj_sol, 0.0, &local, vtsqA);

  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&local, idx);
    const double *vtsq_p = gkyl_array_cfetch(vtsq_ho, linidx);
    const double *vtsqA_p = gkyl_array_cfetch(vtsqA, linidx);
    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(vtsqA_p[m], vtsq_p[m], 1e-10) );
      TEST_MSG("Expected: %.13e in cell (%d)", vtsqA_p[m], idx[0]);
      TEST_MSG("Produced: %.13e", vtsq_p[m]);
    }
  }

  gkyl_array_release(moms_ho); gkyl_array_release(vtsq_ho);
  if (use_gpu) {
    gkyl_array_release(moms); gkyl_array_release(vtsq);
  }
  gkyl_proj_on_basis_release(proj_moms);
  gkyl_vtsq_release(vtsq_up);
  gkyl_array_release(vtsqA);
  gkyl_proj_on_basis_release(proj_sol);
}

void vtsq_vlasov_1x1v_p1() { test_vtsq_vlasov_1x1v(1, false); }

#ifdef GKYL_HAVE_CUDA
void vtsq_vlasov_1x1v_p1_gpu() { test_vtsq_vlasov_1x1v(1, true); }
#endif

TEST_LIST = {
  { "vtsq_vlasov_1x1v_p1", vtsq_vlasov_1x1v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "vtsq_vlasov_1x1v_p1_gpu", vtsq_vlasov_1x1v_p1_gpu },
#endif
  { NULL, NULL },
};
