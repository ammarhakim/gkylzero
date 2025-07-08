// Test updater which modifies the potential at the boundary to account
// for the rarefaction wave that should speed ions up to the sound speed.
//
#include <acutest.h>
#include <math.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_sheath_rarefaction_pot.h>

// allocate array (filled with zeros)
static struct gkyl_array *mkarr(long nc, long size) {
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  gkyl_array_clear(a, 0.);
  return a;
}

// allocate cu_dev array
static struct gkyl_array *mkarr_cu(long nc, long size) {
  struct gkyl_array *a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_sheath_rarefaction_pot(bool use_gpu)
{

  double eV = 1.;
  double qe = -eV, me = 1./(2.014*1836.);
  double qi =  eV, mi = 1.;
  double n0 = 1.;
  double upari = 0.6, upare = 0.6;
  double Tpari = 0.5/3., Tpare = 0.5;
  double phi_p = 1.;

  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {16};
  int ndim = sizeof(lower)/sizeof(lower[0]);

  double cs = sqrt((Tpare+3.*Tpari)/mi);

  enum gkyl_edge_loc edge_lo = GKYL_LOWER_EDGE;
  enum gkyl_edge_loc edge_up = GKYL_UPPER_EDGE;

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create arrays for moments.
  struct gkyl_array *m0e, *m1e, *m2pare, *m0i, *m1i, *m2pari;
  if (use_gpu) {
    m0e = mkarr_cu(basis.num_basis, local_ext.volume);
    m1e = mkarr_cu(basis.num_basis, local_ext.volume);
    m2pare = mkarr_cu(basis.num_basis, local_ext.volume);
    m0i = mkarr_cu(basis.num_basis, local_ext.volume);
    m1i = mkarr_cu(basis.num_basis, local_ext.volume);
    m2pari = mkarr_cu(basis.num_basis, local_ext.volume);

    // Assign moment arrays by setting the zeroth coefficient only.
    // The factor multiplying the value depends on basis dimensionality.
    gkyl_array_shiftc0_cu(m0e, sqrt(2.)*n0);
    gkyl_array_shiftc0_cu(m1e, sqrt(2.)*(n0*upare));
    gkyl_array_shiftc0_cu(m2pare, sqrt(2.)*(n0*pow(upare,2)+n0*Tpare/me));
    gkyl_array_shiftc0_cu(m0i, sqrt(2.)*n0);
    gkyl_array_shiftc0_cu(m1i, sqrt(2.)*(n0*upari));
    gkyl_array_shiftc0_cu(m2pari, sqrt(2.)*(n0*pow(upari,2)+n0*Tpari/mi));
  } else {
    m0e = mkarr(basis.num_basis, local_ext.volume);
    m1e = mkarr(basis.num_basis, local_ext.volume);
    m2pare = mkarr(basis.num_basis, local_ext.volume);
    m0i = mkarr(basis.num_basis, local_ext.volume);
    m1i = mkarr(basis.num_basis, local_ext.volume);
    m2pari = mkarr(basis.num_basis, local_ext.volume);

    // Assign moment arrays by setting the zeroth coefficient only.
    // The factor multiplying the value depends on basis dimensionality.
    gkyl_array_shiftc0(m0e, sqrt(2.)*n0);
    gkyl_array_shiftc0(m1e, sqrt(2.)*(n0*upare));
    gkyl_array_shiftc0(m2pare, sqrt(2.)*(n0*pow(upare,2)+n0*Tpare/me));
    gkyl_array_shiftc0(m0i, sqrt(2.)*n0);
    gkyl_array_shiftc0(m1i, sqrt(2.)*(n0*upari));
    gkyl_array_shiftc0(m2pari, sqrt(2.)*(n0*pow(upari,2)+n0*Tpari/mi));
  };

  // Arrays holding m0 and m1 together.
  struct gkyl_array *momse, *momsi;
  if (use_gpu) {
    momse = mkarr_cu(2*basis.num_basis, local_ext.volume);
    momsi = mkarr_cu(2*basis.num_basis, local_ext.volume);
    gkyl_array_set_offset_cu(momse, 1., m0e, 0*basis.num_basis);
    gkyl_array_set_offset_cu(momse, 1., m1e, 1*basis.num_basis);
    gkyl_array_set_offset_cu(momsi, 1., m0i, 0*basis.num_basis);
    gkyl_array_set_offset_cu(momsi, 1., m1i, 1*basis.num_basis);
  } else {
    momse = mkarr(2*basis.num_basis, local_ext.volume);
    momsi = mkarr(2*basis.num_basis, local_ext.volume);
    gkyl_array_set_offset(momse, 1., m0e, 0*basis.num_basis);
    gkyl_array_set_offset(momse, 1., m1e, 1*basis.num_basis);
    gkyl_array_set_offset(momsi, 1., m0i, 0*basis.num_basis);
    gkyl_array_set_offset(momsi, 1., m1i, 1*basis.num_basis);
  };

  // Define and assign the potential array.
  struct gkyl_array *phi, *phi_wall, *phi_host;
  if (use_gpu) {
    phi = mkarr_cu(basis.num_basis, local_ext.volume);
    phi_wall = mkarr_cu(basis.num_basis, local_ext.volume);
    gkyl_array_shiftc0_cu(phi, sqrtf(2.)*phi_p);
  } else {
    phi = mkarr(basis.num_basis, local_ext.volume);
    phi_wall = mkarr(basis.num_basis, local_ext.volume);
    gkyl_array_shiftc0(phi, sqrtf(2.)*phi_p);
  };
  phi_host = mkarr(basis.num_basis, local_ext.volume);

  // Create and use updaters for each boundary.
  struct gkyl_sheath_rarefaction_pot *sheath_pot_lo = gkyl_sheath_rarefaction_pot_new(edge_lo,
    &local_ext, ghost, &basis, false, &grid, eV, me, mi, -1.0, use_gpu);
  struct gkyl_sheath_rarefaction_pot *sheath_pot_up = gkyl_sheath_rarefaction_pot_new(edge_up,
    &local_ext, ghost, &basis, false, &grid, eV, me, mi, -1.0, use_gpu);

  gkyl_sheath_rarefaction_pot_advance(sheath_pot_lo, momse, m2pare, momsi, m2pari, phi_wall, phi);
  gkyl_sheath_rarefaction_pot_advance(sheath_pot_up, momse, m2pare, momsi, m2pari, phi_wall, phi);

  gkyl_array_copy(phi_host, phi);

  /* Check that the value of the potential at the boundary is phi_p-(Tpare/e)*(1-|upari|/c_s). */
  int idx[] = {1};
  double zlog_lo[] = {-1.};
  long linidx = gkyl_range_idx(&local_ext, idx);
  const double *phi_skin = (const double *)gkyl_array_cfetch(phi_host, linidx);
  double phi_b = basis.eval_expand(zlog_lo, phi_skin);
  TEST_CHECK(gkyl_compare(phi_b, phi_p-(Tpare/eV)*(1.-fabs(upari)/cs), 1e-6));

  idx[0] = cells[0];
  zlog_lo[0] = 1.;
  linidx = gkyl_range_idx(&local_ext, idx);
  phi_skin = (const double *)gkyl_array_cfetch(phi_host, linidx);
  phi_b = basis.eval_expand(zlog_lo, phi_skin);
  TEST_CHECK(gkyl_compare(phi_b, phi_p-(Tpare/eV)*(1.-fabs(upari)/cs), 1e-6));

  gkyl_sheath_rarefaction_pot_release(sheath_pot_lo);
  gkyl_sheath_rarefaction_pot_release(sheath_pot_up);
  gkyl_array_release(phi);
  gkyl_array_release(phi_wall);
  gkyl_array_release(m0e);
  gkyl_array_release(m1e);
  gkyl_array_release(m2pare);
  gkyl_array_release(momse);
  gkyl_array_release(m0i);
  gkyl_array_release(m1i);
  gkyl_array_release(m2pari);
  gkyl_array_release(momsi);
}

void test_sheath_rare_pot(){ test_sheath_rarefaction_pot(false); }

#ifdef GKYL_HAVE_CUDA
void test_sheath_rare_pot_cu(){ test_sheath_rarefaction_pot(true); }
#endif

TEST_LIST = {
  { "test_sheath_rare_pot", test_sheath_rare_pot },
#ifdef GKYL_HAVE_CUDA
  { "test_sheath_rare_pot_cu", test_sheath_rare_pot_cu },
#endif
  { NULL, NULL },
};
