#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_iz.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double B0 = 0.5; 

void eval_m0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_m2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*echarge/emass*1.0e19*1.;  //fabs(x);
}
void eval_m2_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*echarge/emass*1.0e19*3.;  //fabs(x);
}
void eval_bmag(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_jac(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void
test_prim_vars_gk_3x(bool use_gpu)
{
  // use vte = 40 eV for elc grid
  double vmax = 4*sqrt(40*echarge/emass);
  double vmin = -vmax;
  double mumax = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;
  double lower[] = {-2.0,-2.0,-2.0,vmin,0.0}, upper[] = {2.0,2.0,2.0,vmax,mumax};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {16, 16, 16, 8, 4};

  // low d test -- use eval_m2_1v
  /* int cdim = 1, vdim = 2; */
  /* int pdim = cdim + vdim; */
  /* double lower[] = {-2.0,vmin,0}, upper[] = {2.,0,vmax,mumax}; */
  /* int ghost[] = {0, 0, 0}; */
  /* int cells[] = {1, 8, 8}; */
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2_3v, NULL);

  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);

  struct gkyl_dg_prim_vars_type *calc_prim_vars_gk_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(&cbasis, &pbasis, "vtSq", use_gpu);
  
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
 
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m2, 2*cbasis.num_basis);

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &confRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(&confRange, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    double *vtSq_elc_d = gkyl_array_fetch(vtSq_elc, loc);

    calc_prim_vars_gk_vtSq->kernel(calc_prim_vars_gk_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);
    
  }

  gkyl_grid_sub_array_write(&confGrid, &confRange, vtSq_elc, "ctest_prim_vars_vtsq_elc.gkyl");
  
  // check vals
  double p1_vals[] = {1.9898778467478656e+13, -8.1549238327948618e-05,
		      3.8383544370155868e-05, 3.8383544370155868e-05,
		      -8.1549238327948822e-05, -2.1582846978896450e-05,
		      -2.1582846978896450e-05,  3.8383544370155820e-05};
  
  const double *pv = gkyl_array_cfetch(vtSq_elc, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  for (int i=0; i<cbasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }

}

void
test_prim_vars_vlasov_3x(bool use_gpu)
{
  // use vte = 40 eV for elc grid
  double vmax = 4*sqrt(40*echarge/emass);
  double vmin = -vmax;
  int poly_order = 1;
  int cdim = 3, vdim = 3;
  int pdim = cdim + vdim;
  double lower[] = {-2.0,-2.0,-2.0,vmin,-vmin,-vmin}, upper[] = {2.0,2.0,2.0,vmax,vmax,vmax};
  int ghost[] = {0, 0, 0, 0, 0, 0};
  int cells[] = {16, 16, 16, 8, 8, 8};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2_3v, NULL);

  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 5*cbasis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_neut = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);

  struct gkyl_dg_prim_vars_type *calc_prim_vars_vlasov_vtSq = gkyl_dg_prim_vars_vlasov_new(&cbasis, &pbasis, "vtSq", use_gpu);
  
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
 
  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2, (vdim+1)*cbasis.num_basis); //HARDCODED for gk vdim = 2

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &confRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(&confRange, conf_iter.idx);

    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    double *vtSq_neut_d = gkyl_array_fetch(vtSq_neut, loc);

    calc_prim_vars_vlasov_vtSq->kernel(calc_prim_vars_vlasov_vtSq, conf_iter.idx, moms_neut_d, vtSq_neut_d);

  }

  gkyl_grid_sub_array_write(&confGrid, &confRange, vtSq_neut, "ctest_prim_vars_vtsq_neut.gkyl");

  double p1_vals[] = {1.9898778467478656e+13, -8.1549238327948618e-05,
		      3.8383544370155868e-05, 3.8383544370155868e-05,
		      -8.1549238327948822e-05, -2.1582846978896450e-05,
		      -2.1582846978896450e-05,  3.8383544370155820e-05};
  
  const double *pv = gkyl_array_cfetch(vtSq_neut, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  for (int i=0; i<cbasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }

}

void
test_coll_iz(bool use_gpu)
{
  // use vte = 40 eV for elc grid
  double vmax = 4*sqrt(40*echarge/emass);
  double vmin = -vmax;
  double mumax = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;
  double lower[] = {-2.0,-2.0,-2.0,vmin,0.0}, upper[] = {2.0,2.0,2.0,vmax,mumax};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {16, 16, 16, 8, 4};

  /* // low d test */
  /* int cdim = 1, vdim = 1; */
  /* int pdim = cdim + vdim; */
  /* double lower[] = {-2.0,vmin}, upper[] = {2.,0,vmax}; */
  /* int ghost[] = {0, 0}; */
  /* int cells[] = {1, 8}; */
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange, &phaseRange);

  // basis functions
  struct gkyl_basis phaseBasis, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis, cdim, vdim);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v, NULL);
  gkyl_proj_on_basis *projBmag = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *projJac = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_jac, NULL);

  // maxwellian on basis for fdist
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&phaseGrid,
    &basis, &phaseBasis, poly_order+1, use_gpu);

  // coll struct 
  struct gkyl_dg_iz *coll_iz = gkyl_dg_iz_new(&phaseGrid, &basis, &phaseBasis, &confRange, &phaseRange,
    echarge, emass, GKYL_IZ_H, true, use_gpu); // might not need bool is_gk ...


  // create moment arrays
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut_ho = gkyl_array_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc_ho = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *distf_elc_ho = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
  struct gkyl_array *coll_iz_elc_ho = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);	
  // arrays necessary for fmax
  struct gkyl_array *bmag_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i_ho = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag_ho);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot_ho);

  gkyl_array_set_offset(moms_neut_ho, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc_ho, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut_ho, 1.0, m2, 4*basis.num_basis);
  gkyl_array_set_offset(moms_elc_ho, 1.0, m2, 2*basis.num_basis);

  gkyl_array_shiftc(bmag_ho, 0.5*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(jacob_tot_ho, 1.0*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(b_z, 1.0*pow(sqrt(2),cdim), 0);

  // project b_i
  gkyl_array_set_offset(b_i_ho, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i_ho, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i_ho, 1.0, b_z, 2*basis.num_basis);

  struct gkyl_array *moms_neut, *moms_elc, *cflRate, *distf_elc, *coll_iz_elc, *bmag, *jacob_tot, *b_i;
  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_neut = gkyl_array_cu_dev_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
    struct gkyl_array *distf_elc = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
    struct gkyl_array *coll_iz_elc = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);	
    // arrays necessary for fmax
    struct gkyl_array *bmag = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_neut, moms_neut_ho);
    gkyl_array_copy(moms_elc, moms_elc_ho);
    gkyl_array_copy(cflRate, cflRate_ho);
    gkyl_array_copy(distf_elc, distf_elc_ho);
    gkyl_array_copy(coll_iz_elc, coll_iz_elc_ho);
    gkyl_array_copy(bmag, bmag_ho);
    gkyl_array_copy(jacob_tot, jacob_tot_ho);
    gkyl_array_copy(b_i, b_i_ho);
  }
  else {
    moms_neut = moms_neut_ho;
    moms_elc = moms_elc_ho;
    cflRate = cflRate_ho;
    distf_elc = distf_elc_ho;
    coll_iz_elc = coll_iz_elc_ho;
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
    b_i = b_i_ho;
  }
  
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max, &phaseRange, &confRange, moms_elc,
    bmag, jacob_tot, emass, distf_elc);
  // gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, distf_elc, "ctest_distf_elc.gkyl");
  
  struct timespec tm;
  double tm_tot = 0.0;
  int iter = 1;
  for (int t=0; t<iter; ++t) {
    tm = gkyl_wall_clock();
    gkyl_dg_iz_coll(coll_iz, moms_elc, moms_neut, bmag, jacob_tot, b_i, distf_elc, coll_iz_elc, cflRate);
    tm_tot = tm_tot + gkyl_time_diff_now_sec(tm);
  }
  tm_tot = tm_tot/iter;
  printf("Avg time over %d loop(s) is %.e s", iter, tm_tot);									    
  
  //gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, coll_iz_elc, "ctest_coll_iz_elc.gkyl");

  double p1_vals[] = {-2.4014565247326178e+00, -3.8514489897962783e-16, -5.7110482341097731e-17,
		      2.4425487837962348e-16,  3.3375014364519018e-01,  2.0687952656971738e+00,
		      -2.3506351153000692e-17, -4.0308416747049199e-17, -4.0308416747049205e-17,
		      4.0625859545969633e-17, -2.1597670852924058e-17, -6.5294602977688739e-18,
		      -8.9730204140520383e-17, -4.4082448908800082e-17, -1.3946142760251460e-17,
		      -2.8751741937314046e-01,  3.1625898030465171e-18,  2.0621968577456047e-17,
		      -4.8782239729356911e-19,  4.4478436458446171e-18,  9.1089012291912308e-17,
		      2.3503224210675240e-17,  2.3503224210675240e-17, -7.0074177551781357e-18,
		      -5.3010544555376448e-18,  9.7672710613792879e-18, -1.5936792928689385e-18,
		      -1.3945912836727961e-17, -1.5710239442394578e-17, -1.0505618208525673e-17,
		      -1.0505618208525673e-17,  9.7672710613792879e-18,  6.7454560758484550e-02,
		      -6.5759644472814987e-18, -3.3502754511312654e-18, -3.3502754511312651e-18,
		      -5.8110429485091361e-02,  1.3427984002112964e-17,  1.2206926022355283e-17,
		      2.8064629419848202e-18,  3.2284449575695321e-17,  4.2517357754476330e-18,
		      4.2517357754476330e-18, -1.1616526826346376e-17, -8.9436596275453193e-18,
		      -7.8917484916331066e-18, -7.8917484916331066e-18,  6.9430384740724518e-19};
  
  const double *pv = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange, (int[5]) { 2, 2, 2, 5, 3}));

  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }
  
   gkyl_dg_iz_release(coll_iz);
}

void prim_vars_gk_3x() { test_prim_vars_gk_3x(false); }
void prim_vars_vlasov_3x() { test_prim_vars_vlasov_3x(false); }
void coll_iz() { test_coll_iz(false); }

#ifdef GKYL_HAVE_CUDA
void prim_vars_gk_3x_gpu() { test_prim_vars_gk_3x(true); }
void prim_vars_vlasov_3x_gpu() { test_prim_vars_vlasov_3x(true); }
void coll_iz_gpu() { test_coll_iz(true); }
#endif

TEST_LIST = {
  { "prim_vars_gk_3x", prim_vars_gk_3x },
  { "prim_vars_vlasov_3x", prim_vars_vlasov_3x },
  { "coll_iz", coll_iz },
#ifdef GKYL_HAVE_CUDA
  { "prim_vars_gk_3x_gpu", prim_vars_3x_gpu },
  { "coll_iz_gpu", coll_iz_gpu },
#endif  
  { NULL, NULL },
};
