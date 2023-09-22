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

  //gkyl_grid_sub_array_write(&confGrid, &confRange, vtSq_elc, "ctest_prim_vars_vtsq_elc.gkyl");
  
  double p1_vals[] = {1.9898778467478656e+13, -8.1549238327948618e-05,
		      3.8383544370155868e-05, 3.8383544370155868e-05,
		      -8.1549238327948822e-05, -2.1582846978896450e-05,
		      -2.1582846978896450e-05,  3.8383544370155820e-05};
  
  const double *pv = gkyl_array_cfetch(vtSq_elc, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  for (int i=0; i<cbasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2);
  gkyl_array_release(moms_elc); gkyl_array_release(vtSq_elc);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2);
  gkyl_dg_prim_vars_type_release(calc_prim_vars_gk_vtSq);

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
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, (2+vdim)*cbasis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);

  struct gkyl_dg_prim_vars_type *calc_prim_vars_vlasov_vtSq = gkyl_dg_prim_vars_vlasov_new(&cbasis, &pbasis, "vtSq", use_gpu);
  
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
 
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m2, (vdim+1)*cbasis.num_basis);

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &confRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(&confRange, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    double *vtSq_elc_d = gkyl_array_fetch(vtSq_elc, loc);

    calc_prim_vars_vlasov_vtSq->kernel(calc_prim_vars_vlasov_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);

  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, vtSq_neut, "ctest_prim_vars_vtsq_neut.gkyl");

  double p1_vals[] = {1.9898778467478656e+13, -8.1549238327948618e-05,
		      3.8383544370155868e-05, 3.8383544370155868e-05,
		      -8.1549238327948822e-05, -2.1582846978896450e-05,
		      -2.1582846978896450e-05,  3.8383544370155820e-05};
  
  const double *pv = gkyl_array_cfetch(vtSq_elc, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  for (int i=0; i<cbasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }

  gkyl_array_release(m0); gkyl_array_release(m2);
  gkyl_array_release(moms_elc); gkyl_array_release(vtSq_elc);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2);
  gkyl_dg_prim_vars_type_release(calc_prim_vars_vlasov_vtSq);

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
  int charge_state = 1;
  bool all_gk = false;
  
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

  // coll struct. FAILING HERE!
  struct gkyl_dg_iz *coll_iz = gkyl_dg_iz_new(&phaseGrid, &basis, &phaseBasis, &confRange, &phaseRange,
					      echarge, emass, GKYL_IZ_H, charge_state, GKYL_IZ_ELC, all_gk, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
  struct gkyl_array *coll_iz_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);	
  // arrays necessary for fmax
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot);

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2, 4*basis.num_basis); //HARDCODED for gk vdim = 2
  gkyl_array_set_offset(moms_elc, 1.0, m2, 2*basis.num_basis);

  gkyl_array_shiftc(bmag, 0.5*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(jacob_tot, 1.0*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(b_z, 1.0*pow(sqrt(2),cdim), 0);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
    struct gkyl_array *coll_iz_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_cu, cflRate);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_iz_elc_cu, coll_iz_elc);
    gkyl_array_copy(bmag_cu, bmag);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max, &phaseRange, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);

    struct timespec tm;
    double tm_tot = 0.0;
    int iter = 1;
    for (int t=0; t<iter; ++t) {
      tm = gkyl_wall_clock();
      gkyl_dg_iz_coll(coll_iz, moms_elc_cu, moms_neut_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_cu);
      tm_tot = tm_tot + gkyl_time_diff_now_sec(tm);
    }
    tm_tot = tm_tot/iter;
    printf("Avg time over %d loop(s) is %.e s", iter, tm_tot);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_neut_cu);
    gkyl_array_release(cflRate_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_iz_elc_cu); gkyl_array_release(b_i_cu);
  }
  else {
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
									     

  }

  gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, coll_iz_elc, "ctest_coll_iz_elc.gkyl");

  double p1_vals[] = {-2.6709857935892035e+02,  7.7615561179731124e-15, -3.0274206010288621e-14,
      3.1077133685936207e-15, -3.7063309667532565e+01,  2.3005854246601783e+02,
      -1.1324477072409585e-14, -1.7133083509664025e-15, -1.7133083509664027e-15,
      -2.5320494026474838e-16, -4.3157812532942333e-15,  4.4873753019363520e-16,
      3.1318621454795794e-14,  3.2523309778925537e-14, 1.3441936264614712e-14,
      3.1915684068822404e+01,  3.1263520728811039e-15, -8.8144842070210272e-16,
      -2.5986148369981679e-15,  4.1538279223702380e-16,  2.6802523170780072e-15,
      3.2872561552147496e-15,  3.2872561552147488e-15,  1.9610073686573559e-15,
      -5.4610154564928583e-15, -1.0234853592124468e-14,  1.0863003627747311e-15,
      -5.6440969257679846e-15,  2.8470585013498972e-15,  2.2693229075893391e-15,
      2.2693229075893399e-15, -6.9183699693311998e-16,  7.5000073131241232e+00,
      -4.3550153883436369e-16, -2.1102998399444678e-15, -9.1742276554554930e-16,
      -6.4595947989844511e+00,  1.4648335092237129e-15,  1.3619180045528966e-16,
      1.1140172408101608e-15, -2.5308010045766849e-15, -1.1017317798945724e-15,
      9.1145294504346216e-17,  4.2716608578520873e-19, -1.7617649270258086e-15,
      -4.5195669650326138e-16, -4.5195669650326138e-16,  8.5785153401928629e-16};
  
  const double *pv = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange, (int[5]) { 1, 1, 1, 4, 2}));

  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_neut);
  gkyl_array_release(cflRate); gkyl_array_release(distf_elc);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot);
  gkyl_array_release(coll_iz_elc); gkyl_array_release(b_i);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_dg_iz_release(coll_iz);
}

// tests when donor species is gk
void
test_coll_iz_all_gk(bool use_gpu)
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
  int charge_state = 1;
  bool all_gk = true;
  
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

  // coll struct. FAILING HERE!
  struct gkyl_dg_iz *coll_iz = gkyl_dg_iz_new(&phaseGrid, &basis, &phaseBasis, &confRange, &phaseRange,
					      echarge, emass, GKYL_IZ_H, charge_state, GKYL_IZ_ELC, all_gk, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
  struct gkyl_array *coll_iz_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);	
  // arrays necessary for fmax
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot);

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2, 2*basis.num_basis); //HARDCODED for gk vdim = 2
  gkyl_array_set_offset(moms_elc, 1.0, m2, 2*basis.num_basis);

  gkyl_array_shiftc(bmag, 0.5*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(jacob_tot, 1.0*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(b_z, 1.0*pow(sqrt(2),cdim), 0);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
    struct gkyl_array *coll_iz_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange.volume);
    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_cu, cflRate);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_iz_elc_cu, coll_iz_elc);
    gkyl_array_copy(bmag_cu, bmag);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max, &phaseRange, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);

    struct timespec tm;
    double tm_tot = 0.0;
    int iter = 1;
    for (int t=0; t<iter; ++t) {
      tm = gkyl_wall_clock();
      gkyl_dg_iz_coll(coll_iz, moms_elc_cu, moms_neut_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_cu);
      tm_tot = tm_tot + gkyl_time_diff_now_sec(tm);
    }
    tm_tot = tm_tot/iter;
    printf("Avg time over %d loop(s) is %.e s", iter, tm_tot);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_neut_cu);
    gkyl_array_release(cflRate_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_iz_elc_cu); gkyl_array_release(b_i_cu);
  }
  else {
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
									     

  }

  gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, coll_iz_elc, "ctest_coll_iz_elc.gkyl");

  double p1_vals[] = {-2.6709857935892035e+02,  7.7615561179731124e-15, -3.0274206010288621e-14,
      3.1077133685936207e-15, -3.7063309667532565e+01,  2.3005854246601783e+02,
      -1.1324477072409585e-14, -1.7133083509664025e-15, -1.7133083509664027e-15,
      -2.5320494026474838e-16, -4.3157812532942333e-15,  4.4873753019363520e-16,
      3.1318621454795794e-14,  3.2523309778925537e-14, 1.3441936264614712e-14,
      3.1915684068822404e+01,  3.1263520728811039e-15, -8.8144842070210272e-16,
      -2.5986148369981679e-15,  4.1538279223702380e-16,  2.6802523170780072e-15,
      3.2872561552147496e-15,  3.2872561552147488e-15,  1.9610073686573559e-15,
      -5.4610154564928583e-15, -1.0234853592124468e-14,  1.0863003627747311e-15,
      -5.6440969257679846e-15,  2.8470585013498972e-15,  2.2693229075893391e-15,
      2.2693229075893399e-15, -6.9183699693311998e-16,  7.5000073131241232e+00,
      -4.3550153883436369e-16, -2.1102998399444678e-15, -9.1742276554554930e-16,
      -6.4595947989844511e+00,  1.4648335092237129e-15,  1.3619180045528966e-16,
      1.1140172408101608e-15, -2.5308010045766849e-15, -1.1017317798945724e-15,
      9.1145294504346216e-17,  4.2716608578520873e-19, -1.7617649270258086e-15,
      -4.5195669650326138e-16, -4.5195669650326138e-16,  8.5785153401928629e-16};
  
  const double *pv = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange, (int[5]) { 1, 1, 1, 4, 2}));

  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_neut);
  gkyl_array_release(cflRate); gkyl_array_release(distf_elc);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot);
  gkyl_array_release(coll_iz_elc); gkyl_array_release(b_i);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_dg_iz_release(coll_iz);
}

void prim_vars_gk_3x() { test_prim_vars_gk_3x(false); }
void prim_vars_vlasov_3x() { test_prim_vars_vlasov_3x(false); }
void coll_iz() { test_coll_iz(false); }
void coll_iz_all_gk() { test_coll_iz_all_gk(false); }

#ifdef GKYL_HAVE_CUDA
void prim_vars_gk_3x_gpu() { test_prim_vars_gk_3x(true); }
void prim_vars_vlasov_3x_gpu() { test_prim_vars_vlasov_3x(true); }
void coll_iz_gpu() { test_coll_iz(true); }
#endif

TEST_LIST = {
  { "prim_vars_gk_3x", prim_vars_gk_3x },
  { "prim_vars_vlasov_3x", prim_vars_vlasov_3x },
  { "coll_iz", coll_iz },
  { "coll_iz_all_gk", coll_iz_all_gk },
#ifdef GKYL_HAVE_CUDA
  { "coll_iz_gpu", coll_iz_gpu },
#endif  
  { NULL, NULL },
};
