#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_recomb.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double h_ion_mass = GKYL_PROTON_MASS*2.01410177811;
double li_ion_mass = GKYL_PROTON_MASS*6.94;
double ar_ion_mass = GKYL_PROTON_MASS*39.948;
double B0 = 0.5;
double check_fac = 1.e10;

void eval_m0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_m2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/emass*1.0e19*1.;  //fabs(x);
}
void eval_m2_3v_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/emass*1.0e19*3.;  //fabs(x);
}
void eval_m2_3v_h_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/h_ion_mass*1.0e19*3.;  //fabs(x);
}
void eval_m2_3v_li_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/li_ion_mass*1.0e19*3.;  //fabs(x);
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
test_coll_recomb_h(bool use_gpu)
{
  int charge_state = 0;
  bool all_gk = false;
  // use vt = 4 eV for all grids
  double vmax_elc = 4.*sqrt(4.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*4.*echarge/(2.*B0);
  double vmax_ion = 4.*sqrt(4.*echarge/h_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*4.*echarge/(2.*B0);
  int poly_order = 1;
  const int cdim = 3, vdim_gk = 2, vdim_vl = 3;
  int pdim_gk = cdim + vdim_gk, pdim_vl = cdim + vdim_vl;
  char basepath[4000] = ".";

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};

  // for vlasov grid
  double lower_vl[] = {-2.0,-2.0,-2.0,vmin_ion,vmin_ion,vmin_ion}, upper_vl[] = {2.0,2.0,2.0,vmax_ion,vmax_ion,vmax_ion};
  int ghost_vl[] = {0, 0, 0, 0, 0, 0};
  int cells_vl[] = {16, 16, 16, 4, 4, 4};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange_ext, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim_gk, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost_gk, &phaseRange_ext_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim_gk, lower_ion, upper_ion, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost_gk, &phaseRange_ext_ion, &phaseRange_ion);

  // vlasov phase grid
  struct gkyl_rect_grid phaseGrid_vl;
  struct gkyl_range phaseRange_vl, phaseRange_ext_vl;
  gkyl_rect_grid_init(&phaseGrid_vl, pdim_vl, lower_vl, upper_vl, cells_vl);
  gkyl_create_grid_ranges(&phaseGrid_vl, ghost_vl, &phaseRange_ext_vl, &phaseRange_vl);

  // basis functions
  struct gkyl_basis phaseBasis_vl, phaseBasis_gk, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis_gk, cdim, vdim_gk);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // vlasov basis func
  gkyl_cart_modal_hybrid(&phaseBasis_vl, cdim, vdim_vl);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);
  gkyl_proj_on_basis *projM2_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_h_ion, NULL);

  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_self = h_ion_mass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_vl,
    .cbasis = &basis,
    .pbasis = &phaseBasis_vl,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_vl,
    .mass_self = h_ion_mass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_RECVR,
    .all_gk = all_gk,
    .base = basepath,
  };
  
  // coll struct.
  struct gkyl_dg_recomb *coll_recomb_up_elc = gkyl_dg_recomb_new(&rec_inp_elc, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_ion = gkyl_dg_recomb_new(&rec_inp_ion, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_neut = gkyl_dg_recomb_new(&rec_inp_rcvr, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, (vdim_vl+1)*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_recomb = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // arrays necessary for prim_vars
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);
 
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  const double *cv_e; const double *cv_i; const double *cv_r;

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *prim_vars_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim_vl+1)*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_recomb_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(b_i_cu, b_i);
    gkyl_array_copy(prim_vars_cu, prim_vars);
    gkyl_array_copy(coef_recomb_cu, coef_recomb);
	
    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_neut, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(prim_vars_cu);
    gkyl_array_release(coef_recomb_cu);
  }
  else {

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_neut, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    
  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coef_recomb, "ctest_h_coef_recomb.gkyl");
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_i[i]*check_fac, 1e-12) );
  }
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }
  // test against predicted value
  double p1_vals[] = {4.0651315620487753e-19, 0.0000000000000000e+00, 0.0000000000000000e+00,
		      0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
		      0.0000000000000000e+00, 0.0000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }

  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_ion);
  gkyl_array_release(b_i); gkyl_array_release(coef_recomb); gkyl_array_release(prim_vars);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_neut);
}

// tests when donor species is gk
void
test_coll_recomb_all_gk_li(bool use_gpu)
{
  int charge_state = 1;
  bool all_gk = true;
  // use vt = 4 eV for all grids
  double vmax_elc = 4.*sqrt(4.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*4.*echarge/(2.*B0);
  double vmax_ion = 4.*sqrt(4.*echarge/li_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*4.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;

  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
  int ghostc[] = {0, 0, 0};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {2, 2, 2, 16, 8};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghostc, &confRange_ext, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost, &phaseRange_ext_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim, lower_ion, upper_ion, cells);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost, &phaseRange_ext_ion, &phaseRange_ion);

  // basis functions
  struct gkyl_basis phaseBasis, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis, cdim, vdim);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);
  gkyl_proj_on_basis *projM2_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_li_ion, NULL);
  
  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_self = li_ion_mass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_self = li_ion_mass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_RECVR,
    .all_gk = all_gk,
    .base = basepath,
  };
  
  // coll struct.
  struct gkyl_dg_recomb *coll_recomb_up_elc = gkyl_dg_recomb_new(&rec_inp_elc, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_ion = gkyl_dg_recomb_new(&rec_inp_ion, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_recvr = gkyl_dg_recomb_new(&rec_inp_rcvr, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_recomb = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

  // arrays necessary for prim vars
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);

  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  const double *cv_e; const double *cv_i; const double *cv_r;
  
  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *prim_vars_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_recomb_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    
    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(b_i_cu, b_i);
    gkyl_array_copy(prim_vars_cu, prim_vars);
    gkyl_array_copy(coef_recomb_cu, coef_recomb);
    
    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(prim_vars_cu);
    gkyl_array_release(coef_recomb_cu);
  }
  else {

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    
  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coef_recomb, "ctest_li_coef_recomb.gkyl");
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_i[i]*check_fac, 1e-12) );
  }
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }
  // test against predicted value
  double p1_vals[] = {1.4761368114720401e-18, 0.0000000000000000e+00, 0.0000000000000000e+00,
  		      0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
  		      0.0000000000000000e+00, 0.0000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_ion);
  gkyl_array_release(b_i); gkyl_array_release(coef_recomb); gkyl_array_release(prim_vars);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_recvr);
}

void
test_coll_recomb_all_gk_ar(bool use_gpu)
{
  int charge_state = 1;
  bool all_gk = true;
  // use vt = 4 eV for all grids
  double vmax_elc = 4.*sqrt(4.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*4.*echarge/(2.*B0);
  double vmax_ion = 4.*sqrt(4.*echarge/ar_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*4.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;

  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
  int ghostc[] = {0, 0, 0};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {2, 2, 2, 16, 8};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghostc, &confRange_ext, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost, &phaseRange_ext_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim, lower_ion, upper_ion, cells);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost, &phaseRange_ext_ion, &phaseRange_ion);

  // basis functions
  struct gkyl_basis phaseBasis, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis, cdim, vdim);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);
  gkyl_proj_on_basis *projM2_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_li_ion, NULL);

  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_AR,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_self = ar_ion_mass,
    .type_ion = GKYL_ION_AR,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_self = ar_ion_mass,
    .type_ion = GKYL_ION_AR,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_RECVR,
    .all_gk = all_gk,
    .base = basepath,
  };
  
  // coll struct.
  struct gkyl_dg_recomb *coll_recomb_up_elc = gkyl_dg_recomb_new(&rec_inp_elc, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_ion = gkyl_dg_recomb_new(&rec_inp_ion, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_recvr = gkyl_dg_recomb_new(&rec_inp_rcvr, use_gpu);

    struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_recomb = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

  // arrays necessary for prim vars
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);

  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);
  
  const double *cv_e; const double *cv_i; const double *cv_r;

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *prim_vars_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_recomb_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    
    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(b_i_cu, b_i);
    gkyl_array_copy(prim_vars_cu, prim_vars);
    gkyl_array_copy(coef_recomb_cu, coef_recomb);
    
    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc_cu, moms_ion_cu, b_i_cu, prim_vars_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(prim_vars_cu);
    gkyl_array_release(coef_recomb_cu);
  }
  else {

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_e = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_i = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc, moms_ion, b_i, prim_vars, coef_recomb, 0);
    cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));
    
  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coef_recomb, "ctest_ar_coef_recomb.gkyl");
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_i[i]*check_fac, 1e-12) );
  }
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }
  // test against predicted value
  double p1_vals[] = {2.6716460249415115e-18, 0.0000000000000000e+00, 0.0000000000000000e+00,
		      0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
		      0.0000000000000000e+00, 0.0000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_ion);
  gkyl_array_release(b_i); gkyl_array_release(coef_recomb); gkyl_array_release(prim_vars);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_recvr);
  
}

void
test_coll_recomb_init_elem(bool use_gpu)
{
  int charge_state = 0; // charge state of reacting species
  bool all_gk = false;
  // use vt = 4 eV for all grids
  double vmax_elc = 4.*sqrt(4.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*4.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 3, vdim_gk = 2;
  int pdim_gk = cdim + vdim_gk;
  double imass = h_ion_mass;
  
  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  int ghost[] = {0, 0, 0};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim_gk, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost_gk, &phaseRange_ext_elc, &phaseRange_elc);

  // basis functions
  struct gkyl_basis  phaseBasis_gk, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis_gk, cdim, vdim_gk);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  struct gkyl_dg_recomb_inp rec_inp_he = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_HE,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_be = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_BE,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_b = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_B,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_c = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_C,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_n = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_N,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_o = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_ION_O,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
    
  // coll struct.
  struct gkyl_dg_recomb *coll_recomb_up_he = gkyl_dg_recomb_new(&rec_inp_he, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_be = gkyl_dg_recomb_new(&rec_inp_be, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_b = gkyl_dg_recomb_new(&rec_inp_b, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_c = gkyl_dg_recomb_new(&rec_inp_c, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_n = gkyl_dg_recomb_new(&rec_inp_n, use_gpu);
  struct gkyl_dg_recomb *coll_recomb_up_o = gkyl_dg_recomb_new(&rec_inp_o, use_gpu);

  gkyl_dg_recomb_release(coll_recomb_up_he);
  gkyl_dg_recomb_release(coll_recomb_up_be);
  gkyl_dg_recomb_release(coll_recomb_up_b);
  gkyl_dg_recomb_release(coll_recomb_up_c);
  gkyl_dg_recomb_release(coll_recomb_up_n);
  gkyl_dg_recomb_release(coll_recomb_up_o);

}

void coll_recomb_h() { test_coll_recomb_h(false); }
void coll_recomb_all_gk_li() { test_coll_recomb_all_gk_li(false); }
void coll_recomb_all_gk_ar() { test_coll_recomb_all_gk_ar(false); }
void coll_recomb_init_elem() { test_coll_recomb_init_elem(false); }

#ifdef GKYL_HAVE_CUDA
void coll_recomb_h_gpu() { test_coll_recomb_h(true); }
void coll_recomb_li_gpu() { test_coll_recomb_all_gk_li(true); }
void coll_recomb_ar_gpu() { test_coll_recomb_all_gk_ar(true); }
void coll_recomb_init_elem_gpu() { test_coll_recomb_init_elem(true); }
#endif

TEST_LIST = {
#ifdef GKYL_HAVE_ADAS
  { "coll_recomb_h", coll_recomb_h },
  { "coll_recomb_all_gk_li", coll_recomb_all_gk_li },
  { "coll_recomb_all_gk_ar", coll_recomb_all_gk_ar },
  { "coll_recomb_init_elem", coll_recomb_init_elem },
#ifdef GKYL_HAVE_CUDA
  { "coll_recomb_h_gpu", coll_recomb_h_gpu },
  { "coll_recomb_li_gpu", coll_recomb_li_gpu },
  { "coll_recomb_ar_gpu", coll_recomb_ar_gpu },
#endif
#endif
  { NULL, NULL },
};
