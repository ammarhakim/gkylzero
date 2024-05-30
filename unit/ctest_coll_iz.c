#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_transform.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double h_ion_mass = GKYL_PROTON_MASS*2.01410177811;
double li_ion_mass = GKYL_PROTON_MASS*6.94;
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
  fout[0] = 40.*echarge/emass*1.0e19*1.;  //fabs(x);
}
void eval_m2_3v_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/emass*1.0e19*3.;  //fabs(x);
}
void eval_m2_3v_li_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 100.*echarge/emass*1.0e19*3.;  //fabs(x);
}
void eval_m2_3v_h_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/h_ion_mass*1.0e19*3.;  //fabs(x);
}
void eval_m2_3v_li_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/li_ion_mass*1.0e19*3.;  //fabs(x);
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
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;
  double lower[] = {-2.0,-2.0,-2.0}, upper[] = {2.0,2.0,2.0};
  int ghost[] = {0, 0, 0};
  int cells[] = {16, 16, 16};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);

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
  int poly_order = 1;
  int cdim = 3, vdim = 3;
  int pdim = cdim + vdim;
  double lower[] = {-2.0,-2.0,-2.0}, upper[] = {2.0,2.0,2.0};
  int ghost[] = {0, 0, 0};
  int cells[] = {16, 16, 16};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);

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
test_prim_vars_transform_1x(bool use_gpu)
{
  int poly_order = 1;
  int cdim = 1, vdim = 2;
  int pdim = cdim + vdim;
  double lower[] = {-2.0}, upper[] = {2.0};
  int ghost[] = {0};
  int cells[] = {16};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2_3v_elc, NULL);

  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, (2+3)*cbasis.num_basis, confRange.volume);
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, confRange.volume);
  
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);

  struct gkyl_dg_prim_vars_type *calc_prim_vars_trans = gkyl_dg_prim_vars_transform_new(&cbasis, &pbasis, &confRange, "prim_gk", use_gpu);
  
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
 
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m2, (3+1)*cbasis.num_basis);

  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, cbasis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*cbasis.num_basis);
							
  /* gkyl_grid_sub_array_write(&confGrid, &confRange, m0, "ctest_prim_vars_trans_1x_m0.gkyl"); */
  /* gkyl_grid_sub_array_write(&confGrid, &confRange, m2, "ctest_prim_vars_trans_1x_m2.gkyl"); */
  /* gkyl_grid_sub_array_write(&confGrid, &confRange, moms_elc, "ctest_prim_vars_trans_1x_moms_elc.gkyl"); */

  gkyl_dg_prim_vars_transform_set_auxfields(calc_prim_vars_trans, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &confRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(&confRange, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    double *prim_vars_d = gkyl_array_fetch(prim_vars, loc);

    calc_prim_vars_trans->kernel(calc_prim_vars_trans, conf_iter.idx, moms_elc_d, prim_vars_d);

  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, prim_vars, "ctest_prim_vars_trans_1x.gkyl");

  double p1_vals[] = {0.000000000000000e+00, 0.000000000000000e+00, 9.949389233739330e+12,
		      7.676708874031159e-05};
  
  const double *pv = gkyl_array_cfetch(prim_vars, gkyl_range_idx(&confRange, (int[1]) {1}));

  for (int i=0; i<2*cbasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], pv[i], 1e-12) );
  }

  gkyl_array_release(m0); gkyl_array_release(m2);
  gkyl_array_release(moms_elc); gkyl_array_release(prim_vars);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2);
  gkyl_dg_prim_vars_type_release(calc_prim_vars_trans);
  gkyl_array_release(b_x); gkyl_array_release(b_y);
  gkyl_array_release(b_z); gkyl_array_release(b_i);

}

void
test_coll_iz_h_1x(bool use_gpu)
{
  int charge_state = 0; // charge state of reacting species
  bool all_gk = false;
  // use vt = 40 eV for all grids
  double vmax_elc = 4.*sqrt(40.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*40.*echarge/(2.*B0);
  double vmax_ion = 4.*sqrt(40.*echarge/h_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*40.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 1, vdim_gk = 2, vdim_vl = 3;
  int pdim_gk = cdim + vdim_gk, pdim_vl = cdim + vdim_vl;
  char basepath[4000] = ".";

  // for gk grids 
  double lower_elc[] = {-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,vmax_ion,mumax_ion};
  int ghost_gk[] = {0, 0, 0};
  int cells_gk[] = {16, 8, 4};
  int ghost[] = {0};
  
  // for vlasov grid 
  double lower_vl[] = {-2.0,vmin_ion,vmin_ion,vmin_ion}, upper_vl[] = {2.0,vmax_ion,vmax_ion,vmax_ion};
  int ghost_vl[] = {0, 0, 0, 0};
  int cells_vl[] = {16, 4, 4, 4};
  //int cells_vl[] = {2, 2, 2, 4, 4, 4};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

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


  struct gkyl_dg_iz_inp iz_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = h_ion_mass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };

  struct gkyl_dg_iz_inp iz_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_ion = h_ion_mass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ION,
    .all_gk = all_gk,
    .base = basepath,
  };

  struct gkyl_dg_iz_inp iz_inp_donor = {
    .grid = &phaseGrid_vl,
    .cbasis = &basis,
    .pbasis = &phaseBasis_vl,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_vl,
    .mass_ion = h_ion_mass,
    .type_ion = GKYL_ION_H,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_DONOR,
    .all_gk = all_gk,
    .base = basepath,
  };

  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_elc = gkyl_dg_iz_new(&iz_inp_elc, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ion = gkyl_dg_iz_new(&iz_inp_ion, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_neut = gkyl_dg_iz_new(&iz_inp_donor, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, (vdim_vl+2)*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  
  // moment arrays
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_iz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *coef_iz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // arrays necessary for prim vars
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2_ion, (vdim_vl+1)*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  const double *cv_e; const double *cv_i; const double *cv_d;
  
  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim_vl+2)*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *vtSq_iz_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *prim_vars_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_iz_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(vtSq_iz_cu, vtSq_iz);
    gkyl_array_copy(prim_vars_cu, prim_vars_cu);
    gkyl_array_copy(coef_iz_cu, coef_iz);
    gkyl_array_copy(b_i_cu, b_i);
    
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc_cu, moms_neut_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_e = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_elc_cu, moms_neut_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_i = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_neut, moms_elc_cu, moms_neut_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_d = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    
    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_neut_cu);
    gkyl_array_release(prim_vars_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(vtSq_iz_cu);
    gkyl_array_release(coef_iz_cu);
  }
  else {   
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc, moms_neut, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_e = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_elc, moms_neut, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_i = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_neut, moms_elc, moms_neut, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_d = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coef_iz, "ctest_h_coef_iz.gkyl");
  // test that coef are equal for different species
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_i[i]*check_fac, 1e-12) );
  }
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_d[i]*check_fac, 1e-12) );
  }
  // test against predicted value
  double p1_vals[] = {4.1936847897461634e-14, 0.0000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_d[i]*check_fac, 1e-12) );
  }
  
  gkyl_array_release(coef_iz); gkyl_array_release(vtSq_iz); gkyl_array_release(prim_vars);
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_neut);
  gkyl_array_release(b_i);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_dg_iz_release(coll_iz_up_elc);
  gkyl_dg_iz_release(coll_iz_up_ion);
  gkyl_dg_iz_release(coll_iz_up_neut);
 
}

void
test_coll_iz_all_gk_li_1x(bool use_gpu)
{
  int charge_state = 1; // charge state of reacting species
  bool all_gk = true;
  // use vt = 100 eV for all grids
  double vmax_elc = 4.*sqrt(100.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*100.*echarge/(2.*B0);
  double vmax_ion = 4.*sqrt(100.*echarge/li_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*100.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 1, vdim = 2;
  int pdim = cdim + vdim;
  char basepath[4000] = ".";
  
  double lower_elc[] = {-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,vmax_ion,mumax_ion};
  int ghostc[] = {0};
  int ghost[] = {0, 0, 0};
  int cells[] = {2, 16, 8};

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

  struct gkyl_dg_iz_inp iz_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = li_ion_mass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };

  struct gkyl_dg_iz_inp iz_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_ion = li_ion_mass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ION,
    .all_gk = all_gk,
    .base = basepath,
  };

  struct gkyl_dg_iz_inp iz_inp_donor = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_ion = li_ion_mass,
    .type_ion = GKYL_ION_LI,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_DONOR,
    .all_gk = all_gk,
    .base = basepath,
  };
  
  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_li_elc, NULL);
  gkyl_proj_on_basis *projM2_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2_3v_li_ion, NULL);
  gkyl_proj_on_basis *projBmag = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *projJac = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_jac, NULL);

  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_elc = gkyl_dg_iz_new(&iz_inp_elc, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ion = gkyl_dg_iz_new(&iz_inp_ion, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_donor = gkyl_dg_iz_new(&iz_inp_donor, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_donor = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

  // moment arrays
  struct gkyl_array *prim_vars = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_iz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *coef_iz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

  // arrays necessary for prim vars
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);

  gkyl_array_set_offset(moms_donor, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_donor, 1.0, m2_ion, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  const double *cv_e; const double *cv_i; const double *cv_d;
  
  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_donor_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *vtSq_iz_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *prim_vars_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_iz_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_donor_cu, moms_donor);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(vtSq_iz_cu, vtSq_iz);
    gkyl_array_copy(prim_vars_cu, prim_vars_cu);
    gkyl_array_copy(coef_iz_cu, coef_iz);
    gkyl_array_copy(b_i_cu, b_i);
    
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc_cu, moms_donor_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_e = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_elc_cu, moms_donor_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_i = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_donor, moms_elc_cu, moms_donor_cu, b_i_cu, vtSq_iz_cu, prim_vars_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);
    cv_d = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    
    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_donor_cu);
    gkyl_array_release(prim_vars_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(vtSq_iz_cu);
    gkyl_array_release(coef_iz_cu);
  }
  else {
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc, moms_donor, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_e = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_elc, moms_donor, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_i = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
    gkyl_dg_iz_coll(coll_iz_up_donor, moms_elc, moms_donor, b_i, vtSq_iz, prim_vars, coef_iz, 0);
    cv_d = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[1]) { 1}));
  }
    
  //gkyl_grid_sub_array_write(&confGrid, &confRange, coef_iz, "ctest_li_coef_iz.gkyl");
  // test that coef are equal for different species
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_i[i]*check_fac, 1e-12) );
  }
  for (int i=0; i<basis.num_basis; ++i) { 
    TEST_CHECK( gkyl_compare_double(cv_e[i]*check_fac, cv_d[i]*check_fac, 1e-12) );
  }
  //test against predicted value
  double p1_vals[] = {2.3606193318967417e-15, 0.000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_d[i]*check_fac, 1e-12) );
  }
  
  gkyl_array_release(coef_iz); gkyl_array_release(vtSq_iz); gkyl_array_release(prim_vars);
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_donor);
  gkyl_array_release(b_i);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_dg_iz_release(coll_iz_up_elc);
  gkyl_dg_iz_release(coll_iz_up_ion);
  gkyl_dg_iz_release(coll_iz_up_donor);
}

void
test_coll_iz_init_elem(bool use_gpu)
{
  int charge_state = 0; // charge state of reacting species
  bool all_gk = true;
  // use vt = 40 eV for all grids
  double vmax_elc = 4.*sqrt(40.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*40.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 3, vdim_gk = 2;
  int pdim_gk = cdim + vdim_gk;
  double imass = h_ion_mass;
  char basepath[4000] = ".";

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};
  int ghost[] = {0, 0, 0};
  
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


  struct gkyl_dg_iz_inp iz_inp_he = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_HE,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_be = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_BE,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_b = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_B,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_c = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_C,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_n = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_N,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_o = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_O,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_iz_inp iz_inp_ar = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .mass_ion = imass,
    .type_ion = GKYL_ION_AR,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };

  
  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_he = gkyl_dg_iz_new(&iz_inp_he, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_be = gkyl_dg_iz_new(&iz_inp_be, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_b = gkyl_dg_iz_new(&iz_inp_b, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_c = gkyl_dg_iz_new(&iz_inp_c, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_n = gkyl_dg_iz_new(&iz_inp_n, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_o = gkyl_dg_iz_new(&iz_inp_o, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ar = gkyl_dg_iz_new(&iz_inp_ar, use_gpu);

  gkyl_dg_iz_release(coll_iz_up_he);
  gkyl_dg_iz_release(coll_iz_up_be);
  gkyl_dg_iz_release(coll_iz_up_b);
  gkyl_dg_iz_release(coll_iz_up_c);
  gkyl_dg_iz_release(coll_iz_up_n);
  gkyl_dg_iz_release(coll_iz_up_o);
  gkyl_dg_iz_release(coll_iz_up_ar);

}

void prim_vars_gk_3x() { test_prim_vars_gk_3x(false); }
void prim_vars_vlasov_3x() { test_prim_vars_vlasov_3x(false); }
void prim_vars_trans_1x() { test_prim_vars_transform_1x(false); }
void coll_iz_h_1x() { test_coll_iz_h_1x(false); }
void coll_iz_all_gk_li_1x() { test_coll_iz_all_gk_li_1x(false); }
void coll_iz_init_elem() { test_coll_iz_init_elem(false); }

#ifdef GKYL_HAVE_CUDA 
void prim_vars_gk_3x_gpu() { test_prim_vars_gk_3x(true); }
void prim_vars_vlasov_3x_gpu() { test_prim_vars_vlasov_3x(true); }
void prim_vars_trans_1x_gpu() { test_prim_vars_transform_1x(true); }
void coll_iz_h_1x_gpu() { test_coll_iz_h_1x(true); }
void coll_iz_all_gk_li_1x_gpu() { test_coll_iz_all_gk_li_1x(true); }
void coll_iz_init_elem_gpu() { test_coll_iz_init_elem(true); }
#endif

TEST_LIST = {
#ifdef GKYL_HAVE_ADAS
  { "prim_vars_gk_3x", prim_vars_gk_3x },
  { "prim_vars_vlasov_3x", prim_vars_vlasov_3x },
  { "prim_vars_trans_1x", prim_vars_trans_1x },
  { "coll_iz_h_1x", coll_iz_h_1x },
  { "coll_iz_all_gk_li_1x", coll_iz_all_gk_li_1x },
  { "coll_iz_init_elem", coll_iz_init_elem },
#ifdef GKYL_HAVE_CUDA
  { "coll_iz_h_1x_gpu", coll_iz_h_1x_gpu },
  { "coll_iz_all_gk_li_1x_gpu", coll_iz_all_gk_li_1x_gpu },
  { "coll_iz_init_elem_gpu", coll_iz_init_elem_gpu },
#endif
#endif
  { NULL, NULL },
};
