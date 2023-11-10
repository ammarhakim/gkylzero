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
double h_ion_mass = GKYL_PROTON_MASS*2.01410177811;
double li_ion_mass = GKYL_PROTON_MASS*6.94;
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
void eval_m2_3v_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/emass*1.0e19*3.;  //fabs(x);
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
test_coll_iz_h(bool use_gpu)
{
  int charge_state = 0; // charge state of reacting species
  bool all_gk = false;
  // use vt = 40 eV for all grids
  double vmax_elc = 4*sqrt(40*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12*40*echarge/(2*B0);
  double vmax_ion = 4*sqrt(40*echarge/h_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim_gk = 2, vdim_vl = 3;
  int pdim_gk = cdim + vdim_gk, pdim_vl = cdim + vdim_vl;
  char basepath[4000] = ".";

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_ion,mumax_ion};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};

  // for vlasov grid 
  double lower_vl[] = {-2.0,-2.0,-2.0,vmin_ion,vmin_ion,vmin_ion}, upper_vl[] = {2.0,2.0,2.0,vmax_ion,vmax_ion,vmax_ion};
  int ghost_vl[] = {0, 0, 0, 0, 0, 0};
  int cells_vl[] = {16, 16, 16, 4, 4, 4};
  //int cells_vl[] = {2, 2, 2, 4, 4, 4};
  
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
  gkyl_proj_on_basis *projBmag = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *projJac = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_jac, NULL);

  // maxwellian on basis for fdist
  gkyl_proj_maxwellian_on_basis *proj_max_elc = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_elc,
    &basis, &phaseBasis_gk, poly_order+1, use_gpu);
  gkyl_proj_maxwellian_on_basis *proj_max_ion = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_ion,
    &basis, &phaseBasis_gk, poly_order+1, use_gpu);
  gkyl_proj_maxwellian_on_basis *proj_max_neut = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_vl,
    &basis, &phaseBasis_vl, poly_order+1, use_gpu);

  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_elc = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
  echarge, emass, h_ion_mass, GKYL_IZ_H, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ion = gkyl_dg_iz_new(&phaseGrid_ion, &basis, &phaseBasis_gk, &confRange, &phaseRange_ion,
  echarge, emass, h_ion_mass, GKYL_IZ_H, charge_state, GKYL_IZ_ION, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_neut = gkyl_dg_iz_new(&phaseGrid_vl, &basis, &phaseBasis_vl, &confRange, &phaseRange_vl,
 echarge, emass, h_ion_mass, GKYL_IZ_H, charge_state, GKYL_IZ_DONOR, all_gk, basepath, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, (vdim_vl+2)*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_iz_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);

  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_iz_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);

  struct gkyl_array *cflRate_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
  struct gkyl_array *distf_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
  struct gkyl_array *coll_iz_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);	

  // arrays necessary for fmax
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot);

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2_ion, (vdim_vl+1)*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);

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
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_iz_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_iz_elc_cu, coll_iz_elc);
    gkyl_array_copy(bmag_cu, bmag);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);

    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc_cu, moms_neut_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_elc_cu);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_neut_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_iz_elc_cu); gkyl_array_release(b_i_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_elc, &confRange, moms_ion,
      bmag, jacob_tot, h_ion_mass, distf_ion);

    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max_neut, &phaseRange_vl, &confRange, moms_neut, distf_neut);
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_ion, "ctest_distf_ion.gkyl");

    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc, moms_neut, bmag, jacob_tot, b_i, distf_elc, coll_iz_elc, cflRate_elc);
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_elc, moms_neut, bmag, jacob_tot, b_i, distf_ion, coll_iz_ion, cflRate_ion);
    gkyl_dg_iz_coll(coll_iz_up_neut, moms_elc, moms_neut, bmag, jacob_tot, b_i, distf_neut, coll_iz_neut, cflRate_neut);

  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coll_iz_up_ion->prim_vars_donor, "ctest_coll_iz_pvneut_gk.gkyl");
  //gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_iz_elc, "ctest_coll_iz_h_elc.gkyl");
  //gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_iz_ion, "ctest_coll_iz_h_ion.gkyl");
  //gkyl_grid_sub_array_write(&phaseGrid_vl, &phaseRange_vl, coll_iz_neut, "ctest_coll_iz_h_neut.gkyl");
  
  double p1_vals_elc[] = {-2.6709857935892035e+02,  7.7615561179731124e-15, -3.0274206010288621e-14,
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

  double p1_vals_ion[] = {6.0136474363708623e+07,  3.2064301630171979e-10, -1.8565299152428219e-09,
  			  -6.0570764407789609e-10,  8.3576702187810987e+06, -5.1806089123801298e+07,
  			  2.8222875586315712e-09, 4.8287882169437475e-10,  4.8287882169437475e-10,
  			  2.0472687401060010e-09,  1.1570178767292057e-10,  1.9919351944203088e-09,
  			  3.4611839236930172e-09,  5.3367646943141775e-09,  5.3367646943141775e-09,
  			  -7.1999267134092050e+06, -6.0570764407789609e-10,  4.1828177034791418e-10,
  			  2.6699177901041736e-10,  1.9745342021115280e-11, -1.5421051609666856e-09,
  			  -6.0431477565610552e-10, -6.0431477565610563e-10, -5.0164525626716406e-10,
  			  1.3352650343934799e-09,  8.4442763228554184e-11, -1.3154464931638150e-10,
  			  3.3347560965447498e-10, -1.7622242733005421e-10, -4.5889832050750040e-11,
  			  -4.5889832050750047e-11,  8.4442763228554184e-11, -1.6891750936774116e+06,
  			  6.2269932920488623e-11, -7.4936579198133358e-12, -7.4936579198133423e-12,
  			  1.4551826719919532e+06, -1.6543240608186146e-10, -1.3485507070008304e-10,
  			  -2.8531743029089381e-10,  5.6714153091620196e-10, -4.2365841266130257e-10,
  			  -1.0490695482437653e-09, -1.0427773531830465e-10,  1.9674008170038323e-10,
  			  1.7039846159847806e-10,  1.7039846159847808e-10,  1.4405684149657291e-10};

  double p1_vals_neut[] = {-4.6069610426637011e+09,  5.3662813205196892e-07, -5.7742406785171177e-07,
  			   -2.5721356643349079e-07, -2.2048493689564290e+09, -2.2048493689564295e+09,
  			   -2.2048493689564290e+09, -1.0379287078447303e-07, -2.0397967899871426e-08,
  			   -2.0397967899871406e-08, -1.3105576805762424e-07,  1.8892618174858080e-07,
  			   1.8892618174858080e-07, -1.5106892439626306e-07, -1.5129747600827899e-07,
  			   1.6891302540994199e-07, -1.0552207181193717e+09, -1.7108208073490185e-07,
  			   -1.1205381637807305e-08,  3.0900511978041366e-07, -1.0552207181193712e+09,
  			   -1.0552207181193719e+09,  6.2996934984730174e-08, -5.1003142703068993e-08,
  			   -5.1117418509076954e-08, -5.1117418509076960e-08,  9.0363263128474223e-09,
  			   8.9220505068394783e-09,  8.9220505068394701e-09, -1.2628306158864632e-07,
  			   1.6020974308410889e-07,  1.0449237499840138e-10, -2.9115266062146068e-08,
  			   -1.1091105831799326e-08, -1.1091105831799336e-08, -1.6630937426592395e-07,
  			   7.2025259862518362e-08,  2.0117648713637208e-08, -4.6230436234091072e-08,
  			   1.5207788521707361e-07,  1.6020974308410891e-07, -5.0501897300829470e+08,
  			   2.8820931039470327e-08,  8.8077747008315162e-09, -6.2041235568134447e-09,
  			   1.2897411814213646e-08,  2.1029269681248928e-08,  1.0682495003115640e-08,
  			   -6.2041235568134356e-09, -7.1157445244251615e-09, -7.1157445244251623e-09,
  			   -4.6230436234091072e-08, -7.1157445244251648e-09, -7.1157445244251714e-09,
  			   -2.6367758781821597e-09, -6.4405761402105407e-08, -6.4405761402105421e-08,
  			   -8.0273654920368848e-09, -8.0273654920368881e-09, -8.0273654920368766e-09,
  			   -2.6367758781821601e-09, -1.3508112301504971e-08, -1.3508112301504969e-08,
  			   -6.1685178162073123e-10,  5.9448047549161889e+07,  1.7261547132013240e-07,
  			   1.0635782415212232e-07,  2.6305198797567073e-08,  2.8451291189706732e+07,
  			   2.8451291189707011e+07,  1.2510220611021907e-08,  1.1434096001733897e-08,
  			   -1.2644989270422552e-08,  5.4320448464936414e-09,  4.4019980786854008e-09,
  			   4.4019980786854016e-09,  1.3563902713528931e-08, -1.5611158259953413e-08,
  			   -1.5611158259953413e-08,  1.3616527434178023e+07,  1.4423900325963542e-08,
  			   9.1814854510618862e-09,  8.6664620671577663e-09,  6.7917417648736431e-09,
  			   7.3067651487777713e-09,  6.7917417648736455e-09,  8.6664620671577679e-09,
  			   6.4438664528131413e-09,  9.7003361064734845e-09,  9.7003361064734845e-09,
  			   8.1514386832536447e-09, -1.3736437957669291e-08, -3.5627117165062643e-09,
  			   -3.8091971919602151e-09, -5.9756587391970555e-11,  2.2014548973369950e-09,
  			   5.9448047549162090e+07,  1.1257600230421598e-07,  4.6318355136205878e-08,
  			   6.2920424589282563e-09,  2.8451291189706620e+07,  2.8451291189707011e+07,
  			   1.2510220611021912e-08,  1.1434096001733902e-08,  7.3681670682162621e-09,
  			   5.4320448464936381e-09,  1.2533855945720685e-08,  1.2533855945720685e-08,
  			   -6.4492536251098791e-09,  1.4408576248004804e-08,  1.4408576248004802e-08,
  			   1.3616527434177954e+07,  6.2920424589282563e-09,  7.3067651487777605e-09,
  			   6.7917417648736431e-09,  6.7917417648736422e-09, -8.2509271825751117e-10,
  			   -1.3401161021616367e-09, -1.3401161021616373e-09,  6.4438664528131388e-09,
  			   -3.0624206284591980e-10, -3.0624206284591985e-10,  6.2767183809695240e-09,
  			   -1.8551394860657620e-09,  6.4438664528131388e-09,  6.1973809773591897e-09,
  			   -5.9756587391970697e-11,  2.2014548973369954e-09,  5.9448047549162194e+07,
  			   2.5036424417025101e-09,  3.5139322312648938e-10, -4.3740848387668774e-08,
  			   2.8451291189706977e+07,  2.8451291189707000e+07, -7.5029357276168959e-09,
  			   -8.5790603369049075e-09, -8.5790603369049058e-09,  3.5573245442095244e-09,
  			   -3.7499034900876345e-08,  2.5272777764012755e-09,  3.5573245442095278e-09,
  			   2.2540434115040086e-08,  6.2566746792317715e-08,  1.3616527434178412e+07,
  			   -3.7145357103911530e-09, -8.2509271825751334e-10, -1.3401161021616369e-09,
  			   6.7917417648736464e-09, -8.2509271825751189e-10, -1.1346694271481042e-08,
  			   -3.2148364044457605e-09,  1.8325164924416679e-08, -3.0624206284592327e-10,
  			   -3.0624206284592378e-10, -1.8551394860657624e-09,  8.1514386832536447e-09,
  			   1.8672888806199108e-10, -1.9344768896760890e-09, -5.9756587391966535e-11,
  			   -2.1809623651300464e-09};
  
  const double *pv_e = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange_elc, (int[5]) { 1, 1, 1, 4, 2}));
  // compare with output from: "pgkyl ctest_coll_iz_h_elc.gkyl sel --z0 0 --z1 0 --z2 0 --z3 3 --z4 1 pr"

  for (int i=0; i<phaseBasis_gk.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_e[i], 1e-12) );
  }

  const double *pv_i = gkyl_array_cfetch(coll_iz_ion, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 4, 2}));
  // compare with output from: "pgkyl ctest_coll_iz_h_ion.gkyl sel --z0 0 --z1 0 --z2 0 --z3 3 --z4 1 pr"

  for (int i=0; i<phaseBasis_gk.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_i[i], 1e-12) );
  }

  const double *pv_n = gkyl_array_cfetch(coll_iz_neut, gkyl_range_idx(&phaseRange_vl, (int[6]) { 1, 1, 1, 2, 2, 2}));
  // compare with output from: "pgkyl ctest_coll_iz_h_neut.gkyl sel --z0 0 --z1 0 --z2 0 --z3 1 --z4 1 --z5 1 pr"

  for (int i=0; i<phaseBasis_vl.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_neut[i], pv_n[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_neut);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_iz_elc); gkyl_array_release(coll_iz_ion); gkyl_array_release(coll_iz_neut);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_proj_maxwellian_on_basis_release(proj_max_neut);
  gkyl_dg_iz_release(coll_iz_up_elc);
  gkyl_dg_iz_release(coll_iz_up_ion);
  gkyl_dg_iz_release(coll_iz_up_neut);
}

void
test_coll_iz_all_gk_li_1x(bool use_gpu)
{
  int charge_state = 1; // charge state of reacting species
  bool all_gk = true;
  // use vt = 40 eV for all grids
  double vmax_elc = 4*sqrt(40*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12*40*echarge/(2*B0);
  double vmax_ion = 4*sqrt(40*echarge/li_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 1, vdim = 2;
  int pdim = cdim + vdim;
  char basepath[4000] = ".";
  
  double lower_elc[] = {-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,vmax_ion,mumax_ion};
  int ghost[] = {0, 0, 0};
  int cells[] = {2, 16, 8};

  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

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
  gkyl_proj_on_basis *projBmag = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *projJac = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_jac, NULL);

  // maxwellian on basis for fdist
  gkyl_proj_maxwellian_on_basis *proj_max_elc = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_elc,
    &basis, &phaseBasis, poly_order+1, use_gpu);
  gkyl_proj_maxwellian_on_basis *proj_max_ion = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_ion,
    &basis, &phaseBasis, poly_order+1, use_gpu);
  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_elc = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis, &confRange, &phaseRange_elc,
						     echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ion = gkyl_dg_iz_new(&phaseGrid_ion, &basis, &phaseBasis, &confRange, &phaseRange_ion,
  						     echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_ION, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_donor = gkyl_dg_iz_new(&phaseGrid_ion, &basis, &phaseBasis, &confRange, &phaseRange_ion,
  						       echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_DONOR, all_gk, basepath, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_donor = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_iz_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);	
  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_iz_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	
  struct gkyl_array *cflRate_donor = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_donor = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_iz_donor = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	

  // arrays necessary for fmax
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot);

  gkyl_array_set_offset(moms_donor, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_donor, 1.0, m2_ion, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  gkyl_array_shiftc(bmag, 0.5*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(jacob_tot, 1.0*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(b_z, 1.0*pow(sqrt(2),cdim), 0);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_donor_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_iz_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_donor_cu, moms_donor);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_iz_elc_cu, coll_iz_elc);
    gkyl_array_copy(bmag_cu, bmag);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);

    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc_cu, moms_donor_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_elc_cu);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_donor_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_iz_elc_cu); gkyl_array_release(b_i_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_donor,
      bmag, jacob_tot, li_ion_mass, distf_ion);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_donor,
      bmag, jacob_tot, li_ion_mass, distf_donor);

    gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, distf_elc, "ctest_coll_iz_li_distf_elc.gkyl");
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_ion, "ctest_coll_iz_li_distf_ion.gkyl");
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_donor, "ctest_coll_iz_li_distf_donor.gkyl");
    
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc, moms_donor, bmag, jacob_tot, b_i, distf_elc, coll_iz_elc, cflRate_elc);
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_donor, moms_donor, bmag, jacob_tot, b_i, distf_ion, coll_iz_ion, cflRate_ion);
    gkyl_dg_iz_coll(coll_iz_up_donor, moms_donor, moms_donor, bmag, jacob_tot, b_i, distf_donor, coll_iz_donor, cflRate_donor);

  }

  gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_iz_elc, "ctest_coll_iz_li_elc_1x.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_iz_ion, "ctest_coll_iz_li_ion_1x.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_iz_donor, "ctest_coll_iz_li_donor_1x.gkyl");
  
  /* double p1_vals_elc[] = {-3.0211866019063268e+03,  1.1151629206897585e-13,  3.6431890306379876e-13, */
  /* 			  2.9260691682239769e-13, -1.0808072282313658e+02, 1.7261432022756057e+03, */
  /* 			  -1.7533165289662860e-13, -4.8930347399217097e-14, -4.8930347399217103e-14, */
  /* 			  -8.0530046129500689e-16,  2.2975665791155755e-15,  2.2975665791155755e-15, */
  /* 			  -5.7757526206287122e-15, -1.3998079811979073e-13,  3.4431743630114681e-15, */
  /* 			  6.1751500182237166e+01,  5.7589718567933035e-15, -1.6153735565021414e-15, */
  /* 			  3.4109651130671742e-16,  1.1511696065138511e-15, -5.7757526206286972e-15, */
  /* 			  -1.1662891288086177e-15, -1.1662891288086169e-15, -4.9806274555816689e-14, */
  /* 			  -8.0319993619460044e-15, -8.0319993619460028e-15,  3.1076396743227091e-15, */
  /* 			  3.4431743630114685e-15,  1.2131640310202148e-14,  1.1013818754303210e-14, */
  /* 			  1.1013818754303212e-14, -8.0319993619460028e-15,  2.6414781912199704e+01, */
  /* 			  2.6536488176685738e-14,  6.0320281695407716e-15,  6.0320281695407716e-15, */
  /* 			  -1.5091982801911699e+01, -8.5094318488076805e-15, -8.9875904769978933e-15, */
  /* 			  7.3202598929381168e-15,  5.8193829106763044e-15,  6.3657710517723692e-15, */
  /* 			  6.3657710517723692e-15,  8.4622474551621719e-15,  8.2496021962977040e-15, */
  /* 			  8.5227962668457375e-15, 8.5227962668457375e-15, -1.6800794278136792e-16};   */

  /* double p1_vals_ion[] = {1.6670249545888491e+04, -3.2417594348322014e-13, -9.2770344417926952e-13, */
  /* 			  -6.2386751751611901e-14,  5.9636588465764521e+02, -9.5244821738970531e+03, */
  /* 			  2.5270185146855172e-13,  2.3937699859641281e-13,  2.3937699859641281e-13, */
  /* 			  1.0785813614986634e-13, -6.1032275456172119e-15, -6.1032275456172127e-15, */
  /* 			  -7.0135428465832752e-13, -6.8092299440848366e-14, -6.8092299440848353e-14, */
  /* 			  -3.4073132630117146e+02, -6.2386751751611901e-14,  1.5487523743201235e-14, */
  /* 			  4.6921480987920046e-15, -4.7156688232485811e-14, -1.2447648970655591e-13, */
  /* 			  4.7935054164240818e-14,  4.7935054164240812e-14, -6.2752950681203300e-14, */
  /* 			  1.4316773372550237e-15,  1.4316773372550249e-15,  8.6267384861047899e-14, */
  /* 			  -6.8092299440848353e-14, -1.0904114349925496e-14, -4.7362185063352344e-15, */
  /* 			  -4.7362185063352352e-15,  1.4316773372550249e-15, -1.4575101249911003e+02, */
  /* 			  1.0988812099679773e-13,  2.5333708527487868e-15,  2.5333708527487845e-15, */
  /* 			  8.3274273522654596e+01, -1.4070439703451535e-14, -1.1432064553296630e-14, */
  /* 			  2.0155883740287529e-14, -6.8403023299617067e-14, -1.4352759840810866e-13, */
  /* 			  6.9185032983427534e-16, -8.7936894031417234e-15, 2.8434502997949620e-14, */
  /* 			  2.6927077628189564e-14,  2.6927077628189560e-14, -1.0635209926056232e-14}; */

  /* double p1_vals_donor[] = {-1.6670249545888491e+04,  3.2417594348322014e-13,  9.2770344417926952e-13, */
  /* 			    6.2386751751611901e-14, -5.9636588465764521e+02, 9.5244821738970531e+03, */
  /* 			    -2.5270185146855172e-13, -2.3937699859641281e-13, -2.3937699859641281e-13, */
  /* 			    -1.0785813614986634e-13,  6.1032275456172119e-15,  6.1032275456172127e-15, */
  /* 			    7.0135428465832752e-13,  6.8092299440848366e-14,  6.8092299440848353e-14, */
  /* 			    3.4073132630117146e+02,  6.2386751751611901e-14, -1.5487523743201235e-14, */
  /* 			    -4.6921480987920046e-15,  4.7156688232485811e-14,  1.2447648970655591e-13, */
  /* 			    -4.7935054164240818e-14, -4.7935054164240812e-14,  6.2752950681203300e-14, */
  /* 			    -1.4316773372550237e-15, -1.4316773372550249e-15, -8.6267384861047899e-14, */
  /* 			    6.8092299440848353e-14,  1.0904114349925496e-14,  4.7362185063352344e-15, */
  /* 			    4.7362185063352352e-15, -1.4316773372550249e-15,  1.4575101249911003e+02, */
  /* 			    -1.0988812099679773e-13, -2.5333708527487868e-15, -2.5333708527487845e-15, */
  /* 			    -8.3274273522654596e+01,  1.4070439703451535e-14,  1.1432064553296630e-14, */
  /* 			    -2.0155883740287529e-14,  6.8403023299617067e-14,  1.4352759840810866e-13, */
  /* 			    -6.9185032983427534e-16,  8.7936894031417234e-15, -2.8434502997949620e-14, */
  /* 			    -2.6927077628189564e-14, -2.6927077628189560e-14, 1.0635209926056232e-14}; */
    
  /* const double *pv_elc = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange_elc, (int[3]) {1, 8, 1})); */
  /* for (int i=0; i<phaseBasis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_elc[i], 1e-30) ); */
  /* } */

  /* const double *pv_ion = gkyl_array_cfetch(coll_iz_ion, gkyl_range_idx(&phaseRange_ion, (int[3]) { 1, 8, 1})); */
  /* for (int i=0; i<phaseBasis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_ion[i], 1e-30) ); */
  /* } */
  
  /* const double *pv_donor = gkyl_array_cfetch(coll_iz_donor, gkyl_range_idx(&phaseRange_ion, (int[3]) {1, 8, 1})); */
  /* for (int i=0; i<phaseBasis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_donor[i], pv_donor[i], 1e-30) ); */
  /* } */
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion); 
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_donor);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(cflRate_ion); gkyl_array_release(distf_ion);
  gkyl_array_release(cflRate_donor); gkyl_array_release(distf_donor);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_iz_elc); gkyl_array_release(coll_iz_ion); gkyl_array_release(coll_iz_donor); 
  gkyl_proj_on_basis_release(projM0);
  gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_dg_iz_release(coll_iz_up_elc);
  gkyl_dg_iz_release(coll_iz_up_ion);
  gkyl_dg_iz_release(coll_iz_up_donor);
}

// tests when donor species is gk
void
test_coll_iz_all_gk_li_3x(bool use_gpu)
{
  int charge_state = 1; // charge state of reacting species
  bool all_gk = true;
  // use vt = 40 eV for all grids
  double vmax_elc = 4*sqrt(40*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12*40*echarge/(2*B0);
  double vmax_ion = 4*sqrt(40*echarge/li_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;
  char basepath[4000] = ".";
  
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_ion,mumax_ion};
  int ghost[] = {0, 0, 0, 0, 0, 0};
  int cells[] = {2, 2, 2, 16, 8};

  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

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
  gkyl_proj_on_basis *projBmag = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *projJac = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_jac, NULL);

  // maxwellian on basis for fdist
  gkyl_proj_maxwellian_on_basis *proj_max_elc = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_elc,
    &basis, &phaseBasis, poly_order+1, use_gpu);
  gkyl_proj_maxwellian_on_basis *proj_max_ion = gkyl_proj_maxwellian_on_basis_new(&phaseGrid_ion,
    &basis, &phaseBasis, poly_order+1, use_gpu);
  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_elc = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis, &confRange, &phaseRange_elc,
						     echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_ion = gkyl_dg_iz_new(&phaseGrid_ion, &basis, &phaseBasis, &confRange, &phaseRange_ion,
  						     echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_ION, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_donor = gkyl_dg_iz_new(&phaseGrid_ion, &basis, &phaseBasis, &confRange, &phaseRange_ion,
  						       echarge, emass, li_ion_mass, GKYL_IZ_LI, charge_state, GKYL_IZ_DONOR, all_gk, basepath, use_gpu);
  
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_donor = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_iz_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);	
  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_iz_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	
  struct gkyl_array *cflRate_donor = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_donor = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_iz_donor = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	

  // arrays necessary for fmax
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *jacob_tot = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2_elc, 0.0, &confRange, m2_elc);
  gkyl_proj_on_basis_advance(projM2_ion, 0.0, &confRange, m2_ion);
  gkyl_proj_on_basis_advance(projBmag, 0.0, &confRange, bmag);
  gkyl_proj_on_basis_advance(projJac, 0.0, &confRange, jacob_tot);

  gkyl_array_set_offset(moms_donor, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_donor, 1.0, m2_ion, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2_elc, 2*basis.num_basis);

  gkyl_array_shiftc(bmag, 0.5*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(jacob_tot, 1.0*pow(sqrt(2),cdim), 0);
  gkyl_array_shiftc(b_z, 1.0*pow(sqrt(2),cdim), 0);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_donor_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_iz_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_donor_cu, moms_donor);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_iz_elc_cu, coll_iz_elc);
    gkyl_array_copy(bmag_cu, bmag);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);

    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc_cu, moms_donor_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_elc_cu);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_donor_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_iz_elc_cu); gkyl_array_release(b_i_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_donor,
      bmag, jacob_tot, li_ion_mass, distf_ion);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_donor,
      bmag, jacob_tot, li_ion_mass, distf_donor);

    gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, distf_elc, "ctest_coll_iz_li_distf_elc.gkyl");
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_ion, "ctest_coll_iz_li_distf_ion.gkyl");
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_donor, "ctest_coll_iz_li_distf_donor.gkyl");
    
    gkyl_dg_iz_coll(coll_iz_up_elc, moms_elc, moms_donor, bmag, jacob_tot, b_i, distf_elc, coll_iz_elc, cflRate_elc);
    gkyl_dg_iz_coll(coll_iz_up_ion, moms_donor, moms_donor, bmag, jacob_tot, b_i, distf_ion, coll_iz_ion, cflRate_ion);
    gkyl_dg_iz_coll(coll_iz_up_donor, moms_donor, moms_donor, bmag, jacob_tot, b_i, distf_donor, coll_iz_donor, cflRate_donor);

  }

  //gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_iz_elc, "ctest_coll_iz_li_elc_3x.gkyl");
  //gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_iz_ion, "ctest_coll_iz_li_ion_3x.gkyl");
  //gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_iz_donor, "ctest_coll_iz_li_donor_3x.gkyl");
  
  double p1_vals_elc[] = {-3.0211866019063268e+03,  1.1151629206897585e-13,  3.6431890306379876e-13,
			  2.9260691682239769e-13, -1.0808072282313658e+02, 1.7261432022756057e+03,
			  -1.7533165289662860e-13, -4.8930347399217097e-14, -4.8930347399217103e-14,
			  -8.0530046129500689e-16,  2.2975665791155755e-15,  2.2975665791155755e-15,
			  -5.7757526206287122e-15, -1.3998079811979073e-13,  3.4431743630114681e-15,
			  6.1751500182237166e+01,  5.7589718567933035e-15, -1.6153735565021414e-15,
			  3.4109651130671742e-16,  1.1511696065138511e-15, -5.7757526206286972e-15,
			  -1.1662891288086177e-15, -1.1662891288086169e-15, -4.9806274555816689e-14,
			  -8.0319993619460044e-15, -8.0319993619460028e-15,  3.1076396743227091e-15,
			  3.4431743630114685e-15,  1.2131640310202148e-14,  1.1013818754303210e-14,
			  1.1013818754303212e-14, -8.0319993619460028e-15,  2.6414781912199704e+01,
			  2.6536488176685738e-14,  6.0320281695407716e-15,  6.0320281695407716e-15,
			  -1.5091982801911699e+01, -8.5094318488076805e-15, -8.9875904769978933e-15,
			  7.3202598929381168e-15,  5.8193829106763044e-15,  6.3657710517723692e-15,
			  6.3657710517723692e-15,  8.4622474551621719e-15,  8.2496021962977040e-15,
			  8.5227962668457375e-15, 8.5227962668457375e-15, -1.6800794278136792e-16};  

  double p1_vals_ion[] = {1.6670249545888491e+04, -3.2417594348322014e-13, -9.2770344417926952e-13,
			  -6.2386751751611901e-14,  5.9636588465764521e+02, -9.5244821738970531e+03,
			  2.5270185146855172e-13,  2.3937699859641281e-13,  2.3937699859641281e-13,
			  1.0785813614986634e-13, -6.1032275456172119e-15, -6.1032275456172127e-15,
			  -7.0135428465832752e-13, -6.8092299440848366e-14, -6.8092299440848353e-14,
			  -3.4073132630117146e+02, -6.2386751751611901e-14,  1.5487523743201235e-14,
			  4.6921480987920046e-15, -4.7156688232485811e-14, -1.2447648970655591e-13,
			  4.7935054164240818e-14,  4.7935054164240812e-14, -6.2752950681203300e-14,
			  1.4316773372550237e-15,  1.4316773372550249e-15,  8.6267384861047899e-14,
			  -6.8092299440848353e-14, -1.0904114349925496e-14, -4.7362185063352344e-15,
			  -4.7362185063352352e-15,  1.4316773372550249e-15, -1.4575101249911003e+02,
			  1.0988812099679773e-13,  2.5333708527487868e-15,  2.5333708527487845e-15,
			  8.3274273522654596e+01, -1.4070439703451535e-14, -1.1432064553296630e-14,
			  2.0155883740287529e-14, -6.8403023299617067e-14, -1.4352759840810866e-13,
			  6.9185032983427534e-16, -8.7936894031417234e-15, 2.8434502997949620e-14,
			  2.6927077628189564e-14,  2.6927077628189560e-14, -1.0635209926056232e-14};

  double p1_vals_donor[] = {-1.6670249545888491e+04,  3.2417594348322014e-13,  9.2770344417926952e-13,
			    6.2386751751611901e-14, -5.9636588465764521e+02, 9.5244821738970531e+03,
			    -2.5270185146855172e-13, -2.3937699859641281e-13, -2.3937699859641281e-13,
			    -1.0785813614986634e-13,  6.1032275456172119e-15,  6.1032275456172127e-15,
			    7.0135428465832752e-13,  6.8092299440848366e-14,  6.8092299440848353e-14,
			    3.4073132630117146e+02,  6.2386751751611901e-14, -1.5487523743201235e-14,
			    -4.6921480987920046e-15,  4.7156688232485811e-14,  1.2447648970655591e-13,
			    -4.7935054164240818e-14, -4.7935054164240812e-14,  6.2752950681203300e-14,
			    -1.4316773372550237e-15, -1.4316773372550249e-15, -8.6267384861047899e-14,
			    6.8092299440848353e-14,  1.0904114349925496e-14,  4.7362185063352344e-15,
			    4.7362185063352352e-15, -1.4316773372550249e-15,  1.4575101249911003e+02,
			    -1.0988812099679773e-13, -2.5333708527487868e-15, -2.5333708527487845e-15,
			    -8.3274273522654596e+01,  1.4070439703451535e-14,  1.1432064553296630e-14,
			    -2.0155883740287529e-14,  6.8403023299617067e-14,  1.4352759840810866e-13,
			    -6.9185032983427534e-16,  8.7936894031417234e-15, -2.8434502997949620e-14,
			    -2.6927077628189564e-14, -2.6927077628189560e-14, 1.0635209926056232e-14};
    
  const double *pv_elc = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange_elc, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<phaseBasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_elc[i], 1e-30) );
  }

  const double *pv_ion = gkyl_array_cfetch(coll_iz_ion, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<phaseBasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_ion[i], 1e-30) );
  }
  
  const double *pv_donor = gkyl_array_cfetch(coll_iz_donor, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<phaseBasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_donor[i], pv_donor[i], 1e-30) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion); 
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_donor);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(cflRate_ion); gkyl_array_release(distf_ion);
  gkyl_array_release(cflRate_donor); gkyl_array_release(distf_donor);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_iz_elc); gkyl_array_release(coll_iz_ion); gkyl_array_release(coll_iz_donor); 
  gkyl_proj_on_basis_release(projM0);
  gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_dg_iz_release(coll_iz_up_elc);
  gkyl_dg_iz_release(coll_iz_up_ion);
  gkyl_dg_iz_release(coll_iz_up_donor);
}

void
test_coll_iz_init_elem(bool use_gpu)
{
  int charge_state = 0; // charge state of reacting species
  bool all_gk = false;
  // use vt = 40 eV for all grids
  double vmax_elc = 4*sqrt(40*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim_gk = 2;
  int pdim_gk = cdim + vdim_gk;
  double imass = h_ion_mass;
  char basepath[4000] = ".";

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange_ext, &confRange);

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

  // coll struct.
  struct gkyl_dg_iz *coll_iz_up_he = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_HE, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_be = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_BE, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_b = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_B, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_c = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_C, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_n = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_N, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);
  struct gkyl_dg_iz *coll_iz_up_o = gkyl_dg_iz_new(&phaseGrid_elc, &basis, &phaseBasis_gk, &confRange, &phaseRange_elc,
    echarge, emass, imass, GKYL_IZ_O, charge_state, GKYL_IZ_ELC, all_gk, basepath, use_gpu);

  gkyl_dg_iz_release(coll_iz_up_he);
  gkyl_dg_iz_release(coll_iz_up_be);
  gkyl_dg_iz_release(coll_iz_up_b);
  gkyl_dg_iz_release(coll_iz_up_c);
  gkyl_dg_iz_release(coll_iz_up_n);
  gkyl_dg_iz_release(coll_iz_up_o);

}
  
void prim_vars_gk_3x() { test_prim_vars_gk_3x(false); }
void prim_vars_vlasov_3x() { test_prim_vars_vlasov_3x(false); }
void coll_iz_h() { test_coll_iz_h(false); }
void coll_iz_all_gk_li_1x() { test_coll_iz_all_gk_li_1x(false); }
void coll_iz_all_gk_li_3x() { test_coll_iz_all_gk_li_3x(false); }
void coll_iz_init_elem() { test_coll_iz_init_elem(false); }

#ifdef GKYL_HAVE_CUDA
void prim_vars_gk_3x_gpu() { test_prim_vars_gk_3x(true); }
void prim_vars_vlasov_3x_gpu() { test_prim_vars_vlasov_3x(true); }
void coll_iz_h_gpu() { test_coll_iz_h(true); }
void coll_iz_all_gk_li_gpu() { test_coll_iz_all_gk_li(true); }
void coll_iz_init_elem_gpu() { test_coll_iz_init_elem(true); }
#endif

TEST_LIST = {
  { "prim_vars_gk_3x", prim_vars_gk_3x },
  { "prim_vars_vlasov_3x", prim_vars_vlasov_3x },
  //{ "coll_iz_h", coll_iz_h },
  //{ "coll_iz_all_gk_li_1x", coll_iz_all_gk_li_1x },
  //{ "coll_iz_all_gk_li_3x", coll_iz_all_gk_li_3x },
  //{ "coll_iz_init_elem", coll_iz_init_elem },
#ifdef GKYL_HAVE_CUDA
  { "coll_iz_h_gpu", coll_iz_h_gpu },
#endif  
  { NULL, NULL },
};
