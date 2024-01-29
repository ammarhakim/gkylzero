#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis.h>
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
test_coll_recomb_h(bool use_gpu)
{
  int charge_state = 1;
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
  /* double lower_elc[] = {-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,vmax_elc,mumax_elc}; */
  /* double lower_ion[] = {-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,vmax_elc,mumax_ion}; */
  /* int ghost_gk[] = {0, 0, 0}; */
  /* int cells_gk[] = {16, 8, 4}; */

  /* // for vlasov grid  */
  /* double lower_vl[] = {-2.0,vmin_ion,vmin_ion,vmin_ion}, upper_vl[] = {2.0,vmax_ion,vmax_ion,vmax_ion}; */
  /* int ghost_vl[] = {0, 0, 0, 0}; */
  /* int cells_vl[] = {16, 4, 4, 4}; */

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
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
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim_gk, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost_gk, &phaseRange_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim_gk, lower_ion, upper_ion, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost_gk, &phaseRange_ion, &phaseRange_ion);

  // vlasov phase grid
  struct gkyl_rect_grid phaseGrid_vl;
  struct gkyl_range phaseRange_vl, phaseRange_ext_vl;
  gkyl_rect_grid_init(&phaseGrid_vl, pdim_vl, lower_vl, upper_vl, cells_vl);
  gkyl_create_grid_ranges(&phaseGrid_vl, ghost_vl, &phaseRange_vl, &phaseRange_vl);

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


  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_H,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_ion,
    .mass_self = h_ion_mass,
    .type_ion = GKYL_RECOMB_H,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_vl,
    .cbasis = &basis,
    .pbasis = &phaseBasis_vl,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_vl,
    .mass_self = h_ion_mass,
    .type_ion = GKYL_RECOMB_H,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_RECVR,
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
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, (vdim_vl+2)*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_recomb_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);

  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_recomb_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);

  struct gkyl_array *cflRate_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
  struct gkyl_array *distf_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
  struct gkyl_array *coll_recomb_neut = gkyl_array_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);

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
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_recomb_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_elc.volume);

    struct gkyl_array *cflRate_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
    struct gkyl_array *distf_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);
    struct gkyl_array *coll_recomb_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_gk.num_basis, phaseRange_ion.volume);

    struct gkyl_array *cflRate_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
    struct gkyl_array *distf_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);
    struct gkyl_array *coll_recomb_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis_vl.num_basis, phaseRange_vl.volume);

    /* // arrays necessary for fmax */
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_recomb_elc_cu, coll_recomb_elc);
    gkyl_array_copy(cflRate_ion_cu, cflRate_ion);
    gkyl_array_copy(distf_ion_cu, distf_ion);
    gkyl_array_copy(coll_recomb_ion_cu, coll_recomb_ion);
    gkyl_array_copy(cflRate_neut_cu, cflRate_neut);
    gkyl_array_copy(distf_neut_cu, distf_neut);
    gkyl_array_copy(coll_recomb_neut_cu, coll_recomb_neut);
    gkyl_array_copy(bmag_cu, bmag);
    gkyl_array_copy(jacob_tot_cu, jacob_tot);
    gkyl_array_copy(b_i_cu, b_i);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
  					    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion_cu,
  					    bmag_cu, jacob_tot_cu, h_ion_mass, distf_ion_cu);

    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max_neut, &phaseRange_vl, &confRange, moms_neut_cu, distf_neut_cu);

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_recomb_elc_cu, cflRate_elc_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_ion_cu, coll_recomb_ion_cu, cflRate_ion_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_neut, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_neut_cu, coll_recomb_neut_cu, cflRate_neut_cu);

    gkyl_array_copy(coll_recomb_elc, coll_recomb_elc_cu);
    gkyl_array_copy(coll_recomb_ion, coll_recomb_ion_cu);
    gkyl_array_copy(coll_recomb_neut, coll_recomb_neut_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu); gkyl_array_release(moms_neut_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(cflRate_ion_cu); gkyl_array_release(distf_ion_cu);
    gkyl_array_release(cflRate_neut_cu); gkyl_array_release(distf_neut_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(coll_recomb_elc_cu); gkyl_array_release(coll_recomb_ion_cu);
    gkyl_array_release(coll_recomb_neut_cu); gkyl_array_release(b_i_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_elc, &confRange, moms_ion,
      bmag, jacob_tot, h_ion_mass, distf_ion);

    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max_neut, &phaseRange_vl, &confRange, moms_neut, distf_neut);
    gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, distf_ion, "ctest_distf_ion.gkyl");

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_elc, coll_recomb_elc, cflRate_elc);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_ion, coll_recomb_ion, cflRate_ion);
    gkyl_dg_recomb_coll(coll_recomb_up_neut, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_neut, coll_recomb_neut, cflRate_neut);

  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, coll_recomb_up_ion->prim_vars_donor, "ctest_coll_recomb_pvneut_gk.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_recomb_elc, "ctest_coll_recomb_h_elc.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_recomb_ion, "ctest_coll_recomb_h_ion.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_vl, &phaseRange_vl, coll_recomb_neut, "ctest_coll_recomb_h_neut.gkyl");
  
  double p1_vals_elc[] = {-7.5558273382519907e-05,  2.1927534265894878e-21, -8.5546987300309449e-21,
  			  8.8336896621295990e-22, -1.0500966974091547e-05,  6.5091588529467040e-05,
  			  -3.2004281141213150e-21, -4.8438188136532691e-22, -4.8438188136532691e-22,
  			  -7.3025917372293973e-23, -1.2186831692817091e-21,  1.2961221589599148e-22,
  			  8.8504916745170347e-21,  9.1905152181260347e-21,  3.7973336774152327e-21,
  			  9.0463239939151549e-06,  8.8336896621295990e-22, -2.5056339787148950e-22,
  			  -7.3462328357659931e-22,  1.1706188951144658e-22,  7.6071936345083000e-22,
  			  9.3073113525533038e-22,  9.3073113525533019e-22,  5.6175643208369726e-22,
  			  -1.5443937308085060e-21, -2.8926891159862066e-21,  3.0714969639518704e-22,
  			  -1.5958478632955701e-21,  8.0556125945253324e-22,  6.4180530320528906e-22,
  			  6.4180530320528916e-22, -1.9609834563080553e-22,  2.1223584333711394e-06,
  			  -1.2186844989803665e-22, -5.9591522454890456e-22, -2.5884137825447938e-22,
  			  -1.8283594326944808e-06,  4.1380999228898569e-22,  3.8317354561723602e-23,
  			  3.1525585536537977e-22, -7.1461088444487552e-22, -3.1134335606125186e-22,
  			  2.5730490233173281e-23, -1.0143687111327883e-25, -4.9800517384831248e-22,
  			  -1.2783448650928816e-22, -1.2783448650928814e-22,  2.4233620082973625e-22};

  double p1_vals_ion[] = {-3.4798692003357380e-10, -2.1925165945430789e-27,  2.0692633082743279e-26,
  			  1.1927820637455578e-28, -4.6687344508551290e-10,  2.9978214692373033e-10,
  			  -1.2479194032727442e-26, -1.0366191940842621e-27, -1.0366191940842619e-27,
  			  1.9945073094734589e-26,  2.6031981365050588e-27,  2.6031981365050577e-27,
  			  -3.2367995014154995e-27, -3.8033950984148370e-27, -3.8033950984148370e-27,
  			  4.0219995537737419e-10,  1.1927820637455578e-28, -1.4299426067060579e-26,
  			  -5.8481139652777606e-27, -2.4636471180356747e-27,  7.0498779367688624e-27,
  			  1.6232414191770120e-27,  1.6232414191770120e-27,  1.3560912005869649e-26,
  			  -6.9564179719402553e-27, -6.9564179719402553e-27, -4.2990124544372175e-27,
  			  -3.8033950984148363e-27,  2.8178130675450196e-26,  1.0610856351754971e-26,
  			  1.0610856351754971e-26, -6.9564179719402553e-27, -3.1124896339034192e-10,
  			  9.5141396987827718e-27,  3.7269447018246597e-27,  3.7269447018246597e-27,
  			  2.6813330358491620e-10, -1.2818087049436047e-26, -2.0405402624886542e-27,
  			  -3.6661352378806472e-27,  1.5847866788658081e-26,  1.6427081712827953e-26,
  			  -4.1462731635407688e-27, -1.5496709137256227e-27,  7.9995818135617057e-27,
  			  -1.9974881625377186e-27, -1.9974881625377182e-27, -1.7078807004527799e-27};

  double p1_vals_neut[] = {1.3017989126513389e+03, -1.5163616803658902e-13,  1.6316396355577849e-13,
  			   7.2681391920060560e-14,  6.2302903898833654e+02,  6.2302903898833665e+02,
  			   6.2302903898833654e+02,  2.9328975234846853e-14,  5.7638977595947346e-15,
  			   5.7638977595947314e-15,  3.7032710885581366e-14, -5.3385278428460364e-14,
  			   -5.3385278428460358e-14,  4.2687871612813743e-14,  4.2752453934489937e-14,
  			   -4.7730117701227994e-14,  2.9817599296665884e+02,  4.8343032340046113e-14,
  			   3.1663288438633443e-15, -8.7316242791854584e-14,  2.9817599296665873e+02,
  			   2.9817599296665890e+02, -1.7801179715657380e-14,  1.4412067976651881e-14,
  			   1.4444359137489984e-14,  1.4444359137489984e-14, -2.5534142050452252e-15,
  			   -2.5211230442071294e-15, -2.5211230442071266e-15,  3.5684076930533398e-14,
  			   -4.5270812453507405e-14, -2.9526635648437629e-17,  8.2271634924313288e-15,
  			   3.1340376830252393e-15,  3.1340376830252425e-15,  4.6994398384998138e-14,
  			   -2.0352332937949917e-14, -5.6846873628808074e-15,  1.3063434021603912e-14,
  			   -4.2972975846879397e-14, -4.5270812453507411e-14,  1.4270429982849873e+02,
  			   -8.1439926106013987e-15, -2.4888318833690280e-15,  1.7531125671391692e-15,
  			   -3.6444494581730025e-15, -5.9422860648010064e-15, -3.0185756403572110e-15,
  			   1.7531125671391649e-15,  2.0107112690593659e-15,  2.0107112690593667e-15,
  			   1.3063434021603912e-14,  2.0107112690593678e-15,  2.0107112690593690e-15,
  			   7.4507944376671968e-16,  1.8199274833299862e-14,  1.8199274833299865e-14,
  			   2.2683099709795656e-15,  2.2683099709795680e-15,  2.2683099709795645e-15,
  			   7.4507944376671939e-16,  3.8170164113009205e-15,  3.8170164113009197e-15,
  			   1.7430513762638888e-16, -1.6798362942960384e+01, -4.8776325823111631e-14,
  			   -3.0053759636990748e-14, -7.4331167280612656e-15, -8.0395426814531739e+00,
  			   -8.0395426814532538e+00, -3.5350400052526599e-15, -3.2309571546977760e-15,
  			   3.5731218758485955e-15, -1.5349446216610717e-15, -1.2438820860991280e-15,
  			   -1.2438820860991282e-15, -3.8327812282890772e-15,  4.4112786411332433e-15,
  			   4.4112786411332441e-15, -3.8476515090414316e+00, -4.0757926074568973e-15,
  			   -2.5944321356374354e-15, -2.4489008678564632e-15, -1.9191571108682818e-15,
  			   -2.0646883786492559e-15, -1.9191571108682826e-15, -2.4489008678564640e-15,
  			   -1.8208572340547445e-15, -2.7410448837783925e-15, -2.7410448837783925e-15,
  			   -2.3033696000754911e-15,  3.8815348841450622e-15,  1.0067231295614412e-15,
  			   1.0763730616877989e-15,  1.6885547711435368e-17, -6.2206985582566530e-16,
  			   -1.6798362942960441e+01, -3.1810843641414521e-14, -1.3088277455293637e-14,
  			   -1.7779560008288930e-15, -8.0395426814531437e+00, -8.0395426814532538e+00,
  			   -3.5350400052526611e-15, -3.2309571546977772e-15, -2.0820388513837756e-15,
  			   -1.5349446216610713e-15, -3.5417186927271317e-15, -3.5417186927271317e-15,
  			   1.8223794989432931e-15, -4.0714624497153124e-15, -4.0714624497153124e-15,
  			   -3.8476515090414121e+00, -1.7779560008288930e-15, -2.0646883786492527e-15,
  			   -1.9191571108682814e-15, -1.9191571108682814e-15,  2.3314822797874768e-16,
  			   3.7867949575972114e-16,  3.7867949575972129e-16, -1.8208572340547437e-15,
  			   8.6535479837792766e-17,  8.6535479837792779e-17, -1.7736258430873102e-15,
  			   5.2421076354069455e-16, -1.8208572340547437e-15, -1.7512073019283863e-15,
  			   1.6885547711435414e-17, -6.2206985582566540e-16, -1.6798362942960470e+01,
  			   -7.0745964163647679e-16, -9.9293940526708622e-17,  1.2359945817252034e-14,
  			   -8.0395426814532449e+00, -8.0395426814532502e+00,  2.1201207219797089e-15,
  			   2.4242035725345931e-15,  2.4242035725345927e-15, -1.0052008646728916e-15,
  			   1.0596183125353797e-14, -7.1413832911094539e-16, -1.0052008646728926e-15,
  			   -6.3692990563433163e-15, -1.7679620510808059e-14, -3.8476515090415417e+00,
  			   1.0496243627872934e-15,  2.3314822797874827e-16, 3.7867949575972119e-16,
  			   -1.9191571108682826e-15,  2.3314822797874783e-16,  3.2062598593759063e-15,
  			   9.0842325274790317e-16, -5.1781813546591149e-15,  8.6535479837793752e-17,
  			   8.6535479837793876e-17,  5.2421076354069455e-16, -2.3033696000754911e-15,
  			   -5.2764384414925519e-17,  5.4662930469961596e-16,  1.6885547711434234e-17,
  			   6.1627923682597552e-16};
  
  const double *pv_e = gkyl_array_cfetch(coll_recomb_elc, gkyl_range_idx(&phaseRange_elc, (int[5]) { 1, 1, 1, 4, 2}));
  // compare with output from: "pgkyl ctest_coll_recomb_h_elc.gkyl sel --z0 0 --z1 0 --z2 0 --z3 3 --z4 1 pr"

  for (int i=0; i<phaseBasis_gk.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_e[i], 1e-12) );
  }

  const double *pv_i = gkyl_array_cfetch(coll_recomb_ion, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 4, 2}));
  // compare with output from: "pgkyl ctest_coll_recomb_h_ion.gkyl sel --z0 0 --z1 0 --z2 0 --z3 3 --z4 1 pr"

  for (int i=0; i<phaseBasis_gk.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_i[i], 1e-12) );
  }

  const double *pv_n = gkyl_array_cfetch(coll_recomb_neut, gkyl_range_idx(&phaseRange_vl, (int[6]) { 1, 1, 1, 2, 2, 2}));
  // compare with output from: "pgkyl ctest_coll_recomb_h_neut.gkyl sel --z0 0 --z1 0 --z2 0 --z3 1 --z4 1 --z5 1 pr"

  for (int i=0; i<phaseBasis_vl.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_neut[i], pv_n[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion);
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_neut); gkyl_array_release(moms_ion);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(cflRate_ion); gkyl_array_release(distf_ion);
  gkyl_array_release(cflRate_neut); gkyl_array_release(distf_neut);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_recomb_elc); gkyl_array_release(coll_recomb_ion); gkyl_array_release(coll_recomb_neut);
  gkyl_proj_on_basis_release(projM0); gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_proj_maxwellian_on_basis_release(proj_max_neut);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_neut);
}

// tests when donor species is gk
void
test_coll_recomb_all_gk_li(bool use_gpu)
{
  int charge_state = 2;
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

  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {2, 2, 2, 16, 8};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost, &phaseRange_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim, lower_ion, upper_ion, cells);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost, &phaseRange_ion, &phaseRange_ion);

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

  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_LI,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_ion,
    .mass_self = li_ion_mass,
    .type_ion = GKYL_RECOMB_LI,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_ion,
    .mass_self = li_ion_mass,
    .type_ion = GKYL_RECOMB_LI,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_RECVR,
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
  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_recomb_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);	
  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_recomb_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	
  struct gkyl_array *cflRate_recvr = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_recvr = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_recomb_recvr = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	

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

  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);
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
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_recomb_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *cflRate_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *distf_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *coll_recomb_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *cflRate_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *distf_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *coll_recomb_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);


    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_recomb_elc_cu, coll_recomb_elc);
    gkyl_array_copy(cflRate_ion_cu, cflRate_ion);
    gkyl_array_copy(distf_ion_cu, distf_ion);
    gkyl_array_copy(coll_recomb_ion_cu, coll_recomb_ion);
    gkyl_array_copy(cflRate_recvr_cu, cflRate_recvr);
    gkyl_array_copy(distf_recvr_cu, distf_elc);
    gkyl_array_copy(coll_recomb_recvr_cu, coll_recomb_recvr);
    gkyl_array_copy(bmag_cu, bmag);
    gkyl_array_copy(jacob_tot_cu, jacob_tot);
    gkyl_array_copy(b_i_cu, b_i);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion_cu,
    bmag_cu, jacob_tot_cu, li_ion_mass, distf_ion_cu);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion_cu,
    bmag_cu, jacob_tot_cu, li_ion_mass, distf_recvr_cu);

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_recomb_elc_cu, cflRate_elc_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_ion_cu, coll_recomb_ion_cu, cflRate_ion_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_recvr_cu, coll_recomb_recvr_cu, cflRate_recvr_cu);

    gkyl_array_copy(coll_recomb_elc, coll_recomb_elc_cu);
    gkyl_array_copy(coll_recomb_ion, coll_recomb_ion_cu);
    gkyl_array_copy(coll_recomb_recvr, coll_recomb_recvr_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(cflRate_ion_cu); gkyl_array_release(distf_ion_cu);
    gkyl_array_release(cflRate_recvr_cu); gkyl_array_release(distf_recvr_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(coll_recomb_elc_cu);
    gkyl_array_release(coll_recomb_ion_cu); gkyl_array_release(coll_recomb_recvr_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion,
      bmag, jacob_tot, li_ion_mass, distf_ion);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion,
      bmag, jacob_tot, li_ion_mass, distf_recvr);
    
    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_elc, coll_recomb_elc, cflRate_elc);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_ion, coll_recomb_ion, cflRate_ion);
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_recvr, coll_recomb_recvr, cflRate_recvr);

  }

  gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_recomb_elc, "ctest_coll_recomb_li_elc.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_recomb_ion, "ctest_coll_recomb_li_ion.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_recomb_recvr, "ctest_coll_recomb_li_recvr.gkyl");
  
  double p1_vals_elc[] = {-1.1294702866949660e+00,  4.1690353814898518e-17,  1.3620058278830629e-16,
			  1.0939106443269294e-16, -4.0405966621267372e-02,  6.4531845081022388e-01,
			  -6.5547719607554868e-17, -1.8292605120850970e-17, -1.8292605120850973e-17,
			  -3.0106149097859049e-19,  8.5894501888000295e-19,  8.5894501888000295e-19,
			  -2.1592644969974534e-18, -5.2331806345357383e-17,  1.2872303658693029e-18,
			  2.3085791711993757e-02,  2.1529910102395784e-18, -6.0390723063272411e-19,
			  1.2751889412363940e-19,  4.3036463377777294e-19, -2.1592644969974522e-18,
			  -4.3601706556407344e-19, -4.3601706556407325e-19, -1.8620070394284992e-17,
			  -3.0027620989536669e-18, -3.0027620989536665e-18,  1.1617907585341364e-18,
			  1.2872303658693027e-18,  4.5354124272225516e-18,  4.1175149585861100e-18,
			  4.1175149585861108e-18, -3.0027620989536665e-18,  9.8751633813455655e-03,
			  9.9206632552543403e-18,  2.2550730834382822e-18,  2.2550730834382822e-18,
			  -5.6421361498541498e-03, -3.1812501828981977e-18, -3.3600097347004253e-18,
			  2.7366783748946434e-18,  2.1755756762465429e-18,  2.3798428241218515e-18,
			  2.3798428241218519e-18,  3.1636103024006822e-18,  3.0841128952089429e-18,
			  3.1862464691465977e-18,  3.1862464691465973e-18, -6.2809751367415434e-20};

  double p1_vals_ion[] = {-7.6846242945232746e-04, -7.8389955672629886e-21,  4.6163486443475541e-20,
			  -6.1988540948317704e-21, -1.0310005382460250e-03,  4.3905801712689564e-04,
			  -3.4020165836416649e-20, -7.0189248310473829e-21, -7.0189248310473829e-21,
			  -1.8145675374875724e-20,  9.3224172047938343e-21,  9.3224172047938343e-21,
			  2.6168620978263334e-20,  1.0273022834532274e-20,  3.7277302672438594e-21,
			  5.8905814341720308e-04, -6.1988540948317696e-21, -2.8003789361644293e-20,
			  -9.3406860784252285e-21, -2.5663742360810315e-20, -1.2549290890324190e-23,
			  -1.4150557954674405e-21, -1.4150557954674397e-21,  4.0570915168827166e-20,
			  -3.7039667479594889e-21, -3.7039667479594889e-21,  1.9180531191562402e-20,
			  3.7277302672438594e-21,  1.7622216047481774e-20,  6.9591246497611424e-21,
			  6.9591246497611424e-21, -3.7039667479594896e-21, -6.8733369216401657e-04,
			  -4.3672002711301248e-21, -4.1743692221083202e-20,  1.0618648317224111e-20,
			  3.9270542894480212e-04, -1.7417610269871002e-20, -4.9755414143916293e-21,
			  -1.4873830536996144e-20, -1.0029322723470916e-20, -1.1156192785854932e-20,
			  -4.6109002185665186e-21, -5.6240576934890849e-21,  6.4544341022579593e-21,
			  -7.1995860635108777e-21, -7.1995860635108777e-21, -1.2177285274144719e-21};

  double p1_vals_recvr[] = {7.6846242945232746e-04,  7.8389955672629886e-21, -4.6163486443475541e-20,
			    6.1988540948317704e-21,  1.0310005382460250e-03, -4.3905801712689564e-04,
			    3.4020165836416649e-20,  7.0189248310473829e-21,  7.0189248310473829e-21,
			    1.8145675374875724e-20, -9.3224172047938343e-21, -9.3224172047938343e-21,
			    -2.6168620978263334e-20, -1.0273022834532274e-20, -3.7277302672438594e-21,
			    -5.8905814341720308e-04,  6.1988540948317696e-21,  2.8003789361644293e-20,
			    9.3406860784252285e-21,  2.5663742360810315e-20,  1.2549290890324190e-23,
			    1.4150557954674405e-21,  1.4150557954674397e-21, -4.0570915168827166e-20,
			    3.7039667479594889e-21,  3.7039667479594889e-21, -1.9180531191562402e-20,
			    -3.7277302672438594e-21, -1.7622216047481774e-20, -6.9591246497611424e-21,
			    -6.9591246497611424e-21,  3.7039667479594896e-21,  6.8733369216401657e-04,
			    4.3672002711301248e-21,  4.1743692221083202e-20, -1.0618648317224111e-20,
			    -3.9270542894480212e-04,  1.7417610269871002e-20,  4.9755414143916293e-21,
			    1.4873830536996144e-20,  1.0029322723470916e-20,  1.1156192785854932e-20,
			    4.6109002185665186e-21,  5.6240576934890849e-21, -6.4544341022579593e-21,
			    7.1995860635108777e-21,  7.1995860635108777e-21,  1.2177285274144719e-21};
    
  const double *pv_elc = gkyl_array_cfetch(coll_recomb_elc, gkyl_range_idx(&phaseRange_elc, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<phaseBasis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_elc[i], 1e-12) );
  }

  const double *pv_ion = gkyl_array_cfetch(coll_recomb_ion, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_ion[i], 1e-12) );
  }
  
  const double *pv_recvr = gkyl_array_cfetch(coll_recomb_recvr, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1}));
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals_recvr[i], pv_recvr[i], 1e-12) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion); 
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_ion);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(cflRate_ion); gkyl_array_release(distf_ion);
  gkyl_array_release(cflRate_recvr); gkyl_array_release(distf_recvr);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_recomb_elc); gkyl_array_release(coll_recomb_ion); gkyl_array_release(coll_recomb_recvr); 
  gkyl_proj_on_basis_release(projM0);
  gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_recvr);
}

void
test_coll_recomb_all_gk_ar(bool use_gpu)
{
  int charge_state = 2;
  bool all_gk = true;
  // use vt = 40 eV for all grids
  double vmax_elc = 4*sqrt(40*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12*40*echarge/(2*B0);
  double vmax_ion = 4*sqrt(40*echarge/ar_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12*40*echarge/(2*B0);
  int poly_order = 1;
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;

  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  double lower_ion[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_ion[] = {2.0,2.0,2.0,vmax_elc,mumax_ion};
  int ghost[] = {0, 0, 0, 0, 0};
  int cells[] = {2, 2, 2, 16, 8};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim, lower_elc, upper_elc, cells);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost, &phaseRange_elc, &phaseRange_elc);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim, lower_ion, upper_ion, cells);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost, &phaseRange_ion, &phaseRange_ion);

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

  struct gkyl_dg_recomb_inp rec_inp_elc = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_AR,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_ion = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_ion,
    .mass_self = ar_ion_mass,
    .type_ion = GKYL_RECOMB_AR,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ION,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_rcvr = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis = &phaseBasis,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_ion,
    .mass_self = ar_ion_mass,
    .type_ion = GKYL_RECOMB_AR,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_RECVR,
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
  struct gkyl_array *cflRate_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *distf_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
  struct gkyl_array *coll_recomb_elc = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);	
  struct gkyl_array *cflRate_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_recomb_ion = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	
  struct gkyl_array *cflRate_recvr = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *distf_recvr = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
  struct gkyl_array *coll_recomb_recvr = gkyl_array_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);	

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

  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, m2_ion, 2*basis.num_basis);
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
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *cflRate_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *distf_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *coll_recomb_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_elc.volume);
    struct gkyl_array *cflRate_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *distf_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *coll_recomb_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *cflRate_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *distf_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);
    struct gkyl_array *coll_recomb_recvr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, phaseBasis.num_basis, phaseRange_ion.volume);


    // arrays necessary for fmax
    struct gkyl_array *bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *b_i_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_elc_cu, moms_elc);
    gkyl_array_copy(cflRate_elc_cu, cflRate_elc);
    gkyl_array_copy(distf_elc_cu, distf_elc);
    gkyl_array_copy(coll_recomb_elc_cu, coll_recomb_elc);
    gkyl_array_copy(cflRate_ion_cu, cflRate_ion);
    gkyl_array_copy(distf_ion_cu, distf_ion);
    gkyl_array_copy(coll_recomb_ion_cu, coll_recomb_ion);
    gkyl_array_copy(cflRate_recvr_cu, cflRate_recvr);
    gkyl_array_copy(distf_recvr_cu, distf_elc);
    gkyl_array_copy(coll_recomb_recvr_cu, coll_recomb_recvr);
    gkyl_array_copy(bmag_cu, bmag);
    gkyl_array_copy(jacob_tot_cu, jacob_tot);
    gkyl_array_copy(b_i_cu, b_i);

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc_cu,
    bmag_cu, jacob_tot_cu, emass, distf_elc_cu);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion_cu,
    bmag_cu, jacob_tot_cu, ar_ion_mass, distf_ion_cu);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion_cu,
    bmag_cu, jacob_tot_cu, ar_ion_mass, distf_recvr_cu);

    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_elc_cu, coll_recomb_elc_cu, cflRate_elc_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_ion_cu, coll_recomb_ion_cu, cflRate_ion_cu);
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc_cu, moms_ion_cu, bmag_cu, jacob_tot_cu, b_i_cu, distf_recvr_cu, coll_recomb_recvr_cu, cflRate_recvr_cu);

    gkyl_array_copy(coll_recomb_elc, coll_recomb_elc_cu);
    gkyl_array_copy(coll_recomb_ion, coll_recomb_ion_cu);
    gkyl_array_copy(coll_recomb_recvr, coll_recomb_recvr_cu);

    gkyl_array_release(moms_elc_cu); gkyl_array_release(moms_ion_cu);
    gkyl_array_release(cflRate_elc_cu); gkyl_array_release(distf_elc_cu);
    gkyl_array_release(cflRate_ion_cu); gkyl_array_release(distf_ion_cu);
    gkyl_array_release(cflRate_recvr_cu); gkyl_array_release(distf_recvr_cu);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
    gkyl_array_release(b_i_cu); gkyl_array_release(coll_recomb_elc_cu);
    gkyl_array_release(coll_recomb_ion_cu); gkyl_array_release(coll_recomb_recvr_cu);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_elc, &phaseRange_elc, &confRange, moms_elc,
      bmag, jacob_tot, emass, distf_elc);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion,
      bmag, jacob_tot, ar_ion_mass, distf_ion);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max_ion, &phaseRange_ion, &confRange, moms_ion,
      bmag, jacob_tot, ar_ion_mass, distf_recvr);
    
    gkyl_dg_recomb_coll(coll_recomb_up_elc, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_elc, coll_recomb_elc, cflRate_elc);
    gkyl_dg_recomb_coll(coll_recomb_up_ion, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_ion, coll_recomb_ion, cflRate_ion);
    gkyl_dg_recomb_coll(coll_recomb_up_recvr, moms_elc, moms_ion, bmag, jacob_tot, b_i, distf_recvr, coll_recomb_recvr, cflRate_recvr);

  }

  gkyl_grid_sub_array_write(&phaseGrid_elc, &phaseRange_elc, coll_recomb_elc, "ctest_coll_recomb_ar_elc.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_recomb_ion, "ctest_coll_recomb_ar_ion.gkyl");
  gkyl_grid_sub_array_write(&phaseGrid_ion, &phaseRange_ion, coll_recomb_recvr, "ctest_coll_recomb_ar_recvr.gkyl");
  
  /* double p1_vals_elc[] = {-1.1294702866949660e+00,  4.1690353814898518e-17,  1.3620058278830629e-16, */
  /* 			  1.0939106443269294e-16, -4.0405966621267372e-02,  6.4531845081022388e-01, */
  /* 			  -6.5547719607554868e-17, -1.8292605120850970e-17, -1.8292605120850973e-17, */
  /* 			  -3.0106149097859049e-19,  8.5894501888000295e-19,  8.5894501888000295e-19, */
  /* 			  -2.1592644969974534e-18, -5.2331806345357383e-17,  1.2872303658693029e-18, */
  /* 			  2.3085791711993757e-02,  2.1529910102395784e-18, -6.0390723063272411e-19, */
  /* 			  1.2751889412363940e-19,  4.3036463377777294e-19, -2.1592644969974522e-18, */
  /* 			  -4.3601706556407344e-19, -4.3601706556407325e-19, -1.8620070394284992e-17, */
  /* 			  -3.0027620989536669e-18, -3.0027620989536665e-18,  1.1617907585341364e-18, */
  /* 			  1.2872303658693027e-18,  4.5354124272225516e-18,  4.1175149585861100e-18, */
  /* 			  4.1175149585861108e-18, -3.0027620989536665e-18,  9.8751633813455655e-03, */
  /* 			  9.9206632552543403e-18,  2.2550730834382822e-18,  2.2550730834382822e-18, */
  /* 			  -5.6421361498541498e-03, -3.1812501828981977e-18, -3.3600097347004253e-18, */
  /* 			  2.7366783748946434e-18,  2.1755756762465429e-18,  2.3798428241218515e-18, */
  /* 			  2.3798428241218519e-18,  3.1636103024006822e-18,  3.0841128952089429e-18, */
  /* 			  3.1862464691465977e-18,  3.1862464691465973e-18, -6.2809751367415434e-20}; */

  /* double p1_vals_ion[] = {-7.6846242945232746e-04, -7.8389955672629886e-21,  4.6163486443475541e-20, */
  /* 			  -6.1988540948317704e-21, -1.0310005382460250e-03,  4.3905801712689564e-04, */
  /* 			  -3.4020165836416649e-20, -7.0189248310473829e-21, -7.0189248310473829e-21, */
  /* 			  -1.8145675374875724e-20,  9.3224172047938343e-21,  9.3224172047938343e-21, */
  /* 			  2.6168620978263334e-20,  1.0273022834532274e-20,  3.7277302672438594e-21, */
  /* 			  5.8905814341720308e-04, -6.1988540948317696e-21, -2.8003789361644293e-20, */
  /* 			  -9.3406860784252285e-21, -2.5663742360810315e-20, -1.2549290890324190e-23, */
  /* 			  -1.4150557954674405e-21, -1.4150557954674397e-21,  4.0570915168827166e-20, */
  /* 			  -3.7039667479594889e-21, -3.7039667479594889e-21,  1.9180531191562402e-20, */
  /* 			  3.7277302672438594e-21,  1.7622216047481774e-20,  6.9591246497611424e-21, */
  /* 			  6.9591246497611424e-21, -3.7039667479594896e-21, -6.8733369216401657e-04, */
  /* 			  -4.3672002711301248e-21, -4.1743692221083202e-20,  1.0618648317224111e-20, */
  /* 			  3.9270542894480212e-04, -1.7417610269871002e-20, -4.9755414143916293e-21, */
  /* 			  -1.4873830536996144e-20, -1.0029322723470916e-20, -1.1156192785854932e-20, */
  /* 			  -4.6109002185665186e-21, -5.6240576934890849e-21,  6.4544341022579593e-21, */
  /* 			  -7.1995860635108777e-21, -7.1995860635108777e-21, -1.2177285274144719e-21}; */

  /* double p1_vals_recvr[] = {7.6846242945232746e-04,  7.8389955672629886e-21, -4.6163486443475541e-20, */
  /* 			    6.1988540948317704e-21,  1.0310005382460250e-03, -4.3905801712689564e-04, */
  /* 			    3.4020165836416649e-20,  7.0189248310473829e-21,  7.0189248310473829e-21, */
  /* 			    1.8145675374875724e-20, -9.3224172047938343e-21, -9.3224172047938343e-21, */
  /* 			    -2.6168620978263334e-20, -1.0273022834532274e-20, -3.7277302672438594e-21, */
  /* 			    -5.8905814341720308e-04,  6.1988540948317696e-21,  2.8003789361644293e-20, */
  /* 			    9.3406860784252285e-21,  2.5663742360810315e-20,  1.2549290890324190e-23, */
  /* 			    1.4150557954674405e-21,  1.4150557954674397e-21, -4.0570915168827166e-20, */
  /* 			    3.7039667479594889e-21,  3.7039667479594889e-21, -1.9180531191562402e-20, */
  /* 			    -3.7277302672438594e-21, -1.7622216047481774e-20, -6.9591246497611424e-21, */
  /* 			    -6.9591246497611424e-21,  3.7039667479594896e-21,  6.8733369216401657e-04, */
  /* 			    4.3672002711301248e-21,  4.1743692221083202e-20, -1.0618648317224111e-20, */
  /* 			    -3.9270542894480212e-04,  1.7417610269871002e-20,  4.9755414143916293e-21, */
  /* 			    1.4873830536996144e-20,  1.0029322723470916e-20,  1.1156192785854932e-20, */
  /* 			    4.6109002185665186e-21,  5.6240576934890849e-21, -6.4544341022579593e-21, */
  /* 			    7.1995860635108777e-21,  7.1995860635108777e-21,  1.2177285274144719e-21}; */
    
  /* const double *pv_elc = gkyl_array_cfetch(coll_recomb_elc, gkyl_range_idx(&phaseRange_elc, (int[5]) { 1, 1, 1, 8, 1})); */
  /* for (int i=0; i<phaseBasis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_elc[i], pv_elc[i], 1e-12) ); */
  /* } */

  /* const double *pv_ion = gkyl_array_cfetch(coll_recomb_ion, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1})); */
  /* for (int i=0; i<basis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_ion[i], pv_ion[i], 1e-12) ); */
  /* } */
  
  /* const double *pv_recvr = gkyl_array_cfetch(coll_recomb_recvr, gkyl_range_idx(&phaseRange_ion, (int[5]) { 1, 1, 1, 8, 1})); */
  /* for (int i=0; i<basis.num_basis; ++i) { */
  /*   TEST_CHECK( gkyl_compare_double(p1_vals_recvr[i], pv_recvr[i], 1e-12) ); */
  /* } */
  
  gkyl_array_release(m0); gkyl_array_release(m2_elc); gkyl_array_release(m2_ion); 
  gkyl_array_release(b_x); gkyl_array_release(b_y); gkyl_array_release(b_z);
  gkyl_array_release(moms_elc); gkyl_array_release(moms_ion);
  gkyl_array_release(cflRate_elc); gkyl_array_release(distf_elc);
  gkyl_array_release(cflRate_ion); gkyl_array_release(distf_ion);
  gkyl_array_release(cflRate_recvr); gkyl_array_release(distf_recvr);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot); gkyl_array_release(b_i);
  gkyl_array_release(coll_recomb_elc); gkyl_array_release(coll_recomb_ion); gkyl_array_release(coll_recomb_recvr); 
  gkyl_proj_on_basis_release(projM0);
  gkyl_proj_on_basis_release(projM2_elc); gkyl_proj_on_basis_release(projM2_ion);
  gkyl_proj_on_basis_release(projBmag); gkyl_proj_on_basis_release(projJac);
  gkyl_proj_maxwellian_on_basis_release(proj_max_elc);
  gkyl_proj_maxwellian_on_basis_release(proj_max_ion);
  gkyl_dg_recomb_release(coll_recomb_up_elc);
  gkyl_dg_recomb_release(coll_recomb_up_ion);
  gkyl_dg_recomb_release(coll_recomb_up_recvr);
}

void
test_coll_recomb_init_elem(bool use_gpu)
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
  
  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};
  char basepath[4000] = ".";

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim_gk, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost_gk, &phaseRange_elc, &phaseRange_elc);

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
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_HE,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_be = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_BE,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_b = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_B,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_c = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_C,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_n = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_N,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
    .all_gk = all_gk,
    .base = basepath,
  };
  struct gkyl_dg_recomb_inp rec_inp_o = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .phase_rng = &phaseRange_elc,
    .mass_self = emass,
    .type_ion = GKYL_RECOMB_O,
    .charge_state = charge_state,
    .type_self = GKYL_RECOMB_ELC,
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
void coll_recomb_init_elem() { test_coll_recomb_init_elem(false); }

#ifdef GKYL_HAVE_CUDA
void coll_recomb_h_gpu() { test_coll_recomb_h(true); }
void coll_recomb_li_gpu() { test_coll_recomb_all_gk_li(true); }
void coll_recomb_ar_gpu() { test_coll_recomb_all_gk_ar(true); }
void coll_recomb_init_elem_gpu() { test_coll_recomb_init_elem(true); }
#endif

TEST_LIST = {
  /* { "coll_recomb_h", coll_recomb_h }, */
  /* { "coll_recomb_all_gk_li", coll_recomb_all_gk_li }, */
  /* { "coll_recomb_init_elem", coll_recomb_init_elem }, */
#ifdef GKYL_HAVE_CUDA
  { "coll_recomb_h_gpu", coll_recomb_h_gpu },
  { "coll_recomb_li_gpu", coll_recomb_li_gpu },
  { "coll_recomb_ar_gpu", coll_recomb_ar_gpu },
#endif  
  { NULL, NULL },
};
