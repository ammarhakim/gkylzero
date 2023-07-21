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
void eval_m2(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*echarge/emass*1.0e19;  //fabs(x);
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
test_prim_vars_3x(bool use_gpu)
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
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_m2, NULL);

  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 5*cbasis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_neut = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, confRange.volume);

  struct gkyl_dg_prim_vars_type *calc_prim_vars_vlasov_vtSq = gkyl_dg_prim_vars_vlasov_new(&cbasis, &pbasis, "vtSq", use_gpu);
  struct gkyl_dg_prim_vars_type *calc_prim_vars_gk_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(&cbasis, &pbasis, "vtSq", use_gpu);
  
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);
 
  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2, 4*cbasis.num_basis); //HARDCODED for gk vdim = 2
  gkyl_array_set_offset(moms_elc, 1.0, m2, 2*cbasis.num_basis);

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &confRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(&confRange, conf_iter.idx);

    // populate the field variables
    // calculate vtsq
  }

  //gkyl_grid_sub_array_write(&confGrid, &confRange, moms_neut, "ctest_prim_vars_vtsq_neut.gkyl");
  //gkyl_grid_sub_array_write(&confGrid, &confRange, moms_elc, "ctest_prim_vars_vtsq_elc.gkyl");
  // check vals

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
    poly_order+1, 1, eval_m2, NULL);
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
      gkyl_dg_iz_coll(coll_iz, moms_elc_cu, moms_neut_cu, bmag_cu, jacob_tot, b_i_cu, distf_elc_cu, coll_iz_elc_cu, cflRate_cu);
      tm_tot = tm_tot + gkyl_time_diff_now_sec(tm);
    }
    tm_tot = tm_tot/iter;
    printf("Avg time over %d loop(s) is %.e s", iter, tm_tot);

    gkyl_array_copy(coll_iz_elc, coll_iz_elc_cu);
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

  double p1_vals[] = {-8.3463180471577168e-09, -2.5357290020268411e-25,  1.5128738678010901e-24,
		      7.8073497887977595e-25,  3.1614574052200733e-09,  8.3394417363177697e-09,
		      -2.5357290020268416e-25, -1.0248840512211112e-25, -1.0248840512211114e-25,
		      1.9631603927709042e-25,  3.4276915366603188e-26, -5.7240445748561069e-26,
		      -9.0446838856348759e-25, -8.4031898521952369e-25,  6.2395879262310432e-25,
		      -3.1588527640234747e-09,  4.8596089958461935e-26, 1.4873389941337646e-25,
		      9.1505407389989828e-26, -4.3947174976624737e-26,  1.9373994481848335e-25,
		      4.2779924260136905e-26,  4.2779924260136882e-26,  9.6360470203878242e-26,
		      7.5287924529321644e-26,  7.5287924529321644e-26, -9.6583058848471189e-27,
		      -1.0818009629820964e-25, -3.9074761433474810e-26, -7.3410779567240854e-26,
		      -7.3410779567240854e-26,  7.5287924529321644e-26,  2.7031868880483916e-10,
		      -1.3779079306969203e-25, -1.3937080366629122e-26, -1.3937080366629122e-26,
		      -2.7009598038184710e-10,  2.4188048187977241e-26,  1.9294757232832083e-26,
		      6.1412104954585808e-26, -1.9735547477471852e-25, -4.5422335402692920e-27,
		      -4.5422335402692877e-27, -5.4236554558686281e-26, -6.5590691655281548e-26,
		      -6.0701432153221186e-26, -6.0701432153221186e-26,  5.8584528742794487e-26};
  
  const double *cv = gkyl_array_cfetch(coll_iz_elc, gkyl_range_idx(&phaseRange, (int[5]) { 2, 2, 2, 5, 3}));

  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i], cv[i], 1e-12) );
  }
  
   gkyl_dg_iz_release(coll_iz);
}

void prim_vars_3x() { test_prim_vars_3x(false); }
void coll_iz() { test_coll_iz(false); }

#ifdef GKYL_HAVE_CUDA
void prim_vars_3x_gpu() { test_prim_vars_3x(true); }
void coll_iz_gpu() { test_coll_iz(true); }
#endif

TEST_LIST = {
  { "prim_vars_3x", prim_vars_3x },
  { "coll_iz", coll_iz },
#ifdef GKYL_HAVE_CUDA
  { "prim_vars_3x_gpu", prim_vars_3x_gpu },
  { "coll_iz_gpu", coll_iz_gpu },
#endif  
  { NULL, NULL },
};
