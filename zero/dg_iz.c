#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_transform_vlasov_gk.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>

// Functions to extract ADAS data and interpolate

/* Calculate ionization rate coefficient.
 * If ne is out of bounds, use value at nearest ne given
 * If Te>Te_max, use value at highest Te
 * If Te<Te_min, set coefficient to 0
 * Inputs: atomic_data data - struct containing adas data
 *         double* te - array of temperatures to evaluate coefficient at
 *         double* ne - array of densities to evaluate coefficient at
 *         int charge_state - atomic charge state the ion started at
 *         int number - number of (te,ne) points to evaluate
 * 
 */
void calc_ionization_coef(atomic_data data, struct gkyl_array* Teq, struct gkyl_array* M0q, int charge_state, int number, struct gkyl_array* ioniz_data){
  int charge_state_index, *tooSmallTe;
  double* result, *localte, *localne;
  double minLogTe, maxLogTe, minLogNe, maxLogNe;

  //tooSmallTe = (int*)malloc(sizeof(int)*number);
  minLogTe = data.logT[0];
  minLogNe = data.logNe[0];
  maxLogTe = data.logT[data.te_intervals-1];
  maxLogNe = data.logNe[data.ne_intervals-1];

  for(int i=0;i<number;i++){
    double *M0q_d = gkyl_array_fetch(M0q, i);
    double *Teq_d = gkyl_array_fetch(Teq, i);
    //tooSmallTe[i]=1.;
    if(Teq_d[0] > maxLogTe){
      Teq_d[0] = maxLogTe;
    }if(M0q_d[0] < minLogNe){
      M0q_d[0] = minLogNe;
    }else if(M0q_d[0] > maxLogNe){
      M0q_d[0] = maxLogNe;
    }
  }
  
  charge_state_index = data.te_intervals * data.ne_intervals * charge_state;
  bilinear_interp(data.logT,data.logNe,data.logData+charge_state_index,Teq,M0q,
    data.te_intervals,data.ne_intervals,number,ioniz_data);

  for(int i=0;i<number;i++){
    //Multiplies by 0 if Te<Te_min
    double *Teq_d = gkyl_array_fetch(Teq, i);
    double *ioniz_data_d = gkyl_array_fetch(ioniz_data, i);
    if(Teq_d[0] < minLogTe) ioniz_data_d[0] = 0.0;
  }
}

struct gkyl_dg_iz*
gkyl_dg_iz_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_iz_cu_dev_new(grid, cbasis, pbasis, conf_rng, phase_rng, elem_charge, mass_elc, type_ion);
  } 
#endif
  gkyl_dg_iz *up = gkyl_malloc(sizeof(struct gkyl_dg_iz));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cdim = cdim;
  up->use_gpu = use_gpu;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng;
  up->grid = grid;

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;

  // Interpolate ADAS data
  // Resolution, density and temp should be set in input file 
  atomic_data h_ion;
  h_ion = loadadf11("kernels/neutral/scd12_h.dat");

  int resM0=50, resTe=100, qpoints=resM0*resTe;
  double minM0 = 1e15, maxM0 = 1e20;
  double minTe = 1.0, maxTe = 4e3;

  up->minLogM0 = log10(minM0);
  up->minLogTe = log10(minTe);    
  up->maxLogM0 = log10(maxM0);
  up->maxLogTe = log10(maxTe);
  
  double dlogM0 = (up->maxLogM0 - up->minLogM0)/(resM0-1);
  double dlogTe = (up->maxLogTe - up->minLogTe)/(resTe-1);
  up->resM0=resM0;
  up->resTe=resTe;
  up->dlogM0=dlogM0;
  up->dlogTe=dlogTe;

  up->M0q = gkyl_array_new(GKYL_DOUBLE, 1, qpoints);
  up->Teq = gkyl_array_new(GKYL_DOUBLE, 1, qpoints);
  up->ioniz_data = gkyl_array_new(GKYL_DOUBLE, 1, qpoints);

  for(int i=0;i<resM0;i++){
    for(int j=0;j<resTe;j++){
      double *M0q_d = gkyl_array_fetch(up->M0q, i*resTe+j);
      double *Teq_d = gkyl_array_fetch(up->Teq, i*resTe+j);
      M0q_d[0] = up->minLogM0 + dlogM0*i;
      Teq_d[0] = up->minLogTe + dlogTe*j;     
    }
  }

  up->E = 13.6; // need to add all E_iz data for different species and charge states

  calc_ionization_coef(h_ion, up->Teq, up->M0q, h_ion.Z-1, qpoints, up->ioniz_data);

  for(int i=0;i<qpoints;i++){
    double *iz_dat_d = gkyl_array_fetch(up->ioniz_data, i);
    double *M0q_d = gkyl_array_fetch(up->M0q, i);
    double *Teq_d = gkyl_array_fetch(up->Teq, i);
    /* printf("Te = %eeV, M0 = %em^-3, Ioniz. coef = %em^3/s\n", */
    /* 	   pow(10,Teq_d[0]),pow(10,M0q_d[0]),iz_dat_d[0]); */
  }
  
  // allocate fields for prim mom calculation
  up->prim_vars_neut = gkyl_array_new(GKYL_DOUBLE, 2*cbasis->num_basis, up->conf_rng->volume); // elc, ion
  up->vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume); // all
  up->vtSq_iz = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);  // elc
  up->prim_vars_fmax = gkyl_array_new(GKYL_DOUBLE, 2*cbasis->num_basis, up->conf_rng->volume);  //elc
  up->coef_iz = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);  // all 
  up->fmax_iz = gkyl_array_new(GKYL_DOUBLE, pbasis->num_basis, up->phase_rng->volume); // elc, ion

  up->calc_prim_vars_neut_gk = gkyl_dg_prim_vars_transform_vlasov_gk_new(cbasis, pbasis, up->conf_rng, "prim", use_gpu); // elc, ion
  up->calc_prim_vars_elc_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(cbasis, pbasis, "vtSq", use_gpu); // all 
  
  up->proj_max = gkyl_proj_maxwellian_on_basis_new(grid, cbasis, pbasis, poly_order+1, use_gpu); // elc, ion

  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

void gkyl_dg_iz_coll_elc(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coll_iz)) {
    return gkyl_dg_iz_coll_elc_cu(up, moms_elc, moms_neut, bmag, jacob_tot, b_i, f_self, coll_iz, cflrate);
  } 
#endif
  // Set auxiliary variable (b_i) for computation of upar
  gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_neut_gk, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  
  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *m0_elc_d = &moms_elc_d[0];
    const double *m0_neut_d = &moms_neut_d[0];

    double *vtSq_elc_d = gkyl_array_fetch(up->vtSq_elc, loc);
    double *prim_vars_neut_d = gkyl_array_fetch(up->prim_vars_neut, loc);
    double *coef_iz_d = gkyl_array_fetch(up->coef_iz, loc);

    up->calc_prim_vars_elc_vtSq->kernel(up->calc_prim_vars_elc_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);
    up->calc_prim_vars_neut_gk->kernel(up->calc_prim_vars_neut_gk, conf_iter.idx,
				       moms_neut_d, prim_vars_neut_d);

    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    if (log10(m0_elc_av) < up->minLogM0) m0_idx=0;
    else if (log10(m0_elc_av) > up->maxLogM0) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->minLogM0)/(up->dlogM0);
      double *M0q_1 = gkyl_array_fetch(up->M0q, m0_idx*up->resTe);
      double *M0q_2 = gkyl_array_fetch(up->M0q, (1+m0_idx)*up->resTe);
      diff1 = fabs(M0q_1[0]-log10(m0_elc_av));
      diff2 = fabs(M0q_2[0]-log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->minLogTe) t_idx=0;
    else if (log10(temp_elc_av) > up->maxLogTe) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->minLogTe)/(up->dlogTe);
      double *Teq_1 = gkyl_array_fetch(up->Teq, t_idx);
      double *Teq_2 = gkyl_array_fetch(up->Teq, t_idx+1); 
      diff1 = fabs(Teq_1[0]-log10(temp_elc_av));
      diff2 = fabs(Teq_2[0]-log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    double *iz_dat_d = gkyl_array_fetch(up->ioniz_data, q_idx);
    //printf("m0 %g Te %g iz_dat %g\n", m0_elc_av, temp_elc_av, iz_dat_d[0]);
    coef_iz_d[0] = iz_dat_d[0]/cell_av_fac;    
  }

  // Calculate vt_sq_iz 
  gkyl_array_copy_range(up->vtSq_iz, up->vtSq_elc, up->conf_rng);
  gkyl_array_scale_range(up->vtSq_iz, 1/2.0, up->conf_rng);
  gkyl_array_shiftc(up->vtSq_iz, -up->E*up->elem_charge/(3*up->mass_elc)*pow(sqrt(2),up->cdim), 0);

  // Set fmax moments
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->prim_vars_neut, 0, up->conf_rng);
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->vtSq_iz, up->cbasis->num_basis, up->conf_rng);

  // Proj maxwellian on basis
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_elc,
    up->prim_vars_fmax, bmag, jacob_tot, up->mass_elc, up->fmax_iz);

  // copy, scale and accumulate
  gkyl_array_set_range(coll_iz, 2.0, up->fmax_iz, up->phase_rng);
  gkyl_array_accumulate_range(coll_iz, -1.0, f_self, up->phase_rng);
  
  // weak multiply
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_neut, up->conf_rng);
  // coll_iz = n_n*coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
    up->conf_rng, up->phase_rng);
  
  // cfl calculation
  //struct gkyl_range vel_rng;
  /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
  /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
  /* // cfl associated with reaction is a *phase space* cfl */
  /* // Need to loop over velocity space for each configuration space cell */
  /* // to get total cfl rate in each phase space cell */
  /* while (gkyl_range_iter_next(&vel_iter)) { */
  /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
  /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
  /*   cflrate_d[0] += cflr; // frequencies are additive */
  /* } */

}

void gkyl_dg_iz_coll_ion(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coll_iz)) {
    return gkyl_dg_iz_coll_ion_cu(up, moms_elc, moms_neut, bmag, jacob_tot, b_i, coll_iz, cflrate);
  } 
#endif
  // Set auxiliary variable (b_i) for computation of upar
  gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_neut_gk, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  
  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *m0_elc_d = &moms_elc_d[0];
    const double *m0_neut_d = &moms_neut_d[0];

    double *vtSq_elc_d = gkyl_array_fetch(up->vtSq_elc, loc);
    double *prim_vars_neut_d = gkyl_array_fetch(up->prim_vars_neut, loc);
    double *coef_iz_d = gkyl_array_fetch(up->coef_iz, loc);

    up->calc_prim_vars_elc_vtSq->kernel(up->calc_prim_vars_elc_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);
    up->calc_prim_vars_neut_gk->kernel(up->calc_prim_vars_neut_gk, conf_iter.idx,
				       moms_neut_d, prim_vars_neut_d);

    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    if (log10(m0_elc_av) < up->minLogM0) m0_idx=0;
    else if (log10(m0_elc_av) > up->maxLogM0) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->minLogM0)/(up->dlogM0);
      double *M0q_1 = gkyl_array_fetch(up->M0q, m0_idx*up->resTe);
      double *M0q_2 = gkyl_array_fetch(up->M0q, (1+m0_idx)*up->resTe);
      diff1 = fabs(M0q_1[0]-log10(m0_elc_av));
      diff2 = fabs(M0q_2[0]-log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->minLogTe) t_idx=0;
    else if (log10(temp_elc_av) > up->maxLogTe) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->minLogTe)/(up->dlogTe);
      double *Teq_1 = gkyl_array_fetch(up->Teq, t_idx);
      double *Teq_2 = gkyl_array_fetch(up->Teq, t_idx+1); 
      diff1 = fabs(Teq_1[0]-log10(temp_elc_av));
      diff2 = fabs(Teq_2[0]-log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    double *iz_dat_d = gkyl_array_fetch(up->ioniz_data, q_idx);
    //printf("m0 %g Te %g iz_dat %g\n", m0_elc_av, temp_elc_av, iz_dat_d[0]);
    coef_iz_d[0] = iz_dat_d[0]/cell_av_fac;    
  }

  // Proj maxwellian on basis
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_neut,
    up->prim_vars_neut, bmag, jacob_tot, up->mass_elc, coll_iz);

  // weak multiply
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_elc, up->conf_rng);
  // coll_iz = n_e*coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
    up->conf_rng, up->phase_rng);
  
  // cfl calculation
  //struct gkyl_range vel_rng;
  /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
  /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
  /* // cfl associated with reaction is a *phase space* cfl */
  /* // Need to loop over velocity space for each configuration space cell */
  /* // to get total cfl rate in each phase space cell */
  /* while (gkyl_range_iter_next(&vel_iter)) { */
  /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
  /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
  /*   cflrate_d[0] += cflr; // frequencies are additive */
  /* } */

}

void gkyl_dg_iz_coll_neut(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut, const struct gkyl_array *f_self,
  struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coll_iz)) {
    return gkyl_dg_iz_coll_neut_cu(up, moms_elc, moms_neut, f_self, coll_iz, cflrate);
  }
#endif
  
  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *m0_elc_d = &moms_elc_d[0];
    const double *m0_neut_d = &moms_neut_d[0];

    double *vtSq_elc_d = gkyl_array_fetch(up->vtSq_elc, loc);
    double *coef_iz_d = gkyl_array_fetch(up->coef_iz, loc);

    up->calc_prim_vars_elc_vtSq->kernel(up->calc_prim_vars_elc_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);

    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    if (log10(m0_elc_av) < up->minLogM0) m0_idx=0;
    else if (log10(m0_elc_av) > up->maxLogM0) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->minLogM0)/(up->dlogM0);
      double *M0q_1 = gkyl_array_fetch(up->M0q, m0_idx*up->resTe);
      double *M0q_2 = gkyl_array_fetch(up->M0q, (1+m0_idx)*up->resTe);
      diff1 = fabs(M0q_1[0]-log10(m0_elc_av));
      diff2 = fabs(M0q_2[0]-log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->minLogTe) t_idx=0;
    else if (log10(temp_elc_av) > up->maxLogTe) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->minLogTe)/(up->dlogTe);
      double *Teq_1 = gkyl_array_fetch(up->Teq, t_idx);
      double *Teq_2 = gkyl_array_fetch(up->Teq, t_idx+1);
      diff1 = fabs(Teq_1[0]-log10(temp_elc_av));
      diff2 = fabs(Teq_2[0]-log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    double *iz_dat_d = gkyl_array_fetch(up->ioniz_data, q_idx);
    //printf("m0 %g Te %g iz_dat %g\n", m0_elc_av, temp_elc_av, iz_dat_d[0]);
    coef_iz_d[0] = iz_dat_d[0]/cell_av_fac;
  }

  // neut coll_iz = -f_n
  gkyl_array_set_range(coll_iz, -1.0, f_self, up->phase_rng);

  // weak multiply
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_elc, up->conf_rng);
  // coll_iz = n_e*coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
    up->conf_rng, up->phase_rng);
  
  // cfl calculation
  //struct gkyl_range vel_rng;
  /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
  /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
  /* // cfl associated with reaction is a *phase space* cfl */
  /* // Need to loop over velocity space for each configuration space cell */
  /* // to get total cfl rate in each phase space cell */
  /* while (gkyl_range_iter_next(&vel_iter)) { */
  /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
  /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
  /*   cflrate_d[0] += cflr; // frequencies are additive */
  /* } */

}

void
gkyl_dg_iz_release(struct gkyl_dg_iz* up)
{
  gkyl_array_release(up->M0q);
  gkyl_array_release(up->Teq);
  gkyl_array_release(up->ioniz_data);
  gkyl_array_release(up->prim_vars_neut);
  gkyl_array_release(up->vtSq_elc);
  gkyl_array_release(up->vtSq_iz);
  gkyl_array_release(up->prim_vars_fmax);
  gkyl_array_release(up->coef_iz);
  gkyl_array_release(up->fmax_iz);
  gkyl_dg_prim_vars_type_release(up->calc_prim_vars_neut_gk);
  gkyl_dg_prim_vars_type_release(up->calc_prim_vars_elc_vtSq);
  gkyl_proj_maxwellian_on_basis_release(up->proj_max);
  free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_iz*
gkyl_dg_iz_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion)
{
  assert(false);
  return 0;
}

#endif
