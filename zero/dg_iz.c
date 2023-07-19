#include <assert.h>
#include <stdio.h>

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
double* calc_ionization_coef(atomic_data data, double* te, double* ne, int charge_state, int number){
  int charge_state_index, *tooSmallTe;
  double* result, *localte, *localne;
  double minLogTe, maxLogTe, minLogNe, maxLogNe;
  
  localte = (double*)malloc(sizeof(double)*number);
  localne = (double*)malloc(sizeof(double)*number);
  tooSmallTe = (int*)malloc(sizeof(int)*number);
  minLogTe = data.logT[0];
  minLogNe = data.logNe[0];
  maxLogTe = data.logT[data.te_intervals-1];
  maxLogNe = data.logNe[data.ne_intervals-1];

  for(int i=0;i<number;i++){
    localte[i]=te[i];
    localne[i]=ne[i];
    tooSmallTe[i]=1.;
    if(localte[i]<minLogTe){
      tooSmallTe[i]=0;
    }else if(localte[i]>maxLogTe){
      localte[i]=maxLogTe;
    }if(localne[i]<minLogNe){
      localne[i]=minLogNe;
    }else if(localne[i]>maxLogNe){
      localne[i]=maxLogNe;
    }
  }
  charge_state_index = data.te_intervals * data.ne_intervals * charge_state;
  result = bilinear_interp(data.logT,data.logNe,data.logData+charge_state_index,localte,localne,
			 data.te_intervals,data.ne_intervals,number);
  for(int i=0;i<number;i++){
    //Multiplies by 0 if Te<Te_min, 1 otherwise
    result[i]=tooSmallTe[i]*pow(10.0,result[i]);
  }
  return result; 
}

struct gkyl_dg_iz*
gkyl_dg_iz_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool is_gk, bool use_gpu)
{
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

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;

  // Interpolate ADAS data
  // Resolution, density and temp should be set in input file 
  atomic_data h_ion;
  h_ion = loadadf11("scd12_h.dat");

  double *M0q, *Teq, *ioniz_data;
  double minM0, minTe, maxM0, maxTe;
  int resM0=50, resTe=100, qpoints=resM0*resTe;
  minM0 = 1e15;
  maxM0 = 1e20;
  minTe = 1.0;
  maxTe = 4e3;
  double dlogM0 = (log10(maxM0) - log10(minM0))/(resM0-1);
  double dlogTe = (log10(maxTe) - log10(minTe))/(resTe-1);
  up->resM0=resM0;
  up->resTe=resTe;
  up->dlogM0=dlogM0;
  up->dlogTe=dlogTe;

  M0q = (double*) malloc(sizeof(double)*qpoints);
  Teq = (double*) malloc(sizeof(double)*qpoints);

  for(int i=0;i<resM0;i++){
    for(int j=0;j<resTe;j++){
      M0q[i*resTe+j] = log10(minM0) + dlogM0*i;
      Teq[i*resTe+j] = log10(minTe) + dlogTe*j;
    }
  }

  up->E = 13.6;
  up->M0q=M0q;
  up->Teq=Teq;

  up->ioniz_data = calc_ionization_coef(h_ion, up->Teq, up->M0q, h_ion.Z-1, qpoints);

  // Establish vdim vlasov species
  int vdim_vl; 
  int vdim = pdim - up->cdim;
  if (is_gk) {
    if (vdim == 1) {
      vdim_vl = vdim;
    }
    else {
      vdim_vl = vdim+1;
    }
  }
  else  {
    vdim_vl = vdim;
  }
  up->vdim_vl = vdim_vl;

  // allocate fields for prim mom calculation
  // do array declarations require anything for gpus??
  up->upar_neut = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->vtSq_iz = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->prim_vars_fmax = gkyl_array_new(GKYL_DOUBLE, 2*cbasis->num_basis, up->conf_rng->volume);
  up->coef_iz = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->fmax_iz = gkyl_array_new(GKYL_DOUBLE, pbasis->num_basis, up->phase_rng->volume);

  up->calc_prim_vars_neut_upar = gkyl_dg_prim_vars_transform_vlasov_gk_new(cbasis, pbasis, up->conf_rng, "u_par", use_gpu);
  up->calc_prim_vars_elc_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(cbasis, pbasis, "vtSq", use_gpu);

  up->proj_max = gkyl_proj_maxwellian_on_basis_new(grid, cbasis, pbasis, poly_order+1, use_gpu);

  return up;
}

void gkyl_dg_iz_coll(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
  // Set auxiliary variable (b_i) for computation of upar
  gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_neut_upar, 
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
    double *upar_neut_d = gkyl_array_fetch(up->upar_neut, loc);
    double *coef_iz_d = gkyl_array_fetch(up->coef_iz, loc);

    up->calc_prim_vars_elc_vtSq->kernel(up->calc_prim_vars_elc_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);
    up->calc_prim_vars_neut_upar->kernel(up->calc_prim_vars_neut_upar, conf_iter.idx, moms_neut_d, upar_neut_d);

    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    // First consider density in flattened array
    if (log10(m0_elc_av) < up->M0q[0]) m0_idx=0;
    else if (log10(m0_elc_av) > up->M0q[qtotal-1]) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->M0q[0])/(up->dlogM0);
      diff1 = fabs(up->M0q[m0_idx*up->resTe]-log10(m0_elc_av));
      diff2 = fabs(up->M0q[(m0_idx+1)*up->resTe] - log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->Teq[0]) t_idx=0;
    else if (log10(temp_elc_av) > up->Teq[up->resTe-1]) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->Teq[0])/(up->dlogTe);
      diff1 = fabs(up->Teq[t_idx]- log10(temp_elc_av));
      diff2 = fabs(up->Teq[(t_idx+1)] - log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    coef_iz_d[0] = up->ioniz_data[q_idx]/cell_av_fac;
  }

  // Calculate vt_sq_iz 
  gkyl_array_copy_range(up->vtSq_iz, up->vtSq_elc, *up->conf_rng);
  gkyl_array_scale_range(up->vtSq_iz, 1/2.0, *up->conf_rng);
  gkyl_array_shiftc(up->vtSq_iz, -up->E*up->elem_charge/(3*up->mass_elc)*pow(sqrt(2),up->cdim), 0);

  // Set fmax moments
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->upar_neut, 0, *up->conf_rng);
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->vtSq_iz, up->cbasis->num_basis, *up->conf_rng);

  // POM
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_elc,
    up->prim_vars_fmax, bmag, jacob_tot, up->mass_elc, up->fmax_iz);

  // copy, scale and accumulate
  gkyl_array_set_range(coll_iz, 2.0, up->fmax_iz, *up->phase_rng);
  gkyl_array_accumulate_range(coll_iz, -1.0, f_self, *up->phase_rng);
  
  // weak multiply
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_neut, up->conf_rng);
  // coll_iz = coef_iz*coll_iz
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
gkyl_dg_iz_release(gkyl_dg_iz* iz)
{
  free(iz);
}
