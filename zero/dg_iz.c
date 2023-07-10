#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
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

/* Calculate Recombination rate coefficient.
 * If ne is out of bounds, use value at nearest ne given
 * If Te>Te_max, set coefficient to 0
 * If Te<Te_min, use value at highest Te
 * Inputs: atomic_data data - struct containing adas data
 *         double* te - array of temperatures to evaluate coefficient at
 *         double* ne - array of densities to evaluate coefficient at
 *         int charge_state - atomic charge state the ion ends at
 *         int number - number of (te,ne) points to evaluate
 * 
 */
double* calc_recombination_coef(atomic_data data, double* te, double* ne, int charge_state, int number){
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
    tooSmallTe[i]=1;
    if(localte[i]<minLogTe){
      localte[i]=minLogTe;
    }else if(localte[i]>maxLogTe){
      tooSmallTe[i]=0;
    }
    //printf("%f %f% f\n",localne[i],minLogNe,maxLogNe);
    if(localne[i]<minLogNe){
      localne[i]=minLogNe;
    }else if(localne[i]>maxLogNe){
      localne[i]=maxLogNe;
    }
    //printf("%f %f% f\n",localne[i],minLogNe,maxLogNe);
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
gkyl_dg_iz_new(struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool is_gk, bool use_gpu)
{
  gkyl_dg_iz *up = gkyl_malloc(sizeof(struct gkyl_dg_iz));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->cdim = cdim;
  up->basis = cbasis;
  up->use_gpu = use_gpu;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng; 

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;

  // Interpolate ADAS data
  // Resolution, density and temp should be set in input file 
  atomic_data c_recomb, c_ion; 
  c_recomb = loadadf11("acd96_c.dat");
  c_ion = loadadf11("scd96_c.dat");

  double *M0q, *Teq, *recomb_data, *ioniz_data;
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

  up->recomb_data = calc_recombination_coef(c_recomb, up->Teq, up->M0q, c_recomb.Z-1, qpoints);
  up->ioniz_data = calc_ionization_coef(c_ion, up->Teq, up->M0q, c_recomb.Z-1, qpoints);
  for(int i=0;i<qpoints;i++){
    printf("Te = %eeV, M0 = %em^-3, Recomb. coef = %em^3/s, Ioniz. coef = %em^3/s\n",
  	   pow(10,Teq[i]),pow(10,M0q[i]),up->recomb_data[i],up->ioniz_data[i]);
  }

  // Establish vdim for gk and vlasov species
  // will be either 1 & 1 or 2 & 3
  // for gk and vlasov, respectively
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
  up->u_sq_temp = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->m2_temp = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_neut = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_elc = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);

  gkyl_dg_bin_op_mem *mem; // see binop unit test for allocating gpu mem
  up->mem = gkyl_dg_bin_op_mem_new(up->vth_sq_elc->size, up->basis->num_basis);
  
  return up;
}

void gkyl_dg_iz_react_rate(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut, struct gkyl_array *vth_sq_iz,
  struct gkyl_array *cflrate, struct gkyl_array *coef_iz)
{
  // gpu func calls will be the same for bin ops
  // Calculate vtsq from moms using bin ops
  // neutral vtsq
  for (int d=1; d<up->vdim_vl+1; d+=1) {
    gkyl_dg_mul_op_range(*up->basis, 0, up->u_sq_temp, d, moms_neut, d, moms_neut, up->conf_rng);
    gkyl_array_accumulate_range(up->m2_temp, -1.0, up->u_sq_temp, *up->conf_rng); // -M1.M1
  }
  gkyl_dg_div_op_range(up->mem, *up->basis, 0, up->m2_temp, 0, up->m2_temp, 0, moms_neut, up->conf_rng); // -u.u*m0
  gkyl_array_accumulate_offset_range(up->m2_temp, 1.0, moms_neut, (up->vdim_vl+1)*up->basis->num_basis, *up->conf_rng); // m2-u.u*m0
  gkyl_dg_div_op_range(up->mem, *up->basis, 0, up->vth_sq_neut, 0, up->m2_temp, 0, moms_neut, up->conf_rng); // (m2-u.u*m0)/m0
  gkyl_array_scale_range(up->vth_sq_neut, 1/up->vdim_vl, *up->conf_rng); // (m2-u.u*m0)/(vdim*m0)

  // elc vtsq
  gkyl_dg_mul_op_range(*up->basis, 0, up->m2_temp, 1, moms_elc, 1, moms_elc, up->conf_rng); // M1par*M1par
  gkyl_dg_mul_op_range(*up->basis, 0, up->m2_temp, 0, up->m2_temp, 0, moms_elc , up->conf_rng); // uparSq*m0
  gkyl_array_accumulate_offset_range(up->m2_temp, -1.0, moms_elc, (up->vdim_vl+1)*up->basis->num_basis, *up->conf_rng); // uparSq*m0 - m2
  gkyl_dg_div_op_range(up->mem, *up->basis, 0, up->vth_sq_elc, 0, up->m2_temp, 0, moms_elc, up->conf_rng); // (uparSq*m0 - m2)/m0
  gkyl_array_scale_range(up->vth_sq_elc, -1/up->vdim_vl, *up->conf_rng); // (m2 - uparSq*m0) / (vdim*m0)

  //for testing
  /* gkyl_array_accumulate_offset_range(up->vth_sq_neut, 1.0, moms_neut, (up->vdim_vl+1)*up->basis->num_basis, *up->conf_rng); */
  /* gkyl_array_accumulate_offset_range(up->vth_sq_elc, 1.0, moms_elc, (up->vdim_vl+1)*up->basis->num_basis, *up->conf_rng); */
    
  // Calculate vt_sq_iz
  gkyl_array_copy_range(vth_sq_iz, up->vth_sq_elc, *up->conf_rng);
  gkyl_array_scale_range(vth_sq_iz, 1/2.0, *up->conf_rng);
  gkyl_array_shiftc0(vth_sq_iz, -up->E*up->elem_charge/(3*up->mass_elc)*pow(sqrt(2),up->cdim));

  //struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);
   
    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    const double *m0_elc_d = &moms_elc_d[0]; 
    const double *vth_sq_neut_d = gkyl_array_cfetch(up->vth_sq_neut, loc);
    const double *vth_sq_elc_d = gkyl_array_cfetch(up->vth_sq_elc, loc);
    double *coef_iz_d = gkyl_array_fetch(coef_iz, loc);

    // Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac; 
    double temp_elc_av = vth_sq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
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
      printf("%f %f %f %f", diff1, diff2, pow(10.,up->Teq[t_idx]),  pow(10.,up->Teq[t_idx+1]));
      if (diff1 > diff2) t_idx += 1;
    }
    
    q_idx = m0_idx*up->resTe + t_idx;
    coef_iz_d[0] = up->ioniz_data[q_idx]/cell_av_fac; // not sure if this is correct

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
}

void
gkyl_dg_iz_release(gkyl_dg_iz* iz)
{
  free(iz);
}
