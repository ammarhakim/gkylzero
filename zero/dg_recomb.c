#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#include <gkyl_dg_recomb.h>
#include <gkyl_dg_recomb_priv.h>
#include <gkyl_util.h>
#include <gkyl_read_adas_priv.h>

// FIX MASS FOR PMOB ETC
struct gkyl_dg_recomb*
gkyl_dg_recomb_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge,
  double mass_elc, double mass_self, enum gkyl_dg_recomb_type type_ion, int charge_state,
  enum gkyl_dg_recomb_self type_self, bool all_gk, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_recomb_cu_dev_new(grid, cbasis, pbasis, conf_rng, phase_rng, elem_charge, mass_elc, type_ion, charge_state);
  } 
#endif
  gkyl_dg_recomb *up = gkyl_malloc(sizeof(struct gkyl_dg_recomb));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  
  int poly_order = cbasis->poly_order;
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cdim = cdim;
  up->vdim = vdim;
  up->use_gpu = use_gpu;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng;
  up->grid = grid;
  up->all_gk = all_gk;

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;
  up->mass_self = mass_self;
  
  up->type_self = type_self;
  
  // Project ADAS data (H, He, Li)
  struct adas_field data;

  // Conditional to determine element
  if (type_ion == GKYL_RECOMB_H) {
    data.NT = 29;
    data.NN = 24;
    data.logData = fopen("adas-dat/recomb_h.npy", "rb");
    data.logT = fopen("adas-dat/logT_h.npy", "rb");
    data.logN = fopen("adas-dat/logN_h.npy", "rb");
    data.Zmax = 1;
  }
  else if (type_ion == GKYL_RECOMB_HE) {
    data.NT = 30;
    data.NN = 24;
    data.logData = fopen("adas-dat/recomb_he.npy", "rb");
    data.logT = fopen("adas-dat/logT_he.npy", "rb");
    data.logN = fopen("adas-dat/logN_he.npy", "rb");
    data.Zmax = 2;
  }
  else if (type_ion == GKYL_RECOMB_LI) {
    data.NT = 25;
    data.NN = 24;
    data.logData = fopen("adas-dat/recomb_li.npy", "rb");
    data.logT = fopen("adas-dat/logT_li.npy", "rb");
    data.logN = fopen("adas-dat/logN_li.npy", "rb");
    data.Zmax = 3;
  }
  else fprintf(stderr, "Incorrect ion type for recombination.");
  
  long sz = data.NT*data.NN;
  double *minmax;

  if (data.logT == NULL) fprintf(stderr, "Unable to load ADAS 'logT_<elem>.npy' file.");
  if (data.logN == NULL) fprintf(stderr, "Unable to load ADAS 'logN_<elem>.npy' file.");
  if (data.logData == NULL) fprintf(stderr, "Unable to load ADAS 'recomb_<elem>.npy' file.");

  // "duplicate symbol" error
  minmax = minmax_from_numpy(data.logT, data.NT);
  fclose(data.logT);
  double logTmin = minmax[0], logTmax = minmax[1];
  minmax = minmax_from_numpy(data.logN, data.NN);
  fclose(data.logN);
  double logNmin = minmax[0]+6., logNmax = minmax[1]+6.; //adjust for 1/cm^3 to 1/m^3 conversion

  // "duplicate symbol" error
  struct gkyl_array *adas_nodal = array_from_numpy(data.logData, sz, data.Zmax);
  fclose(data.logData);

  if (!adas_nodal) {
    fprintf(stderr, "Unable to read data from adas nodal numpy file!\n");
    return 0;
  }

  struct gkyl_range range_node;
  gkyl_range_init_from_shape(&range_node, 2, (int[]) { data.NT, data.NN } );
  
  // allocate grid and DG array
  struct gkyl_rect_grid tn_grid;
  gkyl_rect_grid_init(&tn_grid, 2,
    (double[]) { logTmin, logNmin},
    (double []) { logTmax, logNmax},
    (int[]) { data.NT-1, data.NN-1 }
  );

  struct gkyl_range adas_rng;
  //int ghost[] = { 0, 0 };
  //gkyl_create_grid_ranges(&tn_grid, ghost, &adas_rng_ext, &adas_rng);
  gkyl_range_init_from_shape(&adas_rng, 2, tn_grid.cells);
  
  struct gkyl_basis adas_basis;
  gkyl_cart_modal_serendip(&adas_basis, 2, 1);
  
  struct gkyl_array *adas_dg =
    gkyl_array_new(GKYL_DOUBLE, adas_basis.num_basis, data.NT*data.NN);

  // "duplicate symbol" error
  create_dg_from_nodal(&tn_grid, &range_node, adas_nodal, adas_dg, charge_state);
  //gkyl_grid_sub_array_write(&tn_grid, &adas_rng, adas_dg, "adas_dg.gkyl");

  // ADAS data pointers
  up->recomb_data = adas_dg;
  up->minLogM0 = logNmin;
  up->minLogTe = logTmin;
  up->maxLogM0 = logNmax;
  up->maxLogTe = logTmax;
  up->dlogTe = tn_grid.dx[0];
  up->dlogM0 = tn_grid.dx[1];
  up->resTe = tn_grid.cells[0];
  up->resM0 = tn_grid.cells[1];
  up->adas_rng = adas_rng;
  up->adas_basis = adas_basis;
  
  // allocate fields for prim mom calculation
  up->vtSq_elc = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume); // all
  up->coef_recomb = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);  // all
  up->calc_prim_vars_elc_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(cbasis, pbasis, "vtSq", use_gpu); // all

  if (type_self == GKYL_RECOMB_RECVR) {
    if (up->all_gk==false) {
      up->calc_prim_vars_ion_udrift = gkyl_dg_prim_vars_transform_vlasov_gk_new(cbasis, pbasis, up->conf_rng, "u_par_i", use_gpu);
      up->calc_prim_vars_ion_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(cbasis, pbasis, "vtSq", use_gpu);

      up->udrift_ion = gkyl_array_new(GKYL_DOUBLE, vdim*cbasis->num_basis, up->conf_rng->volume);
      up->vtSq_ion = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
      up->prim_vars_ion = gkyl_array_new(GKYL_DOUBLE, (vdim+1)*cbasis->num_basis, up->conf_rng->volume);
    }
    
    up->proj_max = gkyl_proj_maxwellian_on_basis_new(grid, cbasis, pbasis, poly_order+1, use_gpu);
  }
  
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

void gkyl_dg_recomb_coll(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *b_i, const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot,
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coll_recomb)) {
    return gkyl_dg_recomb_coll_cu(up, moms_elc, moms_ions, b_i, f_self, coll_recomb, cflrate);
  } 
#endif
  if ((up->all_gk == false) && (up->type_self == GKYL_RECOMB_RECVR)) {
    // Set auxiliary variable (b_i) for computation of udrift_i
    gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_ion_udrift, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  }

  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;
  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);
    const double *moms_elc_d = gkyl_array_cfetch(moms_elc, loc);
    const double *m0_elc_d = &moms_elc_d[0];

    double *vtSq_elc_d = gkyl_array_fetch(up->vtSq_elc, loc);
    double *coef_recomb_d = gkyl_array_fetch(up->coef_recomb, loc);
    up->calc_prim_vars_elc_vtSq->kernel(up->calc_prim_vars_elc_vtSq, conf_iter.idx, moms_elc_d, vtSq_elc_d);

    //Find cell containing value of n,T
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    double cell_val_t;
    double cell_val_m0;
    int m0_idx, t_idx;
    double cell_vals_2d[2];
    double cell_center; 

    if (log_Te_av < up->minLogTe) t_idx=1;
    else if (log_Te_av > up->maxLogTe) t_idx=up->resTe;
    else t_idx = (log_Te_av - up->minLogTe)/(up->dlogTe)+1;
    cell_center = (t_idx - 0.5)*up->dlogTe + up->minLogTe;
    cell_vals_2d[0] = 2.0*(log_Te_av - cell_center)/up->dlogTe; // Te value on cell interval
      
    if (log_m0_av < up->minLogM0) m0_idx=1;
    else if (log_m0_av > up->maxLogM0) m0_idx=up->resM0;
    else m0_idx = (log_m0_av - up->minLogM0)/(up->dlogM0)+1;
    cell_center = (m0_idx - 0.5)*up->dlogM0 + up->minLogM0;
    cell_vals_2d[1] = 2.0*(log_m0_av - cell_center)/up->dlogM0; // M0 value on cell interval
    
    double *recomb_dat_d = gkyl_array_fetch(up->recomb_data, gkyl_range_idx(&up->adas_rng, (int[2]) {t_idx,m0_idx}));
    double adas_eval = up->adas_basis.eval_expand(cell_vals_2d, recomb_dat_d);
    coef_recomb_d[0] = pow(10.0,adas_eval)/cell_av_fac;

    if ((up->all_gk==false) && (up->type_self == GKYL_RECOMB_RECVR)) {
      const double *moms_ion_d = gkyl_array_cfetch(moms_ion, loc);
      double *udrift_ion_d = gkyl_array_fetch(up->udrift_ion, loc);
      double *vtSq_ion_d = gkyl_array_fetch(up->vtSq_ion, loc);
      double *prim_vars_ion_d = gkyl_array_fetch(up->prim_vars_ion, loc);
      
      // condense the following 2 kernels...
      up->calc_prim_vars_ion_udrift->kernel(up->calc_prim_vars_ion_udrift, conf_iter.idx,
					    moms_ion_d, udrift_ion_d);
      up->calc_prim_vars_ion_vtSq->kernel(up->calc_prim_vars_ion_vtSq, conf_iter.idx,
					    moms_ion_d, vtSq_ion_d);
    }
  }

  if (up->type_self == GKYL_RECOMB_RECVR) {
    if (up->all_gk) {
      gkyl_proj_gkmaxwellian_on_basis_lab_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_ion, bmag, jacob_tot, up->mass_self, coll_recomb);
    }
    else {
      // Set fmax moments
      gkyl_array_set_offset_range(up->prim_vars_ion, 1., up->udrift_ion, 0, *up->conf_rng);
      gkyl_array_set_offset_range(up->prim_vars_ion, 1., up->vtSq_ion, (up->vdim)*up->cbasis->num_basis, *up->conf_rng);
      
      // Proj maxwellian on basis
      gkyl_proj_maxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_ion,
					     up->prim_vars_ion, coll_recomb);
    }

  }
  else {
    // copy and weak mult
    gkyl_array_set_range(coll_recomb, -1.0, f_self, *up->phase_rng);
  }

  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_recomb, 0, up->coef_recomb, 0, moms_elc, up->conf_rng);
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_recomb, up->coef_recomb, coll_recomb,
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
gkyl_dg_recomb_release(struct gkyl_dg_recomb* up)
{
  gkyl_array_release(up->recomb_data);
  gkyl_array_release(up->vtSq_elc);
  gkyl_array_release(up->coef_recomb);
  gkyl_dg_prim_vars_type_release(up->calc_prim_vars_elc_vtSq);

  // only used for Vlasov neut coll.
  gkyl_array_release(up->udrift_ion);
  gkyl_array_release(up->vtSq_ion);
  gkyl_array_release(up->prim_vars_ion);
  gkyl_proj_maxwellian_on_basis_release(up->proj_max);
  gkyl_dg_prim_vars_type_release(up->calc_prim_vars_ion_udrift);
  gkyl_dg_prim_vars_type_release(up->calc_prim_vars_ion_vtSq);
  
  free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_recomb*
gkyl_dg_recomb_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge, double mass_elc,
  double mass_self, enum gkyl_dg_recomb_type type_ion, int charge_state, enum gkyl_dg_recomb_self type_self, bool all_gk)
{
  assert(false);
  return 0;
}

#endif
