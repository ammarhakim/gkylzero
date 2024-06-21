#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis_priv.h>
#include <gkyl_range.h>
#include <assert.h>

// create range to loop over quadrature points.
static inline struct gkyl_range
get_qrange(int cdim, int dim, int num_quad, int num_quad_v, bool *is_vdim_p2)
{
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = is_vdim_p2[i-cdim]? num_quad_v : num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(int cdim, const struct gkyl_basis *basis, int num_quad, struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  int vdim = ndim-cdim;
  int num_quad_v = num_quad;
  // hybrid basis have p=2 in velocity space.
  // gkhybrid basis have p=2 in vpar only.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if ((basis->b_type == GKYL_BASIS_MODAL_HYBRID) ||
      (basis->b_type == GKYL_BASIS_MODAL_GKHYBRID)) {
    num_quad_v = num_quad+1;
    is_vdim_p2[0] = true;  // for gkhybrid.
    if (basis->b_type == GKYL_BASIS_MODAL_HYBRID)
      for (int d=0; d<vdim; d++) is_vdim_p2[d] = true;
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];
  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }
  if (num_quad_v <= gkyl_gauss_max) {
    memcpy(ordinates1_v, gkyl_gauss_ordinates[num_quad_v], sizeof(double[num_quad_v]));
    memcpy(weights1_v, gkyl_gauss_weights[num_quad_v], sizeof(double[num_quad_v]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1_v, weights1_v, num_quad_v);
  }

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad, num_quad_v, is_vdim_p2);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *ordinates = gkyl_array_cu_dev_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  }
  else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<cdim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i) {
      ord[i] = is_vdim_p2[i-cdim]? ordinates1_v[iter.idx[i]-qrange.lower[i]] :
        ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    //printf("x=%10.8e, vpar=%10.8e, vperp=%10.8e\n", ord[0], ord[1], ord[2]);

    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i)
      wgt[0] *= is_vdim_p2[i-cdim]? weights1_v[iter.idx[i]-qrange.lower[i]] :
        weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) 
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

gkyl_proj_maxwellian_on_basis*
gkyl_proj_maxwellian_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, const struct gkyl_velocity_map *vel_map, bool use_gpu)
{
  return gkyl_proj_maxwellian_on_basis_inew( &(struct gkyl_proj_maxwellian_on_basis_inp) {
        .grid = grid,
        .conf_basis = conf_basis,
        .phase_basis = phase_basis,
        .num_quad = num_quad,
        .vel_map = vel_map,
        .use_gpu = use_gpu,
      }
    ); 
}

gkyl_proj_maxwellian_on_basis*
gkyl_proj_maxwellian_on_basis_inew(const struct gkyl_proj_maxwellian_on_basis_inp *inp)
{
  gkyl_proj_maxwellian_on_basis *up = gkyl_malloc(sizeof(*up));

  const struct gkyl_rect_grid *grid = inp->grid;
  const struct gkyl_basis *conf_basis = inp->conf_basis;
  const struct gkyl_basis *phase_basis = inp->phase_basis;
  int num_quad = inp->num_quad == 0 ? phase_basis->poly_order+1 : inp->num_quad;

  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->use_gpu = inp->use_gpu;
  up->phase_basis_on_dev = inp->phase_basis_on_dev;
  up->conf_basis_on_dev = inp->conf_basis_on_dev;
  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);

  // MF 2022/08/09: device kernel has arrays hard-coded to 3x, vdim=3, p=2 for now.
  if (up->use_gpu) assert((up->cdim<3 && conf_basis->poly_order<4) || (up->cdim==3 && conf_basis->poly_order<3));

  // initialize data needed for conf-space quadrature 
  up->tot_conf_quad = init_quad_values(up->cdim, conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, false);

  // initialize data needed for phase-space quadrature 
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, false);

  up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad); // Only used in CPU implementation.

  int vdim = up->pdim-up->cdim;
  // To avoid creating iterators over ranges in device kernel, we'll
  // create a map between phase-space and conf-space ordinates.
  int num_quad_v = num_quad;  // Hybrid basis have p=2 in velocity space.
  // hybrid basis have p=2 in velocity space.
  // gkhybrid basis have p=2 in vpar only.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if ((phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) ||
      (phase_basis->b_type == GKYL_BASIS_MODAL_GKHYBRID)) {
    num_quad_v = num_quad+1;
    is_vdim_p2[0] = true;  // for gkhybrid.
    if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID)
      for (int d=0; d<vdim; d++) is_vdim_p2[d] = true;
  }
  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate device copies of arrays needed for quadrature.
    int p2c_qidx_ho[up->phase_qrange.volume];
    up->p2c_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);
    
    up->fm_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_quad, inp->phase_range_ext->volume); // D.L. added 06/06/2024.
    up->den_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 
    up->upar_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 
    up->vtsq_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 
    up->bmag_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 
    up->jactot_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 
    up->expamp_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_conf_quad, inp->conf_range_ext->volume); 

    // Allocate the memory for computing the specific phase nodal to modal calculation
    up->phase_nodal_to_modal_mem = gkyl_cu_mat_mm_array_mem_cu_dev_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0,
      GKYL_NO_TRANS, GKYL_NO_TRANS);
    // Compute the matrix A for the phase nodal to modal memory
    const double *phase_w = (const double*) up->weights->data;
    const double *phaseb_o = (const double*) up->basis_at_ords->data;
    for (int n=0; n<up->tot_quad; ++n){
      for (int k=0; k<up->num_phase_basis; ++k){
        gkyl_mat_set(up->phase_nodal_to_modal_mem->A_ho, k, n, phase_w[n]*phaseb_o[k+up->num_phase_basis*n]);
        double actual = gkyl_mat_get(up->phase_nodal_to_modal_mem->A_ho, k, n);
      }
    }
    // copy to device
    gkyl_mat_copy(up->phase_nodal_to_modal_mem->A_cu, up->phase_nodal_to_modal_mem->A_ho);

    // Create a cuda handle for all cublas operations
    up->cuh = 0;
    cublasCreate_v2(&up->cuh);

    // initialize data needed for conf-space quadrature 
    up->tot_conf_quad = init_quad_values(up->cdim, conf_basis, num_quad,
      &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, up->use_gpu);

    // initialize data needed for phase-space quadrature 
    up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
      &up->ordinates, &up->weights, &up->basis_at_ords, up->use_gpu);

    int pidx[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&up->conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
  }
#endif

  return up;
}

static void
proj_on_basis(const gkyl_proj_maxwellian_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

void
gkyl_proj_gkmaxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *prim_moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu(up, phase_rng, conf_rng, prim_moms,
      bmag, jacob_tot, mass, fmax);
#endif

  double fJacB_floor = 1.e-40;
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;  

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.};
  double expamp_o[tot_conf_quad], upar_o[tot_conf_quad], vtsq_o[tot_conf_quad];
  double bmag_o[tot_conf_quad];
  
  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *prim_moms_d = gkyl_array_cfetch(prim_moms, midx);
    const double *m0_d = prim_moms_d;
    const double *upar_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[2*num_conf_basis];
    const double *bmag_d = gkyl_array_cfetch(bmag, midx);
    const double *jactot_d = gkyl_array_cfetch(jacob_tot, midx);
    
    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      double m0_o = 0., jac_o = 0.;
      upar_o[n] = 0.;
      vtsq_o[n] = 0.;
      bmag_o[n] = 0.;
      for (int k=0; k<num_conf_basis; ++k) {
        m0_o += m0_d[k]*b_ord[k];
        upar_o[n] += upar_d[k]*b_ord[k];
        vtsq_o[n] += vtsq_d[k]*b_ord[k];
        bmag_o[n] += bmag_d[k]*b_ord[k];
        jac_o += jactot_d[k]*b_ord[k];
      }

      // Amplitude of the exponential.
      if ((m0_o > 0.) && (vtsq_o[n]>0.))
        expamp_o[n] = jac_o*m0_o/sqrt(pow(2.0*GKYL_PI*vtsq_o[n], vdim_phys));
      else
        expamp_o[n] = 0.;
    }

    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute Maxwellian at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        const double *xcomp_d = gkyl_array_cfetch(up->ordinates, pqidx);

        // Convert comp velocity coordinate to phys velocity coord.
        const struct gkyl_velocity_map *gvm = up->vel_map;
        long vlinidx = gkyl_range_idx(&gvm->local_ext_vel, vel_iter.idx);
        const double *vmap_d = gkyl_array_cfetch(gvm->vmap, vlinidx);
        double xcomp[1];
        for (int vd=0; vd<vdim; vd++) {
          xcomp[0] = xcomp_d[cdim+vd];
          xmu[cdim+vd] = gvm->vmap_basis->eval_expand(xcomp, vmap_d+vd*gvm->vmap_basis->num_basis);
        }

        double efact = 0.0;        
        // vpar term.
        efact += pow(xmu[cdim]-upar_o[cqidx],2);
        // mu term (only for 2v, vdim_phys=3).
        efact += (vdim_phys-1)*xmu[cdim+1]*bmag_o[cqidx]/mass;

        double *fq = gkyl_array_fetch(up->fun_at_ords, pqidx);
        fq[0] = vtsq_o[cqidx] > 0.0 ? fJacB_floor+expamp_o[cqidx]*exp(-efact/(2.0*vtsq_o[cqidx])) : fJacB_floor;
      }

      // compute expansion coefficients of Maxwellian on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(fmax, lidx));
    }
  }
}

void
gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->p2c_qidx);
    gkyl_array_release(up->fm_quad); // D.L. added 06/06/2024.
    gkyl_array_release(up->den_quad); 
    gkyl_array_release(up->upar_quad); 
    gkyl_array_release(up->vtsq_quad); 
    gkyl_array_release(up->bmag_quad); 
    gkyl_array_release(up->jactot_quad); 
    gkyl_array_release(up->expamp_quad); 
    gkyl_cu_mat_mm_array_mem_release(up->phase_nodal_to_modal_mem);
    cublasDestroy(up->cuh);
  }
#endif
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_array_release(up->fun_at_ords);
  gkyl_velocity_map_release(up->vel_map);
  gkyl_free(up);
}
