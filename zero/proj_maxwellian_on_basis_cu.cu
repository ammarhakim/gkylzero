/* -*- c++ -*- */

extern "C" {
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads and can parallelize over components too
  dimBlock->y = GKYL_MIN2(ncomp, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->y = gkyl_int_div_up(ncomp, dimBlock->y);
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
  printf("dimBlock->y=%d, dimGrid->y=%d, dimBlock->x=%d, dimGrid->x=%d\n", dimBlock->y, dimGrid->y, dimBlock->x,  dimGrid->x);
}

__global__ static void
gkyl_expand_prim_mom_at_quad_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r, 
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  struct gkyl_array* GKYL_RESTRICT den_quad_d,
  struct gkyl_array* GKYL_RESTRICT upar_quad_d,
  struct gkyl_array* GKYL_RESTRICT vtsq_quad_d,
  struct gkyl_array* GKYL_RESTRICT bmag_quad_d, 
  struct gkyl_array* GKYL_RESTRICT jactot_quad_d,
  struct gkyl_array* GKYL_RESTRICT expamp_quad_d,
  //struct gkyl_array* GKYL_RESTRICT fm_quad,
  const struct gkyl_array* GKYL_RESTRICT prim_moms,
  const struct gkyl_array* GKYL_RESTRICT bmag, 
  const struct gkyl_array* GKYL_RESTRICT jacob_tot,
  const struct gkyl_basis* conf_basis_on_dev)
{
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;

  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];
 
  // 2D thread grid.
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_r, tid, cidx);

    // Get conf-space linear index.
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *prim_moms_d = (const double *) gkyl_array_cfetch(prim_moms, lincC);
    const double *den_d = prim_moms_d;
    const double *upar_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[2*num_conf_basis];
    const double *bmag_d = (const double *) gkyl_array_cfetch(bmag, lincC);
    const double *jactot_d = (const double *) gkyl_array_cfetch(jacob_tot, lincC);

    double *den_quad = (double*) gkyl_array_fetch(den_quad_d, lincC);
    double *upar_quad = (double*) gkyl_array_fetch(upar_quad_d, lincC);
    double *vtsq_quad = (double*) gkyl_array_fetch(vtsq_quad_d, lincC);
    double *bmag_quad = (double*) gkyl_array_fetch(bmag_quad_d, lincC);  
    double *jactot_quad = (double*) gkyl_array_fetch(jactot_quad_d, lincC);  
    double *expamp_quad = (double*) gkyl_array_fetch(expamp_quad_d, lincC);  
    printf("den_quad=%d\n", &den_quad_d->size);

    // Sum over basis for given primitive moments.
    int n = linc2;
    //const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);

    // Compute primitive moments at quadrature nodes.
    for (int k=0; k<num_conf_basis; ++k) {
      printf("n=%d, den_quad=%d\n", n, &den_quad_d->size);
      conf_basis_on_dev->modal_to_quad_nodal(den_d, den_quad, n);
      printf("n=%d, upar_quad=%d\n", n, &upar_quad_d->size);
      conf_basis_on_dev->modal_to_quad_nodal(upar_d, upar_quad, n);
      conf_basis_on_dev->modal_to_quad_nodal(vtsq_d, vtsq_quad, n);
      conf_basis_on_dev->modal_to_quad_nodal(bmag_d, bmag_quad, n);
      conf_basis_on_dev->modal_to_quad_nodal(jactot_d, jactot_quad, n);
    }

    // Amplitude of the exponential.
    if ((den_quad[n] > 0.) && (vtsq_quad[n]>0.))
      expamp_quad[n] = jactot_quad[n]*den_quad[n]/sqrt(pow(2.0*GKYL_PI*vtsq_quad[n], vdim_phys));
    else
      expamp_quad[n] = 0.;
  }
}

__global__ static void
gkyl_construct_gkmaxwellian_at_quad_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r, struct gkyl_range vel_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  struct gkyl_array* GKYL_RESTRICT den_quad_d,
  struct gkyl_array* GKYL_RESTRICT upar_quad_d,
  struct gkyl_array* GKYL_RESTRICT vtsq_quad_d,
  struct gkyl_array* GKYL_RESTRICT bmag_quad_d, 
  struct gkyl_array* GKYL_RESTRICT jactot_quad_d,
  struct gkyl_array* GKYL_RESTRICT expamp_quad_d,
  const int *p2c_qidx, struct gkyl_array* GKYL_RESTRICT vmap, struct gkyl_basis* GKYL_RESTRICT vmap_basis,
  double mass, struct gkyl_array* GKYL_RESTRICT fm_quad)
{
  double fJacB_floor = 1.e-40;
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;

  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  // 2D thread grid
  // linc2 = c where c is the component index (from 0 to tot_phase_quad)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get configuration-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_r, cidx);

    double *den_quad = (double*) gkyl_array_fetch(den_quad_d, lincC);
    double *upar_quad = (double*) gkyl_array_fetch(upar_quad_d, lincC);
    double *vtsq_quad = (double*) gkyl_array_fetch(vtsq_quad_d, lincC);
    double *bmag_quad = (double*) gkyl_array_fetch(bmag_quad_d, lincC);  
    double *jactot_quad = (double*) gkyl_array_fetch(jactot_quad_d, lincC);  
    double *expamp_quad = (double*) gkyl_array_fetch(expamp_quad_d, lincC);  

    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lidx = gkyl_range_idx(&phase_r, pidx);

    // Select for a phase space index fq.
    double *fq = (double *) gkyl_array_fetch(fm_quad, lidx);

    int n = linc2;

    int cqidx = p2c_qidx[n];

    for (unsigned int d=cdim; d<pdim; d++) vidx[d-cdim] = pidx[d];
    long vlinidx = gkyl_range_idx(&vel_r, vidx);
    const double *vmap_d = (const double *) gkyl_array_cfetch(vmap, vlinidx);
    const double *xcomp_d = (const double *) gkyl_array_cfetch(phase_ordinates, n);

    // Convert comp velocity coordinate to phys velocity coord.
    double xcomp[1];
    for (int vd=0; vd<vdim; vd++) {
      xcomp[0] = xcomp_d[cdim+vd];
      xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
    }

    double efact = 0.0;
    // vpar term.
    efact += pow(xmu[cdim]-upar_quad[cqidx],2);
    // mu term (only for 2v, vdim_phys=3).
    efact += (vdim_phys-1)*xmu[cdim+1]*bmag_quad[cqidx]/mass;

    fq[n] = vtsq_quad[cqidx] > 0.0 ? fJacB_floor+expamp_quad[cqidx]*exp(-efact/(2.0*vtsq_quad[cqidx])) : fJacB_floor;
  }
}

__global__ static void
gkyl_proj_gkmaxwellian_on_basis_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, 
  const struct gkyl_basis* phase_basis_on_dev, struct gkyl_array* GKYL_RESTRICT fm_quad, struct gkyl_array* GKYL_RESTRICT fmax)
{
  int num_phase_basis = fmax->ncomp;
  int tot_phase_quad = basis_at_ords->size;
  int pidx[GKYL_MAX_DIM];

  // 2D thread grid
  // linc2 = c where c is the component index (from 0 to num_phase_basis)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
    tid < phase_r.volume; tid+= blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    long lidx = gkyl_range_idx(&phase_r, pidx);
    double *fm = (double *) gkyl_array_fetch(fmax, lidx);
    
    // Select for a phase space index fq
    double *fq = (double *) gkyl_array_fetch(fm_quad, lidx);

    // Convert back to modal basis; take the thread id for the basis function 
    // so we cant parallelize over basis functions.
    phase_basis_on_dev->quad_nodal_to_modal(fq, fm, linc2);
     
    // The following is modeled after proj_on_basis in the private header.
    /*
    const double *phase_w = (const double*) phase_weights->data;
    const double *phaseb_o = (const double*) basis_at_ords->data;
    for (int n=0; n<tot_phase_quad; ++n) {
      double tmp = phase_w[n]*fq[n];
      fm[linc2] += tmp*phaseb_o[linc2+num_phase_basis*n];
      if (linc2==0) {
        printf("quad=%d\n", n);
        printf("coeff=%g\n", phase_w[n]*phaseb_o[num_phase_basis*n]);
      }
    }*/
  }
}

void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *prim_moms,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, 
  double mass, struct gkyl_array *fmax)
{
  dim3 dimGrid, dimBlock;
  int num_phase_basis = fmax->ncomp;
  int tot_phase_quad = up->basis_at_ords->size;
  int tot_conf_quad = up->conf_basis_at_ords->size;

  const struct gkyl_velocity_map *gvm = up->vel_map;
  
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *conf_r, tot_conf_quad);
  gkyl_expand_prim_mom_at_quad_cu_ker<<<dimGrid, dimBlock>>>
    (up->grid, *phase_r, *conf_r, 
     up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev, 
     up->ordinates->on_dev,
     up->den_quad->on_dev, 
     up->upar_quad->on_dev, 
     up->vtsq_quad->on_dev, 
     up->bmag_quad->on_dev, 
     up->jactot_quad->on_dev, 
     up->expamp_quad->on_dev, 
     //up->fm_quad->on_dev,
     prim_moms->on_dev, bmag->on_dev, jacob_tot->on_dev,
     up->conf_basis_on_dev);     

  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_r, tot_phase_quad);
  gkyl_construct_gkmaxwellian_at_quad_cu_ker<<<dimGrid, dimBlock>>>
    (up->grid, *phase_r, *conf_r, gvm->local_ext_vel, 
     up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev, 
     up->ordinates->on_dev, 
     up->den_quad->on_dev, 
     up->upar_quad->on_dev, 
     up->vtsq_quad->on_dev, 
     up->bmag_quad->on_dev, 
     up->jactot_quad->on_dev, 
     up->expamp_quad->on_dev, 
     up->p2c_qidx, gvm->vmap->on_dev, gvm->vmap_basis, 
     mass, up->fm_quad->on_dev);
    
  gkyl_array_clear(fmax, 0.0);

  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_r, num_phase_basis);
  gkyl_proj_gkmaxwellian_on_basis_cu_ker<<<dimGrid, dimBlock>>>
    (up->grid, *phase_r, *conf_r, 
     up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, 
     up->weights->on_dev,
     up->phase_basis_on_dev, up->fm_quad->on_dev, fmax->on_dev);
}
