/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

// CUDA kernel to set device pointers to function for normalization coefficient calculation.
__global__ static void
gkyl_bc_emission_spectrum_set_cu_norm_func_ptrs(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  funcs->func = bc_weighted_delta;
  switch (norm_type) {
    case GKYL_SEE_CHUNG_EVERHART:
      funcs->norm = chung_everhart_norm;
      break;
    case GKYL_SEE_GAUSSIAN:
      funcs->norm = gaussian_norm;
      break;
    case GKYL_SEE_MAXWELLIAN:
      funcs->norm = maxwellian_norm;
      break;
    default:
      assert(false);
      break;
  }
};

// CUDA kernel to set device pointers to function for SEY calculation.
__global__ static void
gkyl_bc_emission_spectrum_set_cu_yield_func_ptrs(enum gkyl_bc_emission_spectrum_yield_type yield_type,
  struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  switch (yield_type) {
    case GKYL_SEE_FURMAN_PIVI:
      funcs->yield = furman_pivi_yield;
      break;
    case GKYL_SEE_SCHOU:
      funcs->yield = schou_yield;
      break;
    case GKYL_SEE_CONSTANT:
      funcs->yield = constant_yield;
      break; 
    default:
      assert(false);
      break;
  }
};

void
gkyl_bc_emission_spectrum_choose_func_cu(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  enum gkyl_bc_emission_spectrum_yield_type yield_type, struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  gkyl_bc_emission_spectrum_set_cu_norm_func_ptrs<<<1,1>>>(norm_type, funcs);
  gkyl_bc_emission_spectrum_set_cu_yield_func_ptrs<<<1,1>>>(yield_type, funcs);
}

__global__ static void
gkyl_bc_emission_spectrum_sey_calc_cu_ker(int cdim, int vdim, double *sey_param, struct gkyl_rect_grid grid, const struct gkyl_range ghost_r, struct gkyl_array *gamma, struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < ghost_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange

    gkyl_sub_range_inv_idx(&ghost_r, linc, pidx);

    long loc = gkyl_range_idx(&ghost_r, pidx);

    double *out = (double *) gkyl_array_fetch(gamma, loc);
    
    gkyl_rect_grid_cell_center(&grid, pidx, xc);
    
    funcs->gamma(out, cdim, vdim, xc, sey_param);
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_weight_ker(int cdim, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_array *f_skin, struct gkyl_array *weight, struct gkyl_rect_grid grid, struct gkyl_array *gamma,
  const struct gkyl_range skin_r, const struct gkyl_range ghost_r, const struct gkyl_range conf_r, struct gkyl_bc_emission_spectrum_funcs *funcs,
  double *bc_param)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], fidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < skin_r.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&skin_r, tid, pidx);
    gkyl_copy_int_arr(skin_r.ndim, pidx, fidx);
    fidx[dir] = ghost_r.lower[dir];
    
    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lincS = gkyl_range_idx(&skin_r, pidx);
    long lincG = gkyl_range_idx(&ghost_r, fidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(f_skin, lincS);
    const double* gain = (const double*) gkyl_array_cfetch(gamma, lincG);
    double wLocal[2];
    for (unsigned int k=0; k<weight->ncomp; ++k)
      wLocal[k] = 0.0;

    funcs->func(inp, cdim, dir, edge, xc, gain, &wLocal[0]);

    // get conf-space linear index.
    for (unsigned int i = 0; i < conf_r.ndim; i++)
      cidx[i] = fidx[i];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    double* wptr = (double*) gkyl_array_fetch(weight, lincC);
    for (unsigned int k = 0; k < weight->ncomp; ++k) {
       if (tid < skin_r.volume)
         atomicAdd(&wptr[k], wLocal[k]);
    }
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_accumulate_ker(const struct gkyl_array *f_proj, struct gkyl_array *f_buff, struct gkyl_array *weight, struct gkyl_array *k, const struct gkyl_array *flux, const struct gkyl_range buff_r, const struct gkyl_range conf_r, struct gkyl_bc_emission_spectrum_funcs *funcs, double *bc_param)
{
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < buff_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&buff_r, tid, pidx);
    
    // get conf-space linear index.
    for (unsigned int i = 0; i < conf_r.ndim; i++)
      cidx[i] = pidx[i];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    long lincP = gkyl_range_idx(&buff_r, pidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(f_proj, lincP);
    const double* w = (const double*) gkyl_array_cfetch(weight, lincC);
    const double* bflux = (const double*) gkyl_array_cfetch(flux, lincC);
    double* fac = (double*) gkyl_array_fetch(k, lincC);
    double* out = (double*) gkyl_array_fetch(f_buff, lincP);
    double effective_gamma = w[0]/w[1];
    funcs->norm(fac, bflux, bc_param, effective_gamma);
    
    for (int c=0; c<f_proj->ncomp; ++c)
      out[c] += fac[0]*inp[c];
  }
}

void
gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r, const struct gkyl_range *gamma_r)
{
  int nblocks = ghost_r->nblocks, nthreads = ghost_r->nthreads;

  gkyl_bc_emission_spectrum_sey_calc_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->vdim, up->yield_param_cu, *grid, *ghost_r, gamma->on_dev, up->funcs_cu);
}

void
gkyl_bc_emission_spectrum_advance_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k)
{
  int nblocks = emit_buff_r->nblocks, nthreads = emit_buff_r->nthreads;

  // gkyl_array_clear_range(weight, 0.0, conf_r);

  // // Calculate weighted mean numerator and denominator
  // gkyl_bc_emission_spectrum_advance_cu_weight_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, up->edge, f_skin->on_dev, weight->on_dev, *grid, gamma->on_dev, *skin_r, *ghost_r, *conf_r, up->funcs_cu, up->bc_param_cu);

  // nblocks = buff_r->nblocks;
  // nthreads = buff_r->nthreads;

  // // Finish weighted mean calculation and accumulate to buffer
  // gkyl_bc_emission_spectrum_advance_cu_accumulate_ker<<<nblocks, nthreads>>>(f_proj->on_dev, f_buff->on_dev, weight->on_dev, k->on_dev, flux->on_dev, *buff_r, *conf_r, up->funcs_cu, up->bc_param_cu);
}
