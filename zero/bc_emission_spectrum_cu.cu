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
  struct gkyl_bc_emission_spectrum_funcs *funcs, void *norm_param_cu)
{
  funcs->func = bc_weighted_delta;
  switch (norm_type) {
    case GKYL_SEE_CHUNG_EVERHART:
      funcs->norm = chung_everhart_norm;
      funcs->norm_param = (struct gkyl_bc_emission_spectrum_norm_chung_everhart *) norm_param_cu;
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
  struct gkyl_bc_emission_spectrum_funcs *funcs, void *yield_param_cu)
{
  switch (yield_type) {
    case GKYL_SEE_FURMAN_PIVI:
      funcs->yield = furman_pivi_yield;
      funcs->yield_param = (struct gkyl_bc_emission_spectrum_yield_furman_pivi *) yield_param_cu;
      break;
    case GKYL_SEE_SCHOU:
      funcs->yield = schou_yield;
      funcs->yield_param = (struct gkyl_bc_emission_spectrum_yield_schou *) yield_param_cu;
      break;
    case GKYL_SEE_CONSTANT:
      funcs->yield = constant_yield;
      funcs->yield_param = (struct gkyl_bc_emission_spectrum_yield_constant *) yield_param_cu;
      break; 
    default:
      assert(false);
      break;
  }
};

void
gkyl_bc_emission_spectrum_choose_norm_cu(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  struct gkyl_bc_emission_spectrum_funcs *funcs, void *norm_param)
{
  void *norm_param_cu;
  switch (norm_type) {
    case GKYL_SEE_CHUNG_EVERHART:
      norm_param_cu = (struct gkyl_bc_emission_spectrum_norm_chung_everhart *) gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_norm_chung_everhart));
      gkyl_cu_memcpy(norm_param_cu,
        (struct gkyl_bc_emission_spectrum_norm_chung_everhart *) norm_param,
        sizeof(struct gkyl_bc_emission_spectrum_norm_chung_everhart), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_SEE_GAUSSIAN:
      norm_param_cu = (struct gkyl_bc_emission_spectrum_norm_gaussian *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_norm_gaussian));
      gkyl_cu_memcpy(norm_param_cu,
        (struct gkyl_bc_emission_spectrum_norm_gaussian *) norm_param,
        sizeof(struct gkyl_bc_emission_spectrum_norm_gaussian), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_SEE_MAXWELLIAN:
      norm_param_cu = (struct gkyl_bc_emission_spectrum_norm_maxwellian *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_norm_maxwellian));
      gkyl_cu_memcpy(funcs->norm_param,
        (struct gkyl_bc_emission_spectrum_norm_maxwellian *) norm_param,
        sizeof(struct gkyl_bc_emission_spectrum_norm_maxwellian), GKYL_CU_MEMCPY_H2D);
      break;
    default:
      assert(false);
      break;
  }
  gkyl_bc_emission_spectrum_set_cu_norm_func_ptrs<<<1,1>>>(norm_type, funcs, norm_param_cu);
};

void
gkyl_bc_emission_spectrum_choose_yield_cu(enum gkyl_bc_emission_spectrum_yield_type yield_type, struct gkyl_bc_emission_spectrum_funcs *funcs, void *yield_param)
{
  void *yield_param_cu;
  switch (yield_type) {
    case GKYL_SEE_FURMAN_PIVI: 
      yield_param_cu = (struct gkyl_bc_emission_spectrum_yield_furman_pivi *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_yield_furman_pivi));
      gkyl_cu_memcpy(yield_param_cu,
        (struct gkyl_bc_emission_spectrum_yield_furman_pivi *) yield_param,
        sizeof(struct gkyl_bc_emission_spectrum_yield_furman_pivi), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_SEE_SCHOU:
      yield_param_cu = (struct gkyl_bc_emission_spectrum_yield_schou *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_yield_schou));
      gkyl_cu_memcpy(yield_param_cu,
        (struct gkyl_bc_emission_spectrum_yield_schou *) yield_param,
        sizeof(struct gkyl_bc_emission_spectrum_yield_schou), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_SEE_CONSTANT:
      yield_param_cu = (struct gkyl_bc_emission_spectrum_yield_constant *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_yield_constant));
      gkyl_cu_memcpy(yield_param_cu,
        (struct gkyl_bc_emission_spectrum_yield_constant *) yield_param,
        sizeof(struct gkyl_bc_emission_spectrum_yield_constant), GKYL_CU_MEMCPY_H2D);
      break;
    default:
      assert(false);
      break;
  }
  gkyl_bc_emission_spectrum_set_cu_yield_func_ptrs<<<1,1>>>(yield_type, funcs, yield_param_cu);
};

__global__ static void
gkyl_bc_emission_spectrum_sey_calc_cu_ker(int cdim, int vdim, struct gkyl_rect_grid grid, const struct gkyl_range ghost_r, struct gkyl_array *yield, struct gkyl_bc_emission_spectrum_funcs *funcs)
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
    
    double *out = (double *) gkyl_array_fetch(yield, loc);
    
    gkyl_rect_grid_cell_center(&grid, pidx, xc);
    
    funcs->yield(out, cdim, vdim, xc, funcs->yield_param);
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_weight_ker(int cdim, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_array *bflux, struct gkyl_array *weight, struct gkyl_rect_grid grid,
  struct gkyl_array *yield, const struct gkyl_range impact_buff_r,
  const struct gkyl_range impact_cbuff_r, struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < impact_buff_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&impact_buff_r, tid, pidx);
    
    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lincP = gkyl_range_idx(&impact_buff_r, pidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(bflux, lincP);
    const double* gain = (const double*) gkyl_array_cfetch(yield, lincP);
    double wLocal[2];
    for (unsigned int k=0; k<weight->ncomp; ++k)
      wLocal[k] = 0.0;

    funcs->func(inp, cdim, dir, edge, xc, gain, &wLocal[0]);
    
    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++)
      cidx[i] = pidx[i];
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    double* wptr = (double*) gkyl_array_fetch(weight, lincC);
    for (unsigned int k = 0; k < weight->ncomp; ++k) {
       if (tid < impact_buff_r.volume)
         atomicAdd(&wptr[k], wLocal[k]);
    }
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_accumulate_ker(const struct gkyl_array *spectrum, struct gkyl_array *f_emit, struct gkyl_array *weight, struct gkyl_array *k, const struct gkyl_array *flux, const struct gkyl_range emit_buff_r, const struct gkyl_range impact_cbuff_r, struct gkyl_bc_emission_spectrum_funcs *funcs)
{
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < emit_buff_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&emit_buff_r, tid, pidx);
    
    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++)
      cidx[i] = pidx[i];
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    long lincP = gkyl_range_idx(&emit_buff_r, pidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(spectrum, lincP);
    const double* w = (const double*) gkyl_array_cfetch(weight, lincC);
    const double* boundary_flux = (const double*) gkyl_array_cfetch(flux, lincC);
    double* fac = (double*) gkyl_array_fetch(k, lincC);
    double* out = (double*) gkyl_array_fetch(f_emit, lincP);
    double effective_delta = w[0]/w[1];
    funcs->norm(fac, boundary_flux, funcs->norm_param, effective_delta);

    for (int c=0; c<spectrum->ncomp; ++c)
      out[c] += fac[0]*inp[c];

  }
}

void
gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *gamma_r)
{
  int nblocks = gamma_r->nblocks, nthreads = gamma_r->nthreads;

  gkyl_bc_emission_spectrum_sey_calc_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->vdim, *grid, *gamma_r, yield->on_dev, up->funcs_cu);
}

void
gkyl_bc_emission_spectrum_advance_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k)
{
  int nblocks = impact_buff_r->nblocks, nthreads = impact_buff_r->nthreads;

  // Calculate weighted mean numerator and denominator
  gkyl_bc_emission_spectrum_advance_cu_weight_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, up->edge, bflux->on_dev, weight->on_dev, *up->grid, yield->on_dev, *impact_buff_r, *impact_cbuff_r, up->funcs_cu);

  nblocks = emit_buff_r->nblocks;
  nthreads = emit_buff_r->nthreads;

  // Finish weighted mean calculation and accumulate to buffer
  gkyl_bc_emission_spectrum_advance_cu_accumulate_ker<<<nblocks, nthreads>>>(spectrum->on_dev, f_emit->on_dev, weight->on_dev, k->on_dev, flux->on_dev, *emit_buff_r, *impact_cbuff_r, up->funcs_cu);
}
