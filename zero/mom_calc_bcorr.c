#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc_bcorr_priv.h>
#include <gkyl_mom_bcorr_fpo_vlasov.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_util.h>

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_mom_type *momt, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_calc_bcorr_cu_dev_new(grid, momt);
  } 
#endif
  gkyl_mom_calc_bcorr *up = gkyl_malloc(sizeof(gkyl_mom_calc_bcorr));
  up->grid = *grid;
  up->momt = gkyl_mom_type_acquire(momt);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;
  
  return up;
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

void
gkyl_mom_calc_bcorr_advance(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
  double xc[GKYL_MAX_DIM];
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], viter_idx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  enum gkyl_vel_edge edge;
  
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;
  gkyl_array_clear_range(out, 0.0, conf_rng);

  // outer loop is over configuration space cells; for each
  // config-space cell inner loop walks over the edges of velocity
  // space
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);
    
    for (int d=0; d<conf_rng->ndim; ++d)
      viter_idx[d] = conf_iter.idx[d];
    
    for (int d=0; d<phase_rng->ndim - conf_rng->ndim; ++d) {
      rem_dir[conf_rng->ndim + d] = 1;
      
      // loop over upper edge of velocity space
      edge = d + GKYL_MAX_CDIM;
      viter_idx[conf_rng->ndim + d] = phase_rng->upper[conf_rng->ndim + d];
      gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, viter_idx);
      gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
      
      while (gkyl_range_iter_next(&vel_iter)) {
        copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
        gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
      
        long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
        gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
          gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
        );
      }
      
      // loop over lower edge of velocity space
      edge = d;
      viter_idx[conf_rng->ndim + d] = phase_rng->lower[conf_rng->ndim + d];
      gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, viter_idx);
      gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
      
      while (gkyl_range_iter_next(&vel_iter)) {
        copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
  
        gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
      
        long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
        gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
          gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
        );
      }
      rem_dir[conf_rng->ndim + d] = 0;
      viter_idx[conf_rng->ndim + d] = 0;
    }
  }
}

void
gkyl_mom_calc_bcorr_release(gkyl_mom_calc_bcorr* up)
{
  gkyl_mom_type_release(up->momt);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

// "derived" class constructors
struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, bool use_gpu)
{
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  bcorr_type = gkyl_mom_bcorr_lbo_vlasov_new(cbasis, pbasis, vBoundary, use_gpu);  
  struct gkyl_mom_calc_bcorr* calc = gkyl_mom_calc_bcorr_new(grid, bcorr_type, use_gpu);
  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_mom_calc_bcorr_release method on calculator deallocates specific type data
  gkyl_mom_type_release(bcorr_type);
  return calc;
}

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_pkpm_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, double mass, bool use_gpu)
{
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  bcorr_type = gkyl_mom_bcorr_lbo_pkpm_new(cbasis, pbasis, vBoundary, mass, use_gpu);  
  struct gkyl_mom_calc_bcorr* calc = gkyl_mom_calc_bcorr_new(grid, bcorr_type, use_gpu);
  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_mom_calc_bcorr_release method on calculator deallocates specific type data
  gkyl_mom_type_release(bcorr_type);
  return calc;
}

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, double mass, bool use_gpu)
{
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  bcorr_type = gkyl_mom_bcorr_lbo_gyrokinetic_new(cbasis, pbasis, vBoundary, mass, use_gpu);
  struct gkyl_mom_calc_bcorr* calc = gkyl_mom_calc_bcorr_new(grid, bcorr_type, use_gpu);
  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_mom_calc_bcorr_release method on calculator deallocates specific type data
  gkyl_mom_type_release(bcorr_type);
  return calc;
}

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_fpo_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range,
  const double* vBoundary, const struct gkyl_array* diff_coeff, bool use_gpu)
{
  struct gkyl_mom_type *bcorr_type; // FPO boundary corrections moment type
  bcorr_type = gkyl_mom_bcorr_fpo_vlasov_new(cbasis, pbasis, phase_range, vBoundary, use_gpu);  

  // Set auxfields
  struct gkyl_mom_bcorr_fpo_vlasov_auxfields fpo_bcorr_auxfields = {.D = diff_coeff};
  gkyl_mom_bcorr_fpo_vlasov_set_auxfields(bcorr_type, fpo_bcorr_auxfields);
    
  struct gkyl_mom_calc_bcorr* calc = gkyl_mom_calc_bcorr_new(grid, bcorr_type, use_gpu);

  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_mom_calc_bcorr_release method on calculator deallocates specific type data
  gkyl_mom_type_release(bcorr_type);
  return calc;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_mom_type *momt)
{
  assert(false);
}

void
gkyl_mom_calc_bcorr_advance_cu(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
  assert(false);
}

#endif
