#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

struct gkyl_bc_emission_elastic*
gkyl_bc_emission_elastic_new(enum gkyl_bc_emission_elastic_type elastic_type,
  void *elastic_param, struct gkyl_array *elastic_yield, int dir, enum gkyl_edge_loc edge,
  int cdim, int vdim, struct gkyl_rect_grid *grid, struct gkyl_range *emit_buff_r, int poly_order,
  struct gkyl_basis *basis, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission_elastic *up = gkyl_malloc(sizeof(struct gkyl_bc_emission_elastic));

  up->dir = dir;
  up->cdim = cdim;
  up->vdim = vdim;
  up->edge = edge;
  up->basis = basis;
  up->use_gpu = use_gpu;
  up->elastic_param = elastic_param;

  struct gkyl_array_copy_func *cf = gkyl_malloc(sizeof(*cf));
  cf->func = reflection;
  cf->ctx = up;
  cf->ctx_on_dev = cf->ctx;

  cf->flags = 0;
  cf->on_dev = cf; // CPU function obj points to itself.
  up->reflect_func = cf;

  up->funcs = gkyl_malloc(sizeof(struct gkyl_bc_emission_elastic_funcs));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->funcs_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_elastic_funcs));
    gkyl_bc_emission_elastic_choose_func_cu(elastic_type, up->funcs_cu);
    up->elastic_param_cu = gkyl_cu_malloc(10*sizeof(double));
    gkyl_cu_memcpy(up->elastic_param_cu, up->elastic_param, 10*sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    up->funcs->yield = bc_emission_elastic_choose_yield_func(elastic_type);
    up->funcs_cu = up->funcs;
    up->elastic_param_cu = up->elastic_param;

    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(grid, basis, poly_order + 1, 1,
      up->funcs->yield, up);
    gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, elastic_yield);
    gkyl_proj_on_basis_release(proj);
  }
#else
  up->funcs->yield = bc_emission_elastic_choose_yield_func(elastic_type);
  up->funcs_cu = up->funcs;
  up->elastic_param_cu = up->elastic_param;
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(grid, basis, poly_order + 1, 1,
    up->funcs->yield, up);
  gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, elastic_yield);
  gkyl_proj_on_basis_release(proj);
#endif

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
gkyl_bc_emission_elastic_advance(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_bc_emission_elastic_advance_cu(up, emit_skin_r, buff_arr, f_skin, f_emit, elastic_yield,
      basis);
  }
#endif
  gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_skin, up->dir+up->cdim, emit_skin_r,
    up->reflect_func->on_dev);
  gkyl_dg_mul_op(*basis, 0, f_emit, 0, buff_arr, 0, elastic_yield);
}

void gkyl_bc_emission_elastic_release(struct gkyl_bc_emission_elastic *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->funcs_cu);
    gkyl_cu_free(up->elastic_param_cu);
  }
#endif
  gkyl_free(up->funcs);
  gkyl_free(up->reflect_func);
  gkyl_free(up->elastic_param);
  // Release updater memory.
  gkyl_free(up);
}
