#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

// Increment an int vector by fact*del[d] in each direction d.
static inline void
incr_int_array(int ndim, int fact, const int * GKYL_RESTRICT del,
  const int * GKYL_RESTRICT inp, int *GKYL_RESTRICT out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

// KB - Number of cells in the positive and negative directions
// is assumed to be the same symmetrically, half of the cells in that dimension. Need to
// do this more rigorously.
void
gkyl_bc_emission_pos_neg_ranges(struct gkyl_range *pos, struct gkyl_range *neg,
  int dir, const struct gkyl_range *parent, const int *nghost)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  incr_int_array(ndim, 0, nghost, parent->lower, lo);
  incr_int_array(ndim, 0, nghost, parent->upper, up);

  double nneg = (up[dir]-lo[dir] + 1)/2;
  double npos = (up[dir]-lo[dir] + 1)/2;
  
  up[dir] = lo[dir] + nneg;
  gkyl_sub_range_init(neg, parent, lo, up);

  incr_int_array(ndim, 0, nghost, parent->lower, lo);
  incr_int_array(ndim, 0, nghost, parent->upper, up);
    
  lo[dir] = up[dir] - npos;
  gkyl_sub_range_init(pos, parent, lo, up);
}

void
gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_emission_spectrum_sey_calc_cu(up, yield, grid, ghost_r);
    return;
  }
#endif
  double xc[GKYL_MAX_DIM];
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(ghost_r, iter.idx);

    double *out = gkyl_array_fetch(yield, loc);
    
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    up->funcs->yield(out, up->cdim, up->vdim, xc, up->yield_param);
  }
}

struct gkyl_bc_emission_spectrum*
gkyl_bc_emission_spectrum_new(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  enum gkyl_bc_emission_spectrum_yield_type yield_type, void *norm_param, void *yield_param,
  struct gkyl_array *yield, struct gkyl_array *spectrum, int dir, enum gkyl_edge_loc edge,
  int cdim, int vdim, struct gkyl_range *impact_skin_r, struct gkyl_range *impact_ghost_r,
  struct gkyl_rect_grid *grid, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission_spectrum *up = gkyl_malloc(sizeof(struct gkyl_bc_emission_spectrum));

  up->dir = dir;
  up->cdim = cdim;
  up->vdim = vdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->norm_param = norm_param;
  up->yield_param = yield_param;
  up->impact_skin_r = impact_skin_r;
  up->impact_ghost_r = impact_ghost_r;


  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    ghost[cdim+d] = 0;
  }

  gkyl_bc_emission_pos_neg_ranges(up->pos_r, up->neg_r, dir, impact_skin_r, ghost);

  up->funcs = gkyl_malloc(sizeof(struct gkyl_bc_emission_spectrum_funcs));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->funcs_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_spectrum_funcs));
    gkyl_bc_emission_spectrum_choose_func_cu(norm_type, yield_type, up->funcs_cu);
    up->bc_param_cu = gkyl_cu_malloc(10*sizeof(double));
    gkyl_cu_memcpy(up->bc_param_cu, up->bc_param, 10*sizeof(double), GKYL_CU_MEMCPY_H2D);
    up->sey_param_cu = gkyl_cu_malloc(10*sizeof(double));
    gkyl_cu_memcpy(up->yield_param_cu, up->yield_param, 10*sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    up->funcs->func = bc_weighted_delta;
    up->funcs->norm = bc_emission_spectrum_choose_norm_func(norm_type);
    up->funcs->yield = bc_emission_spectrum_choose_yield_func(yield_type);
    up->funcs_cu = up->funcs;
    up->norm_param_cu = up->norm_param;
    up->yield_param_cu = up->yield_param;
  }
#else
  up->funcs->func = bc_weighted_delta;
  up->funcs->norm = bc_emission_spectrum_choose_norm_func(norm_type);
  up->funcs->yield = bc_emission_spectrum_choose_yield_func(yield_type);
  up->funcs_cu = up->funcs;
  up->norm_param_cu = up->norm_param;
  up->yield_param_cu = up->yield_param;

  gkyl_bc_emission_spectrum_sey_calc(up, yield, grid, impact_ghost_r);
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
gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  const struct gkyl_array *f_skin, const struct gkyl_array *f_proj,
  struct gkyl_array *f_buff, struct gkyl_rect_grid *grid)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_bc_emission_spectrum_advance_cu(up, f_skin, f_proj, f_buff, grid);
  }
#endif
  /* double xc[GKYL_MAX_DIM]; */
  /* int cidx[GKYL_MAX_CDIM], pidx[GKYL_MAX_DIM], fidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 }; */
  /* for (int d=0; d<conf_r->ndim; ++d) rem_dir[d] = 1; */

  /* struct gkyl_range vel_r, vel_ghost_r, vel_buff_r; */
  /* struct gkyl_range_iter conf_iter, vel_iter; */

  /* gkyl_array_clear_range(weight, 0.0, conf_r); */

  /* gkyl_range_iter_init(&conf_iter, conf_r); */
  /* while (gkyl_range_iter_next(&conf_iter)) { */
  /*   long midx = gkyl_range_idx(conf_r, conf_iter.idx); */

  /*   gkyl_copy_int_arr(conf_r->ndim, conf_iter.idx, cidx); */
  /*   cidx[up->dir] = skin_r->lower[up->dir]; */

  /*   gkyl_range_deflate(&vel_r, skin_r, rem_dir, cidx); */
  /*   gkyl_range_deflate(&vel_ghost_r, ghost_r, rem_dir, conf_iter.idx); */
  /*   gkyl_range_deflate(&vel_buff_r, buff_r, rem_dir, conf_iter.idx); */
  /*   gkyl_range_iter_no_split_init(&vel_iter, &vel_r); */

  /*   double *w = gkyl_array_fetch(weight, midx); */

  /*   while (gkyl_range_iter_next(&vel_iter)) { */
  /*     copy_idx_arrays(conf_r->ndim, skin_r->ndim, conf_iter.idx, vel_iter.idx, pidx); */
  /*     pidx[up->dir] = skin_r->lower[up->dir]; */
  /*     gkyl_rect_grid_cell_center(grid, pidx, xc); */

  /*     long loc = gkyl_range_idx(&vel_r, vel_iter.idx); */
  /*     gkyl_copy_int_arr(skin_r->ndim, pidx, fidx); */
  /*     fidx[up->dir] = ghost_r->lower[up->dir]; */
  /*     long loc2 = gkyl_range_idx(ghost_r, fidx); */

  /*     const double *inp = gkyl_array_cfetch(f_skin, loc); */
  /*     const double *gain = gkyl_array_cfetch(gamma, loc2); */
      
  /*     up->funcs->func(inp, up->cdim, up->dir, up->edge, xc, gain, w); */
  /*   } */
  /*   double effective_gamma = w[0]/w[1]; */
  /*   const double *bflux = gkyl_array_cfetch(flux, midx); */
  /*   double *out = gkyl_array_fetch(k, midx); */
    
  /*   up->funcs->norm(out, bflux, up->bc_param, effective_gamma); */
  /*   gkyl_array_accumulate_range(f_buff, out[0], f_proj, &vel_buff_r); */
  /* } */
}

void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->funcs_cu);
    gkyl_cu_free(up->bc_param_cu);
    gkyl_cu_free(up->sey_param_cu);
  }
#endif
  gkyl_free(up->funcs);
  // Release updater memory.
  gkyl_free(up);
}
