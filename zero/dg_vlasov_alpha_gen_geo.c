#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_alpha_gen_geo.h>
#include <gkyl_dg_vlasov_alpha_gen_geo_priv.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,poly_order) lst[0].kernels[poly_order]

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

// calculate alpha gen geo
void
gkyl_dg_alpha_gen_geo(struct gkyl_basis *conf_basis,
  struct gkyl_basis *phase_basis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_rect_grid *grid,const struct gkyl_array* tv_comp,
  const struct gkyl_array* gij, struct gkyl_array* alpha_geo)
{
  // Add GPU capability later...
  /*#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_alpha_gen_geo_cu();
  }
  #endif*/
  
  int num_basis = phase_basis->num_basis;
  int cdim = conf_basis->ndim;
  int vdim = phase_basis->ndim - cdim;
  int poly_order = conf_basis->poly_order;
  int pidx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  // Assert that grid is 3x3v and then select kernel here
  assert(cdim + vdim == 6);
  dg_vlasov_alpha_gen_geof_t alpha_gen_geo = CK(ser_vlasov_alpha_gen_geo_kernels, poly_order); 

  struct gkyl_range_iter phase_iter;

  gkyl_range_iter_init(&phase_iter, phase_rng);
  while (gkyl_range_iter_next(&phase_iter)) {
    long ploc = gkyl_range_idx(phase_rng, phase_iter.idx);
    double *alpha_geo_d = gkyl_array_fetch(alpha_geo, ploc);

    int cidx[3]; 
    for (int d=0; d<cdim; d++) cidx[d] = phase_iter.idx[d];
    long cloc = gkyl_range_idx(conf_rng, cidx);
    const double *gij_d = gkyl_array_cfetch(gij, cloc);
    const double *tv_comp_d = gkyl_array_cfetch(tv_comp, cloc);
    gkyl_rect_grid_cell_center(grid, phase_iter.idx, xc);

    // Call alpha_gen_geo kernel.
    alpha_gen_geo(xc, grid->dx, tv_comp_d, gij_d, alpha_geo_d);

  }
  
}
