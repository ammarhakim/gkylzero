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
gkyl_dg_alpha_gen_geo(struct gkyl_basis basis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_rect_grid *grid,const struct gkyl_array *tv_comp,
  const struct gkyl_array *gij, struct gkyl_array *alpha_geo)
{
  // Add GPU capability later...
  /*#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_alpha_gen_geo_cu();
  }
  #endif*/
  
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  int pidx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  // Assert that grid is 3x3v and then select kernel here
  assert(ndim == 6);
  dg_vlasov_alpha_gen_geof_t alpha_gen_geo = CK(ser_vlasov_alpha_gen_geo_kernels, poly_order); 

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc_c = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *gij_d = gkyl_array_cfetch(gij, loc_c);

    // Loop over velocity space. 
    while (gkyl_range_iter_next(&vel_iter)) {
      long loc_v = gkyl_range_idx(&vel_rng, vel_iter.idx);
      const double *tv_comp_d = gkyl_array_cfetch(tv_comp, loc_v);
      double *alpha_geo_d = gkyl_array_fetch(alpha_geo, loc_v);

      // Get cell center and width.
      copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(grid, pidx, xc);

      // Call alpha_gen_geo kernel.
      alpha_gen_geo(xc, grid->dx, tv_comp_d, gij_d, alpha_geo_d);
    }
  }
  
}
