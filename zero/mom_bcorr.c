#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_bcorr.h>
#include <gkyl_util.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

void
gkyl_mom_bcorr_advance(gkyl_mom_bcorr *bcorr,
  struct gkyl_range phase_rng, struct gkyl_range conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out)
{
  double xc[GKYL_MAX_DIM];
  struct gkyl_range vel_rng, vel_edge;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng.ndim; ++d) rem_dir[d] = 1;

  gkyl_array_clear_range(out, 0.0, conf_rng);

  // the outer loop is over configuration space cells; for each
  // config-space cell the inner loop walks over the velocity space
  // computing the contribution to the moment
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);

    gkyl_range_deflate(&vel_rng, &phase_rng, rem_dir, conf_iter.idx);
    gkyl_sub_range_init(&vel_edge, &vel_rng, vel_rng.upper, vel_rng.upper);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_edge);

    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_rng.ndim, phase_rng.ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
      
      long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      if (bcorr->momt) {
        printf("momt is here!\n");
        fflush(stdout);
      }
      printf("xc: %f\n", xc[0]);
      printf("dx: %f\n", bcorr->grid.dx[0]);
      fflush(stdout);

      gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
        gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx)
      );
    }
  }
}

gkyl_mom_bcorr* gkyl_mom_bcorr_new(const struct gkyl_rect_grid *grid, const struct gkyl_mom_type *momt)
{
  gkyl_mom_bcorr *up = gkyl_malloc(sizeof(gkyl_mom_bcorr));
  up->grid = *grid;
  up->momt = gkyl_mom_type_aquire(momt);
  
  return up;
}

void gkyl_mom_bcorr_release(gkyl_mom_bcorr* bcorr)
{
  gkyl_mom_type_release(bcorr->momt);
  free(bcorr);
}
