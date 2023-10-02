#include <gkyl_calc_bmag.h>
#include <gkyl_calc_bmag_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_bmag*
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, const gkyl_geo_gyrokinetic *app, const struct gkyl_geo_gyrokinetic_geo_inp *ginp, bool use_gpu)
{
  gkyl_calc_bmag *up = gkyl_malloc(sizeof(gkyl_calc_bmag));
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cgrid = cgrid;
  up->pgrid = pgrid;
  up->use_gpu = use_gpu;
  up->kernel = bmag_choose_kernel(up->pbasis->ndim, up->pbasis->b_type, up->pbasis->poly_order);
  up->app = app;
  up->ginp = ginp;
  return up;
}

static inline void bmag_comp(double t, const double *xn, double *fout, void *ctx)
{
  struct bmag_ctx *gc = (struct bmag_ctx*) ctx;
  struct gkyl_range_iter iter;
  double XYZ[gc->cgrid->ndim];

  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&iter, gc->crange);
  for(int i = 0; i < gc->cgrid->ndim; i++){
    citer.idx[i] = fmin(gc->crange->lower[i] + (int) floor((xn[i] - (gc->cgrid->lower[i]) )/gc->cgrid->dx[i]), gc->crange->upper[i]);
  }

  long lidx = gkyl_range_idx(gc->crange, citer.idx);
  const double *mcoeffs = gkyl_array_cfetch(gc->mapc2p, lidx);
  

  double cxc[gc->cgrid->ndim];
  double xyz[gc->cgrid->ndim];
  gkyl_rect_grid_cell_center(gc->cgrid, citer.idx, cxc);
  for(int i = 0; i < gc->cgrid->ndim; i++)
    xyz[i] = (xn[i]-cxc[i])/(gc->cgrid->dx[i]*0.5);
  for(int i = 0; i < gc->cgrid->ndim; i++){
    XYZ[i] = gc->cbasis->eval_expand(xyz, &mcoeffs[i*gc->cbasis->num_basis]);
  }

  double R = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1]);
  double Z = XYZ[2];

  gkyl_range_iter_init(&iter, gc->range);
  iter.idx[0] = fmin(gc->range->lower[0] + (int) floor((R - gc->grid->lower[0])/gc->grid->dx[0]), gc->range->upper[0]);
  iter.idx[1] = fmin(gc->range->lower[1] + (int) floor((Z - gc->grid->lower[1])/gc->grid->dx[1]), gc->range->upper[1]);
  long loc = gkyl_range_idx(gc->range, iter.idx);
  const double *coeffs = gkyl_array_cfetch(gc->bmagdg,loc);




  double xc[2];
  gkyl_rect_grid_cell_center(gc->grid, iter.idx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->grid->dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->grid->dx[1]*0.5);
  fout[0] = gc->basis->eval_expand(xy, coeffs);
}



void
gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, const struct gkyl_range *crange, const struct gkyl_range *crange_ext,
     const struct gkyl_range *prange, const struct gkyl_range *prange_ext, struct gkyl_array *psidg, struct gkyl_array *psibyrdg, struct gkyl_array *psibyr2dg, struct gkyl_array *bphidg, struct gkyl_array* bmag_compdg, struct gkyl_array* mapc2p)
{
  //First Stage is use psidg, psibydg, psibyr2dg to create bmagdg = B(R,Z) continuous
  struct gkyl_array* bmagdg = gkyl_array_new(GKYL_DOUBLE, up->pbasis->num_basis, prange_ext->volume);
  //make psibyr have LIRBT for recovery
  const double **psibyrall = gkyl_malloc(5*sizeof(double*));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, prange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(prange, iter.idx);
    double *bmag_i = gkyl_array_fetch(bmagdg, loc);
    double *psibyr2_i= gkyl_array_fetch(psibyr2dg, loc);
    double *bphi_i= gkyl_array_fetch(bphidg, loc);
    psibyrall[0] = gkyl_array_cfetch(psibyrdg,loc);
    int count = 1;
    int idx_temp[2] = {iter.idx[0], iter.idx[1]};
    for(int i = 0; i<2; i++){
      idx_temp[0] = iter.idx[0];
      idx_temp[1] = iter.idx[1];
      for(int j = -1; j<3; j+=2){
        idx_temp[i] = iter.idx[i] + j;
        loc = gkyl_range_idx(prange, idx_temp);
        psibyrall[count] = gkyl_array_cfetch(psibyrdg, loc);
        count = count+1;
      }
    }
    double scale_factorR = 2.0/(up->pgrid->dx[0]);
    double scale_factorZ = 2.0/(up->pgrid->dx[1]);
    up->kernel(psibyrall,psibyr2_i, bphi_i, bmag_i, scale_factorR, scale_factorZ);
  }

  //Second stage is to convert bmag into computational coordinates
  struct bmag_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->grid = up->pgrid;
  ctx->cgrid = up->cgrid;
  ctx->range = prange;
  ctx->crange = crange;
  ctx->bmagdg = bmagdg;
  ctx->basis = up->pbasis;
  ctx->cbasis = up->cbasis;
  ctx->app = up->app;
  ctx->ginp = up->ginp;
  ctx->mapc2p = mapc2p;
  gkyl_eval_on_nodes *eval_bmag_comp = gkyl_eval_on_nodes_new(up->cgrid, up->cbasis, 1, bmag_comp, ctx);
  gkyl_eval_on_nodes_advance(eval_bmag_comp, 0.0, crange, bmag_compdg); //on ghosts with ext_range
}

void
gkyl_calc_bmag_release(gkyl_calc_bmag* up)
{
  gkyl_free(up);
}
