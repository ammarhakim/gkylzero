#include <gkyl_calc_bmag.h>
#include <gkyl_calc_bmag_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_bmag*
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_basis *fbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, const struct gkyl_rect_grid *fgrid, 
  double sibry, bool use_gpu)
{
  gkyl_calc_bmag *up = gkyl_malloc(sizeof(gkyl_calc_bmag));
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cgrid = cgrid;
  up->pgrid = pgrid;
  up->use_gpu = use_gpu;
  up->kernel = bmag_choose_kernel(up->pbasis->ndim, up->pbasis->b_type, up->pbasis->poly_order);
  up->fbasis = fbasis;
  up->fgrid = fgrid;
  up->sibry = sibry;
  return up;
}

static inline void bmag_comp(double t, const double *xn, double *fout, void *ctx)
{
  struct bmag_ctx *gc = (struct bmag_ctx*) ctx;
  double RZPHI[gc->cgrid->ndim];

  int cidx[GKYL_MAX_CDIM];
  for(int i = 0; i < gc->cgrid->ndim; i++){
    int idxtemp = gc->crange_global->lower[i] + (int) floor((xn[i] - (gc->cgrid->lower[i]) )/gc->cgrid->dx[i]);
    idxtemp = GKYL_MIN2(idxtemp, gc->crange->upper[i]);
    idxtemp = GKYL_MAX2(idxtemp, gc->crange->lower[i]);
    cidx[i] = idxtemp;
  }

  long lidx = gkyl_range_idx(gc->crange, cidx);
  const double *mcoeffs = gkyl_array_cfetch(gc->mapc2p, lidx);
  
  double cxc[gc->cgrid->ndim];
  double xyz[gc->cgrid->ndim];
  gkyl_rect_grid_cell_center(gc->cgrid, cidx, cxc);
  for(int i = 0; i < gc->cgrid->ndim; i++)
    xyz[i] = (xn[i]-cxc[i])/(gc->cgrid->dx[i]*0.5);
  for(int i = 0; i < gc->cgrid->ndim; i++){
    RZPHI[i] = gc->cbasis->eval_expand(xyz, &mcoeffs[i*gc->cbasis->num_basis]);
  }

  double R = RZPHI[0];
  double Z = RZPHI[1];

  int rzidx[2];
  int idxtemp = gc->range->lower[0] + (int) floor((R - gc->grid->lower[0])/gc->grid->dx[0]);
  idxtemp = GKYL_MIN2(idxtemp, gc->range->upper[0]);
  idxtemp = GKYL_MAX2(idxtemp, gc->range->lower[0]);
  rzidx[0] = idxtemp;
  idxtemp = gc->range->lower[1] + (int) floor((Z - gc->grid->lower[1])/gc->grid->dx[1]);
  idxtemp = GKYL_MIN2(idxtemp, gc->range->upper[1]);
  idxtemp = GKYL_MAX2(idxtemp, gc->range->lower[1]);
  rzidx[1] = idxtemp;


  long loc = gkyl_range_idx(gc->range, rzidx);
  const double *coeffs = gkyl_array_cfetch(gc->bmagdg,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(gc->grid, rzidx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->grid->dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->grid->dx[1]*0.5);
  fout[0] = gc->basis->eval_expand(xy, coeffs);
}

static inline void bphi_RZ(double t, const double *xn, double *fout, void *ctx){

  struct fpol_ctx *gc = (struct fpol_ctx*) ctx;
  struct gkyl_range_iter iter;
  double psi;

  // first find the RZ index of where we are
  struct gkyl_range_iter rziter;
  gkyl_range_iter_init(&iter, gc->rzrange);
  for(int i = 0; i < gc->rzgrid->ndim; i++){
    rziter.idx[i] = fmin(gc->rzrange->lower[i] + (int) floor((xn[i] - (gc->rzgrid->lower[i]) )/gc->rzgrid->dx[i]), gc->rzrange->upper[i]);
  }

  // now get the coeffs of psi at that RZ location
  long lidx = gkyl_range_idx(gc->rzrange, rziter.idx);
  const double *psicoeffs = gkyl_array_cfetch(gc->psidg, lidx);
  
  // Now get psi at that location
  double rzxc[gc->rzgrid->ndim];
  double rz[gc->rzgrid->ndim];
  gkyl_rect_grid_cell_center(gc->rzgrid, rziter.idx, rzxc);
  for(int i = 0; i < gc->rzgrid->ndim; i++)
    rz[i] = (xn[i]-rzxc[i])/(gc->rzgrid->dx[i]*0.5);
  psi = gc->rzbasis->eval_expand(rz, &psicoeffs[gc->rzbasis->num_basis]);

  if ( (psi < gc->grid->lower[0]) || (psi > gc->grid->upper[0]) ) // F = F(psi_sep) in the SOL. Works regardless of psi convention
    psi = gc->sibry;

  // now find psi cell this lies in and get coeffs for fpol
  gkyl_range_iter_init(&iter, gc->range);
  iter.idx[0] = fmin(gc->range->lower[0] + (int) floor((psi - gc->grid->lower[0])/gc->grid->dx[0]), gc->range->upper[0]);
  long loc = gkyl_range_idx(gc->range, iter.idx);
  const double *coeffs = gkyl_array_cfetch(gc->fpoldg,loc);

  double fxc[1];
  gkyl_rect_grid_cell_center(gc->grid, iter.idx, fxc);
  double fx[1];
  fx[0] = (psi-fxc[0])/(gc->grid->dx[0]*0.5);
  double R = xn[0];
  fout[0] = gc->basis->eval_expand(fx, coeffs)/R;
}

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, 
  const struct gkyl_range *crange, const struct gkyl_range *crange_ext, const struct gkyl_range *crange_global, 
  const struct gkyl_range *prange, const struct gkyl_range *prange_ext, 
  const struct gkyl_range *frange, const struct gkyl_range* frange_ext, 
  const struct gkyl_array *psidg, const struct gkyl_array *psibyrdg, const struct gkyl_array *psibyr2dg, 
  struct gkyl_array* bmag_compdg, const struct gkyl_array* fpoldg, struct gkyl_array* mapc2p, bool calc_bphi)
{
  // 0th stage is to calculate bphi from fpol on the RZ grid from its representation on the flux grid
  // We will do this with an eval on nodes similar to what is done in bmag comp
  // We will pass an RZ coordinate to the function
  // it will evaluate psi at that coordinate
  // we find which cell that lies in for fpol and the normalized coord
  // then evaluate

  struct fpol_ctx *fctx = gkyl_malloc(sizeof(*fctx));
  fctx->grid = up->fgrid;
  fctx->rzgrid = up->pgrid;
  fctx->range = frange;
  fctx->rzrange = prange;
  fctx->fpoldg= fpoldg;
  fctx->psidg = psidg;
  fctx->basis = up->fbasis;
  fctx->rzbasis = up->pbasis;
  fctx->sibry = up->sibry;

  struct gkyl_array *bphirz = gkyl_array_new(GKYL_DOUBLE, up->pbasis->num_basis, prange_ext->volume);
  if (calc_bphi){
    gkyl_eval_on_nodes *eval_fpol_RZ = gkyl_eval_on_nodes_new(up->pgrid, up->pbasis, 1, bphi_RZ, fctx);
    gkyl_eval_on_nodes_advance(eval_fpol_RZ, 0.0, prange, bphirz); //on ghosts with ext_range
    gkyl_eval_on_nodes_release(eval_fpol_RZ);
  }

  struct gkyl_array *bmagrz = gkyl_array_new(GKYL_DOUBLE, up->pbasis->num_basis, prange_ext->volume);
  //First Stage is use psidg, psibydg, psibyr2dg to create bmagdg = B(R,Z) continuous
  //struct gkyl_array* bmagdg = gkyl_array_new(GKYL_DOUBLE, up->pbasis->num_basis, prange_ext->volume);
  //make psibyr have LIRBT for recovery
  const double **psibyrall = gkyl_malloc(5*sizeof(double*));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, prange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(prange, iter.idx);
    double *bmag_i = gkyl_array_fetch(bmagrz, loc);
    const double *psibyr2_i= gkyl_array_cfetch(psibyr2dg, loc);
    double *bphi_i= gkyl_array_fetch(bphirz, loc);
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
    up->kernel(psibyrall, psibyr2_i, bphi_i, bmag_i, scale_factorR, scale_factorZ);
  }

  // Second stage is to convert bmag into computational coordinates
  struct bmag_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->grid = up->pgrid;
  ctx->cgrid = up->cgrid;
  ctx->range = prange;
  ctx->crange = crange;
  ctx->crange_global = crange_global;
  ctx->bmagdg = bmagrz;
  ctx->basis = up->pbasis;
  ctx->cbasis = up->cbasis;
  ctx->mapc2p = mapc2p;
  gkyl_eval_on_nodes *eval_bmag_comp = gkyl_eval_on_nodes_new(up->cgrid, up->cbasis, 1, bmag_comp, ctx);
  gkyl_eval_on_nodes_advance(eval_bmag_comp, 0.0, crange, bmag_compdg); //on ghosts with ext_range
  gkyl_eval_on_nodes_release(eval_bmag_comp);

  // Free allocated memory
  gkyl_free(fctx);
  gkyl_array_release(bphirz);
  gkyl_array_release(bmagrz);
  gkyl_free(psibyrall);
  gkyl_free(ctx);
}

void
gkyl_calc_bmag_release(gkyl_calc_bmag* up)
{
  gkyl_free(up);
}
