#include <gkyl_calc_cart_bmag.h>
#include <gkyl_calc_cart_bmag_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

gkyl_calc_cart_bmag*
gkyl_calc_cart_bmag_new(const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *local, const struct gkyl_range *local_ext, bool use_gpu)
{
  gkyl_calc_cart_bmag *up = gkyl_malloc(sizeof(gkyl_calc_cart_bmag));
  up->basis = basis;
  up->local = *local;
  up->local_ext = *local_ext;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->kernel = cart_bmag_choose_kernel(up->basis->ndim, up->basis->b_type, up->basis->poly_order);

  up->bmag_rz = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, local_ext->volume);
  up->br_rz = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, local_ext->volume);
  up->bz_rz = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, local_ext->volume);
  return up;
}

void
gkyl_calc_cart_bmag_advance(const gkyl_calc_cart_bmag *up, const struct gkyl_array *psidg, const struct gkyl_array *psibyrdg, const struct gkyl_array *psibyr2dg)
{

  //First Stage is use psidg, psibydg, psibyr2dg to create bmagdg = B(R,Z) continuous
  //make psibyr have LIRBT for recovery
  double scale_factorR = 2.0/(up->grid->dx[0]);
  double scale_factorZ = 2.0/(up->grid->dx[1]);
  const double **psibyrall = gkyl_malloc(5*sizeof(double*));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->local, iter.idx);
    double *bmag_i = gkyl_array_fetch(up->bmag_rz, loc);
    double *br_i = gkyl_array_fetch(up->br_rz, loc);
    double *bz_i = gkyl_array_fetch(up->bz_rz, loc);
    const double *psibyr2_i= gkyl_array_cfetch(psibyr2dg, loc);
    psibyrall[0] = gkyl_array_cfetch(psibyrdg,loc);
    int count = 1;
    int idx_temp[2] = {iter.idx[0], iter.idx[1]};
    for(int i = 0; i<2; i++){
      idx_temp[0] = iter.idx[0];
      idx_temp[1] = iter.idx[1];
      for(int j = -1; j<3; j+=2){
        idx_temp[i] = iter.idx[i] + j;
        loc = gkyl_range_idx(&up->local, idx_temp);
        psibyrall[count] = gkyl_array_cfetch(psibyrdg, loc);
        count = count+1;
      }
    }
    up->kernel(psibyrall, psibyr2_i, bmag_i, br_i, bz_i, scale_factorR, scale_factorZ);
  }



  // Free allocated memory
  gkyl_free(psibyrall);
}


void
gkyl_eval_cart_bmag(const gkyl_calc_cart_bmag *up, double xcart[3], double bcart[3])
{
  double R = sqrt(xcart[0]*xcart[0] + xcart[1]*xcart[1]);
  double Z = xcart[2];
  double phi = atan(xcart[1]/xcart[0]);


  //First find which rz cell we are in to get B, B_R B_Z
  int rz_idx[2];

  int idxtemp = up->local.lower[0] + (int) floor((R - up->grid->lower[0])/up->grid->dx[0]);
  idxtemp = GKYL_MIN2(idxtemp, up->local.upper[0]);
  idxtemp = GKYL_MAX2(idxtemp, up->local.lower[0]);
  rz_idx[0] = idxtemp;

  idxtemp = up->local.lower[1] + (int) floor((Z - up->grid->lower[1])/up->grid->dx[1]);
  idxtemp = GKYL_MIN2(idxtemp, up->local.upper[1]);
  idxtemp = GKYL_MAX2(idxtemp, up->local.lower[1]);
  rz_idx[1] = idxtemp;


  long loc = gkyl_range_idx(&up->local, rz_idx);
  const double *coeffs_bmag = gkyl_array_cfetch(up->bmag_rz,loc);
  const double *coeffs_br = gkyl_array_cfetch(up->br_rz,loc);
  const double *coeffs_bz = gkyl_array_cfetch(up->bz_rz,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(up->grid, rz_idx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(up->grid->dx[0]*0.5);
  xy[1] = (Z-xc[1])/(up->grid->dx[1]*0.5);
  double bmag = up->basis->eval_expand(xy, coeffs_bmag);
  double br = up->basis->eval_expand(xy, coeffs_br);
  double bz = up->basis->eval_expand(xy, coeffs_bz);

  bcart[0] = br*cos(phi);
  bcart[1] = br*sin(phi);
  bcart[2] = bz;
}

void
gkyl_calc_cart_bmag_release(gkyl_calc_cart_bmag* up)
{
  gkyl_array_release(up->bmag_rz);
  gkyl_array_release(up->br_rz);
  gkyl_array_release(up->bz_rz);
  gkyl_free(up);
}
