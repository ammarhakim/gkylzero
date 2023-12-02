#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));
  up->cbasis = cbasis;
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  return up;
}

static inline double calc_metric(double dxdz[3][3], int i, int j) 
{ double sum = 0;   for (int k=0; k<3; ++k) sum += dxdz[k][i-1]*dxdz[k][j-1]; return sum; } 

void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, double *dzc, struct gkyl_array *gFld, const struct gkyl_range *update_range){
  up->gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
      for (int ip=nrange->lower[PH_IDX]; ip<=nrange->upper[PH_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PH_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3];

              if(ip==nrange->lower[PH_IDX]){
                dxdz[0][0] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[6+X_IDX] - mc2p_n[12+X_IDX] )/dzc[0]/2;
                dxdz[1][0] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[6+Y_IDX] - mc2p_n[12+Y_IDX] )/dzc[0]/2;
                dxdz[2][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
              }
              else if(ip==nrange->upper[PH_IDX]){
                dxdz[0][0] = (3*mc2p_n[X_IDX] - 4*mc2p_n[3+X_IDX] + mc2p_n[9+X_IDX] )/dzc[0]/2;
                dxdz[1][0] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[3+Y_IDX] + mc2p_n[9+Y_IDX] )/dzc[0]/2;
                dxdz[2][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
              }
              else{
                dxdz[0][0] = -(mc2p_n[3 +X_IDX] -   mc2p_n[6+X_IDX])/2/dzc[0];
                dxdz[1][0] = -(mc2p_n[3 +Y_IDX] -   mc2p_n[6+Y_IDX])/2/dzc[0];
                dxdz[2][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
              }


              if(ia==nrange->lower[AL_IDX]){
                dxdz[0][1] = (-3*mc2p_n[X_IDX] +  4*mc2p_n[18+X_IDX] -  mc2p_n[24+X_IDX])/dzc[1]/2;
                dxdz[1][1] = (-3*mc2p_n[Y_IDX] +  4*mc2p_n[18+Y_IDX] -  mc2p_n[24+Y_IDX])/dzc[1]/2;
                dxdz[2][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
              }
              else if(ia==nrange->upper[AL_IDX]){
                dxdz[0][1] = (3*mc2p_n[X_IDX] -  4*mc2p_n[15+X_IDX] +  mc2p_n[21+X_IDX] )/dzc[1]/2;
                dxdz[1][1] = (3*mc2p_n[Y_IDX] -  4*mc2p_n[15+Y_IDX] +  mc2p_n[21+Y_IDX] )/dzc[1]/2;
                dxdz[2][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
              }
              else{
                dxdz[0][1] = -(mc2p_n[15 +X_IDX] -  mc2p_n[18 +X_IDX])/2/dzc[1];
                dxdz[1][1] = -(mc2p_n[15 +Y_IDX] -  mc2p_n[18 +Y_IDX])/2/dzc[1];
                dxdz[2][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
              }

              if(it==nrange->lower[TH_IDX]){
                dxdz[0][2] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[30+X_IDX] - mc2p_n[36+X_IDX])/dzc[2]/2;
                dxdz[1][2] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[30+Y_IDX] - mc2p_n[36+Y_IDX])/dzc[2]/2;
                dxdz[2][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
              }
              else if(it==nrange->upper[TH_IDX]){
                dxdz[0][2] = (3*mc2p_n[X_IDX] - 4*mc2p_n[27+X_IDX] + mc2p_n[33+X_IDX] )/dzc[2]/2;
                dxdz[1][2] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[27+Y_IDX] + mc2p_n[33+Y_IDX] )/dzc[2]/2;
                dxdz[2][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
              }
              else{
                dxdz[0][2] = -(mc2p_n[27 +X_IDX] - mc2p_n[30 +X_IDX])/2/dzc[2];
                dxdz[1][2] = -(mc2p_n[27 +Y_IDX] - mc2p_n[30 +Y_IDX])/2/dzc[2];
                dxdz[2][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
              }

              double *gFld_n= gkyl_array_fetch(up->gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[0] = calc_metric(dxdz, 1, 1); 
              gFld_n[1] = calc_metric(dxdz, 1, 2); 
              gFld_n[2] = calc_metric(dxdz, 1, 3); 
              gFld_n[3] = calc_metric(dxdz, 2, 2); 
              gFld_n[4] = calc_metric(dxdz, 2, 3); 
              gFld_n[5] = calc_metric(dxdz, 3, 3); 
      }
    }
  }
  gkyl_nodal_ops_n2m(up->cbasis, up->grid, nrange, update_range, 6, up->gFld_nodal, gFld);
}

void
gkyl_calc_metric_release(gkyl_calc_metric* up)
{
  gkyl_free(up);
}
