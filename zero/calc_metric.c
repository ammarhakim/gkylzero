#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, const int *bcs, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));
  up->cbasis = cbasis;
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  up->bcs = gkyl_malloc((up->cdim)*sizeof(int));
  up->geo_bcs = gkyl_malloc((up->cdim)*sizeof(int));
  for(int i=0; i<up->cdim; i++){
    up->bcs[i] = bcs[i];
    up->geo_bcs[i] = bcs[i];
  }
  up->bcs[1] = 0; // To always use the periodic stencil (same as inner cell stencil when using ghost cells)
  //up->kernels = metric_choose_kernel(up->cdim, up->poly_order, up->bcs);

  return up;
}

static inline double calc_metric(double dxdz[3][3], int i, int j) 
{ double sum = 0;   for (int k=0; k<3; ++k) sum += dxdz[k][i-1]*dxdz[k][j-1]; return sum; } 

void gkyl_calc_metric_advance_nodal(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_xyz_nodal_fd, double *dzc){
  up->gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX }; // arrangement of physical coordinates  
  enum { X_IDX, Y_IDX, Zc_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
      for (int ip=nrange->lower[PH_IDX]; ip<=nrange->upper[PH_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PH_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_xyz_n = gkyl_array_cfetch(mc2p_xyz_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3];
              dxdz[0][0] = -(mc2p_xyz_n[3 +X_IDX] -   mc2p_xyz_n[6+X_IDX])/2/dzc[0];
              dxdz[0][1] = -(mc2p_xyz_n[9 +X_IDX] -  mc2p_xyz_n[12+X_IDX])/2/dzc[1];
              dxdz[0][2] = -(mc2p_xyz_n[15 +X_IDX] - mc2p_xyz_n[18+X_IDX])/2/dzc[2];
              dxdz[1][0] = -(mc2p_xyz_n[3 +Y_IDX] -   mc2p_xyz_n[6+Y_IDX])/2/dzc[0];
              dxdz[1][1] = -(mc2p_xyz_n[9 +Y_IDX] -  mc2p_xyz_n[12+Y_IDX])/2/dzc[1];
              dxdz[1][2] = -(mc2p_xyz_n[15 +Y_IDX] - mc2p_xyz_n[18+Y_IDX])/2/dzc[2];
              dxdz[2][0] = -(mc2p_xyz_n[3 +Zc_IDX] -   mc2p_xyz_n[6+Zc_IDX])/2/dzc[0];
              dxdz[2][1] = -(mc2p_xyz_n[9 +Zc_IDX] -  mc2p_xyz_n[12+Zc_IDX])/2/dzc[1];
              dxdz[2][2] = -(mc2p_xyz_n[15 +Zc_IDX] - mc2p_xyz_n[18+Zc_IDX])/2/dzc[2];
              //printf("mc2p 2 2 = %g %g\n" ,mc2p_xyz_n[15 +Zc_IDX] , mc2p_xyz_n[18+Zc_IDX]);
              //printf("dxdz[2][2] = %g\n",dxdz[2][2]);
              double *gFld_n= gkyl_array_fetch(up->gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[0] = calc_metric(dxdz, 1, 1); 
              gFld_n[1] = calc_metric(dxdz, 1, 2); 
              gFld_n[2] = calc_metric(dxdz, 1, 3); 
              gFld_n[3] = calc_metric(dxdz, 2, 2); 
              gFld_n[4] = calc_metric(dxdz, 2, 3); 
              gFld_n[5] = calc_metric(dxdz, 3, 3); 
              //printf("g11 = %g\n",gFld_n[0]);
              //printf("g12 = %g\n",gFld_n[1]);
              //printf("g13 = %g\n",gFld_n[2]);
              //printf("g22 = %g\n",gFld_n[3]);
              //printf("g23 = %g\n",gFld_n[4]);
              //printf("g33 = %g\n",gFld_n[5]);
      }
    }
  }
}

void gkyl_calc_metric_advance_modal(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *gFld, struct gkyl_range *update_range){
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = 6;
  int num_basis = up->cbasis->num_basis;
  int cpoly_order = up->cbasis->poly_order;
  //initialize the nodes
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, up->cbasis->num_basis);
  up->cbasis->node_list(gkyl_array_fetch(nodes, 0));
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
     gkyl_rect_grid_cell_center(up->grid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {
      const double* temp  = gkyl_array_cfetch(nodes,i);
      for( int j = 0; j < up->grid->ndim; j++){
        if(cpoly_order==1){
          if(j==1 && up->geo_bcs[1]==1) // this conversion is valid if ghost cells are included
            nidx[j] = iter.idx[j] + (temp[j]+1)/2 ;
          else // otherwise this conversion is valid
            nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2)
          nidx[j] = 2*iter.idx[j] + (temp[j] + 1) ; // add another line for if other bcs are used. Done correctly in efit.c
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(gFld, lidx); // pointer to expansion in cell
    double fao[num_basis*num_ret_vals];
  
    for (int i=0; i<num_basis; ++i) {
      double* temp = gkyl_array_fetch(up->gFld_nodal, lin_nidx[i]);
      for (int j=0; j<num_ret_vals; ++j) {
        //printf("temp = %g\n", temp[j]);
        fao[i*num_ret_vals + j] = temp[j];
      }
    }

    for (int i=0; i<num_ret_vals; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fnodal[k] = fao[num_ret_vals*k+i];
      // transform to modal expansion
      up->cbasis->nodal_to_modal(fnodal, &arr_p[num_basis*i]);
    }
  }
}


//void
//gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange, struct gkyl_array *XYZ, struct gkyl_array *gFld)
//{
//  const double **xyz = gkyl_malloc((1+2*up->cdim)*sizeof(double*));
//  struct gkyl_range_iter iter;
//  gkyl_range_iter_init(&iter, crange);
//  while (gkyl_range_iter_next(&iter)) {
//    long loc = gkyl_range_idx(crange, iter.idx);
//    double *gij = gkyl_array_fetch(gFld, loc);
//    xyz[0] = gkyl_array_cfetch(XYZ,loc);
//    int count = 1;
//    int idx_temp[up->cdim];
//    for(int l = 0; l<up->cdim; l++){idx_temp[l] = iter.idx[l]; }
//    for(int i = 0; i<up->cdim; i++){
//      for(int l = 0; l<up->cdim; l++){idx_temp[l] = iter.idx[l]; }
//      for(int j = -1; j<3; j+=2){
//        idx_temp[i] = iter.idx[i] + j;
//        loc = gkyl_range_idx(crange, idx_temp);
//        xyz[count] = gkyl_array_cfetch(XYZ, loc);
//        count = count+1;
//      }
//    }
//    int linker_idx = idx_to_inloup_ker(up->cdim, up->num_cells, iter.idx);
//    //linker_idx = 0; // to always use two sided
//    up->kernels.kernels[linker_idx](xyz,gij);
//  }
//
//  double scale_factor[up->cdim * (up->cdim+1)/2];
//  int count = 0;
//  for(int i=0; i<up->cdim; i++){
//    for(int j=i; j<up->cdim; j++){
//      scale_factor[count] = 4.0/(up->grid->dx[i]*up->grid->dx[j]);
//      count = count+1;
//    }
//  }
//
//  gkyl_range_iter_init(&iter, crange);
//  while (gkyl_range_iter_next(&iter)) {
//    long loc = gkyl_range_idx(crange, iter.idx);
//    double *gij = gkyl_array_fetch(gFld, loc);
//    for(int i=0; i<up->cdim * (up->cdim+1)/2; i++){
//      double *gcomp = &gij[i*(up->cnum_basis)];
//      for(int j=0; j < (up->cnum_basis); j++){
//        gcomp[j] = gcomp[j]*scale_factor[i];
//      }
//    }
//  }
//}

void
gkyl_calc_metric_release(gkyl_calc_metric* up)
{
  gkyl_free(up->bcs);
  gkyl_free(up->geo_bcs);
  gkyl_free(up);
}
