#include "gkyl_mat.h"
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>
#include <gkyl_array.h>

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(int dir, int do_dir, int shift_dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim,
  const struct gkyl_array *yshift, const int *ndonors, const int *cells_do, bool use_gpu) {

  // Allocate space for new updater.
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(struct gkyl_bc_twistshift));

  up->dir = dir;
  up->do_dir = do_dir;
  up->shift_dir = shift_dir;
  up->edge = edge;
  up->basis = basis;
  up->grid = grid;
  up->ndonors = ndonors;
  up->cells_do = cells_do;
  up->use_gpu = use_gpu;
  up->local_range_ext = local_range_ext;

  // Choose the kernel that does the reflection/no reflection/partial
  // reflection.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
    gkyl_bc_twistshift_choose_kernel_cu(basis, cdim, up->kernels_cu);
  } else {
    gkyl_bc_twistshift_choose_kernel(basis, cdim, up->kernels);
    up->kernels_cu = up->kernels;
  }
#else
  gkyl_bc_twistshift_choose_kernels(basis, cdim, up->kernels);
  up->kernels_cu = up->kernels;
#endif

  return up;
}

void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
  double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
  
  size_t linidx =0;
  for(int i = 0; i<cellidx; i++)
    linidx += up->ndonors[i];
  
  linidx += doidx;
  struct gkyl_mat tsmat = gkyl_nmat_get(mats, linidx);
  up->kernels->xlimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, &tsmat);
}

void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
  double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
  size_t linidx =0;
  for(int i = 0; i<cellidx; i++)
    linidx += up->ndonors[i];
  
  linidx += doidx;
  struct gkyl_mat tsmat = gkyl_nmat_get(mats, linidx);
  up->kernels->ylimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, &tsmat);
}

void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
  size_t linidx =0;
  for(int i = 0; i<cellidx; i++)
    linidx += up->ndonors[i];
  
  linidx += doidx;
  struct gkyl_mat tsmat = gkyl_nmat_get(mats, linidx);
  up->kernels->fullcell(dyDo, yOff, ySh, &tsmat);
}



void gkyl_bc_twistshift_mv(struct gkyl_bc_twistshift *up, struct gkyl_nmat *matsdo, struct gkyl_array *fdo, struct gkyl_array *ftar)
{
//  //donor locations an int array
//  int *cells_do[ny*nx*max(ndonors)]
//  // psuedocode
//  loop over y(k):
//    clear vecs
//    loop over x(j):
//      loop over donors(i):
//        // need donor locations which are y indices (same x index)
//        cells_do[k*nx*max(ndonors) + j*max(ndonors) + i] //should have a y coordinate of the donor of this x and this y
//        populate donor vectors -> vecs //
//    do mv(mats, vecs) // vecs is an nmat of total donors
//    loop over x:
//      distf array <- sum_over_donors(vecs)
//
  
  struct gkyl_nmat *vecstar = gkyl_nmat_new(matsdo->num, matsdo->nr, 1);
  struct gkyl_nmat *vecsdo = gkyl_nmat_new(matsdo->num, matsdo->nr, 1);

  // loop over y and call this function
  gkyl_nmat_mv(1.0, 0.0, GKYL_NO_TRANS, matsdo, vecsdo, vecstar);

  // matsdo is an nmat of N_b x N_b matrices with \sum_{i}^{Nx} n_do[i] elements
  // vecsdo has the same number of elements (number of vectors in the nmat)
  // vecstar has that same number of elements
  

  // Create the deflated range only need to loop over the shift dir and extraneous dirs (vpar,mu)
  // only need update on z ghost cells on one edge
  int remDir[up->grid->ndim];
  for(int i=0;i<up->grid->ndim;i++)
    remDir[i]=0; // keep all dirs by default
  remDir[up->dir] = 1; // remove the bc dir (z)
  remDir[up->do_dir] = 1; // remove the dir which the donors come from (x)
  // choose locations in other dirs
  int locDir[up->grid->ndim];
  for(int i=0;i<up->grid->ndim;i++)
    locDir[i]=up->grid->lower[i];
  // set z location to upper or lower
  if(up->edge == GKYL_LOWER_EDGE){
    locDir[up->dir] = up->local_range_ext->lower[up->dir];
  }
  else if(up->edge == GKYL_UPPER_EDGE){
    locDir[up->dir] = up->local_range_ext->upper[up->dir];
  }

  struct gkyl_range update_range;
  gkyl_range_deflate(&update_range, up->local_range_ext, remDir, locDir);

  // Setup for deflated range for looping over x
  int remDir_do[up->local_range_ext->ndim];
  for(int i=0;i<up->local_range_ext->ndim;i++)
    remDir_do[i]=1;
  remDir_do[0] = 0; // keep x only
  int locDir_do[up->local_range_ext->ndim];
  for(int i=0;i<up->local_range_ext->ndim;i++)
    locDir_do[i]=0;
  //decide locations later inside of loop
  if(up->edge == GKYL_LOWER_EDGE){
    locDir_do[up->dir] = up->local_range_ext->lower[up->dir];
  }
  else if(up->edge == GKYL_UPPER_EDGE){
    locDir_do[up->dir] = up->local_range_ext->upper[up->dir];
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  struct gkyl_range xrange;
  struct gkyl_range_iter iterx;
  int lin_idx = 0; //used to access the y donor indices
  int do_idx[up->local_range_ext->ndim];
  int tar_idx[up->local_range_ext->ndim];

  while (gkyl_range_iter_next(&iter)) {
    printf("y index is %d\n", iter.idx[0]);
    locDir_do[up->do_dir] = 0;
    locDir_do[up->shift_dir] = iter.idx[0];
    gkyl_range_deflate(&xrange, up->local_range_ext, remDir_do, locDir_do);
    gkyl_range_iter_init(&iterx, &xrange);
    while (gkyl_range_iter_next(&iterx)) {
      printf("in the x loop, iterx.idx[0] = %d\n", iterx.idx[0]);
      printf("lmax will be = %d\n", up->ndonors[iterx.idx[0]]);
      for(int i = 0; i < up->ndonors[iterx.idx[0]];i++){
        printf(" lin index = %d\n",lin_idx);
        for(int i = 0; i<up->local_range_ext->ndim; i++){ // set the extra indices of vpar and mu
          if(remDir[i] == 0 && i != up->shift_dir){
            do_idx[i] = iter.idx[i - (up->local_range_ext->ndim - update_range.ndim)];
          }
        }
        do_idx[up->do_dir] = iterx.idx[0];
        do_idx[up->shift_dir] = up->cells_do[lin_idx];
        do_idx[up->dir] = locDir_do[up->dir];
        for(int idx = 0; idx<up->local_range_ext->ndim; idx++)
          printf(" do_idx[%d] = %d ", idx, do_idx[idx]);
        printf("\n");

        struct gkyl_mat gkyl_mat_itr = gkyl_nmat_get(vecsdo, lin_idx);
        long loc = gkyl_range_idx(up->local_range_ext, do_idx);
        printf(" do loc = %ld\n",loc);
        // are we at the right index?
        int idx_check[up->local_range_ext->ndim];
        gkyl_range_inv_idx(up->local_range_ext, loc, idx_check);
        for(int ic = 0; ic<up->local_range_ext->ndim; ic++)
          printf(" idx_check[%d] = %d ", ic, idx_check[ic]);
        printf("\n");
        // done checking
        double *fdo_itr = gkyl_array_fetch(fdo, loc);
        for(int ib = 0; ib < vecsdo->nr; ib++){
          gkyl_mat_set(&gkyl_mat_itr, ib, 0, fdo_itr[ib]);
        }
        lin_idx += 1; 
      }
    }
    // now we have populated vecsdo. Do the multiply
    gkyl_nmat_mv(1.0, 0.0, GKYL_NO_TRANS, matsdo, vecsdo, vecstar);
    // now do the second loop over x to populate distf
    gkyl_range_deflate(&xrange, up->local_range_ext, remDir_do, locDir_do);
    gkyl_range_iter_init(&iterx, &xrange);
    int linidx_start = 0;
    while (gkyl_range_iter_next(&iterx)) { // this is a loop over x
      printf("in the 2nd x loop, iterx.idx[0] = %d\n", iterx.idx[0]);
      for(int i = 0; i<up->local_range_ext->ndim; i++){
        if(remDir[i] == 0 && i != up->shift_dir){
          tar_idx[i] = iter.idx[i - (up->local_range_ext->ndim - update_range.ndim)];
        }
      }
      tar_idx[up->do_dir] = iterx.idx[0];
      tar_idx[up->shift_dir] = iter.idx[0]; // yindex
      tar_idx[up->dir] = locDir_do[up->dir];
      linidx_start = 0;
      for(int i = 0; i < iterx.idx[0] - 1 ; i++){
        linidx_start += up->ndonors[i];
      }

      long loc = gkyl_range_idx(up->local_range_ext, tar_idx);
      // are we at the right index?
      for(int idx = 0; idx<up->local_range_ext->ndim; idx++)
        printf(" tar_idx[%d] = %d ", idx, tar_idx[idx]);
      printf("\n");
      int idx_check[up->local_range_ext->ndim];
      gkyl_range_inv_idx(up->local_range_ext, loc, idx_check);
      for(int ic = 0; ic<up->local_range_ext->ndim; ic++)
        printf(" tar idx_check[%d] = %d ", ic, idx_check[ic]);
      printf("\n");
      // done checking
      double *ftar_itr = gkyl_array_fetch(ftar, loc);
      for(int i = 0; i < up->ndonors[iterx.idx[0]]; i++){
        struct gkyl_mat temp = gkyl_nmat_get(vecstar,linidx_start+i);
        for(int n=0; n<matsdo->nr; n++){
          ftar_itr[n] += gkyl_mat_get(&temp, n, 0);
        }
      }
    }
  }

  gkyl_nmat_release(vecstar);
  gkyl_nmat_release(vecsdo);
}

void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up) {
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels_cu);
#endif
  gkyl_free(up->kernels);
  gkyl_free(up);
}
