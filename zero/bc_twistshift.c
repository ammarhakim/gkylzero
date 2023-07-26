#include "gkyl_mat.h"
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim,
  const struct gkyl_array *yshift, const int *ndonors, bool use_gpu) {

  // Allocate space for new updater.
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(struct gkyl_bc_twistshift));

  up->dir = dir;
  up->edge = edge;
  up->basis = basis;
  up->grid = grid;
  up->ndonors = ndonors;
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



void gkyl_bc_twistshift_mv(struct gkyl_bc_twistshift *up, struct gkyl_nmat *matsdo, struct gkyl_nmat *vecsdo, struct gkyl_array *ftar)
{

  
  //allocate target matrix and perform the multiply
  struct gkyl_nmat *vecstar = gkyl_nmat_new(vecsdo->num, vecsdo->nr, vecsdo->nc);
  gkyl_nmat_mv(1.0, 0.0, GKYL_NO_TRANS, matsdo, vecsdo, vecstar);

  // matsdo is an nmat of N_b x N_b matrices with \sum_{i}^{Nx} n_do[i] elements
  // vecsdo has the same number of elements (number of vectors in the nmat)
  // vecstar has that same number of elements
  

  // Create the deflated range (only need update on z ghost cells on one edge)
  //int remDir[] = {0, 0, 1}, locDir[] = {0, 0, lower[2]};
  int remDir[up->grid->ndim];
  for(int i=0;i<up->grid->ndim;i++)
    remDir[i]=0;
  remDir[3] = 1; // z will always be 3rd index. This is hardcoded, any way around this?
  int locDir[up->grid->ndim];
  for(int i=0;i<up->grid->ndim;i++)
    locDir[i]=0;
  if(up->edge == GKYL_LOWER_EDGE)
    locDir[3] = up->local_range_ext->lower[3];
  else if(up->edge == GKYL_UPPER_EDGE)
    locDir[3] = up->local_range_ext->upper[3];

  struct gkyl_range update_range;
  gkyl_range_deflate(&update_range, up->local_range_ext, remDir, locDir);


  // Now loop through the update range and fill ftar with the right values
  // vecstar is already filled at this point. There is a set of ndo[i] mats for each xcell
  // loop through f and 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  //create some iterating ints to avoid allocating inside loop
  size_t linidx_start = 0;
  int cellidx = 0;
  int nb = vecstar->nr; //number of basis elements
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&update_range, iter.idx);
    double *ftar_itr = gkyl_array_fetch(ftar, loc);
    cellidx = iter.idx[0];
    // based on cell index we can find indices of donor matrices
    linidx_start = 0;
    for(int i = 0; i < cellidx-1; i++){
      linidx_start += up->ndonors[i];
    }
    // loop through donors of this cell
    for(int i = 0; i < up->ndonors[cellidx]; i++){
      struct gkyl_mat temp = gkyl_nmat_get(vecstar,linidx_start+i);
      // loop through nb basis coeffs / matrix elements in each mat
      for(int n=0; n<nb; n++){
        ftar_itr[n] += gkyl_mat_get(&temp, n, 0);
      }
    }
  }
  gkyl_nmat_release(vecstar);


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
