#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>
#include <gkyl_util.h>
#include <gkyl_eval_on_nodes.h>

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(const struct gkyl_bc_twistshift_inp *inp)
{

  // Allocate space for new updater.
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(struct gkyl_bc_twistshift));

  up->dir = inp->dir;
  up->shift_dir = inp->shift_dir;
  up->shear_dir = inp->shear_dir;
  up->edge = inp->edge;
  up->basis = inp->basis;
  up->local_ext_r = inp->local_ext_r;
  up->local_r = inp->local_r;
  up->grid = inp->grid;
  up->use_gpu = inp->use_gpu;

  // Assume the poly order of the DG shift is the same as that of the field,
  // unless requested otherwise.
  up->shift_poly_order = inp->basis.poly_order;
  if (inp->shift_poly_order)
    up->shift_poly_order = inp->shift_poly_order;

  double lo1d[1], up1d[1];  int cells1d[1];

  // Create 1D grid and range in the diretion of the shear.
  gkyl_range_init(&up->shear_r, 1, (int[]) {inp->local_r.lower[inp->shear_dir]},
                                   (int[]) {inp->local_r.upper[inp->shear_dir]});
  lo1d[0] = inp->grid.lower[up->shear_dir];
  up1d[0] = inp->grid.upper[up->shear_dir];
  cells1d[0] = inp->grid.cells[up->shear_dir];
  gkyl_rect_grid_init(&up->shear_grid, 1, lo1d, up1d, cells1d);

  // Create 1D grid and range in the diretion of the shift.
  gkyl_range_init(&up->shift_r, 1, &(int) {inp->local_r.lower[inp->shift_dir]},
                                   &(int) {inp->local_r.upper[inp->shift_dir]});
  lo1d[0] = inp->grid.lower[up->shift_dir];
  up1d[0] = inp->grid.upper[up->shift_dir];
  cells1d[0] = inp->grid.cells[up->shift_dir];
  gkyl_rect_grid_init(&up->shift_grid, 1, lo1d, up1d, cells1d);

  // Create 2D grid (and range) the twist-shift takes place in.
  int dimlo, dimup;
  if (up->shift_dir < up->shear_dir) {
    dimlo = up->shift_dir;
    dimup = up->shear_dir;
    up->shift_dir_in_ts_grid = 0;
    up->shear_dir_in_ts_grid = 1;
  }
  else {
    dimlo = up->shear_dir;
    dimup = up->shift_dir;
    up->shift_dir_in_ts_grid = 1;
    up->shear_dir_in_ts_grid = 0;
  }
  gkyl_range_init(&up->ts_r, 2, (int[]) {inp->local_r.lower[dimlo], inp->local_r.lower[dimup]},
                                (int[]) {inp->local_r.upper[dimlo], inp->local_r.upper[dimup]});
  double lo2d[] = {inp->grid.lower[dimlo], inp->grid.lower[dimup]};
  double up2d[] = {inp->grid.upper[dimlo], inp->grid.upper[dimup]};
  int cells2d[] = {inp->grid.cells[dimlo], inp->grid.cells[dimup]};
  gkyl_rect_grid_init(&up->ts_grid, 2, lo2d, up2d, cells2d);

  // Project the shift onto the shift basis.
  gkyl_cart_modal_serendip(&up->shift_b, 1, up->shift_poly_order);
  up->shift = gkyl_array_new(GKYL_DOUBLE, up->shift_b.num_basis, up->shear_r.volume);
  gkyl_eval_on_nodes *evup = gkyl_eval_on_nodes_new(&up->shear_grid, &up->shift_b, 1, inp->shift_func, inp->shift_func_ctx);
  gkyl_eval_on_nodes_advance(evup, 0.0, &up->shear_r, up->shift);
  gkyl_eval_on_nodes_release(evup);

  // Find the donor cells for each target cell.
  ts_find_donors(up);

//  // Choose the kernels that do the subcell and full cell integrals
//  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
//  gkyl_bc_twistshift_choose_kernels(inp->basis, inp->cdim, up->kernels);
//
//  // Allocate matrices and vectors, initialize if needed
//  int total_donor_mats = 0;
//  int ndonors_cum[up->grid->cells[up->shear_dir]+1];
//  for(int i = 0; i<grid->cells[0]; i++){
//    ndonors_cum[i] = total_donor_mats;
//    total_donor_mats += ndonors[i];
//  }
//  ndonors_cum[grid->cells[0]] = total_donor_mats;
//
//#ifdef GKYL_HAVE_CUDA
//  up->ndonors_cum_cu = gkyl_cu_malloc(sizeof(int)*(up->grid->cells[up->shear_dir]+1));
//  cudaMemcpy(up->ndonors_cum_cu, ndonors_cum, (up->grid->cells[up->shear_dir]+1)*sizeof(int),GKYL_CU_MEMCPY_H2D);
//#endif
//
//  up->matsdo_ho = gkyl_nmat_new(total_donor_mats, basis->num_basis, basis->num_basis);
//
//  for (size_t n=0; n<up->matsdo_ho->num; ++n) {
//    struct gkyl_mat m = gkyl_nmat_get(up->matsdo_ho, n);
//    for (size_t j=0; j<up->matsdo_ho->nc; ++j)
//      for (size_t i=0; i<up->matsdo_ho->nr; ++i)
//        gkyl_mat_set(&m, i, j, 0.0);
//  }
//#ifdef GKYL_HAVE_CUDA
//  up->matsdo = gkyl_nmat_cu_dev_new(total_donor_mats, basis->num_basis, basis->num_basis);
//  gkyl_nmat_copy(up->matsdo, up->matsdo_ho);
//  up->vecstar = gkyl_nmat_cu_dev_new(up->matsdo->num, up->matsdo->nr, 1);
//  up->vecsdo = gkyl_nmat_cu_dev_new(up->matsdo->num, up->matsdo->nr, 1);
//#else
//  up->matsdo = up->matsdo_ho;
//  up->vecstar = gkyl_nmat_new(up->matsdo->num, up->matsdo->nr, 1);
//  up->vecsdo = gkyl_nmat_new(up->matsdo->num, up->matsdo->nr, 1);
//#endif
//
//  // try to set up locs to parallelize
//  up->locs = gkyl_malloc(up->vecsdo->num*sizeof(long));
//  up->tar_locs = gkyl_malloc(up->grid->cells[up->shear_dir]*sizeof(long));
//
//#ifdef GKYL_HAVE_CUDA
//  up->locs_cu = gkyl_cu_malloc(up->vecsdo->num*sizeof(long));
//  up->tar_locs_cu = gkyl_cu_malloc(up->grid->cells[up->shear_dir]*sizeof(long));
//#endif
//
//
//  // Create necessary ranges
//  up->yrange = gkyl_malloc(sizeof(struct gkyl_range));
//  up->xrange = gkyl_malloc(sizeof(struct gkyl_range));
//  // Create the deflated range only need to loop over the shift dir and extraneous dirs (vpar,mu)
//  up->remDir = (int*) gkyl_malloc(sizeof(int) * up->local_r->ndim);
//  for(int i=0;i<up->local_r->ndim;i++)
//    up->remDir[i]=1; // remove all dirs by default
//  up->remDir[up->shift_dir] = 0; // keep only the donor dir 
//  for(int i = up->dir + 1; i<up->local_r->ndim; i++) // keep the extaneous dirs (vpar and mu) assumed to be after z
//    up->remDir[i] = 0;
//  // Setup for deflated range for looping over x
//  up->remDir_do = (int*) gkyl_malloc(sizeof(int) * up->local_r->ndim);
//  for(int i=0;i<up->local_r->ndim;i++)
//    up->remDir_do[i]=1; // remove all dirs y default
//  up->remDir_do[up->shear_dir] = 0; // keep x only
//  // Choose locations in other dirs. Some will be reset in loop
//  up->locDir = (int*) gkyl_malloc(sizeof(int) * up->local_r->ndim);
//  up->locDir_do = (int*) gkyl_malloc(sizeof(int) * up->local_r->ndim);
//  for(int i=0;i<up->local_r->ndim;i++){
//    up->locDir[i] = up->local_r->lower[i];
//    up->locDir_do[i]=up->local_r->lower[i];
//  }
//  if(up->local_r->ndim > 2){
//    if(up->edge == GKYL_LOWER_EDGE){
//      up->locDir[up->dir] = up->local_r->lower[up->dir];
//      up->locDir_do[up->dir] = up->local_r->upper[up->dir] - 1 ;
//    }
//    else if(up->edge == GKYL_UPPER_EDGE){
//      up->locDir[up->dir] = up->local_r->upper[up->dir];
//      up->locDir_do[up->dir] = up->local_r->lower[up->dir] + 1 ;
//    }
//  }
//  gkyl_range_deflate(up->yrange, up->local_r, up->remDir, up->locDir);


  return up;
}

//void gkyl_bc_twistshift_copy_matsdo(struct gkyl_bc_twistshift *up)
//{
//#ifdef GKYL_HAVE_CUDA
//  gkyl_nmat_copy(up->matsdo, up->matsdo_ho);
//#endif
//}
//
//void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
//{
//  struct gkyl_range_iter iter;
//  gkyl_range_iter_init(&iter, up->yrange);
//  struct gkyl_range_iter iterx;
//  int lin_idx = 0;
//  int do_idx[up->local_r->ndim];
//  int tar_idx[up->local_r->ndim];
//  int lin_vecdo_idx = 0;
//  int lin_tar_idx = 0;
//  int last_yidx = 0;
//  int dim_diff = up->local_r->ndim - up->yrange->ndim;
//  if(up->local_r->ndim > 2){
//    do_idx[up->dir] = up->locDir_do[up->dir];
//    tar_idx[up->dir] = up->locDir[up->dir];
//  }
//
//  while (gkyl_range_iter_next(&iter)) {
//    if(iter.idx[0] == last_yidx){
//      lin_idx = (iter.idx[0] - 1)*up->matsdo->num;
//    }
//    up->locDir_do[up->shift_dir] = iter.idx[0];
//    if(up->local_r->ndim>3){
//      for(int i = up->dir + 1; i < up->local_r->ndim; i++){
//        up->locDir_do[i] = iter.idx[i - dim_diff];
//        do_idx[i] = iter.idx[i - dim_diff];
//        tar_idx[i] = iter.idx[i - dim_diff];
//      }
//    }
//    gkyl_range_deflate(up->xrange, up->local_r, up->remDir_do, up->locDir_do);
//    gkyl_range_iter_init(&iterx, up->xrange);
//    lin_vecdo_idx = 0;
//    while (gkyl_range_iter_next(&iterx)) {
//      for(int i = 0; i < up->ndonors[iterx.idx[0]-1];i++){
//        do_idx[up->shear_dir] = iterx.idx[0];
//        do_idx[up->shift_dir] = up->cells_do[lin_idx];
//
//        struct gkyl_mat gkyl_mat_itr = gkyl_nmat_get(up->vecsdo, lin_vecdo_idx);
//        long loc = gkyl_range_idx(up->local_r, do_idx);
//        up->locs[lin_vecdo_idx] = loc;
//#ifndef GKYL_HAVE_CUDA
//        const double *fdo_itr = gkyl_array_cfetch(fdo, loc);
//        for(int ib = 0; ib < up->vecsdo->nr; ib++){
//          gkyl_mat_set(&gkyl_mat_itr, ib, 0, fdo_itr[ib]);
//        }
//#endif
//        lin_vecdo_idx += 1;
//        lin_idx += 1; 
//      }
//    }
//#ifdef GKYL_HAVE_CUDA
//    cudaMemcpy(up->locs_cu, up->locs, up->vecsdo->num*sizeof(long),GKYL_CU_MEMCPY_H2D);
//    gkyl_bc_twistshift_set_vecsdo_cu(fdo, up->locs_cu, up->vecsdo);
//#endif
//
//    gkyl_nmat_mv(1.0, 0.0, GKYL_NO_TRANS, up->matsdo, up->vecsdo, up->vecstar);
//
//    gkyl_range_deflate(up->xrange, up->local_r, up->remDir_do, up->locDir_do);
//    gkyl_range_iter_init(&iterx, up->xrange);
//    while (gkyl_range_iter_next(&iterx)) {
//      tar_idx[up->shear_dir] = iterx.idx[0];
//      tar_idx[up->shift_dir] = iter.idx[0];
//
//      long loc = gkyl_range_idx(up->local_r, tar_idx);
//      up->tar_locs[iterx.idx[0]-1] = loc;
//    }
//
//#ifdef GKYL_HAVE_CUDA
//    cudaMemcpy(up->tar_locs_cu, up->tar_locs, up->grid->cells[up->shear_dir]*sizeof(long),GKYL_CU_MEMCPY_H2D);
//    gkyl_bc_twistshift_clear_cu(fdo, up->tar_locs_cu, up->grid->cells[up->shear_dir]);
//#endif
//    gkyl_range_deflate(up->xrange, up->local_r, up->remDir_do, up->locDir_do);
//    gkyl_range_iter_init(&iterx, up->xrange);
//    lin_tar_idx = 0;
//
//#ifdef GKYL_HAVE_CUDA
//    gkyl_bc_twistshift_inc_cu(ftar, up->tar_locs_cu, up->grid->cells[up->shear_dir], up->vecstar, up->ndonors_cum_cu);
//#else
//    while (gkyl_range_iter_next(&iterx)) {
//      double *ftar_itr = gkyl_array_fetch(ftar, up->tar_locs[iterx.idx[0]-1]);
//      for(int n=0; n<ftar->ncomp; n++){
//        ftar_itr[n] = 0;
//      }
//      for(int i = 0; i < up->ndonors[iterx.idx[0]-1]; i++){
//        struct gkyl_mat temp = gkyl_nmat_get(up->vecstar,lin_tar_idx);
//        for(int n=0; n<up->matsdo->nr; n++){
//          ftar_itr[n] += gkyl_mat_get(&temp, n, 0);
//        }
//        lin_tar_idx +=1;
//      }
//    }
//#endif
//    last_yidx = iter.idx[0];
//  }
//}

void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up) {
  // Release memory associated with this updater.
  gkyl_free(up->num_do);
//  gkyl_free(up->shift_dir_idx_do);

  gkyl_free(up);

  // *****   IN PROGRESS   ****** 
//#ifdef GKYL_HAVE_CUDA
//  gkyl_nmat_release(up->matsdo_ho);
//  gkyl_cu_free(up->tar_locs_cu);
//  gkyl_cu_free(up->locs_cu);
//  gkyl_cu_free(up->ndonors_cum_cu);
//  if (up->use_gpu)
//    gkyl_cu_free(up->kernels_cu);
//#endif
//  gkyl_free(up->locs);
//  gkyl_free(up->tar_locs);
//  gkyl_free(up->kernels);
//  gkyl_free(up->yrange);
//  gkyl_free(up->xrange);
//  gkyl_free(up->remDir);
//  gkyl_free(up->locDir);
//  gkyl_free(up->remDir_do);
//  gkyl_free(up->locDir_do);
//  gkyl_nmat_release(up->matsdo);
//  gkyl_nmat_release(up->vecsdo);
//  gkyl_nmat_release(up->vecstar);
}
