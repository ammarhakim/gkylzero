#include <cstdio>
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <gkyl_bc_twistshift.h>
}

#include <cstdio>
#include <cassert>

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)


__global__ void
gkyl_bc_twistshift_inc_cu_kernel(struct gkyl_array* ftar, long* tar_locs, int num_tar_locs, struct gkyl_nmat* vecs_contribution, int* ndonors_cum, const struct gkyl_range local_range_update, int unique_donor_mats, const struct gkyl_rect_grid grid)
{
  int nmu = grid.cells[4];
  int nv = grid.cells[3];
  long lidx = START_ID/ftar->ncomp;
  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&local_range_update, lidx, idx);
  int x_idx = idx[3];
  int mu_idx = idx[2];
  int vpar_idx = idx[1];
  int y_idx = idx[0];
  int set_idx = START_ID%ftar->ncomp;
  int do_start = ndonors_cum[x_idx-1];
  int do_end = ndonors_cum[x_idx+1-1];
  if(lidx<num_tar_locs && set_idx<ftar->ncomp){
    double *ftar_itr = (double*) gkyl_array_fetch(ftar, tar_locs[lidx]);
    int vectar_start = (mu_idx-1)*unique_donor_mats + (vpar_idx-1)*nmu*unique_donor_mats + (y_idx-1)*nv*nmu*unique_donor_mats;
    for(int i = do_start; i < do_end; i++){ // needs to only loop over ndonors[i] elements
      int lin_vectar_idx = i + vectar_start;
      struct gkyl_mat temp = gkyl_nmat_get(vecs_contribution, lin_vectar_idx%unique_donor_mats);
      ftar_itr[set_idx] += gkyl_mat_get(&temp, set_idx, lin_vectar_idx/unique_donor_mats);
    }
  }
}

__global__ void
gkyl_bc_twistshift_clear_cu_kernel(struct gkyl_array* ftar, long* locs, int num_locs)
{
  int set_idx = START_ID%ftar->ncomp;
  int loc_idx = START_ID/ftar->ncomp;
  if(loc_idx < num_locs && set_idx<ftar->ncomp){
    double *ftar_itr = (double*) gkyl_array_fetch(ftar, locs[loc_idx]);
    ftar_itr[set_idx] = 0.0;
  }
}

__global__ void
gkyl_bc_twistshift_set_vecsdo_cu_kernel(const struct gkyl_array* fdo, long* locs, struct gkyl_nmat* vecsdo, int donor_factor, int unique_donor_mats)
{
  int vecdo_idx = START_ID/fdo->ncomp;
  int set_idx = START_ID%fdo->ncomp;
  if(vecdo_idx < vecsdo->num*donor_factor && set_idx<fdo->ncomp){
    const double *fdo_itr = (const double*) gkyl_array_cfetch(fdo, locs[vecdo_idx]);
    struct gkyl_mat temp = gkyl_nmat_get(vecsdo, vecdo_idx%unique_donor_mats);
    gkyl_mat_set(&temp, set_idx, vecdo_idx/unique_donor_mats, fdo_itr[set_idx]);
  }
}


void
gkyl_bc_twistshift_inc_cu(const struct gkyl_array* ftar, long* tar_locs, int num_tar_locs, struct gkyl_nmat* vecs_contribution, int* ndonors_cum, const struct gkyl_range* local_range_update, int unique_donor_mats, const struct gkyl_rect_grid *grid)
{
  int num_threads = ftar->ncomp*num_tar_locs;
  gkyl_bc_twistshift_inc_cu_kernel<<<(num_threads+255)/GKYL_DEFAULT_NUM_THREADS,GKYL_DEFAULT_NUM_THREADS>>>(ftar->on_dev, tar_locs, num_tar_locs, vecs_contribution->on_dev, ndonors_cum, *local_range_update, unique_donor_mats, *grid);
}


void
gkyl_bc_twistshift_clear_cu(struct gkyl_array* ftar, long* locs, int num_locs)
{
  int num_threads = num_locs*ftar->ncomp;
  gkyl_bc_twistshift_clear_cu_kernel<<<(num_threads+255)/GKYL_DEFAULT_NUM_THREADS,GKYL_DEFAULT_NUM_THREADS>>>(ftar->on_dev, locs, num_locs);
}


void
gkyl_bc_twistshift_set_vecsdo_cu(struct gkyl_array* fdo, long* locs, struct gkyl_nmat* vecsdo, int donor_factor, int unique_donor_mats){
  gkyl_bc_twistshift_set_vecsdo_cu_kernel<<<(vecsdo->num*donor_factor*fdo->ncomp+255)/GKYL_DEFAULT_NUM_THREADS, GKYL_DEFAULT_NUM_THREADS>>>(fdo->on_dev, locs, vecsdo->on_dev, donor_factor, unique_donor_mats);
}
