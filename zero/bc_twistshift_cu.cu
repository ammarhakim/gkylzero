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
gkyl_bc_twistshift_inc_cu_kernel(struct gkyl_array* ftar, long* tar_locs, int num_tar_locs, struct gkyl_nmat* vecstar, int* ndonors_cum)
{
  int x_idx = START_ID/ftar->ncomp;
  int set_idx = START_ID%ftar->ncomp;
  int do_start = ndonors_cum[x_idx];
  int do_end = ndonors_cum[x_idx+1];
  if(x_idx<num_tar_locs && set_idx<ftar->ncomp){
    double *ftar_itr = (double*) gkyl_array_fetch(ftar, tar_locs[x_idx]);
    for(int i = do_start; i < do_end; i++){
      struct gkyl_mat temp = gkyl_nmat_get(vecstar, i);
      ftar_itr[set_idx] += gkyl_mat_get(&temp, set_idx, 0);
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
gkyl_bc_twistshift_set_vecsdo_cu_kernel(const struct gkyl_array* fdo, long* locs, struct gkyl_nmat* vecsdo)
{
  int vecdo_idx = START_ID/fdo->ncomp;
  int set_idx = START_ID%fdo->ncomp;
  if(vecdo_idx < vecsdo->num && set_idx<fdo->ncomp){
    const double *fdo_itr = (const double*) gkyl_array_cfetch(fdo, locs[vecdo_idx]);
    struct gkyl_mat temp = gkyl_nmat_get(vecsdo, vecdo_idx);
    gkyl_mat_set(&temp, set_idx, 0, fdo_itr[set_idx]);
  }
}


void
gkyl_bc_twistshift_inc_cu(const struct gkyl_array* ftar, long* tar_locs, int num_tar_locs, struct gkyl_nmat* vecstar, int* ndonors_cum)
{
  int num_threads = ftar->ncomp*num_tar_locs;
  gkyl_bc_twistshift_inc_cu_kernel<<<(num_threads+255)/GKYL_DEFAULT_NUM_THREADS,GKYL_DEFAULT_NUM_THREADS>>>(ftar->on_dev, tar_locs, num_tar_locs, vecstar->on_dev, ndonors_cum);
}


void
gkyl_bc_twistshift_clear_cu(struct gkyl_array* ftar, long* locs, int num_locs)
{
  int num_threads = num_locs*ftar->ncomp;
  gkyl_bc_twistshift_clear_cu_kernel<<<(num_threads+255)/GKYL_DEFAULT_NUM_THREADS,GKYL_DEFAULT_NUM_THREADS>>>(ftar->on_dev, locs, num_locs);
}


void
gkyl_bc_twistshift_set_vecsdo_cu(struct gkyl_array* fdo, long* locs, struct gkyl_nmat* vecsdo){
  gkyl_bc_twistshift_set_vecsdo_cu_kernel<<<(vecsdo->num*fdo->ncomp+255)/GKYL_DEFAULT_NUM_THREADS, GKYL_DEFAULT_NUM_THREADS>>>(fdo->on_dev, locs, vecsdo->on_dev);
}
