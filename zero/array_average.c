#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>

#include <assert.h>

struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_array_average_inp *inp)
{
  // works for p = 1 only
  assert(inp->basis.poly_order <= 2); 

  // Allocate space for new updater.
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = inp->use_gpu;
  up->ndim = inp->basis.ndim;
  up->basis = inp->basis;
  up->basis_avg = inp->basis_avg;

  // copy the total and sub ranges on the updater 
  // (will be used in advance and in the declaration of the integrated weighted array)
  up->local = *inp->local;
  up->local_ext = *inp->local_ext;
  up->local_avg = *inp->local_avg;

  // Set up the array of all dimensions that are conserved after the average (=0 for removed)
  // according to the operation input variable
  up->navg_dim = 0;
  for (unsigned d=0; d < up->ndim; ++d){
    up->avg_dim[d] = inp->avg_dim[d];
    up->navg_dim += inp->avg_dim[d]; 
  }
  assert(up->navg_dim <= up->ndim);

  for (unsigned d=0; d < up->ndim; ++d) 
    up->isdim_sub[d] = 1 - inp->avg_dim[d];

  int k = 0;
  for (unsigned d=0; d < up->ndim; ++d) 
    if(up->isdim_sub[d]){
      up->sub_dir[k] = d;
      k++;
    }

  // Compute the cell sub-dimensional volume
  up->subvol = 1.0;
  for (unsigned d=0; d < up->ndim; ++d){
    if (up->isdim_sub[d] == 0) up->subvol *= 0.5*inp->grid->dx[d];
  }

  // compute the inverse of the volume of the averaging space
  up->vol_avg_inv = 1.;
  for (int d = 0; d < up->ndim; d++)
    up->vol_avg_inv *= up->avg_dim[d]? inp->grid->upper[d] - inp->grid->lower[d] : 1.0;
  up->vol_avg_inv = 1./up->vol_avg_inv;

  // Handle a possible weighted average
  up->isweighted = false;
  up->weight = NULL;
  up->weight_avg = NULL;
  up->div_mem = NULL;
  if (inp->weight) {
    up->isweighted = true;
    // Compute the subdim integral of the weight (for volume division after integration)
    up->weight_avg = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume)
                                  : gkyl_array_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume);
    // create new average routine to integrate the weight
    struct gkyl_array_average_inp inp_integral = {
      .grid = inp->grid,
      .basis = inp->basis,
      .basis_avg = inp->basis_avg,
      .local = inp->local,
      .local_ext = inp->local_ext,
      .local_avg = inp->local_avg,
      .weight = NULL, // No weight -> average only
      .avg_dim = inp->avg_dim,
      .use_gpu = inp->use_gpu
    };
    struct gkyl_array_average *int_w = gkyl_array_average_new(&inp_integral);
    // run the updater to integrate the weight
    gkyl_array_average_advance(int_w, inp->weight, up->weight_avg);
    gkyl_array_average_release(int_w);
    // multiply by the volume to get the integral and not the average
    gkyl_array_scale(up->weight_avg,1./up->vol_avg_inv);
    // Allocate memory to prepare the weak division at the end of the advance routine
    up->div_mem = gkyl_dg_bin_op_mem_new(up->local_avg.volume, up->basis_avg.num_basis);
    // Assign the weight pointer to the input weight array
    up->weight = gkyl_array_acquire(inp->weight);
  } 
  else {
    // Assign the weight pointer to the identity weight
    up->weight = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, 1);
    double *w_0 = gkyl_array_fetch(up->weight, 0);
    w_0[0] = pow(sqrt(2.),up->ndim);
  }

  // Choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_array_average_cu_dev_new(up);
#endif

  return up;
}

void gkyl_array_average_advance(const struct gkyl_array_average *up, 
  const struct gkyl_array * fin, struct gkyl_array *avgout)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, fin, avgout);
    return;
  }
#endif

  // clear the array that will contain the result
  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  struct gkyl_range_iter iter_cmp, iter_avg;
  struct gkyl_range rng_cmp; // this is the complementary range, sub + cmp = full
  
  // We now loop on the range of the averaged array
  gkyl_range_iter_init(&iter_avg, &up->local_avg);
  while (gkyl_range_iter_next(&iter_avg)) {
    long lidx_avg = gkyl_range_idx(&up->local_avg, iter_avg.idx);

    // We need to pass the moving index to the deflate operation as a sub dimensional iterator
    int parent_idx[GKYL_MAX_CDIM] = {0};
    int cnter = 0;
    for (int i = 0; i < up->basis.ndim; i++){
      if (up->isdim_sub[i]) {
        parent_idx[i] = iter_avg.idx[cnter];
        cnter ++;
      }
    }

    // set up the complementary range, which is the range through all averaged dimensions
    gkyl_range_deflate(&rng_cmp, &up->local, up->isdim_sub, parent_idx);
    gkyl_range_iter_no_split_init(&iter_cmp, &rng_cmp);
    while (gkyl_range_iter_next(&iter_cmp)) {
      long lidx_cmp = gkyl_range_idx(&rng_cmp, iter_cmp.idx);

      const double *fin_i = gkyl_array_cfetch(fin, lidx_cmp);
      const double *win_i = up->isweighted? gkyl_array_cfetch(up->weight, lidx_cmp) : 
        gkyl_array_cfetch(up->weight, 0);

      double *avg_i = gkyl_array_fetch(avgout, lidx_avg);

      up->kernel(up->subvol, win_i, fin_i, avg_i);
    }
  }

  // If we provided some weight, we now divide by the integrated weight
  if (up->isweighted){
    gkyl_dg_div_op_range(up->div_mem, up->basis_avg,
      0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg);
  } else{
    // divide by the volume of the averaging domain
    gkyl_array_scale(avgout,up->vol_avg_inv);
  }

}

void gkyl_array_average_release(struct gkyl_array_average *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  // if (up->use_gpu)
  //   gkyl_cu_free(up->on_dev);
#endif

  gkyl_array_release(up->weight);
  if (up->isweighted) {
    gkyl_array_release(up->weight_avg);
    gkyl_dg_bin_op_mem_release(up->div_mem);
  }

  gkyl_free(up);
}
