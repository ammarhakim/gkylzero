#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>
#include <gkyl_dg_bin_ops.h>

#include <assert.h>

struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_array_average_inp *inp)
{
  // works for p = 1 only
  assert(inp->basis.poly_order == 1); 

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
    up->issub_dim[d] = 1 - inp->avg_dim[d];

  int k = 0;
  for (unsigned d=0; d < up->ndim; ++d) 
    if(up->issub_dim[d]){
      up->sub_dir[k] = d;
      k++;
    }

  // Compute the cell sub-dimensional volume
  up->subvol = 1.0;
  for (unsigned d=0; d < up->ndim; ++d){
    if (up->issub_dim[d] == 0) up->subvol *= 0.5*inp->grid->dx[d];
  }

  up->integrant = gkyl_array_new(GKYL_DOUBLE, inp->basis.num_basis, up->local_ext.volume);

  // Handle a possible weighted average
  if (inp->weights) {
    up->isweighted = true;
    up->weights = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    gkyl_array_set(up->weights, 1.0, inp->weights);
    // Compute the subdim integral of the weights (for volume division after integration)
    up->integral_weights = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis_avg.num_basis, up->local_avg.volume)
                                  : gkyl_array_new(GKYL_DOUBLE, up->basis_avg.num_basis, up->local_avg.volume);
    // create new average routine to integrate the weights
    struct gkyl_array_average *int_w;
    // declare a weightless average 
    // (this will simply integrate, not recursive call because here we do not put weights)
    struct gkyl_array_average_inp inp_integral = {
    .grid = inp->grid,
    .basis = inp->basis,
    .basis_avg = inp->basis_avg,
    .local = inp->local,
    .local_ext = inp->local_ext,
    .local_avg = inp->local_avg,
    .weights = NULL, // No weights -> integral only
    .avg_dim = inp->avg_dim,
    .use_gpu = inp->use_gpu
    };
    int_w = gkyl_array_average_new(&inp_integral);
    // run the updater to integrate the weights
    gkyl_array_average_advance(int_w, up->weights, up->integral_weights);
    // release the updater
    gkyl_array_average_release(int_w);
    // Allocate memory to prepare the weak division at the end of the advance routine
    up->div_mem = gkyl_dg_bin_op_mem_new(up->local_avg.volume, up->basis_avg.num_basis);
  } 
  else{
    up->weights = NULL;
    up->integral_weights = NULL;
    up->div_mem = NULL;
    up->isweighted = false;
  }

  // Choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up);

  #ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_array_average_cu_dev_new(up);
  #endif

  return up;
}

void gkyl_array_average_advance(gkyl_array_average *up, 
  const struct gkyl_array * fin, struct gkyl_array *avgout)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, fin, avgout);
    return;
  }
#endif

  // Copy the input array to the integrant local array
  gkyl_array_copy_range(up->integrant, fin, &up->local_ext);

  // If we provided some weights, we integrate a weighted fin
  if(up->isweighted)
    gkyl_dg_mul_op_range(up->basis, 0, up->integrant, 0, up->integrant, 0, up->weights, &up->local);

  // clear the array that will contain the result
  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  struct gkyl_range_iter cmp_iter, sub_iter;
  struct gkyl_range cmp_rng; // this is the complementary range, sub + cmp = full
  
  // We now loop on the range of the averaged array
  gkyl_range_iter_init(&sub_iter, &up->local_avg);
  while (gkyl_range_iter_next(&sub_iter)) {
    long sub_lidx = gkyl_range_idx(&up->local_avg, sub_iter.idx);

    // We need to pass the moving index to the deflate operation as a sub dimensional iterator
    int parent_idx[GKYL_MAX_CDIM] = {0};
    int cnter = 0;
    for (int i = 0; i < up->basis.ndim; i++){
      if(up->issub_dim[i]){
        parent_idx[i] = sub_iter.idx[cnter];
        cnter ++;
      }
    }

    // set up the complementary range, which is the range through all averaged dimensions
    gkyl_range_deflate(&cmp_rng, &up->local, up->issub_dim, parent_idx);
    gkyl_range_iter_no_split_init(&cmp_iter, &cmp_rng);
    while (gkyl_range_iter_next(&cmp_iter)) {
      long cmp_lidx = gkyl_range_idx(&cmp_rng, cmp_iter.idx);

      const double *fin_i = gkyl_array_cfetch(up->integrant, cmp_lidx);
      double *avg_i = gkyl_array_fetch(avgout, sub_lidx);

      up->kernel(up->subvol, NULL, fin_i, avg_i);
    }
  }

  //If we provided some weights, we now divide by the integrated weight
  if(up->isweighted)
    gkyl_dg_div_op_range(up->div_mem, up->basis_avg,
     0, avgout, 0, avgout, 0, up->integral_weights, &up->local_avg);

}

void gkyl_array_average_release(gkyl_array_average *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  if(up->integrant) gkyl_array_release(up->integrant);
  if(up->weights) gkyl_array_release(up->weights);
  if(up->integral_weights) gkyl_array_release(up->integral_weights);
  if(up->div_mem) gkyl_dg_bin_op_mem_release(up->div_mem);
  gkyl_free(up);
}
