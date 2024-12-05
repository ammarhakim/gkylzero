#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>
#include <gkyl_dg_bin_ops.h>

#include <assert.h>

struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis tot_basis, const struct gkyl_basis sub_basis, 
  const struct gkyl_range tot_rng, const struct gkyl_range sub_rng, 
  struct gkyl_array *weights, enum gkyl_array_average_op op, 
  bool use_gpu)
{
  // works for p = 1 only
  assert(tot_basis.poly_order == 1); 

  // Allocate space for new updater.
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = use_gpu;
  up->ndim = tot_basis.ndim;
  up->tot_basis = tot_basis;
  up->sub_basis = sub_basis;
  // copy the total and sub ranges on the updater 
  // (will be used in advance and in the declaration of the integrated weighted array)
  gkyl_sub_range_init(&up->tot_rng, &tot_rng, tot_rng.lower, tot_rng.upper);
  gkyl_sub_range_init(&up->sub_rng, &sub_rng, sub_rng.lower, sub_rng.upper);
  // here there is two "sub":
  // - one from the gkyl_sub_range_init (but it is actually a full copy)
  // - the other from the sub dimensional output array from the averaging operation

  // Set up the array of all dimensions that are conserved after the average (=0 for removed)
  // according to the operation input variable
  for (unsigned d=0; d < up->ndim; ++d) up->isavg_dim[d] = 0;
  switch (op) {
    case GKYL_ARRAY_AVERAGE_OP: // Full integration
      assert(tot_basis.ndim >= 1); // Ensure at least 1 dimension exists
      for (unsigned d=0; d < up->ndim; ++d){ 
        up->sub_dir[d] = 0;
      }
      break;
    case GKYL_ARRAY_AVERAGE_OP_X:
      assert(tot_basis.ndim >= 1); // Ensure at least 1 dimension exists
      up->isavg_dim[0] = 1;
      up->sub_dir[0] = 0;
      break;
    case GKYL_ARRAY_AVERAGE_OP_Y:
      assert(tot_basis.ndim >= 2); // Ensure at least 2 dimensions for Y
      up->isavg_dim[1] = 1;
      up->sub_dir[0] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_Z:
      assert(tot_basis.ndim >= 3); // Ensure at least 3 dimensions for Z
      up->isavg_dim[2] = 1;
      up->sub_dir[0] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XY:
      assert(tot_basis.ndim >= 3); // Ensure at least 3 dimensions for XY reduction
      up->isavg_dim[0] = 1;
      up->isavg_dim[1] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XZ:
      assert(tot_basis.ndim >= 3); // Ensure at least 3 dimensions for XZ
      up->isavg_dim[0] = 1;
      up->isavg_dim[2] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_YZ:
      assert(tot_basis.ndim >= 3); // Ensure at least 3 dimensions for YZ
      up->isavg_dim[1] = 1;
      up->isavg_dim[2] = 1;
      up->sub_dir[0] = 1;
      up->sub_dir[1] = 2;
      break;
    default:
      assert(false && "Invalid operation in switch(op)");
      break;
  }

  
  // Compute the cell sub-dimensional volume
  up->subvol = 1.0;
  for (unsigned d=0; d < up->ndim; ++d){
    if (up->isavg_dim[d] == 0) up->subvol *= 0.5*grid->dx[d];
  }
  // printf("subvolume = %g\n",up->subvol);

  // Handle a possible weighted average
  up->isweighted = false;
  if (weights) {
    up->isweighted = true;
    up->weights = gkyl_array_new(GKYL_DOUBLE, weights->ncomp, weights->size);
    gkyl_array_copy_range(up->weights, weights, &up->tot_rng);
    // I use gkyl_array_copy before and it was not copying on the total range but one less, why?

    // CHECK THE WEIGHT
    // printf("Check the weights in array avg updater\n");
    // struct gkyl_range_iter iter;
    // gkyl_range_iter_init(&iter, &tot_rng);
    // while (gkyl_range_iter_next(&iter)) {
    //   long lidx = gkyl_range_idx(&tot_rng, iter.idx);
    //   const double *w_i = gkyl_array_cfetch(up->weights, lidx);
    //   printf("w[%ld][0]=%g, w[%ld][1]=%g\n",lidx,w_i[0],lidx,w_i[1]);
    // }
    // printf("Check done\n");

    // Compute the subdim integral of the weights (for volume division after integration)
    up->integral_weights = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->sub_basis.num_basis, up->sub_rng.volume)
                                  : gkyl_array_new(GKYL_DOUBLE, up->sub_basis.num_basis, up->sub_rng.volume);
    // create new average routine to integrate the weights
    struct gkyl_array_average *int_w;
    // declare a weightless average 
    // (this will simply integrate, not recursive call because here we do not put weights)
    int_w = gkyl_array_average_new(grid, up->tot_basis, up->sub_basis, up->tot_rng, up->sub_rng, NULL, op, up->use_gpu);
    // run the updater to integrate the weights
    // printf("Integrate weights\n");
    gkyl_array_average_advance(int_w, up->weights, up->integral_weights);
    // printf("Integrate weights done\n");
    // release the updater
    gkyl_array_average_release(int_w);
    // const double *wint = gkyl_array_cfetch(up->integral_weights, 0);
    // printf("integral of the weights = %g \n",wint[0]);

    // Allocate memory to prepare the weak division at the end of the advance routine
    up->div_mem = gkyl_dg_bin_op_mem_new(up->sub_rng.volume, up->sub_basis.num_basis);
  } 

  // Choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up, &up->tot_basis, op);

  #ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_array_average_cu_dev_new(up, grid, up->tot_basis, op);
  #endif

  return up;
}

void gkyl_array_average_advance(gkyl_array_average *up, 
  struct gkyl_array * fin, struct gkyl_array *avgout)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, tot_rng, sub_rng, fin, out);
    return;
  }
#endif

  // If we provided some weights, we integrate a weighted fin
  if(up->isweighted)
    gkyl_dg_mul_op_range(up->tot_basis, 0, fin, 0, fin, 0, up->weights, &up->tot_rng);

  gkyl_array_clear_range(avgout, 0.0, &up->sub_rng);

  if (up->sub_rng.volume > 1){
    struct gkyl_range_iter cmp_iter, sub_iter;
    struct gkyl_range cmp_rng; // this is the complementary range, sub + cmp = full
    // We now loop on the range of the averaged array
    gkyl_range_iter_init(&sub_iter, &up->sub_rng);
    while (gkyl_range_iter_next(&sub_iter)) {
      long sub_lidx = gkyl_range_idx(&up->sub_rng, sub_iter.idx);

      gkyl_range_deflate(&cmp_rng, &up->tot_rng, up->isavg_dim, sub_iter.idx);
      // printf("sub iter loop\n");
      gkyl_range_iter_no_split_init(&cmp_iter, &cmp_rng);

      while (gkyl_range_iter_next(&cmp_iter)) {
        // printf("\tcmp iter loop\n");
        long cmp_lidx = gkyl_range_idx(&cmp_rng, cmp_iter.idx);
        const double *fin_i = gkyl_array_cfetch(fin, cmp_lidx);
        double *avg_i = gkyl_array_fetch(avgout, sub_lidx);

        up->kernel(up->subvol, NULL, fin_i, avg_i);
        // printf("(subdim loop) fin_i[%2.0ld][0] = %6.4g, fin_i[%2.0ld][1] = %4.2g, avg_i[0] = %6.4g\n",cmp_lidx,fin_i[0],cmp_lidx,fin_i[1],avg_i[0]);
      }
    }
  } 
  else // This is the case if we are asking for a full integration
  {
    struct gkyl_range_iter tot_iter;
    // this is the complementary range, sub + cmp = full
    // We now loop on the range of the entire array
    gkyl_range_iter_init(&tot_iter, &up->tot_rng);
    while (gkyl_range_iter_next(&tot_iter)) {
        long tot_lidx = gkyl_range_idx(&up->tot_rng, tot_iter.idx);
        const double *fin_i = gkyl_array_cfetch(fin, tot_lidx);
        const double *win_i = gkyl_array_cfetch(up->weights, tot_lidx);
        double *avg_i = gkyl_array_fetch(avgout, 0);
        up->kernel(up->subvol, NULL, fin_i, avg_i);

        // To check the integration:
        // printf("(full loop) fin_i[%2.0ld][0] = %6.4g, fin_i[%2.0ld][1] = %6.4g, avg_i[0] = %6.4g\n",tot_lidx,fin_i[0],tot_lidx,fin_i[1],avg_i[0]);
    }
  }

  // // If we provided some weights, we now divide by the integrated weight
  if(up->isweighted)
    gkyl_dg_div_op_range(up->div_mem, up->sub_basis,
     0, avgout, 0, avgout, 0, up->integral_weights, &up->sub_rng);

}

void gkyl_array_average_release(gkyl_array_average *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  if(up->weights) gkyl_array_release(up->weights);
  if(up->integral_weights) gkyl_array_release(up->integral_weights);
  if(up->div_mem) gkyl_dg_bin_op_mem_release(up->div_mem);

  gkyl_free(up);
}
