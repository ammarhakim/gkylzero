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
  assert(inp->tot_basis.poly_order == 1); 

  // Allocate space for new updater.
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = inp->use_gpu;
  up->ndim = inp->tot_basis.ndim;
  up->tot_basis = inp->tot_basis;
  up->sub_basis = inp->sub_basis;

  // copy the total and sub ranges on the updater 
  // (will be used in advance and in the declaration of the integrated weighted array)
  up->tot_rng = *inp->tot_rng;
  up->tot_rng_ext = *inp->tot_rng_ext;
  up->sub_rng = *inp->sub_rng;

  // Set up the array of all dimensions that are conserved after the average (=0 for removed)
  // according to the operation input variable
  for (unsigned d=0; d < up->ndim; ++d) up->issub_dim[d] = 0;
  switch (inp->op) {
    case GKYL_ARRAY_AVERAGE_OP: // Full integration
      assert(inp->tot_basis.ndim >= 1); // Ensure at least 1 dimension exists
      for (unsigned d=0; d < up->ndim; ++d){ 
        up->sub_dir[d] = 0;
      }
      break;
    case GKYL_ARRAY_AVERAGE_OP_X: // integration all except x
      assert(inp->tot_basis.ndim >= 1); // Ensure at least 1 dimension exists
      up->issub_dim[0] = 1; // here the first dimension is conserved
      up->sub_dir[0] = 0; // the first dimension in the reduced array is the first dim in the total array
      break;
    case GKYL_ARRAY_AVERAGE_OP_Y: // integration all except y
      assert(inp->tot_basis.ndim >= 2); // Ensure at least 2 dimensions for Y
      up->issub_dim[1] = 1;
      up->sub_dir[0] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_Z: // integration all except z
      assert(inp->tot_basis.ndim >= 3); // Ensure at least 3 dimensions for Z
      up->issub_dim[2] = 1;
      up->sub_dir[0] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XY: // integration all except xy
      assert(inp->tot_basis.ndim >= 3); // Ensure at least 3 dimensions for XY reduction
      up->issub_dim[0] = 1;
      up->issub_dim[1] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XZ: // integration all except xz
      assert(inp->tot_basis.ndim >= 3); // Ensure at least 3 dimensions for XZ
      up->issub_dim[0] = 1;
      up->issub_dim[2] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_YZ: // integration all except yz
      assert(inp->tot_basis.ndim >= 3); // Ensure at least 3 dimensions for YZ
      up->issub_dim[1] = 1;
      up->issub_dim[2] = 1;
      up->sub_dir[0] = 1;
      up->sub_dir[1] = 2;
      break;
    default:
      assert(false && "-array_average: Invalid operation in switch(op)\n");
      break;
  }

  // Compute the cell sub-dimensional volume
  up->subvol = 1.0;
  for (unsigned d=0; d < up->ndim; ++d){
    if (up->issub_dim[d] == 0) up->subvol *= 0.5*inp->grid->dx[d];
  }

  up->integrant = gkyl_array_new(GKYL_DOUBLE, inp->tot_basis.num_basis, up->tot_rng_ext.volume);

  // Handle a possible weighted average
  if (inp->weights) {
    up->isweighted = true;
    up->weights = gkyl_array_new(GKYL_DOUBLE, up->tot_basis.num_basis, up->tot_rng_ext.volume);
    gkyl_array_set(up->weights, 1.0, inp->weights);
    // Compute the subdim integral of the weights (for volume division after integration)
    up->integral_weights = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->sub_basis.num_basis, up->sub_rng.volume)
                                  : gkyl_array_new(GKYL_DOUBLE, up->sub_basis.num_basis, up->sub_rng.volume);
    // create new average routine to integrate the weights
    struct gkyl_array_average *int_w;
    // declare a weightless average 
    // (this will simply integrate, not recursive call because here we do not put weights)
    struct gkyl_array_average_inp inp_integral = {
    .grid = inp->grid,
    .tot_basis = inp->tot_basis,
    .sub_basis = inp->sub_basis,
    .tot_rng = inp->tot_rng,
    .tot_rng_ext = inp->tot_rng_ext,
    .sub_rng = inp->sub_rng,
    .weights = NULL, // No weights -> integral only
    .op = inp->op,
    .use_gpu = inp->use_gpu
    };
    int_w = gkyl_array_average_new(&inp_integral);
    // run the updater to integrate the weights
    gkyl_array_average_advance(int_w, up->weights, up->integral_weights);
    // release the updater
    gkyl_array_average_release(int_w);
    // Allocate memory to prepare the weak division at the end of the advance routine
    up->div_mem = gkyl_dg_bin_op_mem_new(up->sub_rng.volume, up->sub_basis.num_basis);
  } 
  else{
    up->weights = NULL;
    up->integral_weights = NULL;
    up->div_mem = NULL;
    up->isweighted = false;
  }

  // Choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up, &up->tot_basis, inp->op);

  #ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_array_average_cu_dev_new(up, grid, up->tot_basis, inp->op);
  #endif

  return up;
}

void gkyl_array_average_advance(gkyl_array_average *up, 
  const struct gkyl_array * fin, struct gkyl_array *avgout)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, tot_rng, sub_rng, fin, out);
    return;
  }
#endif

  // Copy the input array to the integrant local array
  gkyl_array_copy_range(up->integrant, fin, &up->tot_rng_ext);

  // If we provided some weights, we integrate a weighted fin
  if(up->isweighted)
    gkyl_dg_mul_op_range(up->tot_basis, 0, up->integrant, 0, up->integrant, 0, up->weights, &up->tot_rng);

  // clear the array that will contain the result
  gkyl_array_clear_range(avgout, 0.0, &up->sub_rng);

  struct gkyl_range_iter cmp_iter, sub_iter;
  struct gkyl_range cmp_rng; // this is the complementary range, sub + cmp = full
  
  // We now loop on the range of the averaged array
  // printf("// We now loop on the range of the averaged array\n");
  gkyl_range_iter_init(&sub_iter, &up->sub_rng);
  while (gkyl_range_iter_next(&sub_iter)) {
    long sub_lidx = gkyl_range_idx(&up->sub_rng, sub_iter.idx);

    // We need to pass the moving index to the deflate operation as a sub dimensional iterator
    int parent_idx[GKYL_MAX_CDIM] = {0};
    int cnter = 0;
    for (int i = 0; i < up->tot_basis.ndim; i++){
      if(up->issub_dim[i]){
        parent_idx[i] = sub_iter.idx[cnter];
        cnter ++;
      }
    }

    gkyl_range_deflate(&cmp_rng, &up->tot_rng, up->issub_dim, parent_idx);
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
  if(up->integrant) gkyl_array_release(up->integrant);
  if(up->weights) gkyl_array_release(up->weights);
  if(up->integral_weights) gkyl_array_release(up->integral_weights);
  if(up->div_mem) gkyl_dg_bin_op_mem_release(up->div_mem);
  gkyl_free(up);
}
