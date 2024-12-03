#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>
#include <gkyl_dg_bin_ops.h>

#include <assert.h>

struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis basis, 
  const struct gkyl_array *weights, enum gkyl_array_average_op op, bool use_gpu)
{
  // works for p = 1 only
  assert(basis.poly_order == 1); 

  // Allocate space for new updater.
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = use_gpu;
  up->ndim = basis.ndim;
  up->basis = basis;
  up->isweighted = false;

  if (weights) {
    printf("weighted\n");
    up->isweighted = true;
    up->weights = gkyl_array_new(GKYL_DOUBLE, weights->ncomp, weights->size);
    gkyl_array_copy(up->weights, weights);
  } 

  // Set up the array of all dimensions that are conserved after the average (=0 for removed)
  // according to the operation input variable
  for (unsigned d=0; d < up->ndim; ++d) up->isavg_dim[d] = 0;
  switch (op) {
    case GKYL_ARRAY_AVERAGE_OP: // Full integration
      assert(basis.ndim >= 1); // Ensure at least 1 dimension exists
      for (unsigned d=0; d < up->ndim; ++d){ 
        up->sub_dir[d] = 0;
      }
      break;
    case GKYL_ARRAY_AVERAGE_OP_X:
      assert(basis.ndim >= 1); // Ensure at least 1 dimension exists
      up->isavg_dim[0] = 1;
      up->sub_dir[0] = 0;
      break;
    case GKYL_ARRAY_AVERAGE_OP_Y:
      assert(basis.ndim >= 2); // Ensure at least 2 dimensions for Y
      up->isavg_dim[1] = 1;
      up->sub_dir[0] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_Z:
      assert(basis.ndim >= 3); // Ensure at least 3 dimensions for Z
      up->isavg_dim[2] = 1;
      up->sub_dir[0] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XY:
      assert(basis.ndim >= 3); // Ensure at least 3 dimensions for XY reduction
      up->isavg_dim[0] = 1;
      up->isavg_dim[1] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 1;
      break;
    case GKYL_ARRAY_AVERAGE_OP_XZ:
      assert(basis.ndim >= 3); // Ensure at least 3 dimensions for XZ
      up->isavg_dim[0] = 1;
      up->isavg_dim[2] = 1;
      up->sub_dir[0] = 0;
      up->sub_dir[1] = 2;
      break;
    case GKYL_ARRAY_AVERAGE_OP_YZ:
      assert(basis.ndim >= 3); // Ensure at least 3 dimensions for YZ
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
  printf("subvolume = %g\n",up->subvol);

  #ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_array_average_cu_dev_new(up, grid, basis, op);
  #endif

  // Choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up, &basis, op);

  return up;
}

void gkyl_array_average_advance(gkyl_array_average *up, 
  const struct gkyl_range *full_rng, const struct gkyl_range *sub_rng,
  struct gkyl_array * fin, struct gkyl_array *avgout)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, full_rng, sub_rng, fin, out);
    return;
  }
#endif

  // If we provided some weights, we integrate a weighted fin
  if(up->isweighted)
    gkyl_dg_mul_op(up->basis,1,fin,1,up->weights,1,fin);

  gkyl_array_clear_range(avgout, 0.0, sub_rng);

  if (sub_rng->volume > 1){
    struct gkyl_range_iter cmp_iter, sub_iter;
    struct gkyl_range cmp_rng; // this is the complementary range, sub + cmp = full
    // We now loop on the range of the averaged array
    gkyl_range_iter_init(&sub_iter, sub_rng);
    while (gkyl_range_iter_next(&sub_iter)) {
      long sub_lidx = gkyl_range_idx(sub_rng, sub_iter.idx);

      gkyl_range_deflate(&cmp_rng, full_rng, up->isavg_dim, sub_iter.idx);
      printf("sub iter loop\n");
      gkyl_range_iter_no_split_init(&cmp_iter, &cmp_rng);

      while (gkyl_range_iter_next(&cmp_iter)) {
        printf("\tcmp iter loop\n");
        long cmp_lidx = gkyl_range_idx(&cmp_rng, cmp_iter.idx);
        const double *fin_i = gkyl_array_cfetch(fin, cmp_lidx);
        double *avg_i = gkyl_array_fetch(avgout, sub_lidx);
        printf("\tfin_i[%ld] = %g\n",cmp_lidx,fin_i[0]);

        up->kernel(fin_i, avg_i, up->subvol);
      }
    }
  } 
  else // This is the case if we are asking for a full integration
  {
    struct gkyl_range_iter full_iter;
    // this is the complementary range, sub + cmp = full
    // We now loop on the range of the entire array
    gkyl_range_iter_init(&full_iter, full_rng);
    while (gkyl_range_iter_next(&full_iter)) {
        long full_lidx = gkyl_range_idx(full_rng, full_iter.idx);
        const double *fin_i = gkyl_array_cfetch(fin, full_lidx);
        double *avg_i = gkyl_array_fetch(avgout, 0);
        up->kernel(fin_i, avg_i, up->subvol);
    }
  }

}

void gkyl_array_average_release(gkyl_array_average *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  // if(up->weights)
  //   gkyl_array_release(up->weights);
  // gkyl_free(up->kernel);
  gkyl_free(up);
}
