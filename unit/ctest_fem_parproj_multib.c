// Test the projection onto an FEM basis that is continuous in the parallel
// direction on multiple grid blocks.
//
#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_parproj_multib.h>

// Max number of blocks in this test.
#define CTEST_NUM_BLOCKS_MAX 20

struct test_ctx {
  int num_blocks; // Number of grid blocks.
  int ndim; // Number of dimensions.
  double lower[CTEST_NUM_BLOCKS_MAX*GKYL_MAX_DIM]; // Grid lower limit in each direction for each block.
  double upper[CTEST_NUM_BLOCKS_MAX*GKYL_MAX_DIM]; // Grid upper limit in each direction for each block.
};

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array**
mkarr_mb(int num_blocks, bool use_gpu, long nc, struct gkyl_range *ranges)
{
  // Allocate an array for each block.
  struct gkyl_array** a = gkyl_malloc(num_blocks*sizeof(struct gkyl_array *));
  for (int i=0; i<num_blocks; i++)
    a[i] = mkarr(use_gpu, nc, ranges[i].volume);

  return a;
}

void evalFunc1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  int num_blocks = tctx->num_blocks;
  int ndim = tctx->ndim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++) Lx[d] = upper[ndim*(num_blocks-1)+d] - lower[ndim*0+d];

  fout[0] = sin(2.*M_PI*x);
}

void eval_weight_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  int num_blocks = tctx->num_blocks;
  int ndim = tctx->ndim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++) Lx[d] = upper[ndim*(num_blocks-1)+d] - lower[ndim*0+d];

  fout[0] = 1.0;
}

void
test_1x(int poly_order, const bool isperiodic, bool use_gpu)
{
  int num_blocks = 2;
  // Grid extents and number of cells for each block, listed in flat arrays.
  double lower[] = {-0.5, 0.0}, upper[] = {0.0, 0.5};
  int cells[] = {2, 2};

  int ndim = sizeof(lower)/(num_blocks*sizeof(lower[0]));

  struct test_ctx proj_ctx = {
    .num_blocks = num_blocks, // Number of grid blocks. 
    .ndim = ndim, // Number of position space dimensions.
  };
  for (int bI=0; bI<num_blocks; bI++) {
    for (int d=0; d<ndim; d++) {
      proj_ctx.lower[ndim*bI+d] = lower[ndim*bI+d];
      proj_ctx.upper[ndim*bI+d] = upper[ndim*bI+d];
    }
  }

  // Grids.
  struct gkyl_rect_grid *grid = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_grid));
  for (int bI=0; bI<num_blocks; bI++) {
    gkyl_rect_grid_init(&grid[bI], ndim, &lower[bI*ndim], &upper[bI*ndim], &cells[bI*ndim]);
  }

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Ranges.
  int ghost[] = { 1, 1, 1 };
  struct gkyl_range *localRange = gkyl_malloc(num_blocks*sizeof(struct gkyl_range));
  struct gkyl_range *localRange_ext = gkyl_malloc(num_blocks*sizeof(struct gkyl_range)); 
  for (int i=0; i<num_blocks; i++) {
    gkyl_create_grid_ranges(&grid[i], ghost, &localRange_ext[i], &localRange[i]);
  }

  // create DG field we wish to make continuous.
  struct gkyl_array **rho = mkarr_mb(num_blocks, use_gpu, basis.num_basis, localRange_ext);
  // create array holding continuous field we'll compute.
  struct gkyl_array **phi = mkarr_mb(num_blocks, use_gpu, basis.num_basis, localRange_ext);

  struct gkyl_array **rho_ho, **phi_ho;
  if (use_gpu) {
    rho_ho = mkarr_mb(num_blocks, false, basis.num_basis, localRange_ext);
    phi_ho = mkarr_mb(num_blocks, false, basis.num_basis, localRange_ext);
  }
  else {
    rho_ho = gkyl_malloc(num_blocks*sizeof(struct gkyl_array *));
    phi_ho = gkyl_malloc(num_blocks*sizeof(struct gkyl_array *));
    for (int i=0; i<num_blocks; i++) {
      rho_ho[i] = gkyl_array_acquire(rho[i]);
      phi_ho[i] = gkyl_array_acquire(phi[i]);
    }
  }

  for (int bI=0; bI<num_blocks; bI++) {
    // Project analytic field onto the DG basis.
    gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid[bI], &basis,
      poly_order+1, 1, evalFunc1x, NULL);
    gkyl_proj_on_basis_advance(projob, 0.0, &localRange[bI], rho_ho[bI]);
    gkyl_proj_on_basis_release(projob);
    gkyl_array_copy(rho[bI], rho_ho[bI]);

    printf("b:%d | grid=%g:%g, localRange=%d:%d\n",bI,grid[bI].lower[0],grid[bI].upper[0],localRange[bI].lower[0],localRange[bI].upper[0]);
    char fname0[1024];
    sprintf(fname0, "ctest_fem_parproj_multib_%dx_p%d_rho_b%d.gkyl", ndim, poly_order, bI);
    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, rho_ho[bI], fname0);
  }

  // Create the multiblock grid and range.
  struct gkyl_rect_grid grid_mb;
  int cells_mb[ndim];
  for (int d=0; d<ndim-1; d++) {
    cells_mb[d] = cells[d];
  }
  int tot_cells_par = 0;
  for (int bI=0; bI<num_blocks; bI++)
    tot_cells_par += cells[bI*ndim+ndim-1];
  cells_mb[ndim-1] = tot_cells_par;
  gkyl_rect_grid_init(&grid_mb, ndim, &lower[0*ndim], &upper[(num_blocks-1)*ndim], cells_mb);

  struct gkyl_range localRange_mb, localRange_ext_mb;
  gkyl_create_grid_ranges(&grid_mb, ghost, &localRange_ext_mb, &localRange_mb);

  // Local range of each block in the multiblock range.
  struct gkyl_range *local_in_mb_range = gkyl_malloc(num_blocks*sizeof(struct gkyl_range));
  int rlower[ndim], rupper[ndim];
  for (int d=0; d<ndim-1; d++) {
    rlower[d] = localRange[0].lower[d];
    rupper[d] = localRange[0].upper[d];
  }
  int par_offset = 0;
  for (int bI=0; bI<num_blocks; bI++) {
    rlower[ndim-1] = par_offset+localRange[bI].lower[ndim-1];
    rupper[ndim-1] = par_offset+localRange[bI].upper[ndim-1];
    gkyl_sub_range_init(&local_in_mb_range[bI], &localRange_ext_mb, rlower, rupper);
    par_offset += localRange[bI].upper[ndim-1]-localRange[bI].lower[ndim-1]+1;
  }

  // Allocate the multiblock weight.
  struct gkyl_array *jac_mb = mkarr(use_gpu, basis.num_basis, localRange_ext_mb.volume);
  struct gkyl_array *jac_mb_ho = use_gpu? mkarr(false, basis.num_basis, localRange_ext_mb.volume) : gkyl_array_acquire(jac_mb);

  // Project the weight on the multiblock range.
  gkyl_proj_on_basis *projob_mb = gkyl_proj_on_basis_new(&grid_mb, &basis,
    poly_order+1, 1, eval_weight_1x, NULL);
  gkyl_proj_on_basis_advance(projob_mb, 0.0, &localRange_mb, jac_mb_ho);
  gkyl_proj_on_basis_release(projob_mb);
  gkyl_array_copy(jac_mb, jac_mb_ho);

  // Allocate the multiblock source and output.
  struct gkyl_array *rho_mb = mkarr(use_gpu, basis.num_basis, localRange_ext_mb.volume);
  struct gkyl_array *phi_mb = mkarr(use_gpu, basis.num_basis, localRange_ext_mb.volume);
  struct gkyl_array *rho_mb_ho = use_gpu? mkarr(false, basis.num_basis, localRange_ext_mb.volume) : gkyl_array_acquire(rho_mb);
  struct gkyl_array *phi_mb_ho = use_gpu? mkarr(false, basis.num_basis, localRange_ext_mb.volume) : gkyl_array_acquire(phi_mb);
  
  // Copy the source into the multiblock array.
  for (int bI=0; bI<num_blocks; bI++) {
    gkyl_array_copy_range_to_range(rho_mb, rho[bI], &local_in_mb_range[bI], &localRange[bI]);
  }
  // Write out the multiblock source.
  gkyl_array_copy(rho_mb_ho, rho_mb);
  char fname1[1024];
  sprintf(fname1, "ctest_fem_parproj_multib_%dx_p%d_rho_mb.gkyl", ndim, poly_order);
  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, rho_mb_ho, fname1);

  // parallel FEM projection method.
  struct gkyl_fem_parproj_multib *parproj = gkyl_fem_parproj_multib_new(num_blocks, &localRange_mb,
    localRange, &basis, GKYL_FEM_PARPROJ_NONE, jac_mb, jac_mb, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_multib_set_rhs(parproj, rho_mb, NULL);

  // Solve the problem.
  gkyl_fem_parproj_multib_solve(parproj, phi_mb);
  gkyl_fem_parproj_multib_release(parproj);
  gkyl_array_copy(phi_mb_ho, phi_mb);

  for (int bI=0; bI<num_blocks; bI++) {
    // Copy solution back to local range of each block and write it out.
    gkyl_array_copy_range_to_range(phi[bI], phi_mb, &localRange[bI], &local_in_mb_range[bI]);
    gkyl_array_copy(phi_ho[bI], phi[bI]);

    char fname2[1024];
    sprintf(fname2, "ctest_fem_parproj_multib_%dx_p%d_phi_b%d.gkyl", ndim, poly_order, bI);
    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, phi_ho[bI], fname2);
  }
  // Write out the multiblock source.
  char fname3[1024];
  sprintf(fname3, "ctest_fem_parproj_multib_%dx_p%d_phi_mb.gkyl", ndim, poly_order);
  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, phi_mb_ho, fname3);

//  if (poly_order == 1) {
//    // Solution (checked visually, also checked that phi is actually continuous,
//    // and checked that visually looks like results in g2):
//    const double sol[8] = {-0.9089542445638024, -0.4554124667453318,
//                           -0.8488758876834943,  0.4900987222626481,
//                            0.8488758876834943,  0.490098722262648 ,
//                            0.9089542445638024, -0.4554124667453318};
//    const double *phi_p;
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 1);
//    TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
//    TEST_MSG("Expected: %.13e in cell (%d)", sol[0], 1);
//    TEST_MSG("Produced: %.13e", phi_p[0]);
//    TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 2);
//    TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 3);
//    TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 4);
//    TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
//  } if (poly_order == 2) {
//    // Solution (checked visually against g2):
//    const double sol[12] = {-0.9010465429057769, -0.4272439810948228,  0.0875367707148495,
//                            -0.9039382020247494,  0.4172269800703625,  0.08107082435707  ,
//                             0.9039382020247495,  0.4172269800703625, -0.0810708243570699,
//                             0.9010465429057768, -0.4272439810948229, -0.0875367707148495};
//    const double *phi_p;
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 1);
//    TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 2);
//    TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 3);
//    TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
//    phi_p = gkyl_array_cfetch(phi_mb_ho, 4);
//    TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
//    TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
//  }

  for (int bI=0; bI<num_blocks; bI++) {
    gkyl_array_release(rho[bI]);
    gkyl_array_release(phi[bI]);
    gkyl_array_release(rho_ho[bI]);
    gkyl_array_release(phi_ho[bI]);
  }
  gkyl_free(rho);
  gkyl_free(phi);
  gkyl_free(rho_ho);
  gkyl_free(phi_ho);
  gkyl_free(grid);
  gkyl_free(localRange);
  gkyl_free(localRange_ext);
  gkyl_free(local_in_mb_range);
  gkyl_array_release(rho_mb);
  gkyl_array_release(phi_mb);
  gkyl_array_release(jac_mb);
  gkyl_array_release(rho_mb_ho);
  gkyl_array_release(phi_mb_ho);
  gkyl_array_release(jac_mb_ho);

}


void test_1x_p1_nonperiodic() {test_1x(1, false, false);}

void test_1x_p2_nonperiodic() {test_1x(2, false, false);}

#ifdef GKYL_HAVE_CUDA
// ......... GPU tests ............ //
void gpu_test_1x_p1_nonperiodic() {test_1x(1, false, true);}

void gpu_test_1x_p2_nonperiodic() {test_1x(2, false, true);}
#endif


TEST_LIST = {
  { "test_1x_p1_nonperiodic", test_1x_p1_nonperiodic },
  { "test_1x_p2_nonperiodic", test_1x_p2_nonperiodic },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_1x_p1_nonperiodic", gpu_test_1x_p1_nonperiodic },
  { "gpu_test_1x_p2_nonperiodic", gpu_test_1x_p2_nonperiodic },
#endif
  { NULL, NULL },
};

