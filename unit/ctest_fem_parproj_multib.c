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

static void check_continuity(int num_blocks, struct gkyl_range *ranges,
  struct gkyl_basis basis, struct gkyl_array **fields)
{
  // Check continuity along last dim.
  int ndim = basis.ndim;
  int pardir = ndim-1;
  const int num_nodes_perp_max = 4; // 3x p=1.
  int num_nodes_perp = 1;
  if (ndim == 2)
    num_nodes_perp = 2;
  else if (ndim == 3)
    num_nodes_perp = 4;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));
  for (int bI=0; bI<num_blocks; bI++) {
    // Check continuity within each block.
    int idx_up[ndim];
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ranges[bI]);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[pardir] < ranges[bI].upper[pardir]) {
        int *idx_lo = iter.idx;
        for (int d=0; d<pardir; d++)
          idx_up[d] = idx_lo[d];
        idx_up[pardir] = idx_lo[pardir] + 1;
  
        long lidx_lo = gkyl_range_idx(&ranges[bI], idx_lo);
        long lidx_up = gkyl_range_idx(&ranges[bI], idx_up);
  
        double *arr_lo = gkyl_array_fetch(fields[bI], lidx_lo);
        double *arr_up = gkyl_array_fetch(fields[bI], lidx_up);
  
        double fn_lo[num_nodes_perp_max], fn_up[num_nodes_perp_max];
        for (int i=0; i<num_nodes_perp; i++) {
          const double *node_lo = gkyl_array_cfetch(nodes, num_nodes_perp+i);
          fn_lo[i] = basis.eval_expand(node_lo, arr_lo);
        }
        for (int i=0; i<num_nodes_perp; i++) {
          const double *node_up = gkyl_array_cfetch(nodes, i);
          fn_up[i] = basis.eval_expand(node_up, arr_up);
        }
        for (int i=0; i<num_nodes_perp; i++) {
          TEST_CHECK( gkyl_compare(fn_lo[i], fn_up[i], 1e-12) );
          TEST_MSG( "b%d idx_lo=%d, node %d: lower=%g upper=%g diff=%g\n", bI, idx_lo[0], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
        }
      }
    }

    // Check continuity between blocks.
    if (bI < num_blocks-1) {
      struct gkyl_range perp_range_lo, perp_range_up;
      gkyl_range_shorten_from_below(&perp_range_lo, &ranges[bI], pardir, 1);
      gkyl_range_shorten_from_above(&perp_range_up, &ranges[bI+1], pardir, 1);
      struct gkyl_range_iter perp_iter_lo;
      gkyl_range_iter_init(&perp_iter_lo, &perp_range_lo);
      while (gkyl_range_iter_next(&perp_iter_lo)) {
        int *idx_lo = perp_iter_lo.idx;
        for (int d=0; d<pardir; d++)
          idx_up[d] = idx_lo[d];
        idx_up[pardir] = perp_range_up.lower[pardir];
  
        long lidx_lo = gkyl_range_idx(&ranges[bI], idx_lo);
        long lidx_up = gkyl_range_idx(&ranges[bI+1], idx_up);
  
        double *arr_lo = gkyl_array_fetch(fields[bI], lidx_lo);
        double *arr_up = gkyl_array_fetch(fields[bI+1], lidx_up);
  
        double fn_lo[num_nodes_perp_max], fn_up[num_nodes_perp_max];
        for (int i=0; i<num_nodes_perp; i++) {
          const double *node_lo = gkyl_array_cfetch(nodes, num_nodes_perp+i);
          fn_lo[i] = basis.eval_expand(node_lo, arr_lo);
        }
        for (int i=0; i<num_nodes_perp; i++) {
          const double *node_up = gkyl_array_cfetch(nodes, i);
          fn_up[i] = basis.eval_expand(node_up, arr_up);
        }
        for (int i=0; i<num_nodes_perp; i++) {
          TEST_CHECK( gkyl_compare(fn_lo[i], fn_up[i], 1e-12) );
          TEST_MSG( "b%d-b%d idx_lo=%d, node %d: lower=%g upper=%g diff=%g\n", bI, bI+1, idx_lo[0], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
        }
      }
    }
  }
  gkyl_array_release(nodes);
}

void check_dirichlet_bc(int num_blocks, struct gkyl_range *ranges, struct gkyl_basis basis,
  struct gkyl_array **fields_dg, struct gkyl_array **fields_fem)
{
  // Check that two fields have the same boundary values in last dimension.
  int ndim = basis.ndim;
  int pardir = ndim-1;
  const int num_nodes_perp_max = 4; // 3x p=1.
  int num_nodes_perp = 1;
  if (ndim == 2)
    num_nodes_perp = 2;
  else if (ndim == 3)
    num_nodes_perp = 4;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));

  int lo_up_blocks[] = {0,num_blocks==1? 0 : num_blocks-1};

  for (int e=0; e<2; e++) {
    int bI = lo_up_blocks[e];

    struct gkyl_range perp_range;
    if (e == 0)
      gkyl_range_shorten_from_above(&perp_range, &ranges[bI], pardir, 1);
    else
      gkyl_range_shorten_from_below(&perp_range, &ranges[bI], pardir, 1);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&ranges[bI], iter.idx);
      double *arr_dg = gkyl_array_fetch(fields_dg[bI], lidx);
      double *arr_fem = gkyl_array_fetch(fields_fem[bI], lidx);

      double fn_dg[num_nodes_perp_max], fn_fem[num_nodes_perp_max];

      int off = e==0? 0 : num_nodes_perp;
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node = gkyl_array_cfetch(nodes, off+i);
        fn_dg[i] = basis.eval_expand(node, arr_dg);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node = gkyl_array_cfetch(nodes, off+i);
        fn_fem[i] = basis.eval_expand(node, arr_fem);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        TEST_CHECK( gkyl_compare(fn_dg[i], fn_fem[i], 1e-12) );
        TEST_MSG( "b%d idx=%d, node %d: dg=%g fem=%g diff=%g\n", bI, iter.idx[0], i, fn_dg[i], fn_fem[i], fn_dg[i]-fn_fem[i]);
      }
    }
  }
  gkyl_array_release(nodes);
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
test_1x(int num_blocks, const double *lower, const double *upper, const int *cells,
  int poly_order, bool is_dirichlet, bool use_gpu)
{
  int ndim = 1;

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

//    char fname0[1024];
//    sprintf(fname0, "ctest_fem_parproj_multib_%dx_p%d_rho_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, rho_ho[bI], fname0);
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
//  char fname1[1024];
//  sprintf(fname1, "ctest_fem_parproj_multib_%dx_p%d_rho_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, rho_mb_ho, fname1);

  // parallel FEM projection method.
  struct gkyl_fem_parproj_multib *parproj = gkyl_fem_parproj_multib_new(num_blocks, &localRange_mb,
    localRange, &basis, is_dirichlet? GKYL_FEM_PARPROJ_DIRICHLET : GKYL_FEM_PARPROJ_NONE, jac_mb, jac_mb, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_multib_set_rhs(parproj, rho_mb, rho_mb);

  // Solve the problem.
  gkyl_fem_parproj_multib_solve(parproj, phi_mb);
  gkyl_fem_parproj_multib_release(parproj);
  gkyl_array_copy(phi_mb_ho, phi_mb);

  for (int bI=0; bI<num_blocks; bI++) {
    // Copy solution back to local range of each block and write it out.
    gkyl_array_copy_range_to_range(phi[bI], phi_mb, &localRange[bI], &local_in_mb_range[bI]);
    gkyl_array_copy(phi_ho[bI], phi[bI]);

//    char fname2[1024];
//    sprintf(fname2, "ctest_fem_parproj_multib_%dx_p%d_phi_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, phi_ho[bI], fname2);
  }
  // Write out the multiblock source.
//  char fname3[1024];
//  sprintf(fname3, "ctest_fem_parproj_multib_%dx_p%d_phi_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, phi_mb_ho, fname3);

  bool check_values = false;
  if (num_blocks == 1) {
    if (cells[0] == 4) check_values = true;
  }
  if (num_blocks == 2) {
    if (cells[0] == 2 && cells[1] == 2) check_values = true;
  }

  if (check_values && (!is_dirichlet) ) {
    if (poly_order == 1) {
      // Solution (checked visually, also checked that phi is actually continuous,
      // and checked that visually looks like results in g2):
      const double sol[8] = {-0.9089542445638024, -0.4554124667453318,
                             -0.8488758876834943,  0.4900987222626481,
                              0.8488758876834943,  0.490098722262648 ,
                              0.9089542445638024, -0.4554124667453318};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_mb_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_MSG("Expected: %.13e in cell (%d)", sol[0], 1);
      TEST_MSG("Produced: %.13e", phi_p[0]);
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 2);
      TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 3);
      TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 4);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
    } if (poly_order == 2) {
      // Solution (checked visually against g2):
      const double sol[12] = {-0.9010465429057769, -0.4272439810948228,  0.0875367707148495,
                              -0.9039382020247494,  0.4172269800703625,  0.08107082435707  ,
                               0.9039382020247495,  0.4172269800703625, -0.0810708243570699,
                               0.9010465429057768, -0.4272439810948229, -0.0875367707148495};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_mb_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 2);
      TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 3);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_mb_ho, 4);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
    }
  }
  else {
    check_continuity(num_blocks, localRange, basis, phi_ho);
  }

  if (is_dirichlet)
    check_dirichlet_bc(num_blocks, localRange, basis, rho_ho, phi_ho);

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

void evalFunc2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

  struct test_ctx *tctx = ctx;
  int num_blocks = tctx->num_blocks;
  int ndim = tctx->ndim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++) Lx[d] = upper[ndim*(num_blocks-1)+d] - lower[ndim*0+d];

  double mu = .2;
  double sig = 0.3;

  fout[0] = exp(-(pow(x-mu,2))/(2.0*sig*sig))*sin(2.*M_PI*y);
}

void eval_weight_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

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
test_2x(int num_blocks, const double *lower, const double *upper, const int *cells,
  int poly_order, bool is_dirichlet, bool use_gpu)
{
  int ndim = 2;

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
      poly_order+1, 1, evalFunc2x, NULL);
    gkyl_proj_on_basis_advance(projob, 0.0, &localRange[bI], rho_ho[bI]);
    gkyl_proj_on_basis_release(projob);
    gkyl_array_copy(rho[bI], rho_ho[bI]);

//    char fname0[1024];
//    sprintf(fname0, "ctest_fem_parproj_multib_%dx_p%d_rho_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, rho_ho[bI], fname0);
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
    poly_order+1, 1, eval_weight_2x, NULL);
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
//  char fname1[1024];
//  sprintf(fname1, "ctest_fem_parproj_multib_%dx_p%d_rho_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, rho_mb_ho, fname1);

  // parallel FEM projection method.
  struct gkyl_fem_parproj_multib *parproj = gkyl_fem_parproj_multib_new(num_blocks, &localRange_mb,
    localRange, &basis, is_dirichlet? GKYL_FEM_PARPROJ_DIRICHLET : GKYL_FEM_PARPROJ_NONE, jac_mb, jac_mb, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_multib_set_rhs(parproj, rho_mb, rho_mb);

  // Solve the problem.
  gkyl_fem_parproj_multib_solve(parproj, phi_mb);
  gkyl_fem_parproj_multib_release(parproj);
  gkyl_array_copy(phi_mb_ho, phi_mb);

  for (int bI=0; bI<num_blocks; bI++) {
    // Copy solution back to local range of each block and write it out.
    gkyl_array_copy_range_to_range(phi[bI], phi_mb, &localRange[bI], &local_in_mb_range[bI]);
    gkyl_array_copy(phi_ho[bI], phi[bI]);

//    char fname2[1024];
//    sprintf(fname2, "ctest_fem_parproj_multib_%dx_p%d_phi_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, phi_ho[bI], fname2);
  }
  // Write out the multiblock source.
//  char fname3[1024];
//  sprintf(fname3, "ctest_fem_parproj_multib_%dx_p%d_phi_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, phi_mb_ho, fname3);

  bool check_values = false;
  if (num_blocks == 1) {
    if (cells[0] == 3 && cells[1] == 4) check_values = true;
  }
  if (num_blocks == 2) {
    if ((cells[0] == 3 && cells[1] == 2) && (cells[2] == 3 && cells[3] == 2)) check_values = true;
  }

  if (check_values && (!is_dirichlet) ) {
    if (poly_order == 1) {
      // Solution (checked continuity manually):
      const double sol[48] = {
        // idx = [0,:]
        -4.2253125086607479e-04, -4.2252954845312053e-04,
        -2.1170042428951191e-04, -2.1169957133119468e-04,
        -3.9460357085969891e-04, -3.9460198096965074e-04,
         2.2782447785903476e-04,  2.2782355993558733e-04,
         3.9460357085969891e-04,  3.9460198096965057e-04,
         2.2782447785903476e-04,  2.2782355993558747e-04,
         4.2253125086607485e-04,  4.2252954845312059e-04,
        -2.1170042428951191e-04, -2.1169957133119463e-04,
        // idx = [1,:]
        -6.2761887708076181e-01, -4.3547064325365575e-01,
        -3.1445527945628282e-01, -2.1818343555290270e-01,
        -5.8613569890365669e-01, -4.0668771950060645e-01,
         3.3840560354367566e-01,  2.3480126432979020e-01,
         5.8613569890365669e-01,  4.0668771950060634e-01,
         3.3840560354367566e-01,  2.3480126432979012e-01,
         6.2761887708076181e-01,  4.3547064325365564e-01,
        -3.1445527945628282e-01, -2.1818343555290265e-01,
        // idx = [2,:]
        -2.8612000924778641e-02,  2.8608472382826922e-02,
        -1.4335443172858760e-02,  1.4333675270894658e-02,
        -2.6720858424593201e-02,  2.6717563105605024e-02,
         1.5427294804416765e-02, -1.5425392251111881e-02,
         2.6720858424593201e-02, -2.6717563105605024e-02,
         1.5427294804416765e-02, -1.5425392251111872e-02,
         2.8612000924778638e-02, -2.8608472382826912e-02,
        -1.4335443172858762e-02,  1.4333675270894653e-02,
      };
      for (int k=0; k<cells[1]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {0+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx0); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(0*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(0*cells_mb[1]+k)*basis.num_basis+m], idx0[0], idx0[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx1[] = {1+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx1); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(1*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(1*cells_mb[1]+k)*basis.num_basis+m], idx1[0], idx1[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx2[] = {2+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx2); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(2*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(2*cells_mb[1]+k)*basis.num_basis+m], idx2[0], idx2[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    } if (poly_order == 2) {
      // Solution (checked continuity manually):
      const double sol[96] = {
        // idx = [0,:]
        -1.1330777557967131e-03, -1.5185686084598735e-03,
        -5.3726486726801399e-04, -7.2005081530163346e-04,
        -1.0232451194644426e-03,  1.1007862855941796e-04,
        -5.0054690727467069e-04,  1.4752910728100009e-04,
        -1.1367140547768834e-03, -1.5234420334777174e-03,
         5.2466835809768231e-04,  7.0316877582606832e-04,
        -9.4510889721781985e-04,  1.0194761685324392e-04,
         5.4565887622221867e-04,  1.3663179765785352e-04,
         1.1367140547768836e-03,  1.5234420334777181e-03,
         5.2466835809768199e-04,  7.0316877582606865e-04,
         9.4510889721782039e-04, -1.0194761685324384e-04,
         5.4565887622221878e-04, -1.3663179765785401e-04,
         1.1330777557967133e-03,  1.5185686084598746e-03,
        -5.3726486726801410e-04, -7.2005081530163281e-04,
         1.0232451194644423e-03, -1.1007862855941822e-04,
        -5.0054690727467112e-04, -1.4752910728100011e-04,
        // idx = [1,:]
        -6.7690840335622460e-01, -2.4487341055248341e-01,
        -3.2096570744703012e-01, -1.1611019609633033e-01,
         3.1090049018637228e-01,  6.5761725813252431e-02,
         1.5208504382060301e-01,  2.3789478759409195e-02,
        -6.7908075324501238e-01, -2.4565926388735609e-01,
         3.1344046668863568e-01,  1.1338792028973634e-01,
         2.8715975657749776e-01,  6.0904203791008561e-02,
        -1.6579176276044591e-01,  2.2032257282288872e-02,
         6.7908075324501249e-01,  2.4565926388735596e-01,
         3.1344046668863534e-01,  1.1338792028973639e-01,
        -2.8715975657749782e-01, -6.0904203791008568e-02,
        -1.6579176276044558e-01, -2.2032257282288945e-02,
         6.7690840335622438e-01,  2.4487341055248346e-01,
        -3.2096570744703029e-01, -1.1611019609633018e-01,
        -3.1090049018637222e-01, -6.5761725813252764e-02,
         1.5208504382060276e-01, -2.3789478759409195e-02,
        // idx = [2,:]
        -4.3172337009889768e-02,  5.7316620628891347e-02,
        -2.0470775103125258e-02,  2.7177487526246415e-02,
        -3.8162615108747848e-02,  4.1941972873804139e-03,
        -1.8668233644930851e-02, -5.5683159961537557e-03,
        -4.3310886658612019e-02,  5.7500562435244373e-02,
         1.9990825041209737e-02, -2.6540294417767608e-02,
        -3.5248472134667398e-02,  3.8843908545787469e-03,
         2.0350714875473255e-02, -5.1570096132443387e-03,
         4.3310886658612019e-02, -5.7500562435244360e-02,
         1.9990825041209734e-02, -2.6540294417767608e-02,
         3.5248472134667398e-02, -3.8843908545787400e-03,
         2.0350714875473234e-02,  5.1570096132443335e-03,
         4.3172337009889775e-02, -5.7316620628891354e-02,
        -2.0470775103125262e-02,  2.7177487526246422e-02,
         3.8162615108747848e-02, -4.1941972873804174e-03,
        -1.8668233644930830e-02,  5.5683159961537644e-03,
      };
      for (int k=0; k<cells_mb[1]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {0+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx0); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(0*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(0*cells_mb[1]+k)*basis.num_basis+m], idx0[0], idx0[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx1[] = {1+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx1); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(1*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(1*cells_mb[1]+k)*basis.num_basis+m], idx1[0], idx1[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx2[] = {2+1,k+1};
        linidx = gkyl_range_idx(&localRange_mb, idx2); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(2*cells_mb[1]+k)*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(2*cells_mb[1]+k)*basis.num_basis+m], idx2[0], idx2[1]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    }
  }
  else {
    check_continuity(num_blocks, localRange, basis, phi_ho);
  }

  if (is_dirichlet)
    check_dirichlet_bc(num_blocks, localRange, basis, rho_ho, phi_ho);

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

void evalFunc3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_ctx *tctx = ctx;
  int num_blocks = tctx->num_blocks;
  int ndim = tctx->ndim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++) Lx[d] = upper[ndim*(num_blocks-1)+d] - lower[ndim*0+d];

  double mu[2] = {.2, 0.2};
  double sig = 0.3;

  fout[0] = exp(-(pow(x-mu[0],2)+pow(y-mu[1],2))/(2.0*sig*sig))*sin(2.*M_PI*z);
}

void eval_weight_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

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
test_3x(int num_blocks, const double *lower, const double *upper, const int *cells,
  int poly_order, bool is_dirichlet, bool use_gpu)
{
  int ndim = 3;

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
      poly_order+1, 1, evalFunc3x, NULL);
    gkyl_proj_on_basis_advance(projob, 0.0, &localRange[bI], rho_ho[bI]);
    gkyl_proj_on_basis_release(projob);
    gkyl_array_copy(rho[bI], rho_ho[bI]);

//    char fname0[1024];
//    sprintf(fname0, "ctest_fem_parproj_multib_%dx_p%d_rho_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, rho_ho[bI], fname0);
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
    poly_order+1, 1, eval_weight_3x, NULL);
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
//  char fname1[1024];
//  sprintf(fname1, "ctest_fem_parproj_multib_%dx_p%d_rho_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, rho_mb_ho, fname1);

  // parallel FEM projection method.
  struct gkyl_fem_parproj_multib *parproj = gkyl_fem_parproj_multib_new(num_blocks, &localRange_mb,
    localRange, &basis, is_dirichlet? GKYL_FEM_PARPROJ_DIRICHLET : GKYL_FEM_PARPROJ_NONE, jac_mb, jac_mb, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_multib_set_rhs(parproj, rho_mb, rho_mb);

  // Solve the problem.
  gkyl_fem_parproj_multib_solve(parproj, phi_mb);
  gkyl_fem_parproj_multib_release(parproj);
  gkyl_array_copy(phi_mb_ho, phi_mb);

  for (int bI=0; bI<num_blocks; bI++) {
    // Copy solution back to local range of each block and write it out.
    gkyl_array_copy_range_to_range(phi[bI], phi_mb, &localRange[bI], &local_in_mb_range[bI]);
    gkyl_array_copy(phi_ho[bI], phi[bI]);

//    char fname2[1024];
//    sprintf(fname2, "ctest_fem_parproj_multib_%dx_p%d_phi_b%d.gkyl", ndim, poly_order, bI);
//    gkyl_grid_sub_array_write(&grid[bI], &localRange[bI], 0, phi_ho[bI], fname2);
  }
  // Write out the multiblock source.
//  char fname3[1024];
//  sprintf(fname3, "ctest_fem_parproj_multib_%dx_p%d_phi_mb.gkyl", ndim, poly_order);
//  gkyl_grid_sub_array_write(&grid_mb, &localRange_mb, 0, phi_mb_ho, fname3);

  bool check_values = false;
  if (num_blocks == 1) {
    if (cells[0] == 3 && cells[1] == 3 && cells[2] == 4) check_values = true;
  }
  if (num_blocks == 2) {
    if ((cells[0] == 3 && cells[1] == 3 && cells[2] == 2) &&
        (cells[3] == 3 && cells[4] == 3 && cells[5] == 2)) check_values = true;
  }

  if (check_values && (!is_dirichlet) ) {
    if (poly_order == 1) {
      // Solution (checked visually, also checked that phi is actually continuous,
      // and checked that visually looks like results in g2):
      const double sol[96] = {
        // idx = [0,1,:]
        -2.9175130738000619e-04, -2.9175013189013625e-04,
        -2.0243038272814358e-04, -1.4617587558971766e-04,
        -2.0242956711958172e-04, -1.4617528663516419e-04,
        -1.0142349903065383e-04, -1.0142309038708445e-04,
        -2.7246767537113221e-04, -2.7246657757659220e-04,
        -1.8905051806532438e-04,  1.5730928572099475e-04,
        -1.8904975636533129e-04,  1.5730865190902163e-04,
         1.0914836749545332e-04,  1.0914792772775724e-04,
         2.7246767537113210e-04,  2.7246657757659220e-04,
         1.8905051806532432e-04,  1.5730928572099478e-04,
         1.8904975636533134e-04,  1.5730865190902163e-04,
         1.0914836749545326e-04,  1.0914792772775718e-04,
         2.9175130738000619e-04,  2.9175013189013631e-04,
         2.0243038272814358e-04, -1.4617587558971761e-04,
         2.0242956711958167e-04, -1.4617528663516419e-04,
        -1.0142349903065386e-04, -1.0142309038708456e-04,
        // idx = [1,0,:]
        -2.9175130738000619e-04, -2.0243038272814342e-04,
        -2.9175013189013631e-04, -1.4617587558971766e-04,
        -2.0242956711958167e-04, -1.0142349903065388e-04,
        -1.4617528663516419e-04, -1.0142309038708458e-04,
        -2.7246767537113221e-04, -1.8905051806532432e-04,
        -2.7246657757659220e-04,  1.5730928572099475e-04,
        -1.8904975636533140e-04,  1.0914836749545325e-04,
         1.5730865190902168e-04,  1.0914792772775730e-04,
         2.7246767537113210e-04,  1.8905051806532430e-04,
         2.7246657757659220e-04,  1.5730928572099478e-04,
         1.8904975636533134e-04,  1.0914836749545328e-04,
         1.5730865190902160e-04,  1.0914792772775720e-04,
         2.9175130738000625e-04,  2.0243038272814347e-04,
         2.9175013189013631e-04, -1.4617587558971761e-04,
         2.0242956711958167e-04, -1.0142349903065384e-04,
        -1.4617528663516416e-04, -1.0142309038708461e-04,
        // idx = [1,2,:]
        -1.9756145041229078e-02, -1.3707715786581542e-02,
         1.9753708637470873e-02, -9.8984022577751872e-03,
         1.3706025298371970e-02, -6.8679635934884792e-03,
         9.8971815487547514e-03,  6.8671166098145447e-03,
        -1.8450340332725093e-02, -1.2801688837518459e-02,
         1.8448064965823664e-02,  1.0652308957739044e-02,
         1.2800110084035459e-02,  7.3910584964231056e-03,
        -1.0650995274046001e-02, -7.3901470026747126e-03,
         1.8450340332725097e-02,  1.2801688837518463e-02,
        -1.8448064965823671e-02,  1.0652308957739044e-02,
        -1.2800110084035459e-02,  7.3910584964231134e-03,
        -1.0650995274045996e-02, -7.3901470026747221e-03,
         1.9756145041229078e-02,  1.3707715786581542e-02,
        -1.9753708637470870e-02, -9.8984022577751907e-03,
        -1.3706025298371980e-02, -6.8679635934884914e-03,
         9.8971815487547583e-03,  6.8671166098145560e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,2,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx0); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx1[] = {2,1,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx1); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[32+k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[32+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
        
        int idx2[] = {2,3,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx2); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[64+k*basis.num_basis+m], phi_p[m], 1e-14) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[64+k*basis.num_basis+m], idx2[0], idx0[2], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    } if (poly_order == 2) {
      // Solution (checked visually against g2):
      const double sol[240] = {
        // idx = [0,1,:]
        -8.5122112791349098e-04, -1.1408199279301246e-03, -3.0793150106132720e-04, -4.0361855482954289e-04, -4.1269463518755719e-04, -5.4093592784822922e-04,
        -1.4601008289058379e-04, -7.6870970263647495e-04, 3.9096141311460275e-04, 8.2696226169962582e-05, -1.9568500684258495e-04, -2.7808277408889083e-04,
        5.2397262768387270e-04, -3.7603430197459076e-04, 1.9124892215530208e-04, 1.1083078143344461e-04, 2.9915579185682800e-05, -1.3603140624741703e-04,
        2.5631480990695021e-04, 4.0093329184926548e-05, -8.5395288617411710e-04, -1.1444810732657105e-03, -3.0891972185865852e-04, 3.9415546662673637e-04,
        -4.1401906422464335e-04, 5.2825334837797448e-04, 1.4258678563043659e-04, -7.1001010952258183e-04, 3.6110713159018460e-04, 7.6587829001101642e-05,
        1.9109705007607910e-04, -2.5684804056722943e-04, 4.8396145058753455e-04, 4.0992452786021841e-04, -2.0848529962988769e-04, 1.0264421158757013e-04,
        2.7705850306071661e-05, 1.4829128536231821e-04, -2.7941527377411485e-04, 3.7131815826626999e-05, 8.5395288617411884e-04, 1.1444810732657098e-03,
        3.0891972185865814e-04, 3.9415546662673609e-04, 4.1401906422464346e-04, 5.2825334837797405e-04, 1.4258678563043765e-04, 7.1001010952258291e-04,
        -3.6110713159018493e-04, -7.6587829001100328e-05, 1.9109705007607869e-04, 2.5684804056722867e-04, -4.8396145058753423e-04, 4.0992452786021858e-04,
        -2.0848529962988585e-04, -1.0264421158756809e-04, -2.7705850306071424e-05, 1.4829128536231699e-04, -2.7941527377411458e-04, -3.7131815826627277e-05,
        8.5122112791348968e-04, 1.1408199279301255e-03, 3.0793150106132769e-04, -4.0361855482954419e-04, 4.1269463518755649e-04, -5.4093592784823020e-04,
        -1.4601008289058344e-04, 7.6870970263647463e-04, -3.9096141311460134e-04, -8.2696226169961281e-05, -1.9568500684258492e-04, 2.7808277408888969e-04,
        -5.2397262768387270e-04, -3.7603430197459179e-04, 1.9124892215530127e-04, -1.1083078143344442e-04, -2.9915579185681851e-05, -1.3603140624741602e-04,
        2.5631480990694972e-04, -4.0093329184926751e-05,
        // idx = [1,0,:]
        -8.5122112791349153e-04, -3.0793150106132720e-04, -1.1408199279301255e-03, -4.0361855482954289e-04, -4.1269463518755741e-04, -1.4601008289058415e-04,
        -5.4093592784822966e-04, 3.9096141311460253e-04, -7.6870970263647528e-04, 8.2696226169961824e-05, -1.9568500684258636e-04, 5.2397262768387226e-04,
        -2.7808277408888986e-04, 1.9124892215530219e-04, -3.7603430197459168e-04, 2.9915579185682231e-05, 1.1083078143344385e-04, 2.5631480990695086e-04,
        -1.3603140624741619e-04, 4.0093329184925192e-05, -8.5395288617411776e-04, -3.0891972185865824e-04, -1.1444810732657109e-03, 3.9415546662673685e-04,
        -4.1401906422464319e-04, 1.4258678563043716e-04, 5.2825334837797470e-04, 3.6110713159018574e-04, -7.1001010952258345e-04, 7.6587829001101182e-05,
        1.9109705007608045e-04, 4.8396145058753482e-04, -2.5684804056722797e-04, -2.0848529962988704e-04, 4.0992452786021874e-04, 2.7705850306071156e-05,
        1.0264421158756911e-04, -2.7941527377411507e-04, 1.4829128536231762e-04, 3.7131815826625535e-05, 8.5395288617411971e-04, 3.0891972185865786e-04,
        1.1444810732657100e-03, 3.9415546662673626e-04, 4.1401906422464341e-04, 1.4258678563043749e-04, 5.2825334837797459e-04, -3.6110713159018422e-04,
        7.1001010952258291e-04, -7.6587829001100924e-05, 1.9109705007608010e-04, -4.8396145058753303e-04, 2.5684804056722867e-04, -2.0848529962988666e-04,
        4.0992452786021901e-04, -2.7705850306071342e-05, -1.0264421158756869e-04, -2.7941527377411393e-04, 1.4829128536231675e-04, -3.7131815826626545e-05,
        8.5122112791348979e-04, 3.0793150106132747e-04, 1.1408199279301255e-03, -4.0361855482954446e-04, 4.1269463518755681e-04, -1.4601008289058388e-04,
        -5.4093592784823074e-04, -3.9096141311460242e-04, 7.6870970263647539e-04, -8.2696226169961674e-05, -1.9568500684258601e-04, -5.2397262768387172e-04,
        2.7808277408888893e-04, 1.9124892215530097e-04, -3.7603430197459184e-04, -2.9915579185682275e-05, -1.1083078143344484e-04, 2.5631480990694907e-04,
        -1.3603140624741627e-04, -4.0093329184925863e-05,
        // idx = [1,2,:]
        -3.2433083445698618e-02, -1.1732753971887125e-02, 4.3058932372746953e-02, -1.5378605910671510e-02, 1.5576683008490425e-02, -5.5632514830916188e-03,
        2.0417002688116338e-02, 1.4896345637796328e-02, -2.8669545501892997e-02, 3.1508776228205857e-03, 7.3859048827129778e-03, -1.9776742488664067e-02,
        -1.0371284137159195e-02, 7.2869340853499104e-03, -1.4024452265605638e-02, 1.1398383384044293e-03, -4.1831800144129828e-03, -9.6743068697453770e-03,
        -5.0733828098187637e-03, -1.5132764669568534e-03, -3.2537168436913032e-02, -1.1770406993565723e-02, 4.3197118080720860e-02, 1.5018044924493976e-02,
        1.5626672054915903e-02, 5.4328175898799803e-03, -1.9938313353934981e-02, 1.3758842852513904e-02, -2.6480304687123233e-02, 2.9181365070276312e-03,
        -7.2127377462512807e-03, -1.8266566757537843e-02, -9.5793204650059353e-03, -7.9436716246366681e-03, 1.5288411039333878e-02, 1.0556436223728178e-03,
        -3.8741873778643924e-03, 1.0546207234634739e-02, 5.5306232491248520e-03, -1.4014975610190871e-03, 3.2537168436913053e-02, 1.1770406993565723e-02,
        -4.3197118080720881e-02, 1.5018044924493958e-02, -1.5626672054915910e-02, 5.4328175898799742e-03, -1.9938313353934967e-02, -1.3758842852513912e-02,
        2.6480304687123216e-02, -2.9181365070276317e-03, -7.2127377462513102e-03, 1.8266566757537846e-02, 9.5793204650059301e-03, -7.9436716246366473e-03,
        1.5288411039333913e-02, -1.0556436223727863e-03, 3.8741873778644024e-03, 1.0546207234634748e-02, 5.5306232491248806e-03, 1.4014975610190752e-03,
        3.2433083445698618e-02, 1.1732753971887144e-02, -4.3058932372746966e-02, -1.5378605910671526e-02, -1.5576683008490423e-02, -5.5632514830916543e-03,
        2.0417002688116349e-02, -1.4896345637796302e-02, 2.8669545501893032e-02, -3.1508776228206014e-03, 7.3859048827130082e-03, 1.9776742488664040e-02,
        1.0371284137159235e-02, 7.2869340853499148e-03, -1.4024452265605647e-02, -1.1398383384044384e-03, 4.1831800144129993e-03, -9.6743068697453978e-03,
        -5.0733828098187663e-03, 1.5132764669568367e-03};
      for (int k=0; k<cells[2]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {1,2,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx0); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[k*basis.num_basis+m], idx0[0], idx0[1], idx0[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }

        int idx1[] = {2,1,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx1); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[80+k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[80+k*basis.num_basis+m], idx1[0], idx1[1], idx1[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }

        int idx2[] = {2,3,k+1};
        linidx= gkyl_range_idx(&localRange_mb, idx2); 
        phi_p = gkyl_array_cfetch(phi_mb_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[160+k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[160+k*basis.num_basis+m], idx2[0], idx2[1], idx2[2]);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    }
  }
  else {
    check_continuity(num_blocks, localRange, basis, phi_ho);
  }

  if (is_dirichlet)
    check_dirichlet_bc(num_blocks, localRange, basis, rho_ho, phi_ho);

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

void test_1x_p_bc_dev(int poly_order, bool is_dirichlet, bool use_gpu) {
  // List grid extents and number of cells for each block in flat arrays.
  int num_blocks0 = 1;
  double lower0[] = {-0.5}, upper0[] = {0.5};
  int cells0[] = {4};
  test_1x(num_blocks0, lower0, upper0, cells0, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks1 = 2;
  double lower1[] = {-0.5, 0.0}, upper1[] = {0.0, 0.5};
  int cells1[] = {2, 2};
  test_1x(num_blocks1, lower1, upper1, cells1, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks2 = 2;
  double lower2[] = {-0.5, 0.0}, upper2[] = {0.0, 0.5};
  int cells2[] = {2, 4};
  test_1x(num_blocks2, lower2, upper2, cells2, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks3 = 4;
  double lower3[] = {-0.5, -0.25, 0.0, 0.25}, upper3[] = {-0.25, 0.0, 0.25, 0.5};
  int cells3[] = {2, 6, 2, 4};
  test_1x(num_blocks3, lower3, upper3, cells3, poly_order, is_dirichlet, use_gpu);
}

void test_1x_p1_bcnone() {
  test_1x_p_bc_dev(1, false, false);
}

void test_1x_p1_bcdirichlet() {
  test_1x_p_bc_dev(1, true, false);
}

void test_2x_p_bc_dev(int poly_order, bool is_dirichlet, bool use_gpu) {
  // List grid extents and number of cells for each block in flat arrays.
  int num_blocks0 = 1;
  double lower0[] = {-2.0, -0.5}, upper0[] = {2.0, 0.5};
  int cells0[] = {3, 4};
  test_2x(num_blocks0, lower0, upper0, cells0, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks1 = 2;
  double lower1[] = {-2.0, -0.5, -2.0, 0.0}, upper1[] = {2.0, 0.0, 2.0, 0.5};
  int cells1[] = {3, 2, 3, 2};
  test_2x(num_blocks1, lower1, upper1, cells1, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks2 = 2;
  double lower2[] = {-2.0, -0.5, -2.0, 0.0}, upper2[] = {2.0, 0.0, 2.0, 0.5};
  int cells2[] = {3, 2, 3, 4};
  test_2x(num_blocks2, lower2, upper2, cells2, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks3 = 4;
  double lower3[] = {-2.0, -0.5, -2.0, -0.25, -2.0, 0.0, -2.0, 0.25},
         upper3[] = {2.0, -0.25, 2.0, 0.0, 2.0, 0.25, 2.0, 0.5};
  int cells3[] = {3, 2, 3, 6, 3, 2, 3, 4};
  test_2x(num_blocks3, lower3, upper3, cells3, poly_order, is_dirichlet, use_gpu);
}

void test_2x_p1_bcnone() {
  test_2x_p_bc_dev(1, false, false);
}

void test_2x_p1_bcdirichlet() {
  test_2x_p_bc_dev(1, true, false);
}

void test_3x_p_bc_dev(int poly_order, bool is_dirichlet, bool use_gpu) {
  // List grid extents and number of cells for each block in flat arrays.
  int num_blocks0 = 1;
  double lower0[] = {-2.0, -2.0, -0.5}, upper0[] = {2.0, 2.0, 0.5};
  int cells0[] = {3, 3, 4};
  test_3x(num_blocks0, lower0, upper0, cells0, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks1 = 2;
  double lower1[] = {-2.0, -2.0, -0.5,
                     -2.0, -2.0, 0.0},
         upper1[] = {2.0, 2.0, 0.0,
                     2.0, 2.0, 0.5};
  int cells1[] = {3, 3, 2,
                  3, 3, 2};
  test_3x(num_blocks1, lower1, upper1, cells1, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks2 = 2;
  double lower2[] = {-2.0, -2.0, -0.5,
                     -2.0, -2.0, 0.0},
         upper2[] = {2.0, 2.0, 0.0,
                     2.0, 2.0, 0.5};
  int cells2[] = {3, 3, 2,
                  3, 3, 4};
  test_3x(num_blocks2, lower2, upper2, cells2, poly_order, is_dirichlet, use_gpu);
 
  int num_blocks3 = 4;
  double lower3[] = {-2.0, -2.0, -0.5,
                     -2.0, -2.0, -0.25,
                     -2.0, -2.0, 0.0,
                     -2.0, -2.0, 0.25},
         upper3[] = {2.0, 2.0, -0.25,
                     2.0, 2.0, 0.0,
                     2.0, 2.0, 0.25,
                     2.0, 2.0, 0.5};
  int cells3[] = {3, 3, 2,
                  3, 3, 6,
                  3, 3, 2,
                  3, 3, 4};
  test_3x(num_blocks3, lower3, upper3, cells3, poly_order, is_dirichlet, use_gpu);
}

void test_3x_p1_bcnone() {
  test_3x_p_bc_dev(1, false, false);
}

void test_3x_p1_bcdirichlet() {
  test_3x_p_bc_dev(1, true, false);
}

#ifdef GKYL_HAVE_CUDA
// ......... GPU tests ............ //
void gpu_test_1x_p1_bcnone() {
  test_1x_p_bc_dev(1, false, true);
}

void gpu_test_1x_p1_bcdirchlet() {
  test_1x_p_bc_dev(1, true, true);
}

void gpu_test_2x_p1_bcnone() {
  test_2x_p_bc_dev(1, false, true);
}

void gpu_test_2x_p1_bcdirichlet() {
  test_2x_p_bc_dev(1, true, true);
}

void test_3x_p1_bcnone() {
  test_3x_p_bc_dev(1, false, true);
}

void test_3x_p1_bcdirichlet() {
  test_3x_p_bc_dev(1, true, true);
}
#endif


TEST_LIST = {
  { "test_1x_p1_bcnone", test_1x_p1_bcnone },
  { "test_1x_p1_bcdirichlet", test_1x_p1_bcdirichlet },
  { "test_2x_p1_bcnone", test_2x_p1_bcnone },
  { "test_2x_p1_bcdirichlet", test_2x_p1_bcdirichlet },
  { "test_3x_p1_bcnone", test_3x_p1_bcnone },
  { "test_3x_p1_bcdirichlet", test_3x_p1_bcdirichlet },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_1x_p1_bcnone", gpu_test_1x_p1_bcnone },
  { "gpu_test_1x_p1_bcdirichlet", gpu_test_1x_p1_bcdirichlet },
  { "gpu_test_2x_p1_bcnone", gpu_test_2x_p1_bcnone },
  { "gpu_test_2x_p1_bcdirichlet", gpu_test_2x_p1_bcdirichlet },
  { "gpu_test_3x_p1_bcnone", gpu_test_3x_p1_bcnone },
  { "gpu_test_3x_p1_bcdirichlet", gpu_test_3x_p1_bcdirichlet },
#endif
  { NULL, NULL },
};

