// Test calculation of Vlasov moments of a distribution function.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>

void evalFunc_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = (x*x)*(vx-0.5)*(vx-0.5);
}

void evalFunc_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2];
  fout[0] = (x*x)*(vx-0.5)*(vx-0.5);
}

void evalFunc_2x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double vx = xn[2], vy = xn[3];
  fout[0] = (vx-1)*(vy-2);
}

// allocate array (filled with zeros)
  static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
  static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
  static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
    const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
        d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
        d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void test_bc_wall_1x1v()
{
  // Change these inputs for going to different dimensions
  int cdim = 1;
  int vdim = 1;
  
  int poly_order = 1;
  int ndim = cdim + vdim;
  double lower[ndim], upper[ndim];
  int cells[ndim];
  for (int ix=0; ix<ndim; ix++){
    lower[ix] =-2.0;
    upper[ix] = 2.0;
    cells[ix] = ix<cdim ? 4 : 2; 
  }
  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int ix=0; ix<cdim; ix++){
    confLower[ix] = lower[ix];
    confUpper[ix] = upper[ix];
    confCells[ix] = cells[ix];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[cdim];
  for (int d=0; d<cdim; d++){confGhost[d] = 1;}
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[ndim];
  for (int d=0; d<cdim; d++){ghost[d] = confGhost[d];}
  for (int d=cdim; d<ndim; d++){ghost[d] = 0;}

  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_1x1v, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local_ext, distf);

  // Make a DG equation object
  struct gkyl_dg_eqn* eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &local, GKYL_FIELD_E_B, false);

  //Create the boundary condition
  struct gkyl_array_copy_func* bc;
  for (int bc_dir=0; bc_dir<cdim; bc_dir++){
    bc =  gkyl_vlasov_wall_bc_create(eqn, bc_dir, &basis);
    // struct gkyl_array_copy_func* bc =  gkyl_vlasov_absorb_bc_create(eqn, bc_dir, &basis);
    // Determine the size of the BC buffer
    long buff_sz = 0;
    for (int d=0; d<cdim; ++d) {
      long vol = skin_ghost.lower_skin[d].volume;
      buff_sz = buff_sz > vol ? buff_sz : vol;
    }
    struct gkyl_array *bc_buffer;
    bc_buffer = mkarr(basis.num_basis, buff_sz);

    //Apply BC to the lower ghost cells
    gkyl_array_flip_copy_to_buffer_fn(&bc_buffer, distf, bc_dir, skin_ghost.lower_skin[bc_dir], bc);
    gkyl_array_copy_from_buffer(distf, &bc_buffer, skin_ghost.lower_ghost[bc_dir]);

    //Apply BC to the upper ghost cells
    gkyl_array_flip_copy_to_buffer_fn(&bc_buffer, distf, bc_dir, skin_ghost.upper_skin[bc_dir], bc);
    gkyl_array_copy_from_buffer(distf, &bc_buffer, skin_ghost.upper_ghost[bc_dir]);
  }

  // Check lower ghost cells after applying BC
  struct gkyl_range_iter iter, iter_skin;
  for (int d=0;d<cdim;d++) {
    gkyl_range_iter_init(&iter, &skin_ghost.lower_ghost[d]);
    while (gkyl_range_iter_next(&iter)) {
      // Find the index and value of f at the ghost and adjacent skin cells
      iter_skin = iter; iter_skin.idx[d] = iter.idx[d]+1;
      int linidx_ghost = gkyl_range_idx(skin_ghost.lower_ghost, iter.idx);
      int linidx_skin  = gkyl_range_idx(skin_ghost.lower_skin,  iter_skin.idx);
      const double *val_ghost = gkyl_array_cfetch(distf, linidx_ghost);
      const double *val_skin  = gkyl_array_cfetch(distf, linidx_skin);

      // Flip the skin value to manually apply wall BC to skin cell
      double val_correct[basis.num_basis];
      basis.flip_odd_sign(d,   val_skin,   val_correct);
      basis.flip_odd_sign(d+1, val_correct,val_correct);

      // Check values
      for (int i=0;i<basis.num_basis;i++){
        TEST_CHECK(gkyl_compare(val_ghost[i],val_correct[i],1e-12));
      }}

    // Check upper ghost cells after applying BC
    gkyl_range_iter_init(&iter, &skin_ghost.upper_ghost[d]);
    while (gkyl_range_iter_next(&iter)) {
      // Find the index and value of f at the ghost and adjacent skin cells
      iter_skin = iter; iter_skin.idx[d] = iter.idx[d]-1;
      int linidx_ghost = gkyl_range_idx(skin_ghost.upper_ghost, iter.idx);
      int linidx_skin  = gkyl_range_idx(skin_ghost.upper_skin,  iter_skin.idx);
      const double *val_ghost = gkyl_array_cfetch(distf, linidx_ghost);
      const double *val_skin  = gkyl_array_cfetch(distf, linidx_skin);

      // Flip the skin value to manually apply wall BC to skin cell
      double val_correct[basis.num_basis];
      basis.flip_odd_sign(d,  val_skin,   val_correct);
      basis.flip_odd_sign(d+1,val_correct,val_correct);

      // Check values
      for (int i=0;i<basis.num_basis;i++){
        TEST_CHECK(gkyl_compare(val_ghost[i],val_correct[i],1e-12));
      }}}

  // release memory for moment data object
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_bc_release(bc);
}

void test_bc_wall_1x2v()
{
  // Change these inputs for going to different dimensions
  int cdim = 1; // Make these actual inputs into the function
  int vdim = 2;
  
  int poly_order = 1;
  int ndim = cdim + vdim;
  double lower[ndim], upper[ndim];
  int cells[ndim];
  for (int ix=0; ix<ndim; ix++){
    lower[ix] =-2.0;
    upper[ix] = 2.0;
    cells[ix] = ix<cdim ? 4 : 2; 
  }
  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int ix=0; ix<cdim; ix++){
    confLower[ix] = lower[ix];
    confUpper[ix] = upper[ix];
    confCells[ix] = cells[ix];
  }
  printf("\n%i",confCells[0]);
  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[cdim];
  for (int d=0; d<cdim; d++){confGhost[d] = 1;}
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  printf("\nDo I get here 1?");
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  printf("\nDo I get here 2?");
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[ndim];
  for (int d=0; d<cdim; d++){ghost[d] = confGhost[d];}
  for (int d=cdim; d<ndim; d++){ghost[d] = 0;}

  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_1x1v, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local_ext, distf);

  // Make a DG equation object
  struct gkyl_dg_eqn* eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &local, GKYL_FIELD_E_B, false);

  //Create the boundary condition
  for (int bc_dir=0; bc_dir<cdim; bc_dir++){
    struct gkyl_array_copy_func *bc =  gkyl_vlasov_wall_bc_create(eqn, bc_dir, &basis);
    // struct gkyl_array_copy_func *bc =  gkyl_vlasov_absorb_bc_create(eqn, bc_dir, &basis);

    // Determine the size of the BC buffer
    long buff_sz = 0;
    printf("\nI'm printing");
    for (int d=0; d<cdim; ++d) {
      long vol = skin_ghost.lower_skin[d].volume;
      buff_sz = buff_sz > vol ? buff_sz : vol;
    }
    printf("\n %i %i %li",skin_ghost.lower_skin[0].lower[2],skin_ghost.lower_skin[0].upper[2],skin_ghost.lower_skin[0].volume);
    printf("\nI'm not printing");
    struct gkyl_array *bc_buffer = mkarr(basis.num_basis, buff_sz);//dealocate this

    //Apply BC to the lower ghost cells
//    gkyl_array_flip_copy_to_buffer_fn(&bc_buffer, distf, bc_dir, skin_ghost.lower_skin[bc_dir], bc);
//    gkyl_array_copy_from_buffer(distf, &bc_buffer, skin_ghost.lower_ghost[bc_dir]);
//
//    //Apply BC to the upper ghost cells
//    gkyl_array_flip_copy_to_buffer_fn(&bc_buffer, distf, bc_dir, skin_ghost.upper_skin[bc_dir], bc);
//    gkyl_array_copy_from_buffer(distf, &bc_buffer, skin_ghost.upper_ghost[bc_dir]);
    gkyl_vlasov_bc_release(bc);
    gkyl_array_release(bc_buffer);
  }

//  // Check lower ghost cells after applying BC
//  struct gkyl_range_iter iter, iter_skin;
//  for (int d=0;d<cdim;d++) {
//    gkyl_range_iter_init(&iter, &skin_ghost.lower_ghost[d]);
//    while (gkyl_range_iter_next(&iter)) {
//      // Find the index and value of f at the ghost and adjacent skin cells
//      iter_skin = iter; iter_skin.idx[d] = iter.idx[d]+1;
//      int linidx_ghost = gkyl_range_idx(skin_ghost.lower_ghost, iter.idx);
//      int linidx_skin  = gkyl_range_idx(skin_ghost.lower_skin,  iter_skin.idx);
//      const double *val_ghost = gkyl_array_cfetch(distf, linidx_ghost);
//      const double *val_skin  = gkyl_array_cfetch(distf, linidx_skin);
//
//      // Flip the skin value to manually apply wall BC to skin cell
//      double val_correct[basis.num_basis];
//      basis.flip_odd_sign(d,   val_skin,   val_correct);
//      basis.flip_odd_sign(d+1, val_correct,val_correct);
//
//      // Check values
//      for (int i=0;i<basis.num_basis;i++){
//        TEST_CHECK(gkyl_compare(val_ghost[i],val_correct[i],1e-12));
//      }}
//
//    // Check upper ghost cells after applying BC
//    gkyl_range_iter_init(&iter, &skin_ghost.upper_ghost[d]);
//    while (gkyl_range_iter_next(&iter)) {
//      // Find the index and value of f at the ghost and adjacent skin cells
//      iter_skin = iter; iter_skin.idx[d] = iter.idx[d]-1;
//      int linidx_ghost = gkyl_range_idx(skin_ghost.upper_ghost, iter.idx);
//      int linidx_skin  = gkyl_range_idx(skin_ghost.upper_skin,  iter_skin.idx);
//      const double *val_ghost = gkyl_array_cfetch(distf, linidx_ghost);
//      const double *val_skin  = gkyl_array_cfetch(distf, linidx_skin);
//
//      // Flip the skin value to manually apply wall BC to skin cell
//      double val_correct[basis.num_basis];
//      basis.flip_odd_sign(d,  val_skin,   val_correct);
//      basis.flip_odd_sign(d+1,val_correct,val_correct);
//
//      // Check values
//      for (int i=0;i<basis.num_basis;i++){
//        TEST_CHECK(gkyl_compare(val_ghost[i],val_correct[i],1e-12));
//      }}}
//
//  // release memory for moment data object
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_dg_eqn_release(eqn);
}



TEST_LIST = {
  { "test_bc_wall_1x1v", test_bc_wall_1x1v },
  //{ "test_bc_absorb_1x1v", test_bc_absorb_1x1v },
  { "test_bc_wall_1x2v", test_bc_wall_1x2v },
  //{ "test_bc_absorb_1x2v", test_bc_absorb_1x2v },
  //{ "test_bc_wall_2x2v", test_bc_wall_2x2v },
  //{ "test_bc_absorb_2x2v", test_bc_absorb_2x2v },
  //{ "test_bc_wall_3x2v", test_bc_wall_3x2v },
  //{ "test_bc_absorb_3x2v", test_bc_absorb_3x2v },
  { NULL, NULL },
};
