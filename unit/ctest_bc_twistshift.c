// Test updater that applies the
// twist shift BC
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_mat.h>
#include <gkyl_bc_twistshift.h>

void
test_bc_twistshift()
{
  bool use_gpu = false;

  int poly_order = 1;
  //double lower[] = {-2.0, -2.0, 0., 0.}, upper[] = {2.0, 2.0, 2.0, 5.};
  double lower[] = {-2.0, -2.0, 0.}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {6, 8, 4};
  //int cells[] = {6, 8, 4,2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = 3;
  int dir = 2; //z
  int do_dir = 0; //x
  int shift_dir = 1; //y
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  int vdim = ndim-cdim;
  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  //if (poly_order == 1)  // Force gkhybrid for p=1.
  //  gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  //else
  //  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  //int ghost[] = { 1, 1, 1, 1 };
  int ghost[] = { 1, 1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  
  enum gkyl_basis_type basis_type = basis.b_type;
  assert(basis_type==GKYL_BASIS_MODAL_SERENDIPITY);

  
  int ndonors[grid.cells[0]+2*ghost[0]];
  for(int i = 0; i< grid.cells[0]+2*ghost[0]; i++){
    ndonors[i] = 2;
  }
  int total_donors = 0;
  for(int j = 0; j <grid.cells[1] + 2*ghost[1]; j++){
    for(int i = 0; i< grid.cells[0] + 2*ghost[0]; i++){
      total_donors +=ndonors[i];
    }
  }
  printf("total donors = %d\n", total_donors);
  int cells_do[total_donors];
  for(int i = 0; i<total_donors; i++){
    cells_do[i] = 1;
  }

  struct gkyl_array *yshift = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_bc_twistshift *bcts = gkyl_bc_twistshift_new(dir, do_dir, shift_dir, edge, &local_ext, ghost, &basis, &grid, cdim, yshift, ndonors, cells_do, false);

  //prepare to do the advance
  struct gkyl_nmat *matsdo = gkyl_nmat_new(total_donors, basis.num_basis, basis.num_basis);
  // fill each matrix with its number
  for (size_t n=0; n<matsdo->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(matsdo, n);
    for (size_t j=0; j<matsdo->nc; ++j)
      for (size_t i=0; i<matsdo->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  struct gkyl_array *fdo = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *ftar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  printf("\n trying the advance \n");
  gkyl_bc_twistshift_mv(bcts, matsdo, fdo, ftar);
  printf("\n back in unit test\n");





  gkyl_array_release(yshift);
  gkyl_array_release(fdo);
  gkyl_array_release(ftar);
  gkyl_nmat_release(matsdo);
  gkyl_bc_twistshift_release(bcts);
}


TEST_LIST = {
  { "test_bc_twistshift", test_bc_twistshift },
  { NULL, NULL },
};
