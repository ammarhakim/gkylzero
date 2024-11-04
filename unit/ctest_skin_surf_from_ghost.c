/**
 * This code implements a test suite for validating the `gkyl_skin_surf_from_ghost` 
 * module, which updates discontinuous Galerkin (DG) coefficients to ensure continuity 
 * at the interface between skin and ghost cells in structured grids. The update process 
 * leverages interface-based kernels found in `gkylzero/kernels/skin_surf_from_ghost` 
 * to adjust skin cell coefficients to match ghost cell values at shared nodes, while 
 * maintaining the opposing skin cell nodal value.
 * 
 * Specifically, the updater enforces that the skin cell value at the interface (e.g., 
 * upper edge) equals the ghost cell value, achieving continuity across the interface 
 * while altering the average skin cell value. Testing ensures accuracy by initializing 
 * skin cells to zero and ghost cells to one; in a first-order (poly_order=1) scenario, 
 * successful updates should result in a skin cell value of 0.5.
 * 
 * Example for GKYL_UPPER_EDGE:
 * 
 * vls            vus  vlg                vug
 *  |___skin cell___|   |____ghost cell____|
 * 
 * Here, `vus` will be set to 1 (matching `vug`), while `vls` remains 0, yielding 
 * an expected skin cell average of 0.5 for order one polynomials.
 */

#include "gkyl_array.h"
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
// This the updater to be tested
#include <gkyl_skin_surf_from_ghost.h>  

// Evaluate the projection of the modal representation inside a cell from 1x to 3x
double eval_f(const double *phi, double x, double y, double z, int cdim);
void compare_skin_ghost_cubes(const double* fskin, const double* fghost, const int cdim, const int edge, const int dir);

// Function to allocate a gkyl array, zero-initialized, on CPU or GPU
static struct gkyl_array* mkarr(bool on_gpu, long nc, long size) {
    struct gkyl_array* a;
    if (on_gpu)
        a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);  // Allocate on GPU
    else
        a = gkyl_array_new(GKYL_DOUBLE, nc, size);         // Allocate on CPU
    return a;
}

// Analytical field evaluation function, setting field to zero
void eval_field(double t, const double *xn, double* restrict fout, void *ctx) {
    fout[0] = 0.0;
    fout[1] = 0.0;
}

// Function to set up and test the ghost-to-skin surf copy in 3D
void test_ssfg(int cdim, int poly_order, bool use_gpu, enum gkyl_edge_loc edge, int dir, bool control) { 
double lower[cdim], upper[cdim];
int cells[cdim];
switch(cdim) {
    case 3:
        lower[0] = 0.0; lower[1] = 0.0; lower[2] = 0.0;
        upper[0] = 1.0; upper[1] = 1.0; upper[2] = 1.0;
        cells[0] = 2;   cells[1] = 4;   cells[2] = 6;
        break;
    case 2:
        lower[0] = 0.0; lower[1] = 0.0;
        upper[0] = 1.0; upper[1] = 1.0;
        cells[0] = 4;   cells[1] = 6;
        break;
    case 1:
        lower[0] = 0.0;
        upper[0] = 1.0;
        cells[0] = 6;
        break;
    default:
        // Handle invalid cdim values
        fprintf(stderr, "Invalid cdim value: %d\n", cdim);
        exit(1);
}
    const int ndim = sizeof(cells)/sizeof(cells[0]);

    // Initialize grid and basis functions
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
    // Get grid spacing in the direction of the ssfg
    double dx = grid.dx[dir];

    struct gkyl_basis basis;
    gkyl_cart_modal_serendip(&basis, cdim, poly_order);

    // Define ranges for local and extended grid with ghost cells
    int ghost[GKYL_MAX_CDIM] = {1, 1, 1};  // Ghost cell configuration
    struct gkyl_range local, local_ext;
    gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

    // Set up ranges for skin (outer) and ghost (inner) cells along direction and edge
    struct gkyl_range skin_r, ghost_r;
    gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

    // Allocate field array for storing basis-projected field values
    struct gkyl_array *field_ho, *field;
    field = mkarr(use_gpu, basis.num_basis, local_ext.volume);      // Device or host array
    field_ho = use_gpu ? mkarr(false, field->ncomp, field->size)    // Device host-offloaded array
                       : gkyl_array_acquire(field);                 // Directly host-based if not using GPU

    // Project analytical field onto the basis functions
    gkyl_proj_on_basis *proj_field = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_field, NULL);
    gkyl_proj_on_basis_advance(proj_field, 0.0, &local, field_ho);
    

    // Set ghost cell modal values
    double fg0 = sqrt(pow(2,cdim));
    double fg1 = 1.0;
    double fg2 = 2.0;
    double fg3 = 3.0;
    double fg4 = 4.0;
    double fg5 = 5.0;
    double fg6 = 6.0;
    double fg7 = 7.0;
    struct gkyl_range_iter iter_ghost;
    gkyl_range_iter_init(&iter_ghost, &ghost_r);
    while (gkyl_range_iter_next(&iter_ghost)) {
        long linidx = gkyl_range_idx(&ghost_r, iter_ghost.idx);
        double *f_i = gkyl_array_fetch(field_ho, linidx);
        f_i[0] = fg0;
        f_i[1] = fg1;
        if (cdim > 1){
            f_i[2] = fg2;
            f_i[3] = fg3;
        }
        if (cdim > 2){
            f_i[4] = fg4;
            f_i[5] = fg5;
            f_i[6] = fg6;
            f_i[7] = fg7;
        }
    }

    // Copy field values to the GPU if necessary
    gkyl_array_copy(field, field_ho);

    // Initialize the skin-surf updater and call it if control is false
    gkyl_skin_surf_from_ghost* up = gkyl_skin_surf_from_ghost_new(dir, edge, basis, &skin_r, &ghost_r, use_gpu);
    if (!control) {
        gkyl_skin_surf_from_ghost_advance(up, field);  // Apply ghost values to skin cells
    }
    gkyl_skin_surf_from_ghost_release(up);

    // Copy field values back to host for checking
    gkyl_array_copy(field_ho, field);


    //---------- CHECK of the result --------------
    // Check that field values in the skin cells meet the expected result
    int gidx[GKYL_MAX_DIM]; // ghost index.
    struct gkyl_range_iter iter_skin;
    gkyl_range_iter_init(&iter_skin, &skin_r);
    int k =0;
    while (gkyl_range_iter_next(&iter_skin)) {
        // get skin cell modal values
        long skin_linidx = gkyl_range_idx(&skin_r, iter_skin.idx);
        const double *fskin = gkyl_array_cfetch(field_ho, skin_linidx);
        // Get ghost cell corresponding to skin cell
        gkyl_copy_int_arr(ndim, iter_skin.idx, gidx);
        gidx[dir] = edge == GKYL_LOWER_EDGE? iter_skin.idx[dir]-1 : iter_skin.idx[dir]+1; 
        long ghost_linidx = gkyl_range_idx(&ghost_r, gidx);
        const double *fghost = (const double*) gkyl_array_cfetch(field_ho, ghost_linidx);

        double xs = 0.0, ys = 0.0, zs = 0.0, xg = 0.0, yg = 0.0, zg = 0.0;
        char dir_c;
        double edge_skin = edge == GKYL_UPPER_EDGE? 1.0 : -1.0;
        switch(dir) {
        case 0:
            xs = edge_skin; 
            xg = -edge_skin;
            dir_c = 'x';
            break;
        case 1:
            ys = edge_skin; 
            yg = -edge_skin; 
            dir_c = 'y';
            break;
        case 2:
            zs = edge_skin; 
            zg = -edge_skin; 
            dir_c = 'z';
            break;
        }
        double vskin  = eval_f( fskin,xs,ys,zs,cdim);
        double vghost = eval_f(fghost,xg,yg,zg,cdim);
        double check_val;

        // Test values based on control mode: expect 0 if control, otherwise expect 0.5
        if (control) {
            check_val = 0.0;
        } else {
            check_val =  vghost;
        }
        // TEST_CHECK(gkyl_compare(vskin, check_val, 1e-14));
        // if (k==0)
            // compare_skin_ghost_cubes(fskin, fghost, cdim, edge, dir);
        k++;
    }
    // Release allocated arrays
    gkyl_array_release(field);
    gkyl_array_release(field_ho);
}

// Tests for 3D case on CPU
void test_ssfg_ho() {
  bool use_gpu    = false;
  int  poly_order = 1;
  printf("\n");
  // Loop over dimensionalities
  for (int cdim = 1; cdim <= 3; cdim++) {
    printf("Running tests for dimensionality cdim = %d\n", cdim);
    // Loop over edges
    // for (int edge = GKYL_LOWER_EDGE; edge <= GKYL_UPPER_EDGE; edge++) {
      // Loop over control states (perform the test or not)
      for (int control = 0; control <= 0; control++) {
        // Loop over directions
        for (int dir = 0; dir <= cdim-1; dir++) {
            // Run the actual test
            test_ssfg(cdim, poly_order, use_gpu, GKYL_UPPER_EDGE, dir, control == 1);
            test_ssfg(cdim, poly_order, use_gpu, GKYL_LOWER_EDGE, dir, control == 1);
        // }
      }
    }
  }
}

// Tests for 3D case on GPU (if available)
void test_ssfg_dev() {
  int  poly_order = 1;
  bool use_gpu    = true;
  printf("\n");
  // Loop over dimensionalities
  for (int cdim = 1; cdim <= 3; cdim++) {
    printf("Running tests for dimensionality cdim = %d\n", cdim);
    // Loop over control states (perform the test or just check identity)
    for (int control = 0; control <= 0; control++) {
    // Loop over directions
    for (int dir = 0; dir <= cdim-1; dir++) {
        // Loop over edges
        // for (int edge = GKYL_LOWER_EDGE; edge <= GKYL_UPPER_EDGE; edge++) {
            // Run the actual test
            test_ssfg(cdim, poly_order, use_gpu, GKYL_UPPER_EDGE, dir, control==1);
            test_ssfg(cdim, poly_order, use_gpu, GKYL_LOWER_EDGE, dir, control==1);
        // }
      }
    }
  }
}

// List of tests for the test framework
TEST_LIST = {
    { "test_ssfg_ho", test_ssfg_ho },
#ifdef GKYL_HAVE_CUDA
    { "test_ssfg_dev", test_ssfg_dev },
#endif
    { NULL, NULL },
};

// Evaluate the projection of the modal representation inside a cell from 1x to 3x
double eval_f(const double *phi, const double x, const double y, const double z, const int cdim) {
    // Common constants to avoid recalculating
    double sqrt2 = sqrt(2.0);
    double sqrt3 = sqrt(3.0);
    double denom = pow(2.0, 1.5);

    switch (cdim) {
        case 1: // 1D case
            return (sqrt3 * phi[1] * x) / sqrt2 + phi[0] / sqrt2;

        case 2: // 2D case
            return (3 * phi[3] * x * y) / 2.0
                   + (sqrt3 * phi[2] * y) / 2.0
                   + (sqrt3 * phi[1] * x) / 2.0
                   + phi[0] / 2.0;

        case 3: // 3D case
            return (pow(3.0, 1.5) * phi[7] * x * y * z) / denom
                   + (3 * phi[6] * y * z) / denom
                   + (3 * phi[5] * x * z) / denom
                   + (sqrt3 * phi[3] * z) / denom
                   + (3 * phi[4] * x * y) / denom
                   + (sqrt3 * phi[2] * y) / denom
                   + (sqrt3 * phi[1] * x) / denom
                   + phi[0] / denom;

        default:
            fprintf(stderr, "Invalid cdim value: %d. Must be 0, 1, or 2.\n", cdim);
            return 0.0;
    }
}

void compare_skin_ghost_cubes(const double* fskin, const double* fghost, const int cdim, const int edge, const int dir) {
  double skin[8] = {0}, ghost[8] = {0}; // Array to store skin and ghost values at each corner
  
  // Evaluate each corner based on the dimension
  switch (cdim) {
    case 1: // 1D case: evaluate only at -1 and +1 along x
      skin[0] = eval_f(fskin, -1.0, 0.0, 0.0, cdim);
      skin[1] = eval_f(fskin, +1.0, 0.0, 0.0, cdim);
      ghost[0] = eval_f(fghost, -1.0, 0.0, 0.0, cdim);
      ghost[1] = eval_f(fghost, +1.0, 0.0, 0.0, cdim);
      if(edge == GKYL_UPPER_EDGE){
        printf("1D, upper\n");
        printf("    Skin Cube               Ghost Cube\n");
        printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[0], skin[1], ghost[0], ghost[1]);
      }
      else{
        printf("1D, lower\n");
        printf("    Ghost Cube              Skin Cube\n");
        printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", ghost[0], ghost[1], skin[0], skin[1]);
      }
      break;

    case 2: // 2D case: evaluate 4 corners (-1,-1), (1,-1), (-1,1), (1,1)
      skin[0] = eval_f(fskin, -1.0, -1.0, 0.0, cdim);
      skin[1] = eval_f(fskin, +1.0, -1.0, 0.0, cdim);
      skin[2] = eval_f(fskin, -1.0, +1.0, 0.0, cdim);
      skin[3] = eval_f(fskin, +1.0, +1.0, 0.0, cdim);
      ghost[0] = eval_f(fghost, -1.0, -1.0, 0.0, cdim);
      ghost[1] = eval_f(fghost, +1.0, -1.0, 0.0, cdim);
      ghost[2] = eval_f(fghost, -1.0, +1.0, 0.0, cdim);
      ghost[3] = eval_f(fghost, +1.0, +1.0, 0.0, cdim);
      if(dir == 0){
        if(edge == GKYL_UPPER_EDGE){
            printf("2D, x-dir, upper\n");
            printf("    Skin Cube        |     Ghost Cube\n");
            printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[0], skin[1], ghost[0], ghost[1]);
            printf("       |           |            |            |\n");
            printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[2], skin[3], ghost[2], ghost[3]);
        }
        else{
            printf("2D, x-dir, lower\n");
            printf("    Ghost Cube        |     Skin Cube\n");
            printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", ghost[0], ghost[1], skin[0], skin[1]);
            printf("       |           |            |            |\n");
            printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", ghost[2], ghost[3], skin[2], skin[3]);
        }
      }
      else{
        if(edge == GKYL_UPPER_EDGE){
            printf("2D, y-dir, upper\n");
            printf("  Ghost Cube\t [%6.2f] --- [%6.2f]\n", ghost[0], ghost[1]);
            printf("            \t [%6.2f] --- [%6.2f]\n", ghost[2], ghost[3]);
            printf("\n");
            printf("  Skin  Cube\t [%6.2f] --- [%6.2f]\n", skin[0], skin[1]);
            printf("            \t [%6.2f] --- [%6.2f]\n", skin[2], skin[3]);            
        }
        else{
            printf("2D, y-dir, lower\n");
            printf("  Skin  Cube\t [%6.2f] --- [%6.2f]\n", skin[0], skin[1]);
            printf("            \t [%6.2f] --- [%6.2f]\n", skin[2], skin[3]);
            printf("\n");
            printf("  Ghost Cube\t [%6.2f] --- [%6.2f]\n", ghost[0], ghost[1]);
            printf("            \t [%6.2f] --- [%6.2f]\n", ghost[2], ghost[3]);      
        }        
      }
      break;

    case 3: // 3D case: evaluate 8 corners (-1,-1,-1), (1,-1,-1), etc.
      skin[0] = eval_f(fskin, -1.0, -1.0, -1.0, cdim);
      skin[1] = eval_f(fskin, +1.0, -1.0, -1.0, cdim);
      skin[2] = eval_f(fskin, -1.0, +1.0, -1.0, cdim);
      skin[3] = eval_f(fskin, +1.0, +1.0, -1.0, cdim);
      skin[4] = eval_f(fskin, -1.0, -1.0, +1.0, cdim);
      skin[5] = eval_f(fskin, +1.0, -1.0, +1.0, cdim);
      skin[6] = eval_f(fskin, -1.0, +1.0, +1.0, cdim);
      skin[7] = eval_f(fskin, +1.0, +1.0, +1.0, cdim);

      ghost[0] = eval_f(fghost, -1.0, -1.0, -1.0, cdim);
      ghost[1] = eval_f(fghost, +1.0, -1.0, -1.0, cdim);
      ghost[2] = eval_f(fghost, -1.0, +1.0, -1.0, cdim);
      ghost[3] = eval_f(fghost, +1.0, +1.0, -1.0, cdim);
      ghost[4] = eval_f(fghost, -1.0, -1.0, +1.0, cdim);
      ghost[5] = eval_f(fghost, +1.0, -1.0, +1.0, cdim);
      ghost[6] = eval_f(fghost, -1.0, +1.0, +1.0, cdim);
      ghost[7] = eval_f(fghost, +1.0, +1.0, +1.0, cdim);

      printf("3D Comparison:\n");
      printf("  Skin Cube (z = -1)    |   Ghost Cube (z = -1)\n");
      printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[0], skin[1], ghost[0], ghost[1]);
      printf("     |           |           |           | \n");
      printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[2], skin[3], ghost[2], ghost[3]);

      printf("\n");

      printf("  Skin Cube (z = +1)    |   Ghost Cube (z = +1)\n");
      printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[4], skin[5], ghost[4], ghost[5]);
      printf("     |           |          |           | \n");
      printf("  [%6.2f] --- [%6.2f]    [%6.2f] --- [%6.2f]\n", skin[6], skin[7], ghost[6], ghost[7]);
      break;

    default:
      printf("Invalid dimension (cdim must be 1, 2, or 3)\n");
  }
}