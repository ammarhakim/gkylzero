#include "gkyl_array.h"
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
#include <gkyl_skin_surf_from_ghost.h>  // Include the custom updater

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
void eval_field_3x(double t, const double *xn, double* restrict fout, void *ctx) {
    fout[0] = 0.0;
}

// Function to set up and test the ghost-to-skin surf copy in 3D
void test_3x(int poly_order, bool use_gpu, enum gkyl_edge_loc edge, int dir, bool control) { 
    // Set up 3D grid parameters
    int cdim = 3;
    double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
    int cells[] = {2, 4, 6};  // Cell count in each dimension

    const int ndim = sizeof(cells)/sizeof(cells[0]);

    // Initialize grid and basis functions
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

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
    gkyl_proj_on_basis *proj_field = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_field_3x, NULL);
    gkyl_proj_on_basis_advance(proj_field, 0.0, &local, field_ho);

    // Set ghost cell values to 1.0 for the test case
    struct gkyl_range_iter iter_ghost;
    gkyl_range_iter_init(&iter_ghost, &ghost_r);
    while (gkyl_range_iter_next(&iter_ghost)) {
        long linidx = gkyl_range_idx(&ghost_r, iter_ghost.idx);
        double *f_i = gkyl_array_fetch(field_ho, linidx);
        f_i[0] = 1.0;
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

    // Check that field values in the skin cells meet the expected result
    struct gkyl_range_iter iter_skin;
    gkyl_range_iter_init(&iter_skin, &skin_r);
    while (gkyl_range_iter_next(&iter_skin)) {
        long linidx = gkyl_range_idx(&skin_r, iter_skin.idx);
        double *f_i = gkyl_array_fetch(field_ho, linidx);

        // Test values based on control mode: expect 0 if control, otherwise expect 0.5
        if (control) {
            TEST_CHECK(gkyl_compare(f_i[0], 0.0, 1e-14));
        } else {
            TEST_CHECK(gkyl_compare(f_i[0], 0.5, 1e-14));
        }
    }

    // Release allocated arrays
    gkyl_array_release(field);
    gkyl_array_release(field_ho);
}

// Tests for 3D case on CPU
void test_3x_ho() {
    bool use_gpu = false;
    // Loop over edges
    for (int edge = GKYL_LOWER_EDGE; edge <= GKYL_UPPER_EDGE; edge++) {
        // Loop over control states (perform the test or not)
        for (int control = 0; control <= 1; control++) {
            // Loop over directions
            for (int dir = 1; dir <= 2; dir++) {
                test_3x(1, use_gpu, edge, dir, control==1);
            }
        }
    }
}

// Tests for 3D case on GPU (if available)
void test_3x_dev() {
    bool use_gpu = true;
    // Loop over edges
    for (int edge = GKYL_LOWER_EDGE; edge <= GKYL_UPPER_EDGE; edge++) {
        // Loop over control states (perform the test or not)
        for (int control = 0; control <= 1; control++) {
            // Loop over directions
            for (int dir = 1; dir <= 2; dir++) {
                test_3x(use_gpu, true, edge, dir, control==1);
            }
        }
    }
}

// List of tests for the test framework
TEST_LIST = {
    { "test_3x_ho", test_3x_ho },
#ifdef GKYL_HAVE_CUDA
    { "test_3x_dev", test_3x_dev },
#endif
    { NULL, NULL },
};