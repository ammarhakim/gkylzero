#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>

static void
test_wham(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[] = { 2.0e-6, 0.0, -2.0 };
  double cupper[] = { 3.0e-3, 2*M_PI, 2.0 };
  int cells[] = { 10, 16, 32 };

  const char *fname = "data/unit/wham_hires.geqdsk_psi.gkyl";
  
  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, 3, clower, cupper, cells);

  if (!gkyl_check_file_exists(fname)) {
    fprintf(stderr, "Unable to find file %s!\n", fname);
    goto cleanup;
  }
  
  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, fname);

  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrnodes = psi_grid.cells[0]-1, // cells and not nodes
        .nznodes = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = fl_coord,
        .include_axis = include_axis,
        .write_psi_cubic = false,
      }
    );

  TEST_CHECK( include_axis == gkyl_mirror_grid_gen_is_include_axis(geom) );
  TEST_CHECK( fl_coord == gkyl_mirror_grid_gen_fl_coord(geom) );  
  
  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 2, (int[2]) { cells[0]+1, cells[2]+1 });

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &node_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&node_range, iter.idx);

    const struct gkyl_mirror_grid_gen_geom *g =
      gkyl_array_cfetch(geom->node_geom, loc);

    /* // compute the tangent vectors */
    /* struct gkyl_vec3 e1 = gkyl_vec3_scale(g->Jc, */
    /*   gkyl_vec3_cross(g->e2d, g->e3d)); */

    /* struct gkyl_vec3 e2 = gkyl_vec3_scale(g->Jc, */
    /*   gkyl_vec3_cross(g->e3d, g->e2d)); */
    
    /* struct gkyl_vec3 e3 = gkyl_vec3_scale(g->Jc, */
    /*   gkyl_vec3_cross(g->e3d, g->e2d)); */
  }

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);

  cleanup:

  return;
}

static void
test_wham_no_axis_psi(void)
{
  test_wham(false, GKYL_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_wham_with_axis_psi(void)
{
  test_wham(true, GKYL_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_wham_no_axis_sqrt_psi(void)
{
  test_wham(false, GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

static void
test_wham_with_axis_sqrt_psi(void)
{
  test_wham(true, GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

TEST_LIST = {
  { "wham_no_axis_psi", test_wham_no_axis_psi },
  { "wham_with_axis_psi", test_wham_with_axis_psi },
  /* { "wham_no_axis_sqrt_psi", test_wham_no_axis_sqrt_psi }, */
  /* { "wham_with_axis_sqrt_psi", test_wham_with_axis_sqrt_psi }, */
  { NULL, NULL },
};
