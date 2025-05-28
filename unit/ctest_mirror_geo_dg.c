#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_mirror_geo_dg.h>

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
  struct gkyl_mirror_grid_gen *mirror_grid =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = fl_coord,
        .include_axis = include_axis,
        .write_psi_cubic = false,
      }
    );

  struct gkyl_mirror_geo_gen *mirror_geo = 
    gkyl_mirror_geo_gen_inew(&(struct gkyl_mirror_geo_gen_inp) {
        .comp_grid = &comp_grid,
        .mirror_grid = mirror_grid,
      }
    );

  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  // Create nodal range
  struct gkyl_range nodal_range;
  double dzc[3] = {0.0};
  
  int nodes[3];
  if (poly_order == 1) {
    for (int d=0; d<3; ++d)
      nodes[d] = gkyl_range_shape(&range, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<3; ++d)
      nodes[d] = 2*gkyl_range_shape(&range, d) + 1;
  }
  gkyl_range_init_from_shape(&nodal_range, 3, nodes);

  // Comp grid is 11 x 33. It's the local nodes
  // ext grid is 12 x 18 x 34

  struct gkyl_mirror_geo_dg *mirror_geo_dg = 
    gkyl_mirror_geo_dg_inew(&(struct gkyl_mirror_geo_dg_inp) {
        .comp_grid = &comp_grid,
        .mirror_geo = mirror_geo,
        .local = range,
        .local_ext = ext_range,
        .nodal_range = nodal_range,
        .basis = basis,
      }
    );

  // struct gkyl_range node_range;
  // gkyl_range_init_from_shape(&node_range, 2, (int[2]) { cells[0]+1, cells[2]+1 });

  // struct gkyl_range_iter iter;
  // gkyl_range_iter_init(&iter, &node_range);
  // while (gkyl_range_iter_next(&iter)) {
    // long loc = gkyl_range_idx(&node_range, iter.idx);

    // const struct gkyl_mirror_geo_gen_geom *g =
    //   gkyl_array_cfetch(mirror_geo->nodes_geom, loc);
    
    // const struct gkyl_mirror_grid_gen_geom *grid = 
    //   gkyl_array_cfetch(mirror_grid->nodes_geom, loc);

    // const double *rz = g->rz_coord;
    // // check R,Z coordinates
    // const double *rz_grid = gkyl_array_cfetch(mirror_grid->nodes_rz, loc);
    // TEST_CHECK( gkyl_compare_double(rz[0], rz_grid[0], 1e-14) );
    // TEST_CHECK( gkyl_compare_double(rz[1], rz_grid[1], 1e-14) );

    // // check tangent vectors
    // for (int i=0; i<3; ++i) {
    //   struct gkyl_vec3 tang_cart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, grid->tang[i]);
    //   TEST_CHECK( gkyl_compare_double(tang_cart.x[0], g->tang[i].x[0], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(tang_cart.x[1], g->tang[i].x[1], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(tang_cart.x[2], g->tang[i].x[2], 1e-14) );
    // }

    // // check dual vectors
    // for (int i=0; i<3; ++i) {
    //   struct gkyl_vec3 dual_cart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, grid->dual[i]);
    //   TEST_CHECK( gkyl_compare_double(dual_cart.x[0], g->dual[i].x[0], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(dual_cart.x[1], g->dual[i].x[1], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(dual_cart.x[2], g->dual[i].x[2], 1e-14) );
    // }

    // // check normal vectors
    // for (int i=0; i<3; ++i) {
    //   TEST_CHECK( gkyl_compare_double(g->normal[i].x[0], g->dual[i].x[0] / g->dualmag[i], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->normal[i].x[1], g->dual[i].x[1] / g->dualmag[i], 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->normal[i].x[2], g->dual[i].x[2] / g->dualmag[i], 1e-14) );
    // }

    // // check dual vector magnitudes
    // for (int i=0; i<3; ++i) {
    //   TEST_CHECK( gkyl_compare_double(g->dualmag[i], gkyl_vec3_len(g->dual[i]), 1e-14) );
    // }

    // // Check metric tensors
    // int count = 0;
    // for (int i=0; i<3; ++i) {
    //   for (int j=0; j<3; ++j) {
    //     if (i > j)
    //       continue;
    //     double g_ij = gkyl_vec3_dot(g->tang[i], g->tang[j]);
    //     TEST_CHECK( gkyl_compare_double(g->metric_covar[count], g_ij, 1e-14) );
    //     double gij = gkyl_vec3_dot(g->dual[i], g->dual[j]);
    //     TEST_CHECK( gkyl_compare_double(g->metric_contr[count], gij, 1e-14) );
    //     TEST_CHECK( gkyl_compare_double(g->metric_covar_neut[count], g_ij, 1e-14) );
    //     TEST_CHECK( gkyl_compare_double(g->metric_contr_neut[count], gij, 1e-14) );
    //     count++;
    //   }
    // }

    // // check magnetic field
    // double Bmag = gkyl_vec3_len( grid->B);
    // TEST_CHECK( gkyl_compare_double(g->Bmag, Bmag, 1e-14));
    // TEST_CHECK( gkyl_compare_double(g->Bmag_inv, 1.0/Bmag, 1e-14) );
    // TEST_CHECK( gkyl_compare_double(g->Bmag_inv_sq, 1.0/(Bmag*Bmag), 1e-14) );

    // // check B vector
    // TEST_CHECK( gkyl_compare_double(g->B_covar.x[0], grid->B.x[0], 1e-14) );
    // TEST_CHECK( gkyl_compare_double(g->B_covar.x[1], grid->B.x[1] * rz[0]*rz[0], 1e-14) ); // B_1 = B^1 * r^2
    // TEST_CHECK( gkyl_compare_double(g->B_covar.x[2], grid->B.x[2], 1e-14) );

    // struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, grid->B);
    // TEST_CHECK( gkyl_compare_double(g->B_cart.x[0], B_cart.x[0], 1e-14) );
    // TEST_CHECK( gkyl_compare_double(g->B_cart.x[1], B_cart.x[1], 1e-14) );
    // TEST_CHECK( gkyl_compare_double(g->B_cart.x[2], B_cart.x[2], 1e-14) );


    // // check Jacobian
    // double Jac = gkyl_vec3_triple(g->tang[0], g->tang[1], g->tang[2]);

    // if (g->rz_coord[0] > 0){
    //   TEST_CHECK( gkyl_compare_double(Jac, g->Jc, 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->Jc_inv, 1.0/Jac, 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->JB, g->Jc*Bmag, 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->JB_inv, 1.0/(g->Jc*Bmag), 1e-14) );
    // }

    // // check C = Jc*Bmag/sqrt(g33)
    // double g33 = gkyl_vec3_dot( g->tang[2], g->tang[2]);

    // if (g->rz_coord[0] > 0) {
    //   TEST_CHECK( gkyl_compare_double(g->C, g->JB/sqrt(g33), 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->C, g->Jc*Bmag/sqrt(g33), 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->C, g->Jc*g->Bmag/sqrt(g33), 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->C, g->Jc*g->Bmag/sqrt(g->metric_covar[5]), 1e-14) );
    //   TEST_CHECK( gkyl_compare_double(g->eps2, g->Jc*g->metric_contr[5] - g->JB/g->metric_covar[5], 1e-14) );
    // }
  // }

  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_mirror_geo_gen_release(mirror_geo);
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
  
  { "wham_no_axis_sqrt_psi", test_wham_no_axis_sqrt_psi },
  { "wham_with_axis_sqrt_psi", test_wham_with_axis_sqrt_psi },
  
  { NULL, NULL },
};
