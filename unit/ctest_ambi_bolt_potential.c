// Description: Test for the boltzmann electron model for potential

#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_ambi_bolt_potential_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

void eval_one(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void eval_hat(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 2. - fabs(xn[0]);
}

void
test_ambi_bolt_init_1x()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_ambi_bolt_potential *ambi = gkyl_ambi_bolt_potential_new(&grid, basis, mass_e, charge_e, temp_e, use_gpu);

  TEST_CHECK( ambi->cdim == 1);
  TEST_CHECK( ambi->num_basis == basis->num_basis);
  TEST_CHECK( ambi->use_gpu == use_gpu);
  TEST_CHECK( gkyl_compare_double(ambi->dz, 2./8., 1e-12));
  TEST_CHECK( gkyl_compare_double(ambi->mass_e, mass_e, 1e-12));
  TEST_CHECK( gkyl_compare_double(ambi->charge_e, charge_e, 1e-12));
  TEST_CHECK( gkyl_compare_double(ambi->temp_e, temp_e, 1e-12));

  gkyl_ambi_bolt_potential_release(ambi);
}

void
test_ambi_bolt_sheath_calc_1x()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_ambi_bolt_potential *ambi = gkyl_ambi_bolt_potential_new(&grid, basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2*cdim];
  for (int j=0; j<cdim; ++j) {
    sheath_vals[2*j]   = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    sheath_vals[2*j+1] = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2*j],   0.0);
    gkyl_array_clear(sheath_vals[2*j+1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost); 
    gkyl_skin_ghost_ranges(&upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost);
  }

  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *gamma_i = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_one, NULL);  
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, jacobgeo_inv);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, gamma_i);

  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_LOWER_EDGE, 
    &lower_skin[0], &lower_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[0]);
  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_UPPER_EDGE,
    &upper_skin[0], &upper_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[1]);

  // Serendipity 1x basis is [1/sqrt(2), sqrt(3/2)x].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part. 
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) gamma_i J dz/2))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double *sheath_lower_c = ((double *) gkyl_array_cfetch(sheath_vals[0], 0));
  double *sheath_upper_c = ((double *) gkyl_array_cfetch(sheath_vals[1], 9));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[3], 0, 1e-12));

  // This operation happens after the sheaths are determined in the app, so we should test this
  gkyl_array_accumulate(sheath_vals[0], 1., sheath_vals[1]);

  double *sheath_lower_c_avg = ((double *) gkyl_array_cfetch(sheath_vals[0], 0));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[3], 0, 1e-12));
  double *sheath_upper_c_avg = ((double *) gkyl_array_cfetch(sheath_vals[0], 9));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[3], 0, 1e-12));

  gkyl_ambi_bolt_potential_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_ambi_bolt_phi_calc_1x()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_ambi_bolt_potential *ambi = gkyl_ambi_bolt_potential_new(&grid, basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2*cdim];
  for (int j=0; j<cdim; ++j) {
    sheath_vals[2*j]   = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    sheath_vals[2*j+1] = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2*j],   0.0);
    gkyl_array_clear(sheath_vals[2*j+1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost); 
    gkyl_skin_ghost_ranges(&upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost);
  }

  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *gamma_i = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_one, NULL);  
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, jacobgeo_inv);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, gamma_i);

  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_LOWER_EDGE, 
    &lower_skin[0], &lower_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[0]);
  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_UPPER_EDGE,
    &upper_skin[0], &upper_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[1]);

  gkyl_array_accumulate(sheath_vals[0], 1., sheath_vals[1]);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_ambi_bolt_potential_phi_calc(ambi, &local, &local_ext,
    jacobgeo_inv, M0, sheath_vals[0], phi);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    double *phi_c = ((double *) gkyl_array_cfetch(phi, iter.idx[0]));
    // phi should be the same value as the sheath potential
    TEST_CHECK(gkyl_compare_double(phi_c[0], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
    TEST_CHECK(gkyl_compare_double(phi_c[1], 0.0, 1e-12));
  }

  gkyl_ambi_bolt_potential_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_ambi_bolt_sheath_calc_1x_hat()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_ambi_bolt_potential *ambi = gkyl_ambi_bolt_potential_new(&grid, basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2*cdim];
  for (int j=0; j<cdim; ++j) {
    sheath_vals[2*j]   = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    sheath_vals[2*j+1] = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2*j],   0.0);
    gkyl_array_clear(sheath_vals[2*j+1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost); 
    gkyl_skin_ghost_ranges(&upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost);
  }

  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *gamma_i = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_one, NULL); 
  gkyl_proj_on_basis *proj_hat = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_hat, NULL); 

  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, jacobgeo_inv);
  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, gamma_i);

  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_LOWER_EDGE, 
    &lower_skin[0], &lower_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[0]);
  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_UPPER_EDGE,
    &upper_skin[0], &upper_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[1]);

  // Serendipity 1x basis is [1/sqrt(2), sqrt(3/2)x].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part. 
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) gamma_i J dz/2))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double *sheath_lower_c = ((double *) gkyl_array_cfetch(sheath_vals[0], 0));
  double *sheath_upper_c = ((double *) gkyl_array_cfetch(sheath_vals[1], 9));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[3], 0, 1e-12));

  // This operation happens after the sheaths are determined in the app, so we should test this
  gkyl_array_accumulate(sheath_vals[0], 1., sheath_vals[1]);

  double *sheath_lower_c_avg = ((double *) gkyl_array_cfetch(sheath_vals[0], 0));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[3], 0, 1e-12));
  double *sheath_upper_c_avg = ((double *) gkyl_array_cfetch(sheath_vals[0], 9));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[2], log(1/(sqrt(2*M_PI)*0.5*1/4))*sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c_avg[3], 0, 1e-12));

  gkyl_ambi_bolt_potential_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

void
test_ambi_bolt_phi_calc_1x_hat()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower)/sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_ambi_bolt_potential *ambi = gkyl_ambi_bolt_potential_new(&grid, basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2*cdim];
  for (int j=0; j<cdim; ++j) {
    sheath_vals[2*j]   = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    sheath_vals[2*j+1] = gkyl_array_new(GKYL_DOUBLE, 2*basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2*j],   0.0);
    gkyl_array_clear(sheath_vals[2*j+1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost); 
    gkyl_skin_ghost_ranges(&upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost);
  }

  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);
  struct gkyl_array *gamma_i = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_one, NULL); 
  gkyl_proj_on_basis *proj_hat = gkyl_proj_on_basis_new(&grid, basis, poly_order+1, 1, eval_hat, NULL); 

  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, jacobgeo_inv);
  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, gamma_i);

  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_LOWER_EDGE, 
    &lower_skin[0], &lower_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[0]);
  gkyl_ambi_bolt_potential_sheath_calc(ambi, GKYL_UPPER_EDGE,
    &upper_skin[0], &upper_ghost[0], jacobgeo_inv, gamma_i, M0, sheath_vals[1]);

  gkyl_array_accumulate(sheath_vals[0], 1., sheath_vals[1]);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_ambi_bolt_potential_phi_calc(ambi, &local, &local_ext,
    jacobgeo_inv, M0, sheath_vals[0], phi);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  double phi_sheath = log(1/(sqrt(2*M_PI)*0.5*1/4));
  while (gkyl_range_iter_next(&iter)) {
    double *phi_c = ((double *) gkyl_array_cfetch(phi, iter.idx[0]));
    double *ni_c  = ((double *) gkyl_array_cfetch(M0, iter.idx[0]));
    double ni = ni_c[0]/sqrt(2);
    TEST_CHECK(gkyl_compare_double(phi_c[0]/sqrt(2), phi_sheath + log(ni), 1e-3));
  }

  gkyl_ambi_bolt_potential_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

TEST_LIST = {
  { "test_ambi_bolt_init_1x", test_ambi_bolt_init_1x },
  { "test_ambi_bolt_sheath_calc_1x", test_ambi_bolt_sheath_calc_1x },
  { "test_ambi_bolt_phi_calc_1x", test_ambi_bolt_phi_calc_1x },
  { "test_ambi_bolt_sheath_calc_1x_hat", test_ambi_bolt_sheath_calc_1x_hat },
  { "test_ambi_bolt_phi_calc_1x_hat", test_ambi_bolt_phi_calc_1x_hat },
  { NULL, NULL },
};