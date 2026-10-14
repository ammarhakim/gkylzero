// Exercise filter allocation, runtime resets, modal conversion and diagnostic output.
#include <assert.h>
#include <math.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_array_rio_priv.h>
#include <mpack.h>

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int dir = 0; dir < 3; ++dir) {
    xp[dir] = xc[dir];
  }
}

static void
bfield(double t, const double *xc, double *out, void *ctx)
{
  out[0] = out[1] = 0.0;
  out[2] = 1.0;
}

static void
initial_distribution(double t, const double *xc, double *out, void *ctx)
{
  out[0] = exp(-0.5 * xc[1] * xc[1] - xc[2]);
}

static void
check_filter_state(const struct gk_species *gks)
{
  const struct gkyl_array *arrays[] = {gks->damping.fbar, gks->damping.fbar1, gks->damping.fbarnew};
  for (int stage = 0; stage < 3; ++stage) {
    for (size_t cell = 0; cell < arrays[stage]->size; ++cell) {
      const double *coeff = gkyl_array_cfetch(arrays[stage], cell);
      assert(coeff[0] == 2.0);
      for (size_t mode = 1; mode < arrays[stage]->ncomp; ++mode) {
        assert(coeff[mode] == 0.0);
      }
    }
  }
}

static void
check_metadata(
  const char *prefix, const char *suffix, const char *form, int poly_order, const char *basis
)
{
  char filename[256];
  snprintf(filename, sizeof filename, "%s-ion_%s.gkyl", prefix, suffix);
  FILE *fp = fopen(filename, "rb");
  assert(fp);
  struct gkyl_array_header_info hdr;
  int status = gkyl_header_meta_read_fp(&hdr, fp);
  assert(status == GKYL_ARRAY_RIO_SUCCESS);
  fclose(fp);
  mpack_tree_t tree;
  mpack_tree_init_data(&tree, hdr.meta, hdr.meta_size);
  mpack_tree_parse(&tree);
  mpack_node_t root = mpack_tree_root(&tree);
  char value[128];
  mpack_node_copy_cstr(mpack_node_map_cstr(root, "value_form"), value, sizeof value);
  assert(strcmp(value, form) == 0);
  assert(mpack_node_strlen(mpack_node_map_cstr(root, "Description")) > 0);
  if (basis) {
    assert(mpack_node_u64(mpack_node_map_cstr(root, "poly_order")) == poly_order);
    mpack_node_copy_cstr(mpack_node_map_cstr(root, "basis_type"), value, sizeof value);
    assert(strcmp(value, basis) == 0);
  }
  mpack_error_t error = mpack_tree_destroy(&tree);
  assert(error == mpack_ok);
  gkyl_array_header_info_release(&hdr);
}

int
main(int argc, char **argv)
{
  // Test both initial activation and activation later through the reset API.
  for (int initially_enabled = 0; initially_enabled < 2; ++initially_enabled) {
    struct gkyl_gk inp = {
      .cdim = 1,
      .lower = {-1.0},
      .upper = {1.0},
      .cells = {2},
      .poly_order = 1,
      .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
      .geometry = {.geometry_id = GKYL_GEOMETRY_MAPC2P, .mapc2p = mapc2p, .bfield_func = bfield},
      .num_periodic_dir = 1,
      .periodic_dirs = {0},
      .num_species = 1,
      .species = {{
        .name = "ion",
        .charge = 1.0,
        .mass = 1.0,
        .vdim = 2,
        .lower = {-4.0, 0.0},
        .upper = {4.0, 4.0},
        .cells = {4, 2},
        .polarization_density = 1.0,
        .projection = {.proj_id = GKYL_PROJ_FUNC, .func = initial_distribution},
        .time_rate_diagnostics = true,
      }},
      .field =
        {
          .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
          .electron_mass = 1.0,
          .electron_charge = -1.0,
          .electron_temp = 1.0,
          .zero_init_field = true,
          .is_static = true,
        },
    };
    snprintf(inp.name, sizeof inp.name, "gk_low_pass_filter_%d", initially_enabled);
    struct gkyl_gyrokinetic_damping damping = {
      .type = GKYL_GK_DAMPING_LOW_PASS_FILTER,
      .rate_const = 1.0,
      .write_rate = true,
      .write_fbar = true,
      .cellwise_const = true,
      .do_not_reset_fbar = true,
    };
    if (initially_enabled) {
      inp.species[0].damping = damping;
    }
    gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp);
    gkyl_gyrokinetic_app_apply_ic(app, 0.0);
    struct gk_species *gks = &app->species[0];
    if (!initially_enabled) {
      assert(!gks->damping.fbar && !gks->damping.fbar1 && !gks->damping.fbarnew);
      assert(!gks->fbardot_mom_old && !gks->fbardot_mom_new && !gks->fbardot_integ_diag);
      gkyl_gyrokinetic_app_calc_integrated_mom(app, 0.0);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_reset_species_damping(app, 0.0, "ion", damping);
    }
    assert(gks->fbardot_mom_old && gks->fbardot_mom_new && gks->fbardot_integ_diag);
    assert(gks->damping.fbar_mom_buffer);
    // Exercise the p0 moment embedding and RK evolution before testing resets.
    struct gkyl_update_status p0_status = gkyl_gyrokinetic_update(app, 1e-3);
    assert(p0_status.success);
    gkyl_array_clear(gks->damping.fbar, 2.0);
    gk_species_damping_write(app, gks, 0.0, 0);

    // p0 -> full basis: preserve the constant mode and zero the higher modes.
    damping.cellwise_const = false;
    gkyl_gyrokinetic_app_reset_species_damping(app, 0.0, "ion", damping);
    assert(gks->damping.fbar->ncomp == gks->basis.num_basis);
    assert(!gks->damping.fbar_mom_buffer);
    check_filter_state(gks);
    gk_species_damping_write(app, gks, 0.0, 1);

    // Full basis -> p0: discard higher modes, preserving the constant mode.
    double *coeff = gkyl_array_fetch(gks->damping.fbar, 0);
    coeff[1] = 3.0;
    damping.cellwise_const = true;
    gkyl_gyrokinetic_app_reset_species_damping(app, 0.0, "ion", damping);
    check_filter_state(gks);

    // Disabled state retains only fbar when explicitly requested.
    struct gkyl_array *mom_old = gks->fbardot_mom_old;
    damping.type = GKYL_GK_DAMPING_NONE;
    gkyl_gyrokinetic_app_reset_species_damping(app, 0.0, "ion", damping);
    assert(gks->damping.fbar && !gks->damping.fbar1 && !gks->damping.fbarnew);
    assert(!gks->damping.fbar_mom_buffer);
    size_t nsamples = gkyl_dynvec_size(gks->fbardot_integ_diag);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, 0.0);
    assert(gkyl_dynvec_size(gks->fbardot_integ_diag) == nsamples);

    damping.type = GKYL_GK_DAMPING_LOW_PASS_FILTER;
    damping.cellwise_const = false;
    gkyl_gyrokinetic_app_reset_species_damping(app, 0.0, "ion", damping);
    check_filter_state(gks);
    assert(gks->fbardot_mom_old == mom_old);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, 1e-3);
    assert(status.success);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, status.dt_actual);
    assert(gkyl_dynvec_size(gks->fbardot_integ_diag) == nsamples + 1);
    gkyl_gyrokinetic_app_write_integrated_mom(app);

    check_metadata(inp.name, "damping_fbar_0", "modal", 0, "serendipity");
    check_metadata(inp.name, "damping_fbar_1", "modal", 1, "gkhybrid");
    check_metadata(inp.name, "damping_rate_0", "nodal", 0, "serendipity");
    check_metadata(inp.name, "fdot_integrated_moms", "nodal", 0, 0);
    check_metadata(inp.name, "fbardot_integrated_moms", "nodal", 0, 0);

    // Disabling without preservation releases the remembered distribution too.
    damping.type = GKYL_GK_DAMPING_NONE;
    damping.do_not_reset_fbar = false;
    gkyl_gyrokinetic_app_reset_species_damping(app, status.dt_actual, "ion", damping);
    assert(!gks->damping.fbar && !gks->damping.fbar_initialized);
    gkyl_gyrokinetic_app_release(app);
  }
  return 0;
}
