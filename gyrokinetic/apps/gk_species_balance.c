#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// Only the particle (M0) balance is available: M0 is always component 0. Balances of
// other moments (e.g. toroidal momentum or energy) need more machinery to reuse
// the existing diagnostics without redundancy, or their own diagnostics here.
static bool
gk_species_balance_has_m0(enum gkyl_distribution_moments mom_type)
{
  return mom_type == GKYL_F_MOMENT_M0 || mom_type == GKYL_F_MOMENT_M0M1M2 ||
         mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP || mom_type == GKYL_F_MOMENT_HAMILTONIAN;
}

// Whether this rank reports the integrated bflux through the boundary at
// (dir, edge). The bflux diagnostic sums over the ranks of a boundary plane
// (comm_plane), so only one rank per plane must contribute to the total.
static bool
gk_species_balance_owns_boundary(gkyl_gyrokinetic_app *app, int dir, enum gkyl_edge_loc edge)
{
  bool at_boundary = edge == GKYL_LOWER_EDGE ? app->local.lower[dir] == app->global.lower[dir] :
                                               app->local.upper[dir] == app->global.upper[dir];
  for (int d = 0; d < app->cdim; d++) {
    if (d != dir) {
      at_boundary = at_boundary && app->local.lower[d] == app->global.lower[d];
    }
  }
  return at_boundary;
}

// Particle balance of balance moment m: dN/dt = source - boundary fluxes.
// Combines the last entries of the integrated moment diagnostics of f, df/dt,
// the source and the boundary fluxes, all computed before this is called.
static void
gk_species_balance_calc_particle(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal, int m,
  double tm
)
{
  double row[bal->num_comp];
  for (int k = 0; k < bal->num_comp; k++) {
    row[k] = 0.0;
  }

  if (gks->fdot_dt > 0.0) {
    // Rates are only defined once a step has been taken.
    double f_mom[gks->integ_moms.num_mom];
    gkyl_dynvec_getlast(gks->integ_diag, f_mom);
    double N = f_mom[bal->f_comp[m]];
    gkyl_dynvec_getlast(gks->fdot_integ_diag, f_mom);
    double dNdt = f_mom[bal->f_comp[m]];

    double src = 0.0;
    if (bal->has_source) {
      double src_mom[gks->src.integ_moms.num_mom];
      gkyl_dynvec_getlast(gks->src.integ_diag, src_mom);
      src = src_mom[bal->src_comp[m]];
    }

    // Total flux through the boundaries.
    struct gk_boundary_fluxes *bflux = &gks->bflux;
    int num_int_mom = gks->info.boundary_flux_diagnostics.num_integrated_diag_moments;
    double bflux_tot_local = 0.0, bflux_tot = 0.0;
    for (int b = 0; b < bal->num_boundaries; ++b) {
      if (gk_species_balance_owns_boundary(
            app, bflux->boundaries_dir[b], bflux->boundaries_edge[b]
          )) {
        int mom_idx = bflux->diag_int_mom_idx[bal->bflux_int_mom[m]];
        double bflux_mom[bflux->moms_op[mom_idx].num_mom];
        gkyl_dynvec_getlast(bflux->intmom[b * num_int_mom + bal->bflux_int_mom[m]], bflux_mom);
        bflux_tot_local += bflux_mom[bal->bflux_comp[m]];
      }
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_tot_local, &bflux_tot);

    // Particles gained/lost in the last step, also relative to N (round-off ~1e-16 if conservative).
    double mom_err = src - bflux_tot - dNdt;
    row[0] = dNdt;
    row[1] = src;
    row[2] = bflux_tot;
    row[3] = mom_err;
    row[4] = mom_err * gks->fdot_dt / N;
  }

  gkyl_dynvec_append(bal->diag[m], tm, row);
}

void
gk_species_balance_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  enum gkyl_species_bflux_type *bflux_type
)
{
  struct gkyl_gyrokinetic_species *info = &gks->info;
  // Check that the diagnostics the balance is computed from are requested, and
  // request them otherwise.
  *bal = (struct gk_species_balance){};
  bal->enabled = info->num_balance_moments > 0 && !info->is_static;
  if (!bal->enabled) {
    return;
  }
  bal->num_moms = info->num_balance_moments;
  for (int m = 0; m < bal->num_moms; m++) {
    bal->mom_types[m] = info->balance_moments[m];
    assert(bal->mom_types[m] == GKYL_F_MOMENT_M0); // Only the particle balance is available.
    bal->f_comp[m] = 0; // M0 is always component 0.
    bal->src_comp[m] = 0; // M0 is always component 0.
  }

  // Integrated moments of f and df/dt (the default set contains M0).
  assert(
    info->num_integrated_diag_moments == 0 ||
    gk_species_balance_has_m0(info->integrated_diag_moments[0])
  );
  info->time_rate_diagnostics = true;

  // Integrated moments of the source.
  bal->has_source = info->source.source_id != 0;
  if (bal->has_source) {
    struct gkyl_phase_diagnostics_inp *src_diag = &info->source.diagnostics;
    assert(!src_diag->time_integrated); // The balance needs rates.
    if (src_diag->num_integrated_diag_moments == 0) {
      src_diag->num_integrated_diag_moments = 1;
      src_diag->integrated_diag_moments[0] = GKYL_F_MOMENT_M0M1M2PARM2PERP; // The default set.
    }
    assert(gk_species_balance_has_m0(src_diag->integrated_diag_moments[0]));
  }

  // Integrated moments of the boundary fluxes (same boundaries as gk_species_bflux_init).
  bal->num_boundaries = 0;
  for (int d = 0; d < app->cdim; ++d) {
    if (gks->bc_is_np[d]) {
      if (gks->lower_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        bal->num_boundaries++;
      }
      if (gks->upper_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        bal->num_boundaries++;
      }
    }
  }
  if (bal->num_boundaries > 0) {
    struct gkyl_phase_diagnostics_inp *bflux_diag = &info->boundary_flux_diagnostics;
    assert(!bflux_diag->time_integrated); // The balance needs rates.
    for (int m = 0; m < bal->num_moms; m++) {
      bal->bflux_int_mom[m] = -1;
      for (int i = 0; i < bflux_diag->num_integrated_diag_moments; i++) {
        if (gk_species_balance_has_m0(bflux_diag->integrated_diag_moments[i])) {
          bal->bflux_int_mom[m] = i;
          bal->bflux_comp[m] = 0;
          break;
        }
      }
      if (bal->bflux_int_mom[m] < 0) {
        // Request the integrated moment of the bflux.
        assert(bflux_diag->num_integrated_diag_moments < 12);
        bal->bflux_int_mom[m] = bflux_diag->num_integrated_diag_moments;
        bal->bflux_comp[m] = 0;
        bflux_diag->integrated_diag_moments[bflux_diag->num_integrated_diag_moments++] =
          bal->mom_types[m];
      }
    }
    *bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  }

  bal->num_comp = 5; // mom_dot, src, bflux_tot, mom_err, mom_err_norm=mom_err*dt/N.
  for (int m = 0; m < bal->num_moms; m++) {
    bal->diag[m] = gkyl_dynvec_new(GKYL_DOUBLE, bal->num_comp);
    bal->is_first_write[m] = true;
    // Momentum and energy balances will need their own recipe.
    bal->calc_func[m] = gk_species_balance_calc_particle;
  }
}

void
gk_species_balance_calc(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal, double tm
)
{
  for (int m = 0; m < bal->num_moms; m++) {
    bal->calc_func[m](app, gks, bal, m, tm);
  }
}

void
gk_species_balance_write(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal
)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int m = 0; m < bal->num_moms; m++) {
    if (rank == 0) {
      // File named after the moment, e.g. <name>-<species>_balance_M0.gkyl.
      const char *mom_name = gkyl_distribution_moments_strs[bal->mom_types[m]];
      const char *fmt = "%s-%s_balance_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, mom_name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, mom_name);

      if (bal->is_first_write[m]) {
        struct gkyl_msgpack_map_elem io_meta_bal[] = {
          {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Species moment balance."},
          {
            .key = "Columns",
            .elem_type = GKYL_MP_STRING,
            .cval = "mom_dot, src_mom, bflux_tot_mom, mom_err, mom_err_norm",
          }
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 2};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_bal
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(bal->diag[m], fileNm, mt);
        gkyl_msgpack_data_release(mt);
        bal->is_first_write[m] = false;
      } else {
        gkyl_dynvec_awrite(bal->diag[m], fileNm);
      }
    }
    gkyl_dynvec_clear(bal->diag[m]);
    app->stat.n_diag_io += 1;
  }
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_balance_release(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_balance *bal
)
{
  if (bal->enabled) {
    for (int m = 0; m < bal->num_moms; m++) {
      gkyl_dynvec_release(bal->diag[m]);
    }
  }
}
