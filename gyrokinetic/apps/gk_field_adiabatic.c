#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <time.h>

// Adiabatic species response in the field equation. Methods valid for any adiabatic species
// are named gk_field_adiabatic_*, those specific to the adiabatic electron model
// gk_field_adiab_elc_*.
//
// Adiabatic electron response with a flux-surface average (FSA) in 2x/3x:
//   -div(eps grad_perp phi) + K (phi - <phi>) = rho,  K = (e^2 n0/Te) J,
// with <phi> = int J phi dy dz / int J dy dz. It is solved by the fem_poisson_perp_adiabatic
// updater (Woodbury identity), and the FSA and its extension are computed by flux_surf_avg.

void
gk_field_adiabatic_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->gk_adiab.rho_bg = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->gk_adiab.kJ = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
}

void
gk_field_adiabatic_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_release(gkf->gk_adiab.rho_bg);
  gkyl_array_release(gkf->gk_adiab.kJ);
}

// Smooth a right-hand side along z like the charge density, before each perpendicular solve.
static void
gk_field_adiabatic_rhs_proj(void *ctx, struct gkyl_array *rhs)
{
  gkyl_gyrokinetic_app *app = ctx;
  app->field->fem_projection_par_rho_func(app, app->field, rhs, rhs);
}

void
gk_field_adiab_elc_rhs_phi_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve -div(eps grad_perp phi) + K (phi - <phi>) = rho.
  gkyl_fem_poisson_perp_adiabatic_advance(field->gk_adiab.solver, field->rho_c, field->phi_smooth);

  // Smooth the potential along z.
  field->fem_projection_par_phi_func(app, field, field->phi_smooth, field->phi_smooth);

  field->invert_flr(app, field, field->phi_smooth);
}

void
gk_field_adiab_elc_dphi(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->gk_adiab;
  gkyl_flux_surf_avg_advance(ad->fsa, phi, ad->psi);
  gkyl_flux_surf_avg_inflate(ad->fsa, ad->psi, out);
  gkyl_array_scale_range(out, -1.0, &app->local);
  gkyl_array_accumulate_range(out, 1.0, phi, &app->local);
}

void
gk_field_adiab_elc_energy_accumulate(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, double factor, double *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  gk_field_adiab_elc_dphi(app, gkf, gkf->phi_smooth, ad->phi_buff);
  gkyl_array_integrate_advance(
    ad->calc_energy, ad->phi_buff, 0.5 * factor, ad->kJ, &app->local, &app->local, ad->energy_red
  );

  if (app->use_gpu) {
    double en[2];
    gkyl_cu_memcpy(&en[0], out, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(&en[1], ad->energy_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    en[0] += en[1];
    gkyl_cu_memcpy(out, &en[0], sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    out[0] += ad->energy_red[0];
  }
}

void
gk_field_adiabatic_profiles_calc(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;
  struct gk_species *gks = &app->species[ad->species_idx];
  double q_s = gks->info.charge;

  // Maxwellian moments n_s0*J, upar and T_s/m_s.
  struct gk_species_moment *moms = &gks->lte.moms;
  gk_species_moment_calc(moms, gks->local, app->local, gks->f);

  gkyl_array_set_offset(ad->rho_bg, q_s, moms->marr, 0);
  gkyl_dg_div_op_range(
    moms->mem_geo, &app->basis, 0, ad->kJ, 0, moms->marr, 2, moms->marr, &app->local
  );
  gkyl_array_scale(ad->kJ, q_s * q_s / gks->info.mass);
}

void
gk_field_adiab_elc_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  // 1D (x) basis and ranges of the flux-surface averages.
  gkyl_cart_modal_serendip(&ad->basis_x, 1, app->poly_order);
  gkyl_range_init(&ad->local_x_ext, 1, &app->local_ext.lower[0], &app->local_ext.upper[0]);
  gkyl_sub_range_init(&ad->local_x, &ad->local_x_ext, &app->local.lower[0], &app->local.upper[0]);
  ad->psi = mkarr(app->use_gpu, ad->basis_x.num_basis, ad->local_x_ext.volume);
  ad->phi_buff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Flux-surface average and its extension, reduced over the ranks decomposing z.
  ad->fsa = gkyl_flux_surf_avg_inew(&(struct gkyl_flux_surf_avg_inp){
    .grid = &app->grid,
    .basis = app->basis,
    .basis_x = ad->basis_x,
    .local = &app->local,
    .local_ext = &app->local_ext,
    .local_x = &ad->local_x,
    .local_x_ext = &ad->local_x_ext,
    .jacobgeo = app->gk_geom->geo_int.jacobgeo,
    .comm = app->comm,
    .use_gpu = app->use_gpu,
  });

  // Energy diagnostic, (1/2) int K (phi-<phi>)^2.
  ad->calc_energy = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu
  );
  ad->energy_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  // Field solver with the adiabatic response, smoothing every right-hand side along z.
  struct timespec wst = gkyl_wall_clock();
  ad->solver = gkyl_fem_poisson_perp_adiabatic_inew(&(struct gkyl_fem_poisson_perp_adiabatic_inp){
    .grid = &app->grid,
    .basis = app->basis,
    .local = &app->local,
    .local_ext = &app->local_ext,
    .basis_x = ad->basis_x,
    .local_x = &ad->local_x,
    .local_x_ext = &ad->local_x_ext,
    .fsa = ad->fsa,
    .bcs = &gkf->poisson_bcs,
    .bias_line_list = gkf->info.bias_line_list,
    .epsilon = gkf->epsilon,
    .kJ = ad->kJ,
    .rhs_proj = gk_field_adiabatic_rhs_proj,
    .rhs_proj_ctx = app,
    .use_gpu = app->use_gpu,
  });

  gkyl_gyrokinetic_app_cout(
    app, stdout, "Adiabatic species '%s' with FSA: m = %ld, cond_1(I-G) = %.4e, setup %.3e sec.\n",
    app->species[ad->species_idx].info.name, ad->local_x.volume * ad->basis_x.num_basis,
    gkyl_fem_poisson_perp_adiabatic_cond(ad->solver), gkyl_time_diff_now_sec(wst)
  );
}

void
gk_field_adiab_elc_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  gkyl_fem_poisson_perp_adiabatic_release(ad->solver);
  gkyl_flux_surf_avg_release(ad->fsa);
  gkyl_array_release(ad->psi);
  gkyl_array_release(ad->phi_buff);
  gkyl_array_integrate_release(ad->calc_energy);
  if (app->use_gpu) {
    gkyl_cu_free(ad->energy_red);
  } else {
    gkyl_free(ad->energy_red);
  }
}
