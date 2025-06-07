#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_species_damping_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_damping_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = "tensor",
    }
  );

  // Write out the damping rate.
  const char *fmt = "%s-%s_damping_rate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.rate_host, gks->damping.rate);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->damping.rate_host, fileNm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_damping_write_enabled(app, gks, tm, frame);
  gks->damping.write_func = gk_species_damping_write_disabled;
}

void
gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_damping *damp)
{
  damp->type = gks->info.damping.type;
  damp->evolve = false; // Whether the rate is time dependent.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;

  if (damp->type) {
    // Allocate rate array.
    damp->rate = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    damp->rate_host = damp->rate;
    if (app->use_gpu) {
      damp->rate_host = mkarr(false, damp->rate->ncomp, damp->rate->size); 
    }

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_new(&gks->grid, &gks->basis, 1, 1, 
        gks->info.damping.rate_profile, gks->info.damping.rate_profile_ctx);
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, damp->rate_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(damp->rate, damp->rate_host);
    }

    // Set function pointers chosen at runtime.
    if (damp->evolve) {
      damp->write_func = gk_species_damping_write_enabled;
    }
    else {
      damp->write_func = gk_species_damping_write_init_only;
    }
  }
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp, 
  const struct gkyl_array *fin, struct gkyl_array *f_buffer, struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();
  if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
    gkyl_array_set(f_buffer, 1.0, fin);
    gkyl_array_scale_by_cell(f_buffer, damp->rate);
    gkyl_array_accumulate(rhs, -1.0, f_buffer);

    // Add the frequency to the CFL frequency.
    gkyl_array_accumulate(cflrate, 1.0, damp->rate);
  }
  app->stat.species_damp_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->damping.write_func(app, gks, tm, frame);
}

void
gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp)
{
  if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
    gkyl_array_release(damp->rate);
    if (app->use_gpu) {
      gkyl_array_release(damp->rate_host);
    }
  }
}
