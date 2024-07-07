#include <gkyl_moment_multib.h>
#include <gkyl_moment_multib_priv.h>

#include <mpack.h>

struct gkyl_moment_multib_app *
gkyl_moment_multib_app_new(struct gkyl_moment_multib *mbinp)
{
  struct gkyl_moment_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  
  int num_blocks = mbapp->num_blocks = gkyl_block_geom_num_blocks(mbinp->block_geom);
  mbapp->app = gkyl_malloc(sizeof(struct gkyl_moment_app*[num_blocks]));

  mbapp->stat = (struct gkyl_moment_stat) {
  };
  
  return mbapp;
}

double
gkyl_moment_multib_app_max_dt(gkyl_moment_multib_app *app)
{
  // TODO
  return 0;
}


void
gkyl_moment_multib_app_apply_ic(gkyl_moment_multib_app* app, double t0)
{
  // TODO
}

void
gkyl_moment_multib_app_apply_ic_field(gkyl_moment_multib_app* app, double t0)
{
  // TODO
}

void
gkyl_moment_multib_app_apply_ic_species(gkyl_moment_multib_app* app, int sidx, double t0)
{
  // TODO
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_field(gkyl_moment_multib_app *app,
  int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_species(gkyl_moment_multib_app *app,
  int sidx, int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };  
}

void
gkyl_moment_multib_app_cout(const gkyl_moment_multib_app* app, FILE *fp, const char *fmt, ...)
{
  // TODO
}

void
gkyl_moment_multib_app_write(const gkyl_moment_multib_app* app, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_field(const gkyl_moment_multib_app *app, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_species(const gkyl_moment_multib_app* app, int sidx, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_field_energy(gkyl_moment_multib_app *app)
{
    // TODO
}

void
gkyl_moment_multib_app_write_integrated_mom(gkyl_moment_multib_app *app)
{
  // TODO
}

void
gkyl_moment_multib_app_stat_write(const gkyl_moment_multib_app *app)
{
  // TODO
}

struct gkyl_update_status
gkyl_moment_multib_update(gkyl_moment_multib_app *app, double dt)
{
  // TODO
  return (struct gkyl_update_status) { };
}

void
gkyl_moment_multib_app_calc_field_energy(gkyl_moment_multib_app *app, double tm)
{
  // TODO  
}

void
gkyl_moment_multib_app_get_field_energy(gkyl_moment_multib_app *app, double *vals)
{
  // TODO
}

void
gkyl_moment_multib_app_calc_integrated_mom(gkyl_moment_multib_app *app, double tm)
{
  // TODO
}

struct gkyl_moment_stat
gkyl_moment_multib_app_stat(gkyl_moment_multib_app *app)
{
  return app->stat;
}

void
gkyl_moment_multib_app_release(gkyl_moment_multib_app* mbapp)
{

  
  gkyl_free(mbapp->app);
  gkyl_free(mbapp);
}
