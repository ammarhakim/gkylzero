#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>

#include <mpack.h>

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(struct gkyl_gyrokinetic_multib *mbinp)
{
  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  
  int num_blocks = mbapp->num_blocks = gkyl_block_geom_num_blocks(mbinp->block_geom);
  mbapp->app = gkyl_malloc(sizeof(struct gkyl_gyrokinetic_app*[num_blocks]));

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {
  };
  
  return mbapp;
}


void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  // TO DO
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}


void gkyl_gyrokinetic_multib_app_calc_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_integrated_neut_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}


void gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_coll_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_drag(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_emissivity(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_integrated_moms(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  // TO DO
}


void gkyl_gyrokinetic_multib_app_write_iz_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_iz_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_source_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_max_corr_status(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_read_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  // TO DO
}

struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  // TO DO
  return (struct gkyl_update_status) { };
}

struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  return app->stat;
}

void gkyl_gyrokinetic_multib_app_species_ktm_rhs(gkyl_gyrokinetic_multib_app* app, int update_vol_term)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  gkyl_free(mbapp->app);
  gkyl_free(mbapp);
}
