#include <stdarg.h>

#include <gkyl_alloc.h>

#include <gkyl_gyrokinetic_mb_priv.h>
#include <gkyl_app_priv.h>

gkyl_gyrokinetic_mb_app*
gkyl_gyrokinetic_mb_app_new(struct gkyl_gk_mb *inp)
{
  gkyl_gyrokinetic_mb_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, inp->name);

  mbapp->cdim = inp->cdim;
  mbapp->vdim = inp->vdim;
  mbapp->poly_order = inp->poly_order;

  mbapp->cfl = inp->cfl_frac == 0 ? 1.0 : inp->cfl_frac;

#ifdef GKYL_HAVE_CUDA
  mbapp->use_gpu = inp->use_gpu;
#else
  mbapp->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  mbapp->num_blocks = inp->num_blocks;
  mbapp->num_blocks_local = inp->num_blocks;

  for (int bidx=0; bidx<inp.num_blocks; bidx++) {
    struct gkyl_gk *blinp = &(inp->blocks[bidx]);

    // Block name = <the name of the app>_b#.
    const char *fmt = "%s_b%d";
    int sz = gkyl_calc_strlen(fmt, inp->name, bidx);
    char appNm[sz+1]; // ensures no buffer overflow
    snprintf(appNm, sizeof AppNm, fmt, inp->name, bidx);
    blinp->name = appNm;

    // Fill the input for each block with common values.
    blinp->cdim = inp->cdim;
    blinp->vdim = inp->vdim;
    blinp->poly_order = inp->poly_order;
    blinp->basis_type = inp->basis_type;

    blinp->use_gpu = inp->use_gpu;

    blinp->num_species = inp->num_species;
    for (int s=0; s<inp->num_species; s++)
      blinp->species[s] = inp->species[s]

    blinp->num_neut_species = inp->num_neut_species;
    for (int s=0; s<inp->num_neut_species; s++)
      blinp->neut_species[s] = inp->neut_species[s]

    blinp->skip_field = inp->field;
    blinp->field = inp->field;

    blinp->num_periodic_dir = inp->num_periodic_dir;
    memcpy(blinp->periodic_dirs, inp->periodic_dirs, inp->num_periodic_dir*sizeof(int));

    // Create a new app for each block.
    mbapp->blocks[bidx] = gkyl_gyrokinetic_app_new(blinp);
  }

  return mbapp;
}

void
gkyl_gyrokinetic_mb_app_apply_ic(gkyl_gyrokinetic_mb_app* app, double t0)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_apply_ic(app->blocks[bidx], t0);
  }
}

struct gkyl_update_status
gkyl_gyrokinetic_mb_update(gkyl_gyrokinetic_mb_app* app, double dt)
{
  struct gkyl_update_status status[GKYL_MAX_BLOCKS];
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    status[bidx] = gkyl_gyrokinetic_update(app->blocks[bidx], dt);
  }

  return status[0];
}

void
gkyl_gyrokinetic_mb_app_calc_mom(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_mom(app->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_mb_app_calc_integrated_mom(gkyl_gyrokinetic_mb_app* app, double tm)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app->blocks[bidx], tm);
  }
}

void
gkyl_gyrokinetic_mb_app_write_mom(gkyl_gyrokinetic_mb_app* app, double tm, int frame)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_mom(app->blocks[bidx], tm, frame);
  }
}

void
gkyl_gyrokinetic_mb_app_write_integrated_mom(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_integrated_mom(app->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_mb_app_write(gkyl_gyrokinetic_mb_app* app, double tm, int frame)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write(app->blocks[bidx], tm, frame);
  }
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_mb_app_read_from_frame(gkyl_gyrokinetic_mb_app *app, int frame)
{
  // NYI
}

void
gkyl_gyrokinetic_mb_app_stat_write(gkyl_gyrokinetic_mb_app* app)
{
  // NYI
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_app_stat(gkyl_gyrokinetic_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gk_species_tm(app->blocks[bidx]);
    gk_species_coll_tm(app->blocks[bidx]);
  }
  return app->blocks[0]->stat;
}

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_mb_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_mb_app_cout(const gkyl_gyrokinetic_mb_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gk_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_mb_app_release(gkyl_gyrokinetic_mb_app* mbapp)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_release(mbapp->blocks[bidx])
  }

  gkyl_free(mbapp);
}
