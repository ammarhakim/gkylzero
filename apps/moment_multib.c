#include <gkyl_moment_multib.h>
#include <gkyl_moment_multib_priv.h>
#include <gkyl_rrobin_decomp.h>

#include <mpack.h>

static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

struct gkyl_moment_multib_app *
gkyl_moment_multib_app_new(const struct gkyl_moment_multib *mbinp)
{
  struct gkyl_moment_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);

  mbapp->block_geom = gkyl_block_geom_acquire(mbinp->block_geom);
  int ndim = gkyl_block_geom_ndim(mbapp->block_geom);

  int num_ranks;
  gkyl_comm_get_size(mbapp->comm, &num_ranks);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  // construct round-robin decomposition
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    branks[i] = calc_cuts(ndim, bgi->cuts);
  }
  const struct gkyl_rrobin_decomp *rrd = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = mbapp->num_local_blocks = num_blocks;
  mbapp->app = gkyl_malloc(sizeof(struct gkyl_moment_app*[num_local_blocks]));

  mbapp->stat = (struct gkyl_moment_stat) {
  };

  gkyl_free(branks);
  gkyl_rrobin_decomp_release(rrd);
  
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

  gkyl_comm_release(mbapp->comm);
  gkyl_block_geom_release(mbapp->block_geom);
  
  gkyl_free(mbapp->app);
  gkyl_free(mbapp);
}
