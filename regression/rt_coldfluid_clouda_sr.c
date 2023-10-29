#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_coldfluid.h>
#include <gkyl_wv_cold_sr_fluid.h>
#include <rt_arg_parse.h>

void
evalColdInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double xcell = xn[0];
  double rho, u;

  if((-2<xcell) && (xcell<-1)) {
    rho = 2.0;
    u = 0.01;
  }
  else if((1<xcell) && (xcell<5)) {
    rho = 1.0;
    u = -0.01;
  }
  else {
    rho = 1e-10; // rho = 0 (note small value)
    u = 0.0;
  }  

  fout[0] = rho;
  fout[1] = rho*u; fout[2] = 0.0; fout[3] = 0.0;
}

void
write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 512);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // equation object
  struct gkyl_wv_eqn *coldf = gkyl_wv_cold_sr_fluid_new();

  struct gkyl_moment_species fluid = {
    .name = "cold",

    .equation = coldf,
    .evolve = 1,
    .init = evalColdInit,
    .split_type = GKYL_WAVE_FWAVE,  
    .limiter = GKYL_MIN_MOD,  

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "cold_sr_fluid_clouda",

    .ndim = 1,
    .lower = { -5.0 },
    .upper = { 10.0 }, 
    .cells = { NX },

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 200.0;

  int nframe = 1000;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = (tend-tcurr)/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);
  
  // simulation complete, free resources
  gkyl_wv_eqn_release(coldf);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);  
  
  return 0;
}
