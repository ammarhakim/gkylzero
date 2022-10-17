#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

const double mu_val = .0000002;
const double v_val = .1;

void
eval_fun(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];

  fout[0] = 1.; // rho
  fout[1] = v_val*sin(x)+1; //rho ux
  fout[2] = 0.; //rho uy
  fout[3] = 0.; //rho uz
}

void
D(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sim_ctx *app = ctx;
  double x = xn[0];
  double mu = mu_val;
  fout[1] = mu;
  fout[2] = mu;
  fout[3] = mu;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

 int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 30);

 struct gkyl_vlasov_fluid_species f = {
   .name = "euler_iso",
   .vt = 100, //warning: if vt is small (relative to perturbation), system will have large diffusive errors

   .init = eval_fun,

   .num_eqn = 4,

   .diffusion = {.D = D},

 };

 // VM app
 struct gkyl_vm vm = {
   .name = "euler_iso_diffusive_pert",

   .cdim = 1,
   .vdim = 3,
   .lower = {0.0},
   .upper = {2.0*M_PI},
   .cells = {NX},
   .poly_order = 2, //warning: poly_order = 1 has large diffusive errors
   .basis_type = app_args.basis_type,
   .cfl_frac = .1,
   .num_periodic_dir = 1,
   .periodic_dirs = {0},

   .num_species = 0,
   .species = {},

   .num_fluid_species = 1,
   .fluid_species = {f},

   .skip_field = true,

   .use_gpu = app_args.use_gpu,
 };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = .4;
  double dt = .00001;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);

  int write_counter = 0;
  gkyl_vlasov_app_write(app, tcurr, write_counter);
  write_counter += 1;
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...\n", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    if (!status.success) {
      fprintf(stderr, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;

  }

  gkyl_vlasov_app_write(app, tcurr, write_counter);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  printf("Writing out final state! \n");

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  printf("Done! \n");

  // simulation complete, free objects
  gkyl_vlasov_app_release(app);

  //output statements
  float tempdx = (vm.upper[0]/((float) NX));
  printf("\n");
  printf("Computing Grid Reynolds Number...");
  printf("NX: %f\n",((float) NX));
  printf("dx: %f\n",tempdx);
  printf("mu: %f\n",mu_val);
  printf("v: %f\n",v_val);
  printf("Reynolds Num: %f\n",(tempdx*(v_val+1.0))/mu_val);
  printf("\n");

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);

  return 0;
}
