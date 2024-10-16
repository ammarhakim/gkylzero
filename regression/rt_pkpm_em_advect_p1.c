#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct em_advect_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double B0; // reference magnetic field 
  double omega; // frequency of driver
  double Lx; // Domain size (x-direction).
  int Nx; // Cell count (x-direction).
  double t_end; // Final simulation time.
  double init_dt; // Initial time step guess so first step does not generate NaN
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct em_advect_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = maxwellian(1.0, v, 1.0);
  // T_perp/m = 1.0
  fout[1] = maxwellian(1.0, v, 1.0);
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  // no initial flow
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct em_advect_ctx *app = ctx;

  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct em_advect_ctx *app = ctx;
  double x = xn[0];
  double B_x = app->B0;
  double E_z = cos(app->omega*t);
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = E_z;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct em_advect_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = 1.0e-4;
}

struct em_advect_ctx
create_ctx(void)
{
  double Lx = 2.0*M_PI; 
  int Nx = 2; 
  double dx = Lx/Nx;
  double t_end = 100.0;
  // initial dt guess so first step does not generate NaN
  double init_dt = (Lx/Nx)/(3.0);

  struct em_advect_ctx ctx = {
    .mass = 1.0,
    .charge = -1.0,
    .vt = 1.0,
    .B0 = 1.0, 
    .omega = 0.5, 
    .Lx = Lx,
    .Nx = Nx, 
    .t_end = t_end, 
    .init_dt = init_dt, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_pkpm_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  struct em_advect_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -6.0*ctx.vt},
    .upper = { 6.0*ctx.vt}, 
    .cells = { NV },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
    },
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .is_static = true,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
    .ext_em_evolve = true, 
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= app_args.cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_em_advect_p1",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0, 

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = ctx.init_dt;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };


  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    gkyl_pkpm_app_calc_field_energy(app, tcurr);
    gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
    gkyl_pkpm_app_calc_integrated_mom(app, tcurr);

    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_pkpm_app_write_field_energy(app);
  gkyl_pkpm_app_write_integrated_L2_f(app);
  gkyl_pkpm_app_write_integrated_mom(app);
  gkyl_pkpm_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_comm_release(comm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
