#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct fpo_ctx {
  double charge; // charge
  double mass; // mass
  double vth; // thermal velocity
  double Lx; // size of the box
  double gamma; // collision frequency factor in FPO
  double t_end; // end time of simulation
  int num_frame; // number of frames to write out
};

static inline double sq(double x) { return x*x; }

static inline double
bump_maxwellian(double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA, double bUx, double bUy, double bUz, double bS, double bVt)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  double bv2 = (vx - bUx)*(vx - bUx) + (vy - bUy)*(vy - bUy) + (vz - bUz)*(vz - bUz);
  return n/pow(sqrt(2*M_PI*vt*vt), 3)*exp(-v2/(2*vt*vt)) + n/pow(sqrt(2*M_PI*bVt*bVt), 3)*exp(-bv2/(2*bVt*bVt))*(bA*bA)/(bv2 + bS*bS);
}

void
evalDistFuncSquare(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double width = 1.0; // corresponds to a final vth of 2/3
  if(vx>-width && vx<width && vy>-width && vy<width && vz>-width && vz<width) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
evalDistFuncBump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = bump_maxwellian(1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vth, 
    sqrt(0.15), 4.0*app->vth, 0.0, 0.0, 0.14, 3.0*app->vth);
}

void
evalGamma(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;  
  fout[0] = app->gamma;
}

struct fpo_ctx
create_ctx(void)
{
  double n0 = 1.0;
  double vth = 1.0;
  double gamma = 1.0;

  // Characteristic collision frequency: sqrt(2)*n_b*gamma/(3*sqrt(pi)*vth_a^3) 
  double nu = sqrt(2.0)*n0*gamma/(3.0*sqrt(M_PI)*pow(vth, 3.0));
  double end_time = 0.1;

  printf("End time: %f\n", end_time);

  struct fpo_ctx ctx = {
    .mass = 1.0,
    .charge = 0.0,
    .vth = vth,
    .Lx = 1.0,
    .gamma = 1.0, 
    .t_end = end_time, 
    .num_frame = 10, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct fpo_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 1);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 8); 

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  // Create global range.
  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // Create decomposition.
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &cglobal_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
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
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  struct gkyl_vlasov_species square = {
    .name = "square",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -4.0*ctx.vth, -4.0*ctx.vth, -4.0*ctx.vth },
    .upper = { 4.0*ctx.vth, 4.0*ctx.vth, 4.0*ctx.vth }, 
    .cells = { NV, NV, NV },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncSquare,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_FPO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalGamma,
    },
    
    .num_diag_moments = 1,
    .diag_moments = { "FiveMoments" },
  };

  struct gkyl_vlasov_species bump = {
    .name = "bump",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -5.0*ctx.vth, -5.0*ctx.vth, -5.0*ctx.vth },
    .upper = { 5.0*ctx.vth, 5.0*ctx.vth, 5.0*ctx.vth }, 
    .cells = { NV, NV, NV },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncBump,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_FPO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalGamma,
    },
    
    .num_diag_moments = 1,
    .diag_moments = { "FiveMoments" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "fpo_vlasov_relax_1x3v_p1",

    .cdim = 1, .vdim = 3,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { square, bump },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = tend-tcurr;
  int nframe = ctx.num_frame;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr); 

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_vlasov_app_calc_integrated_mom(app, tcurr);   
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);
  gkyl_comm_release(comm);
  gkyl_vlasov_app_write_integrated_mom(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

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
  printf("Species collisions took %g secs\n", stat.species_coll_mom_tm);  
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
