#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct gk_app_ctx {
  double charge_elc; // electron charge
  double mass_elc; // electron mass
  double charge_ion; // ion charge
  double mass_ion; // ion mass
  double n0; // reference density
  double Te; // electron temperature
  double Ti; // ion temperature
  double B0; // reference magnetic field
  double nu_ion; // ion collision frequency
  double Lz; // Box size in z.
  double wavek; // Wave number of ion acoustic wave.
  double kperp; // perpendicular wave number used in Poisson solve
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  int Nz, Nvpar, Nmu; // Number of in z, vpar and mu.
  double final_time; // end time
  int num_frames; // number of output frames
};

// Initial conditions.
void eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  double wavek = app->wavek;

  double alpha = 0.01;
  double perturb = alpha*cos(wavek*z);

  fout[0] = n0*(1.+perturb);
}

void eval_upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti = app->Ti;
  fout[0] = Ti;
}

// Collision frequency.
void eval_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->B0;
}

struct gk_app_ctx
create_ctx(void)
{
  double mi = 1.0; // ion mass
  double me = mi/1836.16; // electron mass
  double qi = 1.0; // ion charge
  double qe = -1.0; // electron charge

  // Reference density and temperature.
  double n0 = 1.0;
  double Te = 1.0;
  double Ti = 1.0;
  double B0 = 1.0;
  double nu_ion = 0.005;
  double wavek = 0.5;

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);

  // Simulation box size (m).
  double Lz = 2.*M_PI/wavek;

  double kperp = 0.;

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = mi*pow(5.0*vtIon,2)/(2.0*B0);

  int Nz = 64;
  int Nvpar = 96;
  int Nmu = 16;

  double final_time = 50.0; 
  double num_frames = 100;

  struct gk_app_ctx ctx = {
    .charge_elc = qe,  .mass_elc = me,
    .charge_ion = qi,  .mass_ion = mi,
    .n0 = n0,
    .Te = Te,  .Ti = Ti,
    .B0 = B0,
    .nu_ion = nu_ion,
    .Lz = Lz,
    .wavek = wavek,
    .kperp = kperp,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .Nz = Nz,  .Nvpar = Nvpar,  .Nmu = Nmu,
    .final_time = final_time,
    .num_frames = num_frames,
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr-1);
    gkyl_gyrokinetic_app_calc_mom(app); gkyl_gyrokinetic_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nz);
  int NVPAR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nmu);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // Create global range.
  int ccells[] = { NZ };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // Create decomposition.
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d=0; d<cdim; d++) {
    if (app_args.use_mpi)
      cuts[d] = app_args.cuts[d];
    else
      cuts[d] = 1;
  }
#else
  for (int d=0; d<cdim; d++)
    cuts[d] = 1;
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
    assert(false);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  } else {
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

  int my_rank, comm_size;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d=0; d<cdim; d++)
    ncuts *= cuts[d];

  if (ncuts != comm_size) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    goto mpifinalize;
  }

  for (int d=0; d<cdim-1; d++) {
    if (cuts[d] > 1) {
      if (my_rank == 0)
        fprintf(stderr, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts[d], d);
      goto mpifinalize;
    }
  }

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .density = eval_density_ion,
      .upar = eval_upar_ion,
      .temp = eval_temp_ion,      
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = eval_nu_ion,
    },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ADIABATIC,
    .electron_mass = ctx.mass_elc,
    .electron_charge = ctx.charge_elc,
    .electron_density = ctx.n0,
    .electron_temp = ctx.Te,
    .bmag_fac = ctx.B0, 
    .kperpSq = pow(ctx.kperp, 2.),
    .fem_parbc = GKYL_FEM_PARPROJ_PERIODIC, 
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_ion_sound_adiabatic_elc_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { -ctx.Lz/2.0 },
    .upper = {  ctx.Lz/2.0 },
    .cells = { NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .bmag_func = bmag_func, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { ion },
    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.final_time;
  double dt = tend-tcurr;
  int nframe = ctx.num_frames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 10 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
    }
    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
  gkyl_gyrokinetic_app_stat_write(app);
  
  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);

  mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
