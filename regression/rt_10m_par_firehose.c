#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct par_firehose_ctx {
  // Fundamental constants
  double epsilon0; // permittivity of free space
  double mu0; // permeability of free space
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double n0; // reference density
  double B0; // reference magnetic field strength
  double TiPar; // parallel ion temperature
  double TiPerp; // perpendicular ion temperature
  double Te; // electron temperature
  double k0_elc; // closure parameter for electrons
  double k0_ion; // closure parameter for ions
  double noise_amp; // amplitude of perturbation
  int k_init; // initial wavenumber to perturb
  int k_final; // final wavenumber to perturb
  double Lx; // size of the box
  double tend; // end time of the simulation
};

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct par_firehose_ctx *app = ctx;

  double massElc = app->massElc;
  double Te = app->Te;
  double n0 = app->n0;

  double rhoe = massElc*n0;
  double momxe = 0.0;
  double momye = 0.0;
  double momze = 0.0;
  double pxxe = n0*Te + momxe*momxe/rhoe;
  double pxye = momxe*momye/rhoe;
  double pxze = momxe*momze/rhoe;
  double pyye = n0*Te + momye*momye/rhoe;
  double pyze = momye*momye/rhoe;
  double pzze = n0*Te + momze*momze/rhoe;

  fout[0] = rhoe;
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  fout[4] = pxxe; fout[5] = pxye; fout[6] = pxze;  
  fout[7] = pyye; fout[8] = pyze; fout[9] = pzze;
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct par_firehose_ctx *app = ctx;

  double massIon = app->massIon;
  double TiPar = app->TiPar;
  double TiPerp = app->TiPerp;
  double n0 = app->n0;

  double rhoi = massIon*n0;
  double momxi = 0.0;
  double momyi = 0.0;
  double momzi = 0.0;
  double pxxi = n0*TiPar + momxi*momxi/rhoi;
  double pxyi = momxi*momyi/rhoi;
  double pxzi = momxi*momzi/rhoi;
  double pyyi = n0*TiPerp + momyi*momyi/rhoi;
  double pyzi = momyi*momyi/rhoi;
  double pzzi = n0*TiPerp + momzi*momzi/rhoi;

  fout[0] = rhoi;
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  fout[4] = pxxi; fout[5] = pxyi; fout[6] = pxzi;  
  fout[7] = pyyi; fout[8] = pyzi; fout[9] = pzzi;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct par_firehose_ctx *app = ctx;

  pcg64_random_t rng = gkyl_pcg64_init(0);

  double Bx = app->B0;
  double By = 0.0, Bz = 0.0;

  double alpha = app->noise_amp*Bx;
  double Lx = app->Lx;
  double kx = 2.0*M_PI/Lx;
  for (int i = app->k_init; i < app->k_final; ++i) {
    By -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    Bz -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
  }

  // electric field
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

struct par_firehose_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space
  double lightSpeed = 1.0/sqrt(epsilon0*mu0);

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.0; // ion mass
  double chargeIon = 1.0; // ion charge

  // Ion and electron temperature assumed equal, Ti/Te = 1.0
  double n0 = 1.0; // initial number density
  double vAe = 0.0125;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  // Parallel firehose unstable for betaPerp - betaPar + 2 < 0.
  double beta = 300.0/M_PI; // Trace proton plasma beta = 2 mu0 ni Ti/B0^2 = (betaPar + 2 betaPerp)/3
  double dbeta = 100.0; // dbeta = beta_parallel - beta_perp
  double betaPar = beta + 2.0*dbeta/3.0; // parallel proton plasma beta = 2 mu0 ni T_paralleli/B0^2
  double betaPerp = beta - dbeta/3.0; // perp proton plasma beta = 2 mu0 ni T_perpi/B0^2  

  double vtElc = vAe*sqrt(beta);
  double Te = vtElc*vtElc*massElc/2.0;

  double TiPar = vAe*vAe*(betaPar*massElc/2.0);
  double TiPerp = vAe*vAe*(betaPerp*massElc/2.0);

  // frequencies, skin depths, and Debye length
  double omegaCi = chargeIon*B0/massIon;
  double wpe = sqrt(n0*chargeElc*chargeElc/(epsilon0*massElc));
  double de = lightSpeed/wpe;
  double wpi = sqrt(n0*chargeIon*chargeIon/(epsilon0*massIon));
  double di = lightSpeed/wpi;
  double lambdaD = vtElc/wpe;

  // noise levels for perturbation
  double noise_amp = 1.0e-6*B0;
  int k_init = 1;   // first wave mode to perturb with noise, 1 correspond to box size
  int k_final = 48; // last wave mode to perturb with noise

  // k0 for closure
  double k0_elc = 0.1/de; // rho_e ~ 10, k0e = 1.0/rho_e ~ 0.1
  double k0_ion = 0.1/di; // rho_i ~ 420, k0e = 1.0/rho_i ~ 0.002

  // domain size and simulation time
  double Lx = 300.0*di;
  double tend = 100.0/omegaCi;

  struct par_firehose_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .n0 = n0,
    .B0 = B0,
    .noise_amp = noise_amp,
    .k_init = k_init,
    .k_final = k_final,
    .Te = vtElc,
    .TiPar = TiPar,
    .TiPerp = TiPerp, 
    .k0_elc = k0_elc, 
    .k0_ion = k0_ion, 
    .Lx = Lx,
    .tend = tend,
  };
  return ctx;
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

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif  

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 560);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct par_firehose_ctx ctx = create_ctx(); // context for init functions
  
  // electron/ion equations
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_elc);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_ion);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .equation = elc_ten_moment,
    .evolve = 1,
    .init = evalElcInit,
    .ctx = &ctx,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .equation = ion_ten_moment,
    .evolve = 1,
    .init = evalIonInit,
    .ctx = &ctx,
  };  

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int cells[] = { NX };
  struct gkyl_range globalr;
  gkyl_create_global_range(1, cells, &globalr);
  
  // create decomposition
  int cuts[] = { 1 };
#ifdef GKYL_HAVE_MPI  
  if (app_args.use_mpi) {
    cuts[0] = app_args.cuts[0];
  }
#endif 
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(1, cuts, &globalr);

  // construct communcator for use in app
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu        
      }
    );
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu      
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = cuts[0];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }

  // VM app
  struct gkyl_moment app_inp = {
    .name = "10m_par_firehose",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      .ctx = &ctx,
    }, 

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  int nframe = 100;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_moment_app_calc_field_energy(app, tcurr);
  gkyl_moment_app_calc_integrated_mom(app, tcurr);  

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    // Only calculate the integrated moments and field energy every 100 steps
    if (step % 100 == 0) {
      gkyl_moment_app_calc_field_energy(app, tcurr);
      gkyl_moment_app_calc_integrated_mom(app, tcurr);
    }
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_moment_app_calc_field_energy(app, tcurr);
  gkyl_moment_app_calc_integrated_mom(app, tcurr);
  gkyl_moment_app_write_field_energy(app);
  gkyl_moment_app_write_integrated_mom(app);
    
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);
  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(
    app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(
    app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(
    app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(
    app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);

mpifinalize:;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
