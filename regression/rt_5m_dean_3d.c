#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <gkyl_comm.h>

#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct cmfx_ctx {
  double gasGamma; // gas constant
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double beta;
  double vtElc;
  double vtIon;
  double wpe;
  double nuElc;
  double lambdaD; 
  double deltaPhi;
  double sigma2r;
  double rmin;
  double rmax;
  double Lz;
  double tend;
  bool use_gpu;
};

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct cmfx_ctx *app = ctx;

  double gasGamma = app->gasGamma;
  double massElc = app->massElc;
  double Te = app->vtElc*app->vtElc*massElc;
  double n0 = app->n0;

  // Radial profile of electric field; first coordinate is radial coordinate
  double Er = -app->deltaPhi/log(app->rmax/app->rmin)/r;
  double ExB = Er/app->B0;

  double rhoe = n0*massElc*exp(-(r-app->rmin)*(r-app->rmin)/app->sigma2r);
  double rhoeux = rhoe*ExB*sin(theta);
  double rhoeuy = -rhoe*ExB*cos(theta);
  double rhoeuz = 0.0;
  double ere = n0*Te/(gasGamma-1) + 0.5*(rhoeux*rhoeux + rhoeuy*rhoeuy + rhoeuz*rhoeuz)/rhoe;

  fout[0] = rhoe;
  fout[1] = rhoeux; fout[2] = rhoeuy; fout[3] = rhoeuz;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct cmfx_ctx *app = ctx;

  double gasGamma = app->gasGamma;
  double massIon = app->massIon;
  double Ti = app->vtIon*app->vtIon*massIon;
  double n0 = app->n0;

  // Radial profile of electric field; first coordinate is radial coordinate
  double Er = -app->deltaPhi/log(app->rmax/app->rmin)/r;
  double ExB = Er/app->B0;

  double rhoi = n0*massIon*exp(-(r-app->rmin)*(r-app->rmin)/app->sigma2r);
  double rhoiux = rhoi*ExB*sin(theta);
  double rhoiuy = -rhoi*ExB*cos(theta);
  double rhoiuz = 0.0;
  double eri = n0*Ti/(gasGamma-1) + 0.5*(rhoiux*rhoiux + rhoiuy*rhoiuy + rhoiuz*rhoiuz)/rhoi;

  fout[0] = rhoi;
  fout[1] = rhoiux; fout[2] = rhoiuy; fout[3] = rhoiuz;
  fout[4] = eri;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct cmfx_ctx *app = ctx;

  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;

  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;

  // electric field
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // magnetic field
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct cmfx_ctx *app = ctx;

  // Radial profile of electric field; first coordinate is radial coordinate
  double Er = -app->deltaPhi/log(app->rmax/app->rmin)/r;

  double Ex = Er*cos(theta);
  double Ey = Er*sin(theta);
  double Ez = 0.0;

  double Bx = 0.0;
  double By = 0.0;
  double Bz = app->B0;

  // electric field
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // magnetic field
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

// map (r,theta) -> (x,y)
void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1], z = xc[2];
  xp[0] = r*cos(th); xp[1] = r*sin(th); xp[2] = z;
}

struct cmfx_ctx
create_ctx(void)
{
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space
  double lightSpeed = 1.0/sqrt(epsilon0*mu0);

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.153; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  double vAe = 0.5;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  double beta = 0.02;
  double vtElc = vAe*sqrt(beta/2.0);

  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon); //Ti/Te = 1.0

  // frequencies, skin depths, and Debye length
  double omegaCi = chargeIon*B0/massIon;
  double wpe = sqrt(n0*chargeElc*chargeElc/(epsilon0*massElc));
  double de = lightSpeed/wpe;
  double wpi = sqrt(n0*chargeIon*chargeIon/(epsilon0*massIon));
  double di = lightSpeed/wpi;
  double lambdaD = vtElc/wpe;

  // dummy collision frequency (not used). Gives us the denormalized density
  double nuElc = 1.0e-3*omegaCi;

  // domain size and simulation time
  double rmin = 2.0*di;
  double rmax = 18.0*di;
  double L = rmax - rmin;
  double sigma2r = (L/2.0)*(L/2.0);
  double Lz = 16.0*di;
  double tend = 100.0/omegaCi;

  // potential drop and E x B velocity magnitude at inner radius
  double deltaPhi = 0.5;
  
  struct cmfx_ctx ctx = {
    .gasGamma = gasGamma, 
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .wpe = wpe, 
    .nuElc = nuElc, 
    .lambdaD = lambdaD, 
    .deltaPhi = deltaPhi, 
    .sigma2r = sigma2r, 
    .rmin = rmin,
    .rmax = rmax, 
    .Lz = Lz, 
    .tend = tend,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif  
  
  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 16);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 16);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct cmfx_ctx ctx = create_ctx(); // context for init functions

  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gasGamma);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gasGamma);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT }, 
  };  

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif

  // create global range
  int cells[] = { NX, NY, NZ };
  struct gkyl_range globalr;
  gkyl_create_global_range(3, cells, &globalr);
  
  // create decomposition
  int cuts[] = { 1, 1, 1 };
#ifdef GKYL_HAVE_MPI  
  if (app_args.use_mpi) {
    cuts[0] = app_args.cuts[0];
    cuts[1] = app_args.cuts[1];
    cuts[2] = app_args.cuts[2];
  }
#endif  
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(3, cuts, &globalr);
  
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
        .decomp = decomp
      }
    );
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = cuts[0]*cuts[1]*cuts[2];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }

  struct gkyl_moment app_inp = {
    .name = "5m_cmfx",

    .ndim = 3,
    .lower = { ctx.rmin, 0.0, -ctx.Lz },
    .upper = { ctx.rmax, 2.0*M_PI, ctx.Lz }, 
    .cells = { NX, NY, NZ },

    .mapc2p = mapc2p, // mapping of computational to physical space

    .num_periodic_dir = 2,
    .periodic_dirs = { 1,2 },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      .ctx = &ctx,

      .is_ext_em_static = true,
      .ext_em_func = evalExtEmInit, 
      
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    },

    .has_low_inp = true,
    .low_inp = {
      .comm = comm,
      .local_range = decomp->ranges[my_rank]
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Compute thermal velocities and inner radius ExB velocity from speed of light normalization
  double vtElc_denorm = ctx.vtElc*299792458.0;
  double vtIon_denorm = ctx.vtIon*299792458.0;
  double Er_max = -ctx.deltaPhi/log(ctx.rmax/ctx.rmin)/ctx.rmin;
  double ExB_max = Er_max/ctx.B0*299792458.0;
  // Compute electron temperature in ev from thermal velocity
  double Telc_denorm = vtElc_denorm*vtElc_denorm*9.109383632e-31;
  // Compute density in 1/m^3 from collisionality (and plasma parameter)
  // Note plasma parameter is computed as Lambda = wpe/nu_ee, which is missing a factor of 2pi/ln(Lambda) (1/2-1/3)
  double nElc_denorm = pow(8.85418781e-12*Telc_denorm, 3)/(ctx.wpe/ctx.nuElc*ctx.wpe/ctx.nuElc*pow(1.60217663e-19, 6));
  // Compute Debye length in meters from density and temperature
  double lambdaD_denorm = sqrt(8.85418781e-12*Telc_denorm/(1.60217663e-19*1.60217663e-19*nElc_denorm));
  // Compute background magnetic field in Tesla from beta, density, and temperature
  double B0_denorm = sqrt(2.0*1.256637062e-6*nElc_denorm*Telc_denorm/ctx.beta);
  // Compute axial extent in meters
  double axial_denorm = ctx.Lz/ctx.lambdaD*lambdaD_denorm;
  // Compute radial extent in meters
  double radial_denorm = (ctx.rmax-ctx.rmin)/ctx.lambdaD*lambdaD_denorm;

  gkyl_moment_app_cout(app, stdout, "Axial direction is %g Debye lengths\n", ctx.Lz/ctx.lambdaD);
  gkyl_moment_app_cout(app, stdout, "Electron Plasma parameter (wpe/nu_ee) = %g\n", ctx.wpe/ctx.nuElc);
  gkyl_moment_app_cout(app, stdout, "Electron thermal velocity in m/s = %g\n", vtElc_denorm);
  gkyl_moment_app_cout(app, stdout, "Ion thermal velocity in m/s = %g\n", vtIon_denorm);
  gkyl_moment_app_cout(app, stdout, "ExB velocity at inner radius normalized to electron thermal velocity = %g\n", ExB_max/vtElc_denorm);
  gkyl_moment_app_cout(app, stdout, "ExB velocity at inner radius normalized to ion thermal velocity = %g\n", ExB_max/vtIon_denorm);
  gkyl_moment_app_cout(app, stdout, "Electron temperature in ev = %g\n", Telc_denorm/1.60217663e-19);
  gkyl_moment_app_cout(app, stdout, "Electron density in 1/m^3 = %g\n", nElc_denorm);
  gkyl_moment_app_cout(app, stdout, "Electron Debye length in meters = %g\n", lambdaD_denorm);
  gkyl_moment_app_cout(app, stdout, "Background magnetic field in Tesla = %g\n", B0_denorm);
  gkyl_moment_app_cout(app, stdout, "Axial extent in meters = %g\n", 2.0*axial_denorm);
  gkyl_moment_app_cout(app, stdout, "Radial extent in meters = %g\n", radial_denorm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);
  gkyl_moment_app_calc_integrated_mom(app, tcurr);  

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_moment_app_calc_integrated_mom(app, tcurr);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_write_integrated_mom(app);
    
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app);  
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif    
  
  return 0;
}
