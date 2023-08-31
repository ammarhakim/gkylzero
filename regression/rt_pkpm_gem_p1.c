#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct pkpm_gem_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double nbOverN0;
  double vAe;
  double B0;
  double guide;
  double w0;
  double psi0;
  double noise_amp;
  double noise_index;
  int k_init;
  int k_final;
  double T_e;
  double T_i;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double Lx;
  double Ly;
  double tend;
  double min_dt;
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double temp, double mass)
{
  double v2 = v*v;
  return n/sqrt(2.0*M_PI*temp/mass)*exp(-v2/(2.0*temp/mass));
}

static inline double
sech2(double x)
{
  return 1.0/(cosh(x)*cosh(x));
}

static inline void
noise_init(double noise_amp, double noise_index, int k_init, int k_final, double Lx, double Ly, double x, double y, double noise[3])
{
  pcg64_random_t rng = gkyl_pcg64_init(0);
  double kindex = (noise_index + 1.0) / 2.0;
  double B_amp = 0.0; 
  double B_phase = 0.0;
  for (int i = k_init; i < k_final; ++i) {
    B_amp = gkyl_pcg64_rand_double(&rng);
    B_phase = gkyl_pcg64_rand_double(&rng);

    noise[0] -= 2.0*(2.0*M_PI/Ly)*(Lx/(i*2.0*M_PI))*B_amp*sin(2.0*M_PI*y/Ly)*(cos(2.0*M_PI*y/Ly)+1)*cos(i*2.0*M_PI*x/Lx +  2.0*M_PI*B_phase)*pow(i,kindex);
    noise[1] += B_amp*(cos(2.0*M_PI*y/Ly) + 1.0)*(cos(2.0*M_PI*y/Ly) + 1.0)*sin(i*2.0*M_PI*x/Lx + 2.0*M_PI*B_phase)*pow(i,kindex);
    noise[2] += (2.0*M_PI*i/Lx)*B_amp*(cos(2.0*M_PI*y/Ly) + 1.0)*(cos(2.0*M_PI*y/Ly) + 1.0)*cos(i*2.0*M_PI*x/Lx + 2*M_PI*B_phase)*pow(i,kindex) + 
                 2.0*(2.0*M_PI/Ly)*(2.0*M_PI/Ly)*(Lx/(i*2.0*M_PI))*B_amp*(sin(2.0*M_PI*y/Ly)*sin(2.0*M_PI*y/Ly) - cos(2.0*M_PI*y/Ly)*(cos(2.0*M_PI*y/Ly)+1.0))*cos(i*2.0*M_PI*x/Lx +  2.*M_PI*B_phase)*pow(i,kindex);
  }
  double kdiff = floor(k_final) - floor(k_init) + 1.0;
  noise[0] = noise_amp*noise[0]/sqrt(2.0*kdiff*kdiff/3.0);
  noise[1] = noise_amp*noise[1]/sqrt(2.0*kdiff*kdiff/3.0);
  noise[2] = noise_amp*noise[2]/sqrt(2.0*kdiff*kdiff/3.0);
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double nbOverN0 = app->nbOverN0;
  double B0 = app->B0;
  double w0 = app->w0;

  double b1x = B0*(tanh(y/w0)-tanh((y+Ly)/w0) - tanh((y-Ly)/w0));
  double n = (0.5*(B0*B0 - b1x*b1x) + nbOverN0*(T_i+T_e))/(T_i+T_e); // only correct for mu0,n0 = 1
  
  double fv = maxwellian(n, vx, T_e, me);
    
  fout[0] = fv;
  fout[1] = T_e/me*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double nbOverN0 = app->nbOverN0;
  double B0 = app->B0;
  double w0 = app->w0;

  double b1x = B0*(tanh(y/w0)-tanh((y+Ly)/w0) - tanh((y-Ly)/w0));
  double n = (0.5*(B0*B0 - b1x*b1x) + nbOverN0*(T_i+T_e))/(T_i+T_e); // only correct for mu0,n0 = 1
  
  double fv = maxwellian(n, vx, T_i, mi);
    
  fout[0] = fv;
  fout[1] = T_i/mi*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double T_e = app->T_e;
  double T_i = app->T_i;
  double noise_amp = app->noise_amp;
  double noise_index = app->noise_index;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double TeFrac = T_e/(T_i+T_e);

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double noise[3] = {0.0};
  noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y, noise);

  double Jx = 0.0;
  double Jy = 0.0;
  double Jz = B0/w0*(sech2(y/w0) - sech2((y+Ly)/w0) - sech2((y-Ly)/w0)) 
              -psi0*cos(pi_2*x/Lx)*cos(M_PI*y/Ly)*((pi_2/Lx)*(pi_2/Lx) + (M_PI/Ly)*(M_PI/Ly)) + noise[2];

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = me*Jz/qe*TeFrac;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double T_e = app->T_e;
  double T_i = app->T_i;
  double noise_amp = app->noise_amp;
  double noise_index = app->noise_index;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double TiFrac = T_i/(T_i+T_e);

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double noise[3] = {0.0};
  noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y, noise);

  double Jx = 0.0;
  double Jy = 0.0;
  double Jz = B0/w0*(sech2(y/w0) - sech2((y+Ly)/w0) - sech2((y-Ly)/w0)) 
              -psi0*cos(pi_2*x/Lx)*cos(M_PI*y/Ly)*((pi_2/Lx)*(pi_2/Lx) + (M_PI/Ly)*(M_PI/Ly)) + noise[2];

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = mi*Jz/qi*TiFrac;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;

  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double guide = app->guide;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double noise_amp = app->noise_amp;
  double noise_index = app->noise_index;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double noise[3] = {0.0};
  noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y, noise);

  double b1x = -B0*(tanh(y/w0) - tanh((y+Ly)/w0) - tanh((y-Ly)/w0));
  double b1y = 0.0;
  double b1z = guide;

  double E_x = 0.0;
  double E_y = 0.0;
  double E_z = 0.0;
  double B_x = b1x + psi0 * (M_PI / Ly) * cos(2 * M_PI * x / Lx) * sin(M_PI * y / Ly) + noise[0];
  double B_y = b1y - psi0 * (2 * M_PI / Lx) * sin(2 * M_PI * x / Lx) * cos(M_PI * y / Ly) + noise[1];
  double B_z = b1z;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_gem_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 25.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 0.2; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density
  double nbOverN0 = 0.2; // background number density

  double vAe = 1.0/4.0;
  double vAi = vAe/sqrt(massIon);
  double beta_elc = 1.0/6.0;
  // B0 has 1/sqrt(2) factor because B is equally subdivided between
  // guide field and in-plane field
  double B0 = vAe*sqrt(mu0*n0*massElc);  
  double guide = 0.1*B0;
  double tot_B = sqrt(B0*B0 + guide*guide);

  double T_e = beta_elc*tot_B*tot_B/(2.0*mu0*n0);
  double T_i = T_e/Te_Ti;
  double vtElc = sqrt(2.0*T_e/massElc);
  double vtIon = sqrt(2.0*T_i/massIon);

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;

  // Layer width and perturbation
  double w0 = 0.5*di;
  double psi0 = 0.1*B0*di;

  // noise levels for perturbation
  double noise_amp = 0.01*B0;
  int k_init = 1;            // first wave mode to perturb with noise, 1.0 correspond to box size
  int k_final = 20;          // last wave mode to perturb with noise
  double noise_index = -1.0; // spectral index of the noise

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double Lx = 8.0*M_PI*di;
  double Ly = 4.0*M_PI*di;
  double tend = 20.0/omegaCi;
  
  struct pkpm_gem_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .T_e = T_e,
    .T_i = T_i,
    .n0 = n0,
    .nbOverN0 = nbOverN0,
    .vAe = vAe,
    .B0 = B0,
    .guide = guide,
    .w0 = w0,
    .psi0 = psi0,
    .noise_amp = noise_amp,
    .noise_index = noise_index,
    .k_init = k_init,
    .k_final = k_final,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .Lx = Lx,
    .Ly = Ly,
    .tend = tend,
    .min_dt = 1.0e-2, 
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
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 32);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_gem_ctx ctx = create_ctx(); // context for init functions

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int cells[] = { NX, NY };
  struct gkyl_range globalr;
  gkyl_create_global_range(2, cells, &globalr);
  
  // create decomposition
  int cuts[] = { 1, 1 };
#ifdef GKYL_HAVE_MPI  
  if (app_args.use_mpi) {
    cuts[0] = app_args.cuts[0];
    cuts[1] = app_args.cuts[1];
  }
#endif 
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(2, cuts, &globalr);

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

  int ncuts = cuts[0]*cuts[1];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }

  // electron momentum                                                                                              
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .diffusion = {.D = 1.0e-3, .order=4},
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -8.0 * ctx.vtElc},
    .upper = { 8.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .num_diag_moments = 0,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  // ion momentum                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .diffusion = {.D = 1.0e-3, .order=4},
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -8.0 * ctx.vtIon},
    .upper = { 8.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .num_diag_moments = 0,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_gem_p1",

    .cdim = 2, .vdim = 1,
    .lower = { -ctx.Lx/2.0, -ctx.Ly/2.0 },
    .upper = { ctx.Lx/2.0, ctx.Ly/2.0 },
    .cells = { NX, NY },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    
    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { fluid_elc, fluid_ion },
    .field = field,

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
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 20;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_field_energy(app, tcurr);
  gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
      gkyl_vlasov_app_calc_field_energy(app, tcurr);
      gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
      gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
    }
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      gkyl_vlasov_app_cout(app, stdout, "** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_vlasov_app_write(app, tcurr, 1000);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_vlasov_app_calc_field_energy(app, tcurr);
  gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
  gkyl_vlasov_app_write_field_energy(app);
  gkyl_vlasov_app_write_integrated_L2_f(app);
  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);

  gkyl_vlasov_app_cout(app, stdout, "Species BCs took %g secs\n", stat.species_bc_tm);
  gkyl_vlasov_app_cout(app, stdout, "Fluid Species BCs took %g secs\n", stat.fluid_species_bc_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field BCs took %g secs\n", stat.field_bc_tm);

  gkyl_vlasov_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  
  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
