#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct ot_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double n0;
  double sigma;
  double Telc;
  double Tion;
  double B0;
  double delta_u0;
  double delta_B0;
  double L;
  double tend;
  bool use_gpu;
};

static inline double
maxwelljuttner3D(double n, double px, double py, double pz, double ux, double uy, double uz, double T, double K_2)
{
  double gamma = 1.0/sqrt(1 - ux*ux - uy*uy - uz*uz);
  // T is in units of m_e c^2
  return n/(4*M_PI*T*K_2)*exp(-(gamma/T)*(sqrt(1.0 + px*px + py*py + pz*pz) - ux*px - uy*py - uz*pz));
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3], vz = xn[4];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double Lx = app->L;
  double Ly = app->L;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double n = app->n0;
  // Half the current to the electrons
  double Jz = 0.5*(B0y*(4.0*M_PI/Lx)*cos(4.0*M_PI*x/Lx) + B0x*(2.0*M_PI/Ly)*cos(2.0*M_PI*y/Ly)) / app->mu0;

  double vdrift_x = -u0x*sin(2.0*M_PI*y/Ly);
  double vdrift_y = u0y*sin(2.0*M_PI*x/Lx);
  double vdrift_z = Jz / (n*qe);

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  // double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  //double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  
  double fv = maxwelljuttner3D(n, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, app->Telc, K_2);
    
  fout[0] = fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3], vz = xn[4];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double Lx = app->L;
  double Ly = app->L;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double n = app->n0;
  // Half the current to the positrons
  double Jz = 0.5*(B0y*(4.0*M_PI/Lx)*cos(4.0*M_PI*x/Lx) + B0x*(2.0*M_PI/Ly)*cos(2.0*M_PI*y/Ly)) / app->mu0;

  double vdrift_x = -u0x*sin(2.0*M_PI*y/Ly);
  double vdrift_y = u0y*sin(2.0*M_PI*x/Lx);
  double vdrift_z = Jz / (n*qi);

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  //double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  //double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  
  double fv = maxwelljuttner3D(n, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, app->Tion, K_2);
    
  fout[0] = fv;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ot_ctx *app = ctx;

  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double Lx = app->L;
  double Ly = app->L;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double n = app->n0;
  double Jz = 0.5*(B0y*(4.0*M_PI/Lx)*cos(4.0*M_PI*x/Lx) + B0x*(2.0*M_PI/Ly)*cos(2.0*M_PI*y/Ly)) / app->mu0;

  double B_x = -B0x*sin(2.0*M_PI*y/Ly);
  double B_y = B0y*sin(4.0*M_PI*x/Lx);
  double B_z = app->B0;

  double E_x = 0.0;
  double E_y = 0.0;
  double E_z = 0.0;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct ot_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double n0 = 1.0; // initial number density
  
  double sigma = 1.0; // B^2/(2*mu0*n_i m_i c^2) = 1.0
  // T_i = T_e = 0.04*m_e c^2 ~ 0.02 MeV
  double Telc = 0.04; // T_e/m_e c^2 = 0.04
  double Tion = 0.04; 

  double B0 = sqrt(sigma*massIon*2.0*mu0);
  // ion Alfvén velocity
  double vAi = sqrt(sigma/(sigma+1));

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;

  // OT initial conditions
  double delta_u0 = 0.5*vAi;
  double delta_B0 = 0.5*B0;

  // domain size and simulation time
  double L = 20.0*M_PI;
  double tend = 10.0/omegaCi;
  
  struct ot_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .n0 = n0,
    .sigma = sigma,
    .Telc = Telc,
    .Tion = Tion,
    .B0 = B0,
    .delta_u0 = delta_u0,
    .delta_B0 = delta_B0,
    .L = L,
    .tend = tend,
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct ot_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0, -1.0, -1.0 },
    .upper = { 1.0, 1.0, 1.0 }, 
    .cells = { 8, 8, 8 },

    .ctx = &ctx,
    .init = evalDistFuncElc,
    
    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -1.0, -1.0, -1.0 },
    .upper = { 1.0, 1.0, 1.0}, 
    .cells = { 8, 8, 8 },

    .ctx = &ctx,
    .init = evalDistFuncIon,
    
    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "ot_vlasov_sr_2x3v",

    .cdim = 2, .vdim = 3,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.L, ctx.L },
    .cells = { 16, 16 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 10;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
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

  gkyl_vlasov_app_stat_write(app);

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
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}