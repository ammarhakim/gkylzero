#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_app_ctx {
  // Geometry and magnetic field.
  double R_axis;    // Magnetic axis major radius [m].
  double R0;        // Major radius of the simulation box [m].
  double a_mid;     // Minor radius at outboard midplane [m].
  double r0;        // Minor radius of the simulation box [m].
  double B0;        // Magnetic field magnitude in the simulation box [T].
  double qSep;      // Safety factor at the separatrix.
  double sSep;      // Magnetic shear at the separatrix.
  double kappa;     // Elongation (=1 for no elongation).
  double delta;     // Triangularity (=0 for no triangularity).
  double q0;        // Magnetic safety factor in the center of domain.
  double Lz;        // Domain size along magnetic field.
  double z_min;  double z_max;

  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Collisions.
  double nuFrac;  double nuIon;

  // Source parameters.
  double n_src;  double Ti_src;

  // Grid parameters.
  int num_cell_z;
  int num_cell_vpar;
  int num_cell_mu;
  int poly_order;
  double vpar_max_ion;  double mu_max_ion;
  double final_time;   int num_frames;
};

double qprofile(double r, double a_mid, double qSep, double sSep) 
{
  // Magnetic safety factor as a function of minor radius r.
  return qSep/(1.-sSep*((r-a_mid)/a_mid));
}

double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return R_axis + r*cos(theta + asin(delta)*sin(theta));
}

double Z_rtheta(double r, double theta, void *ctx)
{
  // Z (height) as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*sin(theta);
}

// Partial derivatives of R(r,theta) and Z(r,theta)
double dRdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return cos(theta + asin(delta)*sin(theta));
}
double dRdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return -r*sin(theta + asin(delta)*sin(theta))*(1.+asin(delta)*cos(theta));
}
double dZdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*sin(theta);
}
double dZdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*cos(theta);
}

double Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r,theta,ctx)*( dRdr(r,theta,ctx) * dZdtheta(r,theta,ctx)
		                -dRdtheta(r,theta,ctx) * dZdr(r,theta,ctx) );
}

struct integrand_ctx {
  struct gk_app_ctx *app_ctx;
  double r;
};
double integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}
double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  double B0 = app->B0;
  double R_axis = app->R_axis;
  double a_mid = app->a_mid;
  double qSep = app->qSep;
  double sSep = app->sSep;

  return ( B0*R_axis/(2.*M_PI*qprofile(r,a_mid,qSep,sSep)) )*integral.res;
}

double alpha(double r, double theta, double phi, void *ctx)
{
  double twrap = theta;
  while (twrap < -M_PI) twrap = twrap+2.*M_PI;
  while (M_PI < twrap) twrap = twrap-2.*M_PI;

  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., twrap, 7, 1e-10);
  } else {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, twrap, 0., 7, 1e-10);
    integral.res = -integral.res;
  }

  double B0 = app->B0;
  double R_axis = app->R_axis;

  return phi - B0*R_axis*integral.res/dPsidr(r,theta,ctx);
}

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}
double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

// Ion source profiles.
void density_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Lz = app->Lz;
  double n_src = app->n_src;
  double z = xn[0];

  if (fabs(z) < Lz/4.) {
    fout[0] = app->n_src;
  } else {
    fout[0] = 1e-40;
  }
}
void upar_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti_src = app->Ti_src;
  double z = xn[0];
  fout[0] = Ti_src;
}

// Ion initial conditions
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Lz = app->Lz;
  double mi = app->mi;

  double z = xn[0];
  double Ls = Lz/4.;
  double xn0[] = {0.};
  double n_src0[1];
  density_src(t,xn0,n_src0,ctx);
  double effSrc  = n_src0[0];
  double eV = GKYL_ELEMENTARY_CHARGE;
  double Te_src  = 410*eV;
  double c_ss    = sqrt(5./3.*Te_src/mi);
  double nPeak   = 89.*sqrt(5.)/3./c_ss*Ls*effSrc/2.;
  double perturb = 0.;
  if (fabs(z) <= Ls) {
    fout[0] = nPeak*(1.+sqrt(1.-pow(z/Ls,2)))/2.*(1.+perturb);
  } else {
    fout[0] = nPeak/2.*(1.+perturb);
  }
}
void upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double z = xn[0];
  fout[0] = 0.0;
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;

  double z = xn[0];
  fout[0] = Ti0;
}

// Collision frequencies.
void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;
  double q0 = app->q0;
  double kappa = app->kappa;

  double x = xc[0], y = xc[1], z = xc[2];

  double r = x+r0;

  // Map to cylindrical (R, Z, phi) coordinates.
  double R   = R_rtheta(r, z, ctx);
  double Z   = kappa*r*sin(z);
  double phi = -q0/r0*y - alpha(r, z, 0, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);

  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

void bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;

  double x = 0., y = 0., z = xc[2];
  double r = x+r0;
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);
  fout[0] = sqrt(Bt*Bt + Bp*Bp);
}

struct gk_app_ctx
create_ctx(void)
{
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double R_axis    = 0.4;                // Magnetic axis major radius [m].
  double B_axis    = 0.224374548;        // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 0.5948;             // Major radius of the LCFS at the outboard midplane [m].
  double R0        = R_LCFSmid+0.025;    // Major radius of the simulation box [m].
  double a_mid     = R_LCFSmid-R_axis;   // Minor radius at outboard midplane [m].
  double r0        = R0-R_axis;          // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double qSep      = 3.69546081;         // Safety factor at the separatrix.
  double sSep      = 2.27976219;         // Magnetic shear at the separatrix.
  double kappa     = 1.57;               // Elongation (=1 for no elongation).
  double delta     = 0.6;                // Triangularity (=0 for no triangularity).
  double Lz        = 0.62*2.*M_PI;    // Domain size along magnetic field.

  double q0       = qprofile(0.+r0,a_mid,qSep,sSep);    // Magnetic safety factor in the center of domain.
  double epsilon0 = r0/R0;              // Inverse aspect ratio in the center of the domain.

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double mi  = mp;   // Hydrogen ions.
  double Te0 = 178*eV;
  double Ti0 = 70*eV;
  double n0  = 1.78e18;   // [1/m^3]

  double nuFrac = 1.0;
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.

  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Source parameters
  double n_src = 1.95e22;
  double Ti_src = 40*eV;

  // Grid parameters
  int num_cell_z = 64;
  int num_cell_vpar = 16;
  int num_cell_mu = 45;
  int poly_order = 1;

  double vpar_max_ion = 4.*vti;
  double mu_max_ion = mi*pow(1.5*4*vti,2)/(2*B0);

  double final_time = 1e-3;
  int num_frames = 200;

  struct gk_app_ctx ctx = {
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .r0     = r0    ,
    .B0     = B0    ,
    .qSep   = qSep  ,
    .sSep   = sSep  ,
    .kappa  = kappa ,
    .delta  = delta ,
    .q0     = q0    ,
    .Lz     = Lz    ,
    .z_min = -Lz/2.,  .z_max = Lz/2.,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .nuFrac = nuFrac,  .nuIon = nuIon,
  
    .n_src = n_src,  .Ti_src = Ti_src,
  
    .num_cell_z = num_cell_z,
    .num_cell_vpar = num_cell_vpar,
    .num_cell_mu = num_cell_mu,
    .poly_order = poly_order,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,

    .final_time = final_time,  .num_frames = num_frames,
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = {  ctx.num_cell_vpar, ctx.num_cell_mu },
    .polarization_density = ctx.n0,
    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_ion,
      .ctx_temp = &ctx,
      .temp = temp_ion,      
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN, 
        .ctx_density = &ctx,
        .density = density_src,
        .ctx_upar = &ctx,
        .upar= upar_ion_src,
        .ctx_temp = &ctx,
        .temp = temp_ion_src,  
      }, 
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ADIABATIC,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .bmag_fac = ctx.B0,
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_ltx_adiabatic_elc_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { ctx.z_min },
    .upper = { ctx.z_max },
    .cells = { ctx.num_cell_z },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0., 0.},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // magnetic field magnitude
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 1,
    .species = { ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  printf("Creating app object ...\n");
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.final_time;
  double dt = tend-tcurr;
  int nframe = ctx.num_frames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  printf("Applying initial conditions ...\n");
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  printf("Computing initial diagnostics ...\n");
  write_data(&io_trig, app, tcurr);
  printf("Computing initial field energy ...\n");
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  printf("Starting main loop ...\n");
  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
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
  printf(" ... finished\n");
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
  
  return 0;
}