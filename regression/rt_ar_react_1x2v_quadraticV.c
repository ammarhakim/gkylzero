#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

/*
This is a test to see if our ionization/recombination operators produce coronal equilibrium
Need at least 1e21 argon density otr distribution of charge states may be slightly different
Need ne*time > 10^20

*/

struct gk_app_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massAr; // Li mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac;
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // argon reference density
  double n0Ar1; // argon1 reference density
  double n0Ar2; // argon2 reference density
  double n0Ar3; // argon3 reference density
  double n0Ar4; // argon4 reference density
  double n0Ar5; // argon5 reference density
  double n0Ar6; // argon6 reference density
  double n0Ar7; // argon7 reference density
  double n0Ar8; // argon8 reference density
  double n0Ar9; // argon9 reference density
  double n0Ar10; // argon10 reference density
  double n0Ar11; // argon11 reference density
  double n0Ar12; // argon12 reference density
  double n0Ar13; // argon13 reference density
  double n0Ar14; // argon14 reference density
  double n0Ar15; // argon15 reference density
  double n0Ar16; // argon16 reference density
  double n0Ar17; // argon17 reference density
  double n0Ar18; // argon18 reference density

  // Simulation parameters
  double Lz; // Box size in z.
  double kperp; // perpendicular wave number used in Poisson solve
  double n_src; // Source density.
  double T_src; // Temp density.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Li ions
  double mu_max_Ar; // Velocity space extents in mu for Li ions  
  double finalTime; // end time
  int numFrames; // number of output frames
};

// Source profiles.
void eval_source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_src;
  double Ls = app->Lz/4.;

  if (fabs(z) < Ls) {
    fout[0] = 1.;
  } else {
    fout[0] = 1.e-40;
  }
  fout[0] = n_src*fout[0];
}
void eval_source_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_src;
  double Ls = app->Lz/4.;

  if (fabs(z) < Ls) {
    fout[0] = 1.;
  } else {
    fout[0] = 1.e-40;
  }
  fout[0] = n_src*fout[0];
}

void eval_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
mapc2p_vel_power(const double* vc, double* GKYL_RESTRICT vp, double vmax, double mumax, double vpow, double mupow)
{
  double cvpar = vc[0];
  double cmu = vc[1];

  if ( cvpar < 0.0 )
    vp[0] = -fabs(vmax * pow(fabs(cvpar), vpow));
  else
    vp[0] = fabs(vmax * pow(cvpar, vpow));
  
  vp[1] = mumax * pow(cmu, mupow);
}

void
mapc2p_vel_elc(double t, const double* vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double vmax = app->vpar_max_elc;
  double mumax = app->mu_max_elc;
  mapc2p_vel_power(vc, vp, vmax, mumax, 2.0, 2.0);
}

void eval_source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double T_src = app->T_src;
  fout[0] = T_src;
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

// Initial conditions.
void eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  double n0Ararr[9] = {0};

  n0Ararr[1] = app->n0Ar1;
  n0Ararr[2] = app->n0Ar2;
      
  // fout[0] = n0 - 10*n0Ar;
  for (int i=1; i<=2; i++)
    n0 -= n0Ararr[i] * i;
  fout[0] = n0;
}

void eval_density_ar0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar;
  fout[0] = n0Ar;
}

void eval_density_ar1(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar1;
  fout[0] = n0Ar;
}

void eval_density_ar2(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar2;
  fout[0] = n0Ar;
}

void eval_density_ar3(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar3;
  fout[0] = n0Ar;
}

void eval_density_ar4(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar4;
  fout[0] = n0Ar;
}

void eval_density_ar5(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar5;
  fout[0] = n0Ar;
}

void eval_density_ar6(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar6;
  fout[0] = n0Ar;
}

void eval_density_ar7(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar7;
  fout[0] = n0Ar;
}

void eval_density_ar8(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar8;
  fout[0] = n0Ar;
}

void eval_density_ar9(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar9;
  fout[0] = n0Ar;
}

void eval_density_ar10(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar10;
  fout[0] = n0Ar;
}

void eval_density_ar11(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar11;
  fout[0] = n0Ar;
}

void eval_density_ar12(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar12;
  fout[0] = n0Ar;
}

void eval_density_ar13(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar13;
  fout[0] = n0Ar;
}

void eval_density_ar14(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar14;
  fout[0] = n0Ar;
}

void eval_density_ar15(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar15;
  fout[0] = n0Ar;
}

void eval_density_ar16(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar17;
  fout[0] = n0Ar;
}

void eval_density_ar17(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar17;
  fout[0] = n0Ar;
}

void eval_density_ar18(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar18;
  fout[0] = n0Ar;
}

void eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Te = app->Te;
  fout[0] = Te;
}

void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti = app->Ti;
  fout[0] = Ti;
}

void eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
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

double plasma_frequency(double n, double m)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  return sqrt(n*eV*eV/m/eps0);
}
double coulomb_log(double ns, double nr, double ms, double mr, double Ts, double Tr, double qs, double qr)
{

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double hbar = GKYL_PLANCKS_CONSTANT_H/2/M_PI;
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  double wps = plasma_frequency(ns,ms);
  double wpr = plasma_frequency(nr,mr);
  double inner1 = wps*wps/(Ts/ms + 3*Ts/ms) + wpr*wpr/(Tr/mr + 3*Ts/ms);
  double u = 3*(vts*vts + vtr*vtr);
  double msr = ms*mr/(ms+mr);
  double inner2 = fmax(fabs(qs*qr)/(4*M_PI*eps0*msr*u*u), hbar/(2*sqrt(eV)*msr*u));
  double inner = (1/inner1)*(1/inner2/inner2) + 1;
  return 0.5*log(inner);
}

double norm_nu_func(double nuFrac, double ns, double nr, double ms, double mr, double qs, double qr, double Ts, double Tr)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double clog = coulomb_log(ns,nr,ms,mr,Ts, Tr, qs, qr);
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  return nuFrac/ms*(1/mr+1/ms)*qs*qs*qr*qr*clog/(6*pow(M_PI,1.5)*eps0*eps0);
}

struct gk_app_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // Proton mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.

  double mi = 2.014*mp; // D ion mass
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Reference density and temperature.
  double Te = 200*eV;
  double Ti = 200*eV;
  double TAr = 200*eV;
  double n0 = 1.0e21;
  //double n0Ar = n0 * 1e-18;
  double n0Ar = n0 * 1e-20;
  double n0Ar1 = n0 * 1e-1;
  double n0Ar2 = n0 * 1e-2;
  double n0Ar3 = n0Ar * 1e5;
  double n0Ar4 = n0Ar * 1e7;
  double n0Ar5 = n0Ar * 1e7;
  double n0Ar6 = n0Ar * 1e10;
  double n0Ar7 = n0Ar * 1e10;
  double n0Ar8 = n0Ar * 1e11;
  double n0Ar9 = n0Ar * 1e11;
  double n0Ar10 = n0Ar * 1e11;
  double n0Ar11 = n0Ar * 1e11;
  double n0Ar12 = n0Ar * 1e11;
  double n0Ar13 = n0Ar * 1e10;
  double n0Ar14 = n0Ar * 1e10;
  double n0Ar15 = n0Ar * 1e10;
  double n0Ar16 = n0Ar * 1e10;
  double n0Ar17 = n0Ar * 1e10;
  double n0Ar18 = n0Ar * 1e10;

  // Geometry and magnetic field.
  double B_axis = 0.5;
  double R0     = 0.85;
  double a0     = 0.15;
  double R      = R0 + a0;
  double B0     = B_axis*(R0/R);

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(Ti/mAr);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq
  
  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lz = 4.;

  // Perpendicular wavenumber in SI units:
  double kperpRhos = 0.3;
  double kperp = kperpRhos / rho_s;

  // Source parameters.
  //double n_src = 2.870523e+21;
  //double n_src = 2.870523e2*n0;
  double n_src = 0.0;
  double T_src = 2.*Te;

  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = (3./2.)*0.5*me*pow(4.0*vtElc,2)/(2.0*B0);

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = (3./2.)*0.5*mi*pow(4.0*vtIon,2)/(2.0*B0);

  double vpar_max_Ar = 6.0*vtAr;
  double mu_max_Ar = (3./2.)*0.5*mAr*pow(4.0*vtAr,2)/(2.0*B0);
  
  double finalTime = 1e-7; 
  double numFrames = 1;

  struct gk_app_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti,
    .TAr = TAr,
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac, 
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar,
    .n0Ar1 = n0Ar1,
    .n0Ar2 = n0Ar2,
    .n0Ar3 = n0Ar3,
    .n0Ar4 = n0Ar4,
    .n0Ar5 = n0Ar5,
    .n0Ar6 = n0Ar6,
    .n0Ar7 = n0Ar7,
    .n0Ar8 = n0Ar8,
    .n0Ar9 = n0Ar9,
    .n0Ar10 = n0Ar10,
    .n0Ar11 = n0Ar11,
    .n0Ar12 = n0Ar12,
    .n0Ar13 = n0Ar13,
    .n0Ar14 = n0Ar14,
    .n0Ar15 = n0Ar15,
    .n0Ar16 = n0Ar16,
    .n0Ar17 = n0Ar17,
    .n0Ar18 = n0Ar18,
    .Lz = Lz, 
    .kperp = kperp, 
    .n_src = n_src,
    .T_src = T_src,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
    .finalTime = finalTime, 
    .numFrames = numFrames,
  };
  return ctx;
}

struct all_of_1species_reactions {
  struct gkyl_gyrokinetic_react_type recombination[GKYL_MAX_SPECIES];
  struct gkyl_gyrokinetic_react_type ionization[GKYL_MAX_SPECIES];
};

struct all_of_1species_reactions*
create_isonuclear_react(double ion_mass, double electron_mass, int max_z, enum gkyl_ion_type ion_id, char* electron_name, char* ion_name)
{
  static struct all_of_1species_reactions reactions[4];  // index = gkyl_react_self_type
  // 0 = elc, 1=ion, 2 = donor, 3 = recvr
  char ion_nm[128];
  char donor_nm[128];  // Same as recvr_nm
  for (int j=0; j<4; j++) {
    for (int i=0; i<max_z; i++) {
      snprintf(ion_nm, sizeof(ion_nm), "%s%d", ion_name, i+1);
      snprintf(donor_nm, sizeof(donor_nm), "%s%d", ion_name, i);
      reactions[j].recombination[i].react_id = GKYL_REACT_RECOMB;
      reactions[j].recombination[i].ion_id = ion_id;
      strcpy(reactions[j].recombination[i].elc_nm, electron_name);
      reactions[j].recombination[i].charge_state = i;
      reactions[j].recombination[i].ion_mass = ion_mass;
      reactions[j].recombination[i].elc_mass = electron_mass;
      strcpy(reactions[j].recombination[i].ion_nm, ion_nm);
      strcpy(reactions[j].recombination[i].recvr_nm, donor_nm);
      reactions[j].recombination[i].type_self = j;
      
      reactions[j].ionization[i].react_id = GKYL_REACT_IZ;
      reactions[j].ionization[i].ion_id = ion_id;
      strcpy(reactions[j].ionization[i].elc_nm, electron_name);
      reactions[j].ionization[i].charge_state = i;
      reactions[j].ionization[i].ion_mass = ion_mass;
      reactions[j].ionization[i].elc_mass = electron_mass;
      strcpy(reactions[j].ionization[i].ion_nm, ion_nm);
      strcpy(reactions[j].ionization[i].donor_nm, donor_nm);
      reactions[j].ionization[i].type_self = j;
    }
  }
  return reactions;
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

  // Reactions for argon isonuclear sequence
  struct all_of_1species_reactions *reactions = create_isonuclear_react(ctx.massAr, ctx.massElc, 18, GKYL_ION_AR, "elc", "Ar");

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 2);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 6);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 4);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0, 0.0},
    .upper = { 1.0, 1.0}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
    },

     .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Te,      
      .self_nu = evalNuElc,
      .ctx = &ctx,
      .num_cross_collisions = 3,
      .collide_with = { "ion", "Ar1", "Ar2" },
      },

    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
      },
    
    /*.bcx = {
          .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
          .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
	  },*/
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .n_ref = ctx.n0,
      .T_ref = ctx.Ti,
      .normNu = true,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "Ar1", "Ar2" },
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar1+ ions
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar1,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_ar1,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0Ar1,
      .T_ref = ctx.TAr,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "ion", "Ar2" },
    },
	
    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },


    /*    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar2+ ions
  struct gkyl_gyrokinetic_species Ar2 = {
    .name = "Ar2",
    .charge = 2*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar2,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar2,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0Ar2,
      .T_ref = ctx.TAr,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "ion", "Ar1" },
    },

        .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
      },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
    },
    */
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };
  
  // field
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_PERIODIC, 
    .kperpSq = pow(ctx.kperp, 2.),
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "ar_react_1x2v_quad_wJacobvel",

    //.cfl_frac = 0.1,

    .cdim = 1, .vdim = 2,
    .lower = { -ctx.Lz/2.0 },
    .upper = { ctx.Lz/2.0 },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // mapping of computational to physical space
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 4,
    .species = { elc, ion, Ar1, Ar2},

    .field = field,



    .use_gpu = app_args.use_gpu,
    .skip_field=true,
  };
  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.finalTime;
  double dt = tend-tcurr;
  int nframe = ctx.numFrames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    //gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    //gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 1 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
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
  
  return 0;
}
