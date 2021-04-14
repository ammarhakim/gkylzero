#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

void
evalElcInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double charge = 1.0;
  double gasGamma = 5./3.;
  double elcCharge = -charge;
  double ionCharge = charge;
  double elcMass = 1/1836.2;
  double ionMass = 1.0;

  double rho, pr;
  if (xn[0] < 0.5) {
    rho = 1.0*elcMass/ionMass;
    pr = 5.0e-5/(gasGamma-1);
  }
  else {
    rho = 0.125*elcMass/ionMass;
    pr = 5.0e-6/(gasGamma-1);
  }
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gasGamma-1);  
}

void
evalIonInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double gasGamma = 5./3.;

  double rho, pr;
  if (xn[0] < 0.5) {
    rho = 1.0;
    pr = 5.0e-5/(gasGamma-1);
  }
  else {
    rho = 0.125;
    pr = 5.0e-6/(gasGamma-1);
  }
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gasGamma-1);    
}

void
evalFieldInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double Bx = 0.75e-2;
  double Bz = xn[0] < 0.5 ? 1.0e-2 : -1.0e-2;

  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx, fout[4] = 0.0; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

int
main(int argc, char **argv)
{
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(5.0/3.0);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(5.0/3.0);

  struct gkyl_moment_species elc = {
    .name = "elc",

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
  };  

  struct gkyl_moment_field maxwell = {
    .epsilon0 = 1.0, .mu0 = 1.0,

    .evolve = 1,
    .init = evalFieldInit,
  };

  printf("field -> %g\n", maxwell.epsilon0);
  
  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_riem",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { 1024  },

    .num_species = 2,
    .species = { elc, ion },

    .field = maxwell,
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.8;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  
  gkyl_moment_app_write(app, tcurr, 0);

  gkyl_moment_app_write(app, tcurr, 1);

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app); 
  
  return 0;
}
