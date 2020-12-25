#include <math.h>
#include <time.h>
#include <stdio.h>

#include <gkyl_const.h>
#include <gkyl_vlasov.h>

struct maxwell_ctx {
    double L, kwave, lwave;
};

const double L = 1.0; // length of domain

void evalFieldFunc(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct maxwell_ctx *dat = ctx;

  double x = xn[0], y = xn[1];
  double kwave = dat->kwave, lwave = dat->lwave, L = dat->L;
  double freq = 2*M_PI/L*sqrt(kwave*kwave+lwave*lwave)*GKYL_SPEED_OF_LIGHT;
  double tperiod = 2*M_PI/freq;
  double phi = 2*M_PI/L*(kwave*x + lwave*y);
  double knorm = sqrt(kwave*kwave + lwave*lwave);
  double kxn = kwave/knorm;
  double kyn = lwave/knorm;
  
  // n = (-1, 1, 1), n_hat = 1/math.sqrt(3)
  double E0 = 1.0/sqrt(3.0);
  double Ex = -E0*cos(phi);
  double Ey = E0*cos(phi);
  double Ez = E0*cos(phi);
  
  double Bx = E0/GKYL_SPEED_OF_LIGHT*cos(phi)*2*M_PI/L*kyn;
  double By = -E0/GKYL_SPEED_OF_LIGHT*cos(phi)*2*M_PI/L*kxn;
  double Bz = E0/GKYL_SPEED_OF_LIGHT*cos(phi)*2*M_PI/L*(-kxn - kyn);
  
  fout[0] = Ex; fout[1] = Ey, fout[2] = Ez;
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct maxwell_ctx
create_ctx(void)
{
  struct maxwell_ctx ctx = {
    .L = L,
    .kwave = 2.0,
    .lwave = 2.0
  };

  return ctx;  
}

int
main(int argc, char *argv[])
{
  struct maxwell_ctx ctx = create_ctx();

  struct gkyl_em_field field = {
    .epsilon0 = GKYL_EPSILON0,
    .mu0 = GKYL_MU0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalFieldFunc
  };

  struct gkyl_vm vm = {
    .name = "maxwell_2d",

    .cdim = 2, .vdim = 0,
    .lower = { 0.0, 0.0 },
    .upper = { L, L },
    .cells = { 16, 16 },
    .poly_order = 2,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1},

    .field = field,
    .num_species = 0, // pure EM simulation
  };

  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // initialize simulation
  gkyl_vlasov_app_init_sim(app);
  gkyl_vlasov_app_write(app, 0.0, 0);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);  
  
  return 0;
}
