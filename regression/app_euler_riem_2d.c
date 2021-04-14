#include <math.h>
#include <stdio.h>

#include <gkyl_util.h>
#include <gkyl_moment.h>

struct euler_ctx {
    double gas_gamma; // gas constant
};

void
evalEulerInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double sloc = 0.8;
  double rho, u, v, pr;

  double x = xn[0], y = xn[1];

  double upLeft[] = {0.0, 0.3, 0.5323, 1.206, 0.0};
  double upRight[] = {0.0, 1.5, 1.5, 0.0, 0.0};
  double loLeft[] = {0.0, 0.029, 0.138, 1.206, 1.206};
  double loRight[] = {0.0, 0.3, 0.5323, 0.0, 1.206};

  if (y>sloc) {
    if (x<sloc) {
      pr = upLeft[1];
      rho = upLeft[2];
      u = upLeft[3];
      v = upLeft[4];
    }
    else {
      pr = upRight[1];
      rho = upRight[2];
      u = upRight[3];
      v = upRight[4];
    }
  }
  else {
    if (x<sloc) {
        pr = loLeft[1];
        rho = loLeft[2];
        u = loLeft[3];
        v = loLeft[4];
    }
    else {
      pr = loRight[1];
      rho = loRight[2];
      u = loRight[3];
      v = loRight[4];
    }
  }
  fout[0] = rho;
  fout[1] = rho*u; fout[2] = rho*v; fout[3] = 0.0;
  fout[4] = 0.5*rho*(u*u+v*v) + pr/(gas_gamma-1);
}

struct euler_ctx
euler_ctx(void)
{
  return (struct euler_ctx) { .gas_gamma = 1.4 };
}

int
main(int argc, char **argv)
{
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  struct gkyl_moment_species elc = {
    .name = "euler",
    .charge = 0.0, .mass = 1.0,

    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,
  };

  // VM app
  struct gkyl_moment moments = {
    .name = "euler_riem_2d",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { 128, 128 },

    .num_species = 1,
    .species = { elc },
  };

  
  return 0;
}
