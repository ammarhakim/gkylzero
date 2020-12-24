#include <math.h>
#include <time.h>
#include <stdio.h>

#define SHOW_TIME(msg, tdiff) printf("%s %g\n", msg, 1.0*(tdiff)/CLOCKS_PER_SEC)

#include <gkyl_vlasov.h>

struct weibel_ctx {
    // parameters for plasma streams
    double nElc10, nElc20;
    double vthElc10, vthElc20;
    double uxElc10, uxElc20;
    double uyElc10, uyElc20;

    // perturbation parameters
    double kx, ky;
    double alpha; // ratio of E_y/E_x
    double perturb_n;
};

inline double
maxwellian2D(double n, double vx, double vy, double ux, double uy, double vth)
{
  double v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy);
  return n/(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double nElc10 = app->nElc10, nElc20 = app->nElc20;
  double uxElc10 = app->uxElc10, uxElc20 = app->uxElc20;
  double uyElc10 = app->uyElc10, uyElc20 = app->uyElc20;
  double vthElc10 = app->vthElc10, vthElc20 = app->vthElc20;
  double kx = app->kx, ky = app->ky, perturb_n = app->perturb_n;  
  
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];
  
  double fv = maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20);
    
  fout[0] = (1.0+app->perturb_n*cos(kx*x+ky*y))*fv;
}

void
evalFieldFunc(double t, const double* restrict xn, double* restrict fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double perturb_n = app->perturb_n, alpha = app->alpha;
  double kx = app->kx, ky = app->ky;
  
  double x = xn[0], y = xn[1];
  
  double E_x = -perturb_n*sin(kx*x+ky*y)/(kx+ky*alpha);
  double E_y = alpha*E_x;
  double B_z = kx*E_y-ky*E_x;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct weibel_ctx
create_ctx(void)
{
  double ud = 0.3;
  double k0 = 1.0, theta = 45.0/180.0*M_PI;
  double kx = k0*cos(theta), ky = k0*sin(theta);

  double massElc = 1.0, R = 0.333333333333333;
  double TElc10 = massElc*R*ud*R*ud;
  double TElc20 = massElc*R*ud*R*ud;
  double vthElc10 = sqrt(TElc10/massElc);
  double vthElc20 = sqrt(TElc20/massElc);  
  
  struct weibel_ctx ctx = {
    .nElc10 = 0.5,
    .nElc20 = 0.5,
    .uxElc10 = 0.0,
    .uxElc20 = 0.0,
    .uyElc10 = ud,
    .uyElc20 = -ud,
    .vthElc10 = vthElc10,
    .vthElc20 = vthElc20,

    .kx = kx, .ky = ky,
    .alpha = 1.18281106421231,
    .perturb_n = 1e-8,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct weibel_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .lower = { -0.9, -0.9 },
    .upper = { 0.9, 0.9 }, 
    .cells = { 8, 8 },

    .evolve = 1,
    .ctx = &ctx,
    .init = evalDistFunc,

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_em_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "weibel_4d",

    .cdim = 2, .vdim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 2*M_PI/ctx.kx, 2*M_PI/ctx.ky },
    .cells = { 16, 16 },
    .poly_order = 2,

    .num_species = 1,
    .species = { elc },
    .field = field
  };

  clock_t tstart, tend;

  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // initialize simulation
  gkyl_vlasov_app_init_sim(app);
  gkyl_vlasov_app_write(app, 0.0, 0);

  gkyl_vlasov_app_calc_mom(app);
  gkyl_vlasov_app_mom_write(app, 0.0, 0);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);
  
  return 0;
}
