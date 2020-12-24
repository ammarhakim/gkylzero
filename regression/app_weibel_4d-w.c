#include <math.h>
#include <time.h>

#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_mom.h>

#define SHOW_TIME(msg, tdiff) printf("%s %g\n", msg, 1.0*(tdiff)/CLOCKS_PER_SEC)

struct weibel_app {
    int polyOrder, cells[4];
    double lower[4], upper[4];

    // physical parameters
    double massElc;
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
  struct weibel_app *app = ctx;

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

void evalFieldFunc(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct weibel_app *app = ctx;

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

void
init_phase_ranges(int cdim, int pdim, const int *cells,
  struct gkyl_range *phase_local, struct gkyl_range *phase_local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  for (int i=cdim; i<pdim; ++i) {
    lower_ext[i] = 0;
    upper_ext[i] = cells[i]-1;

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(phase_local, pdim, lower, upper);
  gkyl_range_init(phase_local_ext, pdim, lower_ext, upper_ext);
}

void
init_conf_ranges(int cdim, const int *cells,
  struct gkyl_range *conf_local, struct gkyl_range *conf_local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(conf_local, cdim, lower, upper);    
  gkyl_range_init(conf_local_ext, cdim, lower_ext, upper_ext);
}

struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(sizeof(double)*nc, size);
}

struct weibel_app
parse_args(int argc, char **argv)
{
  double ud = 0.3;
  double k0 = 1.0, theta = 45.0/180.0*M_PI;
  double kx = k0*cos(theta), ky = k0*sin(theta);

  double massElc = 1.0, R = 0.333333333333333;
  double TElc10 = massElc*R*ud*R*ud;
  double TElc20 = massElc*R*ud*R*ud;
  double vthElc10 = sqrt(TElc10/massElc);
  double vthElc20 = sqrt(TElc20/massElc);  
  
  struct weibel_app app = {
    .polyOrder = 2,
    
    .cells = { 32, 32, 32, 32 }, 
    .lower = {0.0, 0.0, -0.9, -0.9},
    .upper = {2*M_PI/kx, 2*M_PI/ky, 0.9, 0.9},
    
    .massElc = 1.0,
    
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
  return app;
}

int
main(int argc, char **argv)
{
  struct weibel_app app = parse_args(argc, argv);
  
  int cdim = 2, vdim = 2, pdim = cdim+vdim;
  int polyOrder = app.polyOrder;
  clock_t tstart, tend;

  // grid
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, app.lower, app.upper, app.cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, app.lower, app.upper, app.cells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, pdim, polyOrder);
  gkyl_cart_modal_serendip(&confBasis, cdim, polyOrder);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalDistFunc, &app);
  gkyl_proj_on_basis *projField = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    polyOrder+1, 8, evalFieldFunc, &app);

  // moment objects
  struct gkyl_mom_type *m0t = gkyl_vlasov_mom_new(&confBasis, &basis, "M0");
  struct gkyl_mom_type *m1it = gkyl_vlasov_mom_new(&confBasis, &basis, "M1i");
  struct gkyl_mom_type *m2ijt = gkyl_vlasov_mom_new(&confBasis, &basis, "M2ij");

  // moment updaters
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, m0t);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, m1it);

  struct gkyl_range phase_local, phase_local_ext, conf_local, conf_local_ext;
  init_phase_ranges(cdim, pdim, app.cells, &phase_local, &phase_local_ext);
  init_conf_ranges(cdim, app.cells, &conf_local, &conf_local_ext);

  // create distribution function, fields
  struct gkyl_array *distf = mkarr(basis.numBasis, phase_local_ext.volume);
  struct gkyl_array *em = mkarr(confBasis.numBasis*8, conf_local_ext.volume);

  // fields to store moments
  struct gkyl_array *m0 = mkarr(confBasis.numBasis*m0t->num_mom, conf_local_ext.volume);
  struct gkyl_array *m1i = mkarr(confBasis.numBasis*m1it->num_mom, conf_local_ext.volume);

  tstart = clock();
  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &phase_local, distf);
  tend = clock();  
  SHOW_TIME("Project on basis took", tend-tstart);

  // project field on basis
  gkyl_proj_on_basis_advance(projField, 0.0, &conf_local, em);

  // compute moments
  tstart = clock();
  gkyl_mom_calc_advance(m0calc, &phase_local, &conf_local, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &phase_local, &conf_local, distf, m1i);
  tend = clock();
  SHOW_TIME("Moment calc took", tend-tstart);

  tstart = clock();
  // write out to file
  gkyl_grid_array_write(&grid, &phase_local, distf, "distf_0.gkyl");
  tend = clock();
  SHOW_TIME("Output took", tend-tstart);

  gkyl_grid_array_write(&confGrid, &conf_local, em, "field_0.gkyl");
  gkyl_grid_array_write(&confGrid, &conf_local, m0, "m0_0.gkyl");
  gkyl_grid_array_write(&confGrid, &conf_local, m1i, "m1i_0.gkyl");

  // release resources
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projField);
  
  gkyl_array_release(distf);
  gkyl_array_release(em);
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  
  gkyl_mom_type_release(m0t);
  gkyl_mom_type_release(m1it);
  gkyl_mom_type_release(m2ijt);
  gkyl_mom_calc_release(m0calc);
  gkyl_mom_calc_release(m1icalc);

  return 0;
}
