#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

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
  gkyl_range_init(conf_local_ext, cdim, lower_ext, upper_ext);
  gkyl_sub_range_init(conf_local, conf_local_ext, lower, upper);
}

struct gkyl_array*
mkarr(long numComp, long size)
{
  return gkyl_array_new(sizeof(double)*numComp, size);
}

struct maxwell_app {
    int ndim, polyOrder;
    int cells[3];
    double L, kwave, lwave;
};

void evalFieldFunc(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct maxwell_app *dat = ctx;

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

struct maxwell_app
parse_args(int argc, char **argv)
{
  // replace with something read from CLI or a config file
  struct maxwell_app app = {
    .ndim = 2,
    .polyOrder = 2,
    .cells = { 16, 16 },
    .L = 1.0,
    .kwave = 2.0,
    .lwave = 2.0
  };

  return app;
}

int
main(int argc, char *argv[])
{
  struct maxwell_app app = parse_args(argc, argv);  
  int cdim = app.ndim;
  int polyOrder = app.polyOrder;

  // grid
  double lower[] = {0.0, 0.0};
  double upper[] = {app.L, app.L};
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, app.cells);

  // basis functions
  struct gkyl_basis confBasis;
  gkyl_cart_modal_serendip(&confBasis, cdim, polyOrder);

  // Maxwell equation object
  struct gkyl_dg_eqn* maxwell_eqn = gkyl_dg_maxwell_new(
    &confBasis, GKYL_SPEED_OF_LIGHT, 0, 0);

  // updaters
  gkyl_proj_on_basis *projField = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    confBasis.polyOrder+1, 8, evalFieldFunc, &app);

  // solver for Maxwell equation
  int max_up_dirs[] = {0, 1}, zero_flux_flags[] = {0, 0};
  gkyl_hyper_dg *maxwellSlvr = gkyl_hyper_dg_new(&confGrid, &confBasis, maxwell_eqn,
    cdim, max_up_dirs, zero_flux_flags, 1);
  
  // ranges
  struct gkyl_range conf_local, conf_local_ext;
  init_conf_ranges(cdim, app.cells, &conf_local, &conf_local_ext);

  // fields
  struct gkyl_array *em = mkarr(confBasis.numBasis*8, conf_local_ext.volume);
  struct gkyl_array *emNew = mkarr(confBasis.numBasis*8, conf_local_ext.volume);

  // run simulation
  gkyl_proj_on_basis_advance(projField, 0.0, &conf_local, em);

  // write fields
  gkyl_grid_array_write(&confGrid, &conf_local, em, "field_0.gkyl");

  // cleanup objects
  gkyl_proj_on_basis_release(projField);
  gkyl_hyper_dg_release(maxwellSlvr);
  gkyl_dg_eqn_release(maxwell_eqn);
  gkyl_array_release(em);
  gkyl_array_release(emNew);
  
  return 0;

}
