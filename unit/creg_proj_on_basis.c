#include <math.h>
#include <time.h>

#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_mom.h>

#define SHOW_TIME(msg, tdiff) printf("%s %g\n", msg, 1.0*(tdiff)/CLOCKS_PER_SEC)

#define epsilon0 1.0
#define mu0 1.0
#define lightSpeed (1.0/sqrt(epsilon0*mu0))
#define elcErrorSpeedFactor 0.0
#define mgnErrorSpeedFactor 0.0

#define massElc 1.0

#define ud 0.3
#define R 0.333333333333333

#define nElc10 0.5
#define nElc20 0.5
#define uxElc10 0.0
#define uyElc10 ud
#define uxElc20 0.0
#define uyElc20 -ud
#define TElc10 (massElc*(R*ud*R*ud))
#define TElc20 (massElc*(R*ud*R*ud))

#define k0 1.0
#define theta (45.0/180.0*M_PI)
#define kx (k0*cos(theta))
#define ky (k0*sin(theta))
#define perturb_n 1e-8
#define alpha 1.18281106421231 // ratio of E_y/E_x 

#define vthElc10 sqrt(TElc10/massElc)
#define vthElc20 sqrt(TElc20/massElc)

inline double
maxwellian2D(double n, double vx, double vy, double ux, double uy, double vth)
{
  double v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy);
  return n/(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void evalDistFunc(double t, const double * restrict xn, double* restrict fout)
{
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];
  double fv = maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20);
  fout[0] = (1.0+perturb_n*cos(kx*x+ky*y))*fv;
}

void evalFieldFunc(double t, const double * restrict xn, double* restrict fout)
{
  double x = xn[0], y = xn[1];
  
  double E_x = -perturb_n*sin(kx*x+ky*y)/(kx+ky*alpha);
  double E_y = alpha*E_x;
  double B_z = kx*E_y-ky*E_x;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

int
main(void)
{
  const int cdim = 2, vdim = 2, pdim = cdim+vdim;
  const int polyOrder = 2;
  clock_t tstart, tend;

  // grid
  double lower[] = {0.0, 0.0, -0.9, -0.9};
  double upper[] = {2*M_PI/kx, 2*M_PI/ky, 0.9, 0.9};
  //int cells[] = {32, 32, 32, 32};
  int cells[] = {4, 4, 8, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, pdim, polyOrder);
  struct gkyl_basis confBasis;
  gkyl_cart_modal_serendip(&confBasis, cdim, polyOrder);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalDistFunc);
  gkyl_proj_on_basis *projField = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    polyOrder+1, 8, evalFieldFunc);

  // moment objects
  struct gkyl_mom_type *m0 = gkyl_vlasov_mom_new(&confBasis, &basis, "M0");
  struct gkyl_mom_type *m1i = gkyl_vlasov_mom_new(&confBasis, &basis, "M1i");
  struct gkyl_mom_type *m2ij = gkyl_vlasov_mom_new(&confBasis, &basis, "M2ij");

  // global index range on phase-space grid (no ghost cells in
  // velocity space)
  int phase_lower[] = { -1, -1, 0, 0 };
  int phase_upper[] = { cells[0], cells[1], cells[2]-1, cells[3]-1 }; 

  // create range for distribution function
  struct gkyl_range dist_range;
  gkyl_range_init(&dist_range, pdim, phase_lower, phase_upper);

  struct gkyl_range em_range;
  gkyl_range_init(&em_range, cdim, phase_lower, phase_upper);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(sizeof(double)*basis.numBasis,
    dist_range.volume);
  // create EM field
  struct gkyl_array *em = gkyl_array_new(sizeof(double)*confBasis.numBasis*8,
    em_range.volume); 

  // create update region
  int local_lower[] = {
    dist_range.lower[0]+1, dist_range.lower[1]+1,
    dist_range.lower[2], dist_range.lower[3]
  };
  int local_upper[] = {
    dist_range.upper[0]-1, dist_range.upper[1]-1,
    dist_range.upper[2], dist_range.upper[3]
  };
  struct gkyl_range dist_local_range;
  gkyl_sub_range_init(&dist_local_range, &dist_range, local_lower, local_upper);

  struct gkyl_range em_local_range;
  gkyl_sub_range_init(&em_local_range, &em_range, local_lower, local_upper);

  tstart = clock();
  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &dist_local_range, distf);
  tend = clock();
  SHOW_TIME("Project on basis took", tend-tstart);

  // project field on basis
  gkyl_proj_on_basis_advance(projField, 0.0, &em_local_range, em);

  tstart = clock();
  // write out to file
  FILE *fp = fopen("distf_0.gkyl", "wb");
  gkyl_rect_grid_write(&grid, fp);
  gkyl_sub_array_write(&dist_local_range, distf, fp);
  fclose(fp);
  tend = clock();
  SHOW_TIME("Output took", tend-tstart);

  // write out to file
  fp = fopen("field_0.gkyl", "wb");
  gkyl_rect_grid_write(&confGrid, fp);
  gkyl_sub_array_write(&em_local_range, em, fp);
  fclose(fp);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projField);
  gkyl_array_release(distf);
  gkyl_array_release(em);
  gkyl_mom_type_release(m0);
  gkyl_mom_type_release(m1i);
  gkyl_mom_type_release(m2ij);

  return 0;
}
