#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat.h>
#include <gkyl_null_comm.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_ten_moment_grad_closure_priv.h>
#include <gkyl_moment_non_ideal_priv.h>

// Makes indexing cleaner
static const unsigned T11 = 0;
static const unsigned T12 = 1;
static const unsigned T13 = 2;
static const unsigned T22 = 3;
static const unsigned T23 = 4;
static const unsigned T33 = 5;

// 1D stencil locations (L: lower, U: upper)
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper)
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

// 3D stencil locations (L: lower, U: upper)
enum loc_3d {
  LLL_3D, LLU_3D,
  LUL_3D, LUU_3D,
  ULL_3D, ULU_3D,
  UUL_3D, UUU_3D
};

struct gkyl_ten_moment_grad_closure {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  double k0; // damping coefficient
  double cfl; // CFL number
  double mass;

  struct gkyl_comm *comm;
};

void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

void
inv_setup(double w, double Bx, double By, double Bz, struct gkyl_mat *qinv)
{
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  double Bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
  if (Bmag > 0.0) {
    bx = w*Bx/Bmag;
    by = w*By/Bmag;
    bz = w*Bz/Bmag;
  }

  gkyl_mat_set(qinv, 0, 0, 1.0);
  gkyl_mat_set(qinv, 0, 1, -3.0*bz);
  gkyl_mat_set(qinv, 0, 2, 3.0*by);
  gkyl_mat_set(qinv, 1, 0, bz);
  gkyl_mat_set(qinv, 1, 1, 1.0);
  gkyl_mat_set(qinv, 1, 2, -bx);
  gkyl_mat_set(qinv, 1, 3, -2.0*bz);
  gkyl_mat_set(qinv, 1, 4, 2.0*by);
  gkyl_mat_set(qinv, 2, 0, -by);
  gkyl_mat_set(qinv, 2, 1, bx);
  gkyl_mat_set(qinv, 2, 2, 1.0);
  gkyl_mat_set(qinv, 2, 4, -2.0*bz);
  gkyl_mat_set(qinv, 2, 5, 2.0*by);
  gkyl_mat_set(qinv, 3, 1, 2.0*bz);
  gkyl_mat_set(qinv, 3, 3, 1.0);
  gkyl_mat_set(qinv, 3, 4, -2.0*bx);
  gkyl_mat_set(qinv, 3, 6, -bz);
  gkyl_mat_set(qinv, 3, 7, by);
  gkyl_mat_set(qinv, 4, 1, -by);
  gkyl_mat_set(qinv, 4, 2, bz);
  gkyl_mat_set(qinv, 4, 3, bx);
  gkyl_mat_set(qinv, 4, 4, 1.0);
  gkyl_mat_set(qinv, 4, 5, -bx);
  gkyl_mat_set(qinv, 4, 7, -bz);
  gkyl_mat_set(qinv, 4, 8, by);
  gkyl_mat_set(qinv, 5, 2, -2.0*by);
  gkyl_mat_set(qinv, 5, 4, 2.0*bx);
  gkyl_mat_set(qinv, 5, 5, 1.0);
  gkyl_mat_set(qinv, 5, 8, -bz);
  gkyl_mat_set(qinv, 5, 9, by);
  gkyl_mat_set(qinv, 6, 3, 3.0*bz);
  gkyl_mat_set(qinv, 6, 6, 1.0);
  gkyl_mat_set(qinv, 6, 7, -3.0*bx);
  gkyl_mat_set(qinv, 7, 3, -by);
  gkyl_mat_set(qinv, 7, 4, 2.0*bz);
  gkyl_mat_set(qinv, 7, 6, bx);
  gkyl_mat_set(qinv, 7, 7, 1.0);
  gkyl_mat_set(qinv, 7, 8, -2.0*bx);
  gkyl_mat_set(qinv, 8, 4, -2.0*by);
  gkyl_mat_set(qinv, 8, 5, bz);
  gkyl_mat_set(qinv, 8, 7, 2.0*bx);
  gkyl_mat_set(qinv, 8, 8, 1.0);
  gkyl_mat_set(qinv, 8, 9, -bx);
  gkyl_mat_set(qinv, 9, 5, -3.0*by);
  gkyl_mat_set(qinv, 9, 8, 3.0*bx);
  gkyl_mat_set(qinv, 9, 9, 1.0);
}

void
var_setup(const gkyl_ten_moment_grad_closure *gces,
  int start, int end,
  const double *fluid_d[],
  double rho[], double p[], double Tij[][6])
{
  for (int j = start; j <= end; ++j) {
    rho[j] = fluid_d[j][RHO];
    p[j] = (fluid_d[j][P11] - fluid_d[j][MX] * fluid_d[j][MX] / fluid_d[j][RHO]
          + fluid_d[j][P22] - fluid_d[j][MY] * fluid_d[j][MY] / fluid_d[j][RHO]
          + fluid_d[j][P33] - fluid_d[j][MZ] * fluid_d[j][MZ] / fluid_d[j][RHO])/3.0;
    Tij[j][T11] = (fluid_d[j][P11] - fluid_d[j][MX] * fluid_d[j][MX] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T12] = (fluid_d[j][P12] - fluid_d[j][MX] * fluid_d[j][MY] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T13] = (fluid_d[j][P13] - fluid_d[j][MX] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T22] = (fluid_d[j][P22] - fluid_d[j][MY] * fluid_d[j][MY] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T23] = (fluid_d[j][P23] - fluid_d[j][MY] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T33] = (fluid_d[j][P33] - fluid_d[j][MZ] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
  }
}

double
calc_unmag_heat_flux(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double alpha = 1.0/gces->k0;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double vth_avg = 0.0;
  double n_avg = 0.0;
  double limit = 0.75;
  double cfla = 0.0;
  
  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double dTdx[6] = {0.0};
    double Tij[2][6] = {0.0};
    double rho[2] = {0.0};
    double p[2] = {0.0};
    double q[10] = {0.0};
    var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
    p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);
    vth_avg = sqrt(p_avg/rho_avg);
    n_avg = rho_avg/gces->mass;
    
    for (int k = T11; k <= T33; ++k) {
      dTdx[k] = calc_sym_grad_1D(dx, Tij[L_1D][k], Tij[U_1D][k]);
    }

    q[Q111] = alpha*n_avg*vth_avg*(dTdx[T11] + dTdx[T11] + dTdx[T11])/3.0;
    q[Q112] = alpha*n_avg*vth_avg*(dTdx[T12] + dTdx[T12])/3.0;
    q[Q113] = alpha*n_avg*vth_avg*(dTdx[T13] + dTdx[T13])/3.0;
    q[Q122] = alpha*n_avg*vth_avg*dTdx[T22]/3.0;
    q[Q123] = alpha*n_avg*vth_avg*dTdx[T23]/3.0;
    q[Q133] = alpha*n_avg*vth_avg*dTdx[T33]/3.0;

    int signx[2] = { -1.0, 1.0 };
    for (int j = L_1D; j <= U_1D; ++j) {
      rhs_d[j][P11] += signx[j]*q[Q111]/dx;
      rhs_d[j][P12] += signx[j]*q[Q112]/dx;
      rhs_d[j][P13] += signx[j]*q[Q113]/dx;
      rhs_d[j][P22] += signx[j]*q[Q122]/dx;
      rhs_d[j][P23] += signx[j]*q[Q123]/dx;
      rhs_d[j][P33] += signx[j]*q[Q133]/dx;
    }

    cfla = dt/(dx*dx);
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTx[2][6] = {0.0}; 
    double dTy[2][6] = {0.0};
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    double q[4][10] = {0.0};
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);
    vth_avg = sqrt(p_avg/rho_avg);
    // n_avg = rho_avg/gces->mass;
    n_avg = rho_avg;

    for (int k = T11; k <= T33; ++k) {
      dTx[L_1D][k] = calc_sym_grad_1D(dx, Tij[LL_2D][k], Tij[UL_2D][k]);
      dTx[U_1D][k] = calc_sym_grad_1D(dx, Tij[LU_2D][k], Tij[UU_2D][k]);
      dTdx[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][k], dTx[U_1D][k]);
      dTdx[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][k], dTx[L_1D][k]);

      dTy[L_1D][k] = calc_sym_grad_1D(dy, Tij[LL_2D][k], Tij[LU_2D][k]);
      dTy[U_1D][k] = calc_sym_grad_1D(dy, Tij[UL_2D][k], Tij[UU_2D][k]);
      dTdy[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][k], dTy[U_1D][k]);
      dTdy[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][k], dTy[L_1D][k]);
    }

    int compx[4] = { L_1D, U_1D, L_1D, U_1D };
    int compy[4] = { L_1D, L_1D, U_1D, U_1D };

    int signx[4] = { 1.0, 1.0, -1.0, -1.0 };
    int signy[4] = { 1.0, -1.0, 1.0, -1.0 };

    double chi = alpha*n_avg*vth_avg;

    for (int j = LL_2D; j <= UU_2D; ++j) {
      q[j][Q111] = chi*dTdx[compx[j]][T11];
      q[j][Q112] = chi*(2.0*dTdx[compx[j]][T12] + dTdy[compy[j]][T11])/3.0;
      q[j][Q113] = chi*2.0*dTdx[compx[j]][T13]/3.0;
      q[j][Q122] = chi*(dTdx[compx[j]][T22] + 2.0*dTdy[compy[j]][T12])/3.0;
      q[j][Q123] = chi*(dTdx[compx[j]][T23] + dTdy[compy[j]][T13])/3.0;
      q[j][Q133] = chi*dTdx[compx[j]][T33]/3.0;
      q[j][Q222] = chi*dTdy[compy[j]][T22];
      q[j][Q223] = chi*2.0*dTdy[compy[j]][T23]/3.0;
      q[j][Q233] = chi*dTdy[compy[j]][T33]/3.0;

      rhs_d[j][P11] += signx[j]*q[j][Q111]/(2.0*dx) + signy[j]*q[j][Q112]/(2.0*dy);
      rhs_d[j][P12] += signx[j]*q[j][Q112]/(2.0*dx) + signy[j]*q[j][Q122]/(2.0*dy);
      rhs_d[j][P13] += signx[j]*q[j][Q113]/(2.0*dx) + signy[j]*q[j][Q123]/(2.0*dy);
      rhs_d[j][P22] += signx[j]*q[j][Q122]/(2.0*dx) + signy[j]*q[j][Q222]/(2.0*dy);
      rhs_d[j][P23] += signx[j]*q[j][Q123]/(2.0*dx) + signy[j]*q[j][Q223]/(2.0*dy);
      rhs_d[j][P33] += signx[j]*q[j][Q133]/(2.0*dx) + signy[j]*q[j][Q233]/(2.0*dy);
    }

    double dmin = fmin(dx, dy);
    cfla = dt/(dmin*dmin);
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];

    double dTx[4][6] = {0.0}; 
    double dTy[4][6] = {0.0};
    double dTz[4][6] = {0.0};
    double dTdx[4][6] = {0.0};
    double dTdy[4][6] = {0.0};
    double dTdz[4][6] = {0.0};
    double Tij[8][6] = {0.0};
    double rho[8] = {0.0};
    double p[8] = {0.0};
    double q[8][10] = {0.0};
    var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                   rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
    p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                                 p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);
    vth_avg = sqrt(p_avg/rho_avg);
    // n_avg = rho_avg/gces->mass;
    n_avg = rho_avg;

    for (int k = T11; k <= T33; ++k) {
      dTx[LL_2D][k] = calc_sym_grad_1D(dx, Tij[LLL_3D][k], Tij[ULL_3D][k]);
      dTx[LU_2D][k] = calc_sym_grad_1D(dx, Tij[LUL_3D][k], Tij[UUL_3D][k]);
      dTx[UL_2D][k] = calc_sym_grad_1D(dx, Tij[LLU_3D][k], Tij[ULU_3D][k]);
      dTx[UU_2D][k] = calc_sym_grad_1D(dx, Tij[LUU_3D][k], Tij[UUU_3D][k]);
      dTdx[LL_2D][k] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][k], dTx[LU_2D][k], dTx[UL_2D][k], dTx[UU_2D][k]);
      dTdx[LU_2D][k] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][k], dTx[LL_2D][k], dTx[UL_2D][k], dTx[UU_2D][k]);
      dTdx[UL_2D][k] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][k], dTx[LL_2D][k], dTx[LU_2D][k], dTx[UU_2D][k]);
      dTdx[UU_2D][k] = calc_sym_grad_limiter_3D(limit, dTx[UU_2D][k], dTx[LL_2D][k], dTx[LU_2D][k], dTx[UL_2D][k]);

      dTy[LL_2D][k] = calc_sym_grad_1D(dy, Tij[LLL_3D][k], Tij[LUL_3D][k]);
      dTy[LU_2D][k] = calc_sym_grad_1D(dy, Tij[ULL_3D][k], Tij[UUL_3D][k]);
      dTy[UL_2D][k] = calc_sym_grad_1D(dy, Tij[LLU_3D][k], Tij[LUU_3D][k]);
      dTy[UU_2D][k] = calc_sym_grad_1D(dy, Tij[ULU_3D][k], Tij[UUU_3D][k]);
      dTdy[LL_2D][k] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][k], dTy[LU_2D][k], dTy[UL_2D][k], dTy[UU_2D][k]);
      dTdy[LU_2D][k] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][k], dTy[LL_2D][k], dTy[UL_2D][k], dTy[UU_2D][k]);
      dTdy[UL_2D][k] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][k], dTy[LL_2D][k], dTy[LU_2D][k], dTy[UU_2D][k]);
      dTdy[UU_2D][k] = calc_sym_grad_limiter_3D(limit, dTy[UU_2D][k], dTy[LL_2D][k], dTy[LU_2D][k], dTy[UL_2D][k]);

      dTz[LL_2D][k] = calc_sym_grad_1D(dz, Tij[LLL_3D][k], Tij[LLU_3D][k]);
      dTz[LU_2D][k] = calc_sym_grad_1D(dz, Tij[ULL_3D][k], Tij[ULU_3D][k]);
      dTz[UL_2D][k] = calc_sym_grad_1D(dz, Tij[LUL_3D][k], Tij[LUU_3D][k]);
      dTz[UU_2D][k] = calc_sym_grad_1D(dz, Tij[UUL_3D][k], Tij[UUU_3D][k]);
      dTdz[LL_2D][k] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][k], dTz[LU_2D][k], dTz[UL_2D][k], dTz[UU_2D][k]);
      dTdz[LU_2D][k] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][k], dTz[LL_2D][k], dTz[UL_2D][k], dTz[UU_2D][k]);
      dTdz[UL_2D][k] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][k], dTz[LL_2D][k], dTz[LU_2D][k], dTz[UU_2D][k]);
      dTdz[UU_2D][k] = calc_sym_grad_limiter_3D(limit, dTz[UU_2D][k], dTz[LL_2D][k], dTz[LU_2D][k], dTz[UL_2D][k]);
    }

    int compx[8] = { LL_2D, UL_2D, LU_2D, UU_2D, LL_2D, UL_2D, LU_2D, UU_2D };
    int compy[8] = { LL_2D, UL_2D, LL_2D, UL_2D, LU_2D, UU_2D, LU_2D, UU_2D };
    int compz[8] = { LL_2D, LL_2D, UL_2D, UL_2D, LU_2D, LU_2D, UU_2D, UU_2D };

    int signx[8] = { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 };
    int signy[8] = { 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0 };
    int signz[8] = { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };

    double chi = alpha*n_avg*vth_avg;

    for (int j = LLL_3D; j <= UUU_3D; ++j) {
      q[j][Q111] = chi*(dTdx[compx[j]][T11] + dTdx[compx[j]][T11] + dTdx[compx[j]][T11])/3.0;
      q[j][Q112] = chi*(dTdx[compx[j]][T12] + dTdx[compx[j]][T12] + dTdy[compy[j]][T11])/3.0;
      q[j][Q113] = chi*(dTdx[compx[j]][T13] + dTdx[compx[j]][T13] + dTdz[compz[j]][T11])/3.0;
      q[j][Q122] = chi*(dTdx[compx[j]][T22] + dTdy[compy[j]][T12] + dTdy[compy[j]][T12])/3.0;
      q[j][Q123] = chi*(dTdx[compx[j]][T23] + dTdy[compy[j]][T13] + dTdz[compz[j]][T12])/3.0;
      q[j][Q133] = chi*(dTdx[compx[j]][T33] + dTdz[compz[j]][T13] + dTdz[compz[j]][T13])/3.0;
      q[j][Q222] = chi*(dTdy[compy[j]][T22] + dTdy[compy[j]][T22] + dTdy[compy[j]][T22])/3.0;
      q[j][Q223] = chi*(dTdy[compy[j]][T23] + dTdy[compy[j]][T23] + dTdz[compz[j]][T22])/3.0;
      q[j][Q233] = chi*(dTdy[compy[j]][T33] + dTdz[compz[j]][T23] + dTdz[compz[j]][T23])/3.0;
      q[j][Q333] = chi*(dTdz[compz[j]][T33] + dTdz[compz[j]][T33] + dTdz[compz[j]][T33])/3.0;

      rhs_d[j][P11] += signx[j]*q[j][Q111]/(4.0*dx) + signy[j]*q[j][Q112]/(4.0*dy) + signz[j]*q[j][Q113]/(4.0*dz);
      rhs_d[j][P12] += signx[j]*q[j][Q112]/(4.0*dx) + signy[j]*q[j][Q122]/(4.0*dy) + signz[j]*q[j][Q123]/(4.0*dz);
      rhs_d[j][P13] += signx[j]*q[j][Q113]/(4.0*dx) + signy[j]*q[j][Q123]/(4.0*dy) + signz[j]*q[j][Q133]/(4.0*dz);
      rhs_d[j][P22] += signx[j]*q[j][Q122]/(4.0*dx) + signy[j]*q[j][Q222]/(4.0*dy) + signz[j]*q[j][Q223]/(4.0*dz);
      rhs_d[j][P23] += signx[j]*q[j][Q123]/(4.0*dx) + signy[j]*q[j][Q223]/(4.0*dy) + signz[j]*q[j][Q233]/(4.0*dz);
      rhs_d[j][P33] += signx[j]*q[j][Q133]/(4.0*dx) + signy[j]*q[j][Q233]/(4.0*dy) + signz[j]*q[j][Q333]/(4.0*dz);
    }
    
    double dmin = fmin(fmin(dx, dy), dz);
    cfla = dt/(dmin*dmin);
  }

  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_mag_heat_flux_linsolve(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double alpha = 1.0/gces->k0;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double vth_avg = 0.0;
  double n_avg = 0.0;
  double limit = 0.75;
  double cfla = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  double w = 10.0;
  if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTx[2][6] = {0.0}; 
    double dTy[2][6] = {0.0};
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    struct gkyl_mat *q = gkyl_mat_new(10, 4, 0.0);
    struct gkyl_mat *q_inv = gkyl_mat_new(10, 10, 0.0);
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);
    vth_avg = sqrt(p_avg/rho_avg);
    n_avg = rho_avg;

    for (int k = T11; k <= T33; ++k) {
      dTx[L_1D][k] = calc_sym_grad_1D(dx, Tij[LL_2D][k], Tij[UL_2D][k]);
      dTx[U_1D][k] = calc_sym_grad_1D(dx, Tij[LU_2D][k], Tij[UU_2D][k]);
      dTdx[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][k], dTx[U_1D][k]);
      dTdx[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][k], dTx[L_1D][k]);

      dTy[L_1D][k] = calc_sym_grad_1D(dy, Tij[LL_2D][k], Tij[LU_2D][k]);
      dTy[U_1D][k] = calc_sym_grad_1D(dy, Tij[UL_2D][k], Tij[UU_2D][k]);
      dTdy[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][k], dTy[U_1D][k]);
      dTdy[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][k], dTy[L_1D][k]);
    }
    
    int compx[4] = { L_1D, U_1D, L_1D, U_1D };
    int compy[4] = { L_1D, L_1D, U_1D, U_1D };

    int signx[4] = { 1.0, 1.0, -1.0, -1.0 };
    int signy[4] = { 1.0, -1.0, 1.0, -1.0 };

    for (int j = LL_2D; j <= UU_2D; ++j) {
      Bx += em_tot_d[j][BX]/4.0;
      By += em_tot_d[j][BY]/4.0;
      Bz += em_tot_d[j][BZ]/4.0;
    }

    inv_setup(w, Bx, By, Bz, q_inv);

    double chi = alpha*n_avg*vth_avg;

    gkyl_mat_set(q, Q111, LL_2D, chi*dTdx[compx[LL_2D]][T11]);
    gkyl_mat_set(q, Q112, LL_2D, chi*(2.0*dTdx[compx[LL_2D]][T12]
      + dTdy[compy[LL_2D]][T11])/3.0);
    gkyl_mat_set(q, Q113, LL_2D, chi*2.0*dTdx[compx[LL_2D]][T13]/3.0);
    gkyl_mat_set(q, Q122, LL_2D, chi*(dTdx[compx[LL_2D]][T22]
      + 2.0*dTdy[compy[LL_2D]][T12])/3.0);
    gkyl_mat_set(q, Q123, LL_2D, chi*(dTdx[compx[LL_2D]][T23]
      + dTdy[compy[LL_2D]][T13])/3.0);
    gkyl_mat_set(q, Q133, LL_2D, chi*dTdx[compx[LL_2D]][T33]/3.0);
    gkyl_mat_set(q, Q222, LL_2D, chi*dTdy[compy[LL_2D]][T22]);
    gkyl_mat_set(q, Q223, LL_2D, chi*2.0*dTdy[compy[LL_2D]][T23]/3.0);
    gkyl_mat_set(q, Q233, LL_2D, chi*dTdy[compy[LL_2D]][T33]/3.0);
    gkyl_mat_set(q, Q333, LL_2D, 0.0);

    gkyl_mat_set(q, Q111, LU_2D, chi*dTdx[compx[LU_2D]][T11]);
    gkyl_mat_set(q, Q112, LU_2D, chi*(2.0*dTdx[compx[LU_2D]][T12]
      + dTdy[compy[LU_2D]][T11])/3.0);
    gkyl_mat_set(q, Q113, LU_2D, chi*2.0*dTdx[compx[LU_2D]][T13]/3.0);
    gkyl_mat_set(q, Q122, LU_2D, chi*(dTdx[compx[LU_2D]][T22]
      + 2.0*dTdy[compy[LU_2D]][T12])/3.0);
    gkyl_mat_set(q, Q123, LU_2D, chi*(dTdx[compx[LU_2D]][T23]
      + dTdy[compy[LU_2D]][T13])/3.0);
    gkyl_mat_set(q, Q133, LU_2D, chi*dTdx[compx[LU_2D]][T33]/3.0);
    gkyl_mat_set(q, Q222, LU_2D, chi*dTdy[compy[LU_2D]][T22]);
    gkyl_mat_set(q, Q223, LU_2D, chi*2.0*dTdy[compy[LU_2D]][T23]/3.0);
    gkyl_mat_set(q, Q233, LU_2D, chi*dTdy[compy[LU_2D]][T33]/3.0);
    gkyl_mat_set(q, Q333, LU_2D, 0.0);

    gkyl_mat_set(q, Q111, UL_2D, chi*dTdx[compx[UL_2D]][T11]);
    gkyl_mat_set(q, Q112, UL_2D, chi*(2.0*dTdx[compx[UL_2D]][T12]
      + dTdy[compy[UL_2D]][T11])/3.0);
    gkyl_mat_set(q, Q113, UL_2D, chi*2.0*dTdx[compx[UL_2D]][T13]/3.0);
    gkyl_mat_set(q, Q122, UL_2D, chi*(dTdx[compx[UL_2D]][T22]
      + 2.0*dTdy[compy[UL_2D]][T12])/3.0);
    gkyl_mat_set(q, Q123, UL_2D, chi*(dTdx[compx[UL_2D]][T23]
      + dTdy[compy[UL_2D]][T13])/3.0);
    gkyl_mat_set(q, Q133, UL_2D, chi*dTdx[compx[UL_2D]][T33]/3.0);
    gkyl_mat_set(q, Q222, UL_2D, chi*dTdy[compy[UL_2D]][T22]);
    gkyl_mat_set(q, Q223, UL_2D, chi*2.0*dTdy[compy[UL_2D]][T23]/3.0);
    gkyl_mat_set(q, Q233, UL_2D, chi*dTdy[compy[UL_2D]][T33]/3.0);
    gkyl_mat_set(q, Q333, UL_2D, 0.0);

    gkyl_mat_set(q, Q111, UU_2D, chi*dTdx[compx[UU_2D]][T11]);
    gkyl_mat_set(q, Q112, UU_2D, chi*(2.0*dTdx[compx[UU_2D]][T12]
      + dTdy[compy[UU_2D]][T11])/3.0);
    gkyl_mat_set(q, Q113, UU_2D, chi*2.0*dTdx[compx[UU_2D]][T13]/3.0);
    gkyl_mat_set(q, Q122, UU_2D, chi*(dTdx[compx[UU_2D]][T22]
      + 2.0*dTdy[compy[UU_2D]][T12])/3.0);
    gkyl_mat_set(q, Q123, UU_2D, chi*(dTdx[compx[UU_2D]][T23]
      + dTdy[compy[UU_2D]][T13])/3.0);
    gkyl_mat_set(q, Q133, UU_2D, chi*dTdx[compx[UU_2D]][T33]/3.0);
    gkyl_mat_set(q, Q222, UU_2D, chi*dTdy[compy[UU_2D]][T22]);
    gkyl_mat_set(q, Q223, UU_2D, chi*2.0*dTdy[compy[UU_2D]][T23]/3.0);
    gkyl_mat_set(q, Q233, UU_2D, chi*dTdy[compy[UU_2D]][T33]/3.0);
    gkyl_mat_set(q, Q333, UU_2D, 0.0);

    gkyl_mem_buff ipiv = gkyl_mem_buff_new(sizeof(long[10]));
    gkyl_mat_linsolve_lu(q_inv, q, gkyl_mem_buff_data(ipiv));
      
    rhs_d[LL_2D][P11] += signx[LL_2D]*gkyl_mat_get(q, Q111, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q112, LL_2D)/(2.0*dy);
    rhs_d[LL_2D][P12] += signx[LL_2D]*gkyl_mat_get(q, Q112, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q122, LL_2D)/(2.0*dy);
    rhs_d[LL_2D][P13] += signx[LL_2D]*gkyl_mat_get(q, Q113, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q123, LL_2D)/(2.0*dy);
    rhs_d[LL_2D][P22] += signx[LL_2D]*gkyl_mat_get(q, Q122, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q222, LL_2D)/(2.0*dy);
    rhs_d[LL_2D][P23] += signx[LL_2D]*gkyl_mat_get(q, Q123, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q223, LL_2D)/(2.0*dy);
    rhs_d[LL_2D][P33] += signx[LL_2D]*gkyl_mat_get(q, Q133, LL_2D)/(2.0*dx)
      + signy[LL_2D]*gkyl_mat_get(q, Q233, LL_2D)/(2.0*dy);

    rhs_d[LU_2D][P11] += signx[LU_2D]*gkyl_mat_get(q, Q111, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q112, LU_2D)/(2.0*dy);
    rhs_d[LU_2D][P12] += signx[LU_2D]*gkyl_mat_get(q, Q112, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q122, LU_2D)/(2.0*dy);
    rhs_d[LU_2D][P13] += signx[LU_2D]*gkyl_mat_get(q, Q113, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q123, LU_2D)/(2.0*dy);
    rhs_d[LU_2D][P22] += signx[LU_2D]*gkyl_mat_get(q, Q122, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q222, LU_2D)/(2.0*dy);
    rhs_d[LU_2D][P23] += signx[LU_2D]*gkyl_mat_get(q, Q123, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q223, LU_2D)/(2.0*dy);
    rhs_d[LU_2D][P33] += signx[LU_2D]*gkyl_mat_get(q, Q133, LU_2D)/(2.0*dx)
      + signy[LU_2D]*gkyl_mat_get(q, Q233, LU_2D)/(2.0*dy);

    rhs_d[UL_2D][P11] += signx[UL_2D]*gkyl_mat_get(q, Q111, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q112, UL_2D)/(2.0*dy);
    rhs_d[UL_2D][P12] += signx[UL_2D]*gkyl_mat_get(q, Q112, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q122, UL_2D)/(2.0*dy);
    rhs_d[UL_2D][P13] += signx[UL_2D]*gkyl_mat_get(q, Q113, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q123, UL_2D)/(2.0*dy);
    rhs_d[UL_2D][P22] += signx[UL_2D]*gkyl_mat_get(q, Q122, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q222, UL_2D)/(2.0*dy);
    rhs_d[UL_2D][P23] += signx[UL_2D]*gkyl_mat_get(q, Q123, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q223, UL_2D)/(2.0*dy);
    rhs_d[UL_2D][P33] += signx[UL_2D]*gkyl_mat_get(q, Q133, UL_2D)/(2.0*dx)
      + signy[UL_2D]*gkyl_mat_get(q, Q233, UL_2D)/(2.0*dy);

    rhs_d[UU_2D][P11] += signx[UU_2D]*gkyl_mat_get(q, Q111, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q112, UU_2D)/(2.0*dy);
    rhs_d[UU_2D][P12] += signx[UU_2D]*gkyl_mat_get(q, Q112, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q122, UU_2D)/(2.0*dy);
    rhs_d[UU_2D][P13] += signx[UU_2D]*gkyl_mat_get(q, Q113, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q123, UU_2D)/(2.0*dy);
    rhs_d[UU_2D][P22] += signx[UU_2D]*gkyl_mat_get(q, Q122, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q222, UU_2D)/(2.0*dy);
    rhs_d[UU_2D][P23] += signx[UU_2D]*gkyl_mat_get(q, Q123, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q223, UU_2D)/(2.0*dy);
    rhs_d[UU_2D][P33] += signx[UU_2D]*gkyl_mat_get(q, Q133, UU_2D)/(2.0*dx)
      + signy[UU_2D]*gkyl_mat_get(q, Q233, UU_2D)/(2.0*dy);
   
    double dmin = fmin(dx, dy);
    cfla = dt/(dmin*dmin);
    gkyl_mem_buff_release(ipiv);
    gkyl_mat_release(q);
    gkyl_mat_release(q_inv);
  }
    
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_mag_heat_flux_kernel(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double alpha = 1.0/gces->k0;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double vth_avg = 0.0;
  double n_avg = 0.0;
  double limit = 0.75;
  double cfla = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  double w = 10.0;
  if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTx[2][6] = {0.0}; 
    double dTy[2][6] = {0.0};
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    double q[4][10] = {0.0};
    double q_src[4][10] = {0.0};
    double q_out[4][10] = {0.0};
    struct gkyl_mat *qsrc = gkyl_mat_new(10, 4, 0.0);
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);
    vth_avg = sqrt(p_avg/rho_avg);
    // n_avg = rho_avg/gces->mass;
    n_avg = rho_avg;

    for (int k = T11; k <= T33; ++k) {
      dTx[L_1D][k] = calc_sym_grad_1D(dx, Tij[LL_2D][k], Tij[UL_2D][k]);
      dTx[U_1D][k] = calc_sym_grad_1D(dx, Tij[LU_2D][k], Tij[UU_2D][k]);
      dTdx[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][k], dTx[U_1D][k]);
      dTdx[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][k], dTx[L_1D][k]);

      dTy[L_1D][k] = calc_sym_grad_1D(dy, Tij[LL_2D][k], Tij[LU_2D][k]);
      dTy[U_1D][k] = calc_sym_grad_1D(dy, Tij[UL_2D][k], Tij[UU_2D][k]);
      dTdy[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][k], dTy[U_1D][k]);
      dTdy[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][k], dTy[L_1D][k]);
    }
    
    int compx[4] = { L_1D, U_1D, L_1D, U_1D };
    int compy[4] = { L_1D, L_1D, U_1D, U_1D };

    int signx[4] = { 1.0, 1.0, -1.0, -1.0 };
    int signy[4] = { 1.0, -1.0, 1.0, -1.0 };

    for (int j = LL_2D; j <= UU_2D; ++j) {
      Bx += em_tot_d[j][BX]/4.0;
      By += em_tot_d[j][BY]/4.0;
      Bz += em_tot_d[j][BZ]/4.0;
    }

    double chi = alpha*n_avg*vth_avg;

    q[LL_2D][Q111] = chi*dTdx[compx[LL_2D]][T11];
    q[LL_2D][Q112] = chi*(2.0*dTdx[compx[LL_2D]][T12]
      + dTdy[compy[LL_2D]][T11])/3.0;
    q[LL_2D][Q113] = chi*2.0*dTdx[compx[LL_2D]][T13]/3.0;
    q[LL_2D][Q122] = chi*(dTdx[compx[LL_2D]][T22]
      + 2.0*dTdy[compy[LL_2D]][T12])/3.0;
    q[LL_2D][Q123] = chi*(dTdx[compx[LL_2D]][T23] + dTdy[compy[LL_2D]][T13])/3.0;
    q[LL_2D][Q133] = chi*dTdx[compx[LL_2D]][T33]/3.0;
    q[LL_2D][Q222] = chi*dTdy[compy[LL_2D]][T22];
    q[LL_2D][Q223] = chi*2.0*dTdy[compy[LL_2D]][T23]/3.0;
    q[LL_2D][Q233] = chi*dTdy[compy[LL_2D]][T33]/3.0;
    q[LL_2D][Q333] = 0.0;
    
    q[LU_2D][Q111] = chi*dTdx[compx[LU_2D]][T11];
    q[LU_2D][Q112] = chi*(2.0*dTdx[compx[LU_2D]][T12]
      + dTdy[compy[LU_2D]][T11])/3.0;
    q[LU_2D][Q113] = chi*2.0*dTdx[compx[LU_2D]][T13]/3.0;
    q[LU_2D][Q122] = chi*(dTdx[compx[LU_2D]][T22]
      + 2.0*dTdy[compy[LU_2D]][T12])/3.0;
    q[LU_2D][Q123] = chi*(dTdx[compx[LU_2D]][T23] + dTdy[compy[LU_2D]][T13])/3.0;
    q[LU_2D][Q133] = chi*dTdx[compx[LU_2D]][T33]/3.0;
    q[LU_2D][Q222] = chi*dTdy[compy[LU_2D]][T22];
    q[LU_2D][Q223] = chi*2.0*dTdy[compy[LU_2D]][T23]/3.0;
    q[LU_2D][Q233] = chi*dTdy[compy[LU_2D]][T33]/3.0;
    q[LU_2D][Q333] = 0.0;

    q[UL_2D][Q111] = chi*dTdx[compx[UL_2D]][T11];
    q[UL_2D][Q112] = chi*(2.0*dTdx[compx[UL_2D]][T12]
      + dTdy[compy[UL_2D]][T11])/3.0;
    q[UL_2D][Q113] = chi*2.0*dTdx[compx[UL_2D]][T13]/3.0;
    q[UL_2D][Q122] = chi*(dTdx[compx[UL_2D]][T22]
      + 2.0*dTdy[compy[UL_2D]][T12])/3.0;
    q[UL_2D][Q123] = chi*(dTdx[compx[UL_2D]][T23] + dTdy[compy[UL_2D]][T13])/3.0;
    q[UL_2D][Q133] = chi*dTdx[compx[UL_2D]][T33]/3.0;
    q[UL_2D][Q222] = chi*dTdy[compy[UL_2D]][T22];
    q[UL_2D][Q223] = chi*2.0*dTdy[compy[UL_2D]][T23]/3.0;
    q[UL_2D][Q233] = chi*dTdy[compy[UL_2D]][T33]/3.0;
    q[UL_2D][Q333] = 0.0;

    q[UU_2D][Q111] = chi*dTdx[compx[UU_2D]][T11];
    q[UU_2D][Q112] = chi*(2.0*dTdx[compx[UU_2D]][T12]
      + dTdy[compy[UU_2D]][T11])/3.0;
    q[UU_2D][Q113] = chi*2.0*dTdx[compx[UU_2D]][T13]/3.0;
    q[UU_2D][Q122] = chi*(dTdx[compx[UU_2D]][T22]
      + 2.0*dTdy[compy[UU_2D]][T12])/3.0;
    q[UU_2D][Q123] = chi*(dTdx[compx[UU_2D]][T23] + dTdy[compy[UU_2D]][T13])/3.0;
    q[UU_2D][Q133] = chi*dTdx[compx[UU_2D]][T33]/3.0;
    q[UU_2D][Q222] = chi*dTdy[compy[UU_2D]][T22];
    q[UU_2D][Q223] = chi*2.0*dTdy[compy[UU_2D]][T23]/3.0;
    q[UU_2D][Q233] = chi*dTdy[compy[UU_2D]][T33]/3.0;
    q[UU_2D][Q333] = 0.0;

    anisotropic_diffusion_kernel(w, Bx, By, Bz, q[LL_2D], q_out[LL_2D]);
    anisotropic_diffusion_kernel(w, Bx, By, Bz, q[LU_2D], q_out[LU_2D]);
    anisotropic_diffusion_kernel(w, Bx, By, Bz, q[UL_2D], q_out[UL_2D]);
    anisotropic_diffusion_kernel(w, Bx, By, Bz, q[UU_2D], q_out[UU_2D]);
    
    rhs_d[LL_2D][P11] += signx[LL_2D]*q_out[LL_2D][Q111]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q112]/(2.0*dy);
    rhs_d[LL_2D][P12] += signx[LL_2D]*q_out[LL_2D][Q112]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q122]/(2.0*dy);
    rhs_d[LL_2D][P13] += signx[LL_2D]*q_out[LL_2D][Q113]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q123]/(2.0*dy);
    rhs_d[LL_2D][P22] += signx[LL_2D]*q_out[LL_2D][Q122]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q222]/(2.0*dy);
    rhs_d[LL_2D][P23] += signx[LL_2D]*q_out[LL_2D][Q123]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q223]/(2.0*dy);
    rhs_d[LL_2D][P33] += signx[LL_2D]*q_out[LL_2D][Q133]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q233]/(2.0*dy);

    rhs_d[LU_2D][P11] += signx[LU_2D]*q_out[LU_2D][Q111]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q112]/(2.0*dy);
    rhs_d[LU_2D][P12] += signx[LU_2D]*q_out[LU_2D][Q112]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q122]/(2.0*dy);
    rhs_d[LU_2D][P13] += signx[LU_2D]*q_out[LU_2D][Q113]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q123]/(2.0*dy);
    rhs_d[LU_2D][P22] += signx[LU_2D]*q_out[LU_2D][Q122]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q222]/(2.0*dy);
    rhs_d[LU_2D][P23] += signx[LU_2D]*q_out[LU_2D][Q123]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q223]/(2.0*dy);
    rhs_d[LU_2D][P33] += signx[LU_2D]*q_out[LU_2D][Q133]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q233]/(2.0*dy);

    rhs_d[UL_2D][P11] += signx[UL_2D]*q_out[UL_2D][Q111]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q112]/(2.0*dy);
    rhs_d[UL_2D][P12] += signx[UL_2D]*q_out[UL_2D][Q112]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q122]/(2.0*dy);
    rhs_d[UL_2D][P13] += signx[UL_2D]*q_out[UL_2D][Q113]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q123]/(2.0*dy);
    rhs_d[UL_2D][P22] += signx[UL_2D]*q_out[UL_2D][Q122]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q222]/(2.0*dy);
    rhs_d[UL_2D][P23] += signx[UL_2D]*q_out[UL_2D][Q123]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q223]/(2.0*dy);
    rhs_d[UL_2D][P33] += signx[UL_2D]*q_out[UL_2D][Q133]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q233]/(2.0*dy);

    rhs_d[UU_2D][P11] += signx[UU_2D]*q_out[UU_2D][Q111]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q112]/(2.0*dy);
    rhs_d[UU_2D][P12] += signx[UU_2D]*q_out[UU_2D][Q112]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q122]/(2.0*dy);
    rhs_d[UU_2D][P13] += signx[UU_2D]*q_out[UU_2D][Q113]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q123]/(2.0*dy);
    rhs_d[UU_2D][P22] += signx[UU_2D]*q_out[UU_2D][Q122]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q222]/(2.0*dy);
    rhs_d[UU_2D][P23] += signx[UU_2D]*q_out[UU_2D][Q123]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q223]/(2.0*dy);
    rhs_d[UU_2D][P33] += signx[UU_2D]*q_out[UU_2D][Q133]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q233]/(2.0*dy);
   
    double dmin = fmin(dx, dy);
    cfla = dt/(dmin*dmin);
  }
}

double
calc_mag_heat_flux_rotate(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double alpha = 1.0/gces->k0;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double vth_avg = 0.0;
  double n_avg = 0.0;
  double limit = 0.75;
  double cfla = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  double w = 10.0;
  if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTx[2][6] = {0.0}; 
    double dTy[2][6] = {0.0};
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    double q[4][10] = {0.0};
    double q_src[4][10] = {0.0};
    double q_out[4][10] = {0.0};
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);
    vth_avg = sqrt(p_avg/rho_avg);
    // n_avg = rho_avg/gces->mass;
    n_avg = rho_avg;

    for (int k = T11; k <= T33; ++k) {
      dTx[L_1D][k] = calc_sym_grad_1D(dx, Tij[LL_2D][k], Tij[UL_2D][k]);
      dTx[U_1D][k] = calc_sym_grad_1D(dx, Tij[LU_2D][k], Tij[UU_2D][k]);
      dTdx[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][k], dTx[U_1D][k]);
      dTdx[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][k], dTx[L_1D][k]);

      dTy[L_1D][k] = calc_sym_grad_1D(dy, Tij[LL_2D][k], Tij[LU_2D][k]);
      dTy[U_1D][k] = calc_sym_grad_1D(dy, Tij[UL_2D][k], Tij[UU_2D][k]);
      dTdy[L_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][k], dTy[U_1D][k]);
      dTdy[U_1D][k] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][k], dTy[L_1D][k]);
    }
    
    int compx[4] = { L_1D, U_1D, L_1D, U_1D };
    int compy[4] = { L_1D, L_1D, U_1D, U_1D };

    int signx[4] = { 1.0, 1.0, -1.0, -1.0 };
    int signy[4] = { 1.0, -1.0, 1.0, -1.0 };

    for (int j = LL_2D; j <= UU_2D; ++j) {
      Bx += em_tot_d[j][BX]/4.0;
      By += em_tot_d[j][BY]/4.0;
      Bz += em_tot_d[j][BZ]/4.0;
    }

    double chi = alpha*n_avg*vth_avg;

    q[LL_2D][Q111] = chi*dTdx[compx[LL_2D]][T11];
    q[LL_2D][Q112] = chi*(2.0*dTdx[compx[LL_2D]][T12]
      + dTdy[compy[LL_2D]][T11])/3.0;
    q[LL_2D][Q113] = chi*2.0*dTdx[compx[LL_2D]][T13]/3.0;
    q[LL_2D][Q122] = chi*(dTdx[compx[LL_2D]][T22]
      + 2.0*dTdy[compy[LL_2D]][T12])/3.0;
    q[LL_2D][Q123] = chi*(dTdx[compx[LL_2D]][T23] + dTdy[compy[LL_2D]][T13])/3.0;
    q[LL_2D][Q133] = chi*dTdx[compx[LL_2D]][T33]/3.0;
    q[LL_2D][Q222] = chi*dTdy[compy[LL_2D]][T22];
    q[LL_2D][Q223] = chi*2.0*dTdy[compy[LL_2D]][T23]/3.0;
    q[LL_2D][Q233] = chi*dTdy[compy[LL_2D]][T33]/3.0;
    q[LL_2D][Q333] = 0.0;

    q[LU_2D][Q111] = chi*dTdx[compx[LU_2D]][T11];
    q[LU_2D][Q112] = chi*(2.0*dTdx[compx[LU_2D]][T12]
      + dTdy[compy[LU_2D]][T11])/3.0;
    q[LU_2D][Q113] = chi*2.0*dTdx[compx[LU_2D]][T13]/3.0;
    q[LU_2D][Q122] = chi*(dTdx[compx[LU_2D]][T22]
      + 2.0*dTdy[compy[LU_2D]][T12])/3.0;
    q[LU_2D][Q123] = chi*(dTdx[compx[LU_2D]][T23] + dTdy[compy[LU_2D]][T13])/3.0;
    q[LU_2D][Q133] = chi*dTdx[compx[LU_2D]][T33]/3.0;
    q[LU_2D][Q222] = chi*dTdy[compy[LU_2D]][T22];
    q[LU_2D][Q223] = chi*2.0*dTdy[compy[LU_2D]][T23]/3.0;
    q[LU_2D][Q233] = chi*dTdy[compy[LU_2D]][T33]/3.0;
    q[LU_2D][Q333] = 0.0;

    q[UL_2D][Q111] = chi*dTdx[compx[UL_2D]][T11];
    q[UL_2D][Q112] = chi*(2.0*dTdx[compx[UL_2D]][T12]
      + dTdy[compy[UL_2D]][T11])/3.0;
    q[UL_2D][Q113] = chi*2.0*dTdx[compx[UL_2D]][T13]/3.0;
    q[UL_2D][Q122] = chi*(dTdx[compx[UL_2D]][T22]
      + 2.0*dTdy[compy[UL_2D]][T12])/3.0;
    q[UL_2D][Q123] = chi*(dTdx[compx[UL_2D]][T23] + dTdy[compy[UL_2D]][T13])/3.0;
    q[UL_2D][Q133] = chi*dTdx[compx[UL_2D]][T33]/3.0;
    q[UL_2D][Q222] = chi*dTdy[compy[UL_2D]][T22];
    q[UL_2D][Q223] = chi*2.0*dTdy[compy[UL_2D]][T23]/3.0;
    q[UL_2D][Q233] = chi*dTdy[compy[UL_2D]][T33]/3.0;
    q[UL_2D][Q333] = 0.0;

    q[UU_2D][Q111] = chi*dTdx[compx[UU_2D]][T11];
    q[UU_2D][Q112] = chi*(2.0*dTdx[compx[UU_2D]][T12]
      + dTdy[compy[UU_2D]][T11])/3.0;
    q[UU_2D][Q113] = chi*2.0*dTdx[compx[UU_2D]][T13]/3.0;
    q[UU_2D][Q122] = chi*(dTdx[compx[UU_2D]][T22]
      + 2.0*dTdy[compy[UU_2D]][T12])/3.0;
    q[UU_2D][Q123] = chi*(dTdx[compx[UU_2D]][T23] + dTdy[compy[UU_2D]][T13])/3.0;
    q[UU_2D][Q133] = chi*dTdx[compx[UU_2D]][T33]/3.0;
    q[UU_2D][Q222] = chi*dTdy[compy[UU_2D]][T22];
    q[UU_2D][Q223] = chi*2.0*dTdy[compy[UU_2D]][T23]/3.0;
    q[UU_2D][Q233] = chi*dTdy[compy[UU_2D]][T33]/3.0;
    q[UU_2D][Q333] = 0.0;
    
    double e1[3] = { 1.0, 0.0, 0.0 };
    double e2[3] = { 0.0, 1.0, 0.0 };
    double e3[3] = { 0.0, 0.0, 1.0 };
     
    double tau1[3] = { 0.0 };
    double tau2[3] = { 0.0 };
    double b[3] = { 0.0 };

    b_unit_vec(Bx, By, Bz, 1.0e-15, tau1, tau2, b);
    rotate_source(q[LL_2D], q_src[LL_2D], e1, e2, e3, tau1, tau2, b);
    anisotropic_diffusion(w, q_src[LL_2D], q[LL_2D]);
    rotate_source(q[LL_2D], q_out[LL_2D], tau1, tau2, b, e1, e2, e3);

    b_unit_vec(Bx, By, Bz, 1.0e-15, tau1, tau2, b);
    rotate_source(q[LU_2D], q_src[LU_2D], e1, e2, e3, tau1, tau2, b);
    anisotropic_diffusion(w, q_src[LU_2D], q[LU_2D]);
    rotate_source(q[LU_2D], q_out[LU_2D], tau1, tau2, b, e1, e2, e3);

    b_unit_vec(Bx, By, Bz, 1.0e-15, tau1, tau2, b);
    rotate_source(q[UL_2D], q_src[UL_2D], e1, e2, e3, tau1, tau2, b);
    anisotropic_diffusion(w, q_src[UL_2D], q[UL_2D]);
    rotate_source(q[UL_2D], q_out[UL_2D], tau1, tau2, b, e1, e2, e3);

    b_unit_vec(Bx, By, Bz, 1.0e-15, tau1, tau2, b);
    rotate_source(q[UU_2D], q_src[UU_2D], e1, e2, e3, tau1, tau2, b);
    anisotropic_diffusion(w, q_src[UU_2D], q[UU_2D]);
    rotate_source(q[UU_2D], q_out[UU_2D], tau1, tau2, b, e1, e2, e3);

    rhs_d[LL_2D][P11] += signx[LL_2D]*q_out[LL_2D][Q111]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q112]/(2.0*dy);
    rhs_d[LL_2D][P12] += signx[LL_2D]*q_out[LL_2D][Q112]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q122]/(2.0*dy);
    rhs_d[LL_2D][P13] += signx[LL_2D]*q_out[LL_2D][Q113]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q123]/(2.0*dy);
    rhs_d[LL_2D][P22] += signx[LL_2D]*q_out[LL_2D][Q122]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q222]/(2.0*dy);
    rhs_d[LL_2D][P23] += signx[LL_2D]*q_out[LL_2D][Q123]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q223]/(2.0*dy);
    rhs_d[LL_2D][P33] += signx[LL_2D]*q_out[LL_2D][Q133]/(2.0*dx)
      + signy[LL_2D]*q_out[LL_2D][Q233]/(2.0*dy);

    rhs_d[LU_2D][P11] += signx[LU_2D]*q_out[LU_2D][Q111]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q112]/(2.0*dy);
    rhs_d[LU_2D][P12] += signx[LU_2D]*q_out[LU_2D][Q112]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q122]/(2.0*dy);
    rhs_d[LU_2D][P13] += signx[LU_2D]*q_out[LU_2D][Q113]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q123]/(2.0*dy);
    rhs_d[LU_2D][P22] += signx[LU_2D]*q_out[LU_2D][Q122]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q222]/(2.0*dy);
    rhs_d[LU_2D][P23] += signx[LU_2D]*q_out[LU_2D][Q123]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q223]/(2.0*dy);
    rhs_d[LU_2D][P33] += signx[LU_2D]*q_out[LU_2D][Q133]/(2.0*dx)
      + signy[LU_2D]*q_out[LU_2D][Q233]/(2.0*dy);

    rhs_d[UL_2D][P11] += signx[UL_2D]*q_out[UL_2D][Q111]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q112]/(2.0*dy);
    rhs_d[UL_2D][P12] += signx[UL_2D]*q_out[UL_2D][Q112]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q122]/(2.0*dy);
    rhs_d[UL_2D][P13] += signx[UL_2D]*q_out[UL_2D][Q113]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q123]/(2.0*dy);
    rhs_d[UL_2D][P22] += signx[UL_2D]*q_out[UL_2D][Q122]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q222]/(2.0*dy);
    rhs_d[UL_2D][P23] += signx[UL_2D]*q_out[UL_2D][Q123]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q223]/(2.0*dy);
    rhs_d[UL_2D][P33] += signx[UL_2D]*q_out[UL_2D][Q133]/(2.0*dx)
      + signy[UL_2D]*q_out[UL_2D][Q233]/(2.0*dy);

    rhs_d[UU_2D][P11] += signx[UU_2D]*q_out[UU_2D][Q111]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q112]/(2.0*dy);
    rhs_d[UU_2D][P12] += signx[UU_2D]*q_out[UU_2D][Q112]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q122]/(2.0*dy);
    rhs_d[UU_2D][P13] += signx[UU_2D]*q_out[UU_2D][Q113]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q123]/(2.0*dy);
    rhs_d[UU_2D][P22] += signx[UU_2D]*q_out[UU_2D][Q122]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q222]/(2.0*dy);
    rhs_d[UU_2D][P23] += signx[UU_2D]*q_out[UU_2D][Q123]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q223]/(2.0*dy);
    rhs_d[UU_2D][P33] += signx[UU_2D]*q_out[UU_2D][Q133]/(2.0*dx)
      + signy[UU_2D]*q_out[UU_2D][Q233]/(2.0*dy);
   
    double dmin = fmin(dx, dy);
    cfla = dt/(dmin*dmin);
  }
    
  return fmax(alpha*vth_avg*cfla, cfl);
}

gkyl_ten_moment_grad_closure*
gkyl_ten_moment_grad_closure_new(struct gkyl_ten_moment_grad_closure_inp inp)
{
  gkyl_ten_moment_grad_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_grad_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;
  up->cfl = inp.cfl;
  up->mass = inp.mass;

  if (inp.comm)
    up->comm = gkyl_comm_acquire(inp.comm);
  else
    up->comm = gkyl_null_comm_new();

  return up;
}

struct gkyl_ten_moment_grad_closure_status
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces,
  const struct gkyl_range *heat_flux_range, const struct gkyl_range *update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, double dt, struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };

  double cfla = 0.0, cfl = gces->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // delibrately a double

  gkyl_array_clear(rhs, 0.0);

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(update_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(heat_flux_range, offsets_centers);

  const double* fluid_d[sz[ndim-1]];
  const double* em_tot_d[sz[ndim-1]];
  double *rhs_d[sz[ndim-1]];

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_range);
  int counter = 0;
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_range, iter_vertex.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
      rhs_d[i] = gkyl_array_fetch(rhs, linc_center + offsets_vertices[i]);
    }

    /* cfla = calc_unmag_heat_flux(gces, fluid_d, gkyl_array_fetch(cflrate, linc_center), */
    /*   cfla, dt, rhs_d); */
    cfla = calc_mag_heat_flux_rotate(gces, fluid_d, em_tot_d,
      gkyl_array_fetch(cflrate, linc_center), cfla, dt, rhs_d);
    /* cfla = calc_mag_heat_flux_linsolve(gces, fluid_d, em_tot_d, */
    /*   gkyl_array_fetch(cflrate, linc_center), cfla, dt, rhs_d); */
    /* cfla = calc_mag_heat_flux_kernel(gces, fluid_d, em_tot_d, */
    /*   gkyl_array_fetch(cflrate, linc_center), cfla, dt, rhs_d); */
    counter = counter + 1;
  }

  if (cfla > cflm)
    is_cfl_violated = 1.0;

  // compute actual CFL, status & max-speed across all domains
  double red_vars[2] = { cfla, is_cfl_violated };
  double red_vars_global[2] = { 0.0, 0.0 };
  gkyl_comm_allreduce(gces->comm, GKYL_DOUBLE, GKYL_MAX, 2, &red_vars, &red_vars_global);

  cfla = red_vars_global[0];
  is_cfl_violated = red_vars_global[1];

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);
  
  if (is_cfl_violated > 0.0)
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_ten_moment_grad_closure_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
    };
  
  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).
  
  return (struct gkyl_ten_moment_grad_closure_status) {
    .success = is_cfl_violated > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
  };
}

void
gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure* up)
{
  free(up);
}
