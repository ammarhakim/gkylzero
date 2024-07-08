#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_moment_non_ideal_priv.h>

// Makes indexing cleaner
static const unsigned T11 = 0;
static const unsigned T12 = 1;
static const unsigned T13 = 2;
static const unsigned T22 = 3;
static const unsigned T23 = 4;
static const unsigned T33 = 5;

static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

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

static void
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

static void
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

static void
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
  const double *fluid_d[], double *cflrate, double *temp_grad_d,
  double cfl, double dt)
{
  const int ndim = gces->ndim;
  double alpha = 1.0/gces->k0;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double vth_avg = 0.0;
  double n_avg = 0.0;
  double Ax[6] = {0.0};
  double Ay[6] = {0.0};
  double Az[6] = {0.0};
  double Bx[6] = {0.0};
  double By[6] = {0.0};
  double Bz[6] = {0.0};
  double Cx[6] = {0.0};
  double Cy[6] = {0.0};
  double Cz[6] = {0.0};
  double Dx[6] = {0.0};
  double Dy[6] = {0.0};
  double Dz[6] = {0.0};
  double limit = 0.75;
  double cfla = 0.0;
  
  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double dTdx[6] = {0.0};
    double Tij[2][6] = {0.0};
    double rho[2] = {0.0};
    double p[2] = {0.0};
    var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
    p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);
    vth_avg = sqrt(p_avg/rho_avg);
    n_avg = rho_avg/gces->mass;
    
    dTdx[T11] = calc_sym_grad_1D(dx, Tij[L_1D][T11], Tij[U_1D][T11]);
    dTdx[T12] = calc_sym_grad_1D(dx, Tij[L_1D][T12], Tij[U_1D][T12]);
    dTdx[T13] = calc_sym_grad_1D(dx, Tij[L_1D][T13], Tij[U_1D][T13]);
    dTdx[T22] = calc_sym_grad_1D(dx, Tij[L_1D][T22], Tij[U_1D][T22]);
    dTdx[T23] = calc_sym_grad_1D(dx, Tij[L_1D][T23], Tij[U_1D][T23]);
    dTdx[T33] = calc_sym_grad_1D(dx, Tij[L_1D][T33], Tij[U_1D][T33]);

    for (int k = 0; k < 6; ++k) {
      temp_grad_d[k] = alpha*n_avg*vth_avg*dTdx[k];
    }

    cfla = dt/(dx*dx);
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);
    vth_avg = sqrt(p_avg/rho_avg);
    n_avg = rho_avg/gces->mass;

    Ax[T11] = calc_sym_grad_1D(dx, Tij[LL_2D][T11], Tij[UL_2D][T11]);
    Ax[T12] = calc_sym_grad_1D(dx, Tij[LL_2D][T12], Tij[UL_2D][T12]);
    Ax[T13] = calc_sym_grad_1D(dx, Tij[LL_2D][T13], Tij[UL_2D][T13]);
    Ax[T22] = calc_sym_grad_1D(dx, Tij[LL_2D][T22], Tij[UL_2D][T22]);
    Ax[T23] = calc_sym_grad_1D(dx, Tij[LL_2D][T23], Tij[UL_2D][T23]);
    Ax[T33] = calc_sym_grad_1D(dx, Tij[LL_2D][T33], Tij[UL_2D][T33]);

    Bx[T11] = calc_sym_grad_1D(dx, Tij[LU_2D][T11], Tij[UU_2D][T11]);
    Bx[T12] = calc_sym_grad_1D(dx, Tij[LU_2D][T12], Tij[UU_2D][T12]);
    Bx[T13] = calc_sym_grad_1D(dx, Tij[LU_2D][T13], Tij[UU_2D][T13]);
    Bx[T22] = calc_sym_grad_1D(dx, Tij[LU_2D][T22], Tij[UU_2D][T22]);
    Bx[T23] = calc_sym_grad_1D(dx, Tij[LU_2D][T23], Tij[UU_2D][T23]);
    Bx[T33] = calc_sym_grad_1D(dx, Tij[LU_2D][T33], Tij[UU_2D][T33]);

    dTdx[L_1D][T11] = calc_sym_grad_limiter_2D(limit, Ax[T11], Bx[T11]);
    dTdx[L_1D][T12] = calc_sym_grad_limiter_2D(limit, Ax[T12], Bx[T12]);
    dTdx[L_1D][T13] = calc_sym_grad_limiter_2D(limit, Ax[T13], Bx[T13]);
    dTdx[L_1D][T22] = calc_sym_grad_limiter_2D(limit, Ax[T22], Bx[T22]);
    dTdx[L_1D][T23] = calc_sym_grad_limiter_2D(limit, Ax[T23], Bx[T23]);
    dTdx[L_1D][T33] = calc_sym_grad_limiter_2D(limit, Ax[T33], Bx[T33]);

    dTdx[U_1D][T11] = calc_sym_grad_limiter_2D(limit, Bx[T11], Ax[T11]);
    dTdx[U_1D][T12] = calc_sym_grad_limiter_2D(limit, Bx[T12], Ax[T12]);
    dTdx[U_1D][T13] = calc_sym_grad_limiter_2D(limit, Bx[T13], Ax[T13]);
    dTdx[U_1D][T22] = calc_sym_grad_limiter_2D(limit, Bx[T22], Ax[T22]);
    dTdx[U_1D][T23] = calc_sym_grad_limiter_2D(limit, Bx[T23], Ax[T23]);
    dTdx[U_1D][T33] = calc_sym_grad_limiter_2D(limit, Bx[T33], Ax[T33]);

    Ay[T11] = calc_sym_grad_1D(dy, Tij[LL_2D][T11], Tij[LU_2D][T11]);
    Ay[T12] = calc_sym_grad_1D(dy, Tij[LL_2D][T12], Tij[LU_2D][T12]);
    Ay[T13] = calc_sym_grad_1D(dy, Tij[LL_2D][T13], Tij[LU_2D][T13]);
    Ay[T22] = calc_sym_grad_1D(dy, Tij[LL_2D][T22], Tij[LU_2D][T22]);
    Ay[T23] = calc_sym_grad_1D(dy, Tij[LL_2D][T23], Tij[LU_2D][T23]);
    Ay[T33] = calc_sym_grad_1D(dy, Tij[LL_2D][T33], Tij[LU_2D][T33]);

    By[T11] = calc_sym_grad_1D(dy, Tij[UL_2D][T11], Tij[UU_2D][T11]);
    By[T12] = calc_sym_grad_1D(dy, Tij[UL_2D][T12], Tij[UU_2D][T12]);
    By[T13] = calc_sym_grad_1D(dy, Tij[UL_2D][T13], Tij[UU_2D][T13]);
    By[T22] = calc_sym_grad_1D(dy, Tij[UL_2D][T22], Tij[UU_2D][T22]);
    By[T23] = calc_sym_grad_1D(dy, Tij[UL_2D][T23], Tij[UU_2D][T23]);
    By[T33] = calc_sym_grad_1D(dy, Tij[UL_2D][T33], Tij[UU_2D][T33]);

    dTdy[L_1D][T11] = calc_sym_grad_limiter_2D(limit, Ay[T11], By[T11]);
    dTdy[L_1D][T12] = calc_sym_grad_limiter_2D(limit, Ay[T12], By[T12]);
    dTdy[L_1D][T13] = calc_sym_grad_limiter_2D(limit, Ay[T13], By[T13]);
    dTdy[L_1D][T22] = calc_sym_grad_limiter_2D(limit, Ay[T22], By[T22]);
    dTdy[L_1D][T23] = calc_sym_grad_limiter_2D(limit, Ay[T23], By[T23]);
    dTdy[L_1D][T33] = calc_sym_grad_limiter_2D(limit, Ay[T33], By[T33]);

    dTdy[U_1D][T11] = calc_sym_grad_limiter_2D(limit, By[T11], Ay[T11]);
    dTdy[U_1D][T12] = calc_sym_grad_limiter_2D(limit, By[T12], Ay[T12]);
    dTdy[U_1D][T13] = calc_sym_grad_limiter_2D(limit, By[T13], Ay[T13]);
    dTdy[U_1D][T22] = calc_sym_grad_limiter_2D(limit, By[T22], Ay[T22]);
    dTdy[U_1D][T23] = calc_sym_grad_limiter_2D(limit, By[T23], Ay[T23]);
    dTdy[U_1D][T33] = calc_sym_grad_limiter_2D(limit, By[T33], Ay[T33]);

    for (int j = L_1D; j <= U_1D; ++j) {
      for (int k = 0; k < 6; ++k) {
        temp_grad_d[6*j*ndim + k] = alpha*n_avg*vth_avg*dTdx[j][k];
        temp_grad_d[6*(j*ndim + 1) + k] = alpha*n_avg*vth_avg*dTdy[j][k];
      }
    }

    double dmin = fmin(dx, dy);
    cfla = dt/(dmin*dmin);
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];

    double dTdx[4][6] = {0.0};
    double dTdy[4][6] = {0.0};
    double dTdz[4][6] = {0.0};
    double Tij[8][6] = {0.0};
    double rho[8] = {0.0};
    double p[8] = {0.0};
    var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                   rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
    p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                                 p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);
    vth_avg = sqrt(p_avg/rho_avg);
    n_avg = rho_avg/gces->mass;

    Ax[T11] = calc_sym_grad_1D(dx, Tij[LLL_3D][T11], Tij[ULL_3D][T11]);
    Ax[T12] = calc_sym_grad_1D(dx, Tij[LLL_3D][T12], Tij[ULL_3D][T12]);
    Ax[T13] = calc_sym_grad_1D(dx, Tij[LLL_3D][T13], Tij[ULL_3D][T13]);
    Ax[T22] = calc_sym_grad_1D(dx, Tij[LLL_3D][T22], Tij[ULL_3D][T22]);
    Ax[T23] = calc_sym_grad_1D(dx, Tij[LLL_3D][T23], Tij[ULL_3D][T23]);
    Ax[T33] = calc_sym_grad_1D(dx, Tij[LLL_3D][T33], Tij[ULL_3D][T33]);

    Bx[T11] = calc_sym_grad_1D(dx, Tij[LUL_3D][T11], Tij[UUL_3D][T11]);
    Bx[T12] = calc_sym_grad_1D(dx, Tij[LUL_3D][T12], Tij[UUL_3D][T12]);
    Bx[T13] = calc_sym_grad_1D(dx, Tij[LUL_3D][T13], Tij[UUL_3D][T13]);
    Bx[T22] = calc_sym_grad_1D(dx, Tij[LUL_3D][T22], Tij[UUL_3D][T22]);
    Bx[T23] = calc_sym_grad_1D(dx, Tij[LUL_3D][T23], Tij[UUL_3D][T23]);
    Bx[T33] = calc_sym_grad_1D(dx, Tij[LUL_3D][T33], Tij[UUL_3D][T33]);

    Cx[T11] = calc_sym_grad_1D(dx, Tij[LLU_3D][T11], Tij[ULU_3D][T11]);
    Cx[T12] = calc_sym_grad_1D(dx, Tij[LLU_3D][T12], Tij[ULU_3D][T12]);
    Cx[T13] = calc_sym_grad_1D(dx, Tij[LLU_3D][T13], Tij[ULU_3D][T13]);
    Cx[T22] = calc_sym_grad_1D(dx, Tij[LLU_3D][T22], Tij[ULU_3D][T22]);
    Cx[T23] = calc_sym_grad_1D(dx, Tij[LLU_3D][T23], Tij[ULU_3D][T23]);
    Cx[T33] = calc_sym_grad_1D(dx, Tij[LLU_3D][T33], Tij[ULU_3D][T33]);

    Dx[T11] = calc_sym_grad_1D(dx, Tij[LUU_3D][T11], Tij[UUU_3D][T11]);
    Dx[T12] = calc_sym_grad_1D(dx, Tij[LUU_3D][T12], Tij[UUU_3D][T12]);
    Dx[T13] = calc_sym_grad_1D(dx, Tij[LUU_3D][T13], Tij[UUU_3D][T13]);
    Dx[T22] = calc_sym_grad_1D(dx, Tij[LUU_3D][T22], Tij[UUU_3D][T22]);
    Dx[T23] = calc_sym_grad_1D(dx, Tij[LUU_3D][T23], Tij[UUU_3D][T23]);
    Dx[T33] = calc_sym_grad_1D(dx, Tij[LUU_3D][T33], Tij[UUU_3D][T33]);

    dTdx[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, Ax[T11], Bx[T11], Cx[T11], Dx[T11]);
    dTdx[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, Ax[T12], Bx[T12], Cx[T12], Dx[T12]);
    dTdx[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, Ax[T13], Bx[T13], Cx[T13], Dx[T13]); 
    dTdx[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, Ax[T22], Bx[T22], Cx[T22], Dx[T22]);
    dTdx[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, Ax[T23], Bx[T23], Cx[T23], Dx[T23]);
    dTdx[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, Ax[T33], Bx[T33], Cx[T33], Dx[T33]);

    dTdx[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, Bx[T11], Ax[T11], Cx[T11], Dx[T11]);
    dTdx[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, Bx[T12], Ax[T12], Cx[T12], Dx[T12]);
    dTdx[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, Bx[T13], Ax[T13], Cx[T13], Dx[T13]); 
    dTdx[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, Bx[T22], Ax[T22], Cx[T22], Dx[T22]);
    dTdx[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, Bx[T23], Ax[T23], Cx[T23], Dx[T23]);
    dTdx[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, Bx[T33], Ax[T33], Cx[T33], Dx[T33]);

    dTdx[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, Cx[T11], Ax[T11], Bx[T11], Dx[T11]);
    dTdx[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, Cx[T12], Ax[T12], Bx[T12], Dx[T12]);
    dTdx[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, Cx[T13], Ax[T13], Bx[T13], Dx[T13]); 
    dTdx[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, Cx[T22], Ax[T22], Bx[T22], Dx[T22]);
    dTdx[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, Cx[T23], Ax[T23], Bx[T23], Dx[T23]);
    dTdx[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, Cx[T33], Ax[T33], Bx[T33], Dx[T33]);

    dTdx[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, Dx[T11], Ax[T11], Bx[T11], Cx[T11]);
    dTdx[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, Dx[T12], Ax[T12], Bx[T12], Cx[T12]);
    dTdx[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, Dx[T13], Ax[T13], Bx[T13], Cx[T13]); 
    dTdx[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, Dx[T22], Ax[T22], Bx[T22], Cx[T22]);
    dTdx[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, Dx[T23], Ax[T23], Bx[T23], Cx[T23]);
    dTdx[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, Dx[T33], Ax[T33], Bx[T33], Cx[T33]);

    Ay[T11] = calc_sym_grad_1D(dy, Tij[LLL_3D][T11], Tij[LUL_3D][T11]);
    Ay[T12] = calc_sym_grad_1D(dy, Tij[LLL_3D][T12], Tij[LUL_3D][T12]);
    Ay[T13] = calc_sym_grad_1D(dy, Tij[LLL_3D][T13], Tij[LUL_3D][T13]);
    Ay[T22] = calc_sym_grad_1D(dy, Tij[LLL_3D][T22], Tij[LUL_3D][T22]);
    Ay[T23] = calc_sym_grad_1D(dy, Tij[LLL_3D][T23], Tij[LUL_3D][T23]);
    Ay[T33] = calc_sym_grad_1D(dy, Tij[LLL_3D][T33], Tij[LUL_3D][T33]);

    By[T11] = calc_sym_grad_1D(dy, Tij[ULL_3D][T11], Tij[UUL_3D][T11]);
    By[T12] = calc_sym_grad_1D(dy, Tij[ULL_3D][T12], Tij[UUL_3D][T12]);
    By[T13] = calc_sym_grad_1D(dy, Tij[ULL_3D][T13], Tij[UUL_3D][T13]);
    By[T22] = calc_sym_grad_1D(dy, Tij[ULL_3D][T22], Tij[UUL_3D][T22]);
    By[T23] = calc_sym_grad_1D(dy, Tij[ULL_3D][T23], Tij[UUL_3D][T23]);
    By[T33] = calc_sym_grad_1D(dy, Tij[ULL_3D][T33], Tij[UUL_3D][T33]);

    Cy[T11] = calc_sym_grad_1D(dy, Tij[LLU_3D][T11], Tij[LUU_3D][T11]);
    Cy[T12] = calc_sym_grad_1D(dy, Tij[LLU_3D][T12], Tij[LUU_3D][T12]);
    Cy[T13] = calc_sym_grad_1D(dy, Tij[LLU_3D][T13], Tij[LUU_3D][T13]);
    Cy[T22] = calc_sym_grad_1D(dy, Tij[LLU_3D][T22], Tij[LUU_3D][T22]);
    Cy[T23] = calc_sym_grad_1D(dy, Tij[LLU_3D][T23], Tij[LUU_3D][T23]);
    Cy[T33] = calc_sym_grad_1D(dy, Tij[LLU_3D][T33], Tij[LUU_3D][T33]);

    Dy[T11] = calc_sym_grad_1D(dy, Tij[ULU_3D][T11], Tij[UUU_3D][T11]);
    Dy[T12] = calc_sym_grad_1D(dy, Tij[ULU_3D][T12], Tij[UUU_3D][T12]);
    Dy[T13] = calc_sym_grad_1D(dy, Tij[ULU_3D][T13], Tij[UUU_3D][T13]);
    Dy[T22] = calc_sym_grad_1D(dy, Tij[ULU_3D][T22], Tij[UUU_3D][T22]);
    Dy[T23] = calc_sym_grad_1D(dy, Tij[ULU_3D][T23], Tij[UUU_3D][T23]);
    Dy[T33] = calc_sym_grad_1D(dy, Tij[ULU_3D][T33], Tij[UUU_3D][T33]);

    dTdy[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, Ay[T11], By[T11], Cy[T11], Dy[T11]);
    dTdy[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, Ay[T12], By[T12], Cy[T12], Dy[T12]);
    dTdy[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, Ay[T13], By[T13], Cy[T13], Dy[T13]); 
    dTdy[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, Ay[T22], By[T22], Cy[T22], Dy[T22]);
    dTdy[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, Ay[T23], By[T23], Cy[T23], Dy[T23]);
    dTdy[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, Ay[T33], By[T33], Cy[T33], Dy[T33]);

    dTdy[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, By[T11], Ay[T11], Cy[T11], Dy[T11]);
    dTdy[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, By[T12], Ay[T12], Cy[T12], Dy[T12]);
    dTdy[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, By[T13], Ay[T13], Cy[T13], Dy[T13]); 
    dTdy[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, By[T22], Ay[T22], Cy[T22], Dy[T22]);
    dTdy[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, By[T23], Ay[T23], Cy[T23], Dy[T23]);
    dTdy[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, By[T33], Ay[T33], Cy[T33], Dy[T33]);

    dTdy[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, Cy[T11], Ay[T11], By[T11], Dy[T11]);
    dTdy[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, Cy[T12], Ay[T12], By[T12], Dy[T12]);
    dTdy[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, Cy[T13], Ay[T13], By[T13], Dy[T13]); 
    dTdy[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, Cy[T22], Ay[T22], By[T22], Dy[T22]);
    dTdy[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, Cy[T23], Ay[T23], By[T23], Dy[T23]);
    dTdy[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, Cy[T33], Ay[T33], By[T33], Dy[T33]);

    dTdy[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, Dy[T11], Ay[T11], By[T11], Cy[T11]);
    dTdy[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, Dy[T12], Ay[T12], By[T12], Cy[T12]);
    dTdy[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, Dy[T13], Ay[T13], By[T13], Cy[T13]); 
    dTdy[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, Dy[T22], Ay[T22], By[T22], Cy[T22]);
    dTdy[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, Dy[T23], Ay[T23], By[T23], Cy[T23]);
    dTdy[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, Dy[T33], Ay[T33], By[T33], Cy[T33]);

    Az[T11] = calc_sym_grad_1D(dz, Tij[LLL_3D][T11], Tij[LLU_3D][T11]);
    Az[T12] = calc_sym_grad_1D(dz, Tij[LLL_3D][T12], Tij[LLU_3D][T12]);
    Az[T13] = calc_sym_grad_1D(dz, Tij[LLL_3D][T13], Tij[LLU_3D][T13]);
    Az[T22] = calc_sym_grad_1D(dz, Tij[LLL_3D][T22], Tij[LLU_3D][T22]);
    Az[T23] = calc_sym_grad_1D(dz, Tij[LLL_3D][T23], Tij[LLU_3D][T23]);
    Az[T33] = calc_sym_grad_1D(dz, Tij[LLL_3D][T33], Tij[LLU_3D][T33]);

    Bz[T11] = calc_sym_grad_1D(dz, Tij[ULL_3D][T11], Tij[ULU_3D][T11]);
    Bz[T12] = calc_sym_grad_1D(dz, Tij[ULL_3D][T12], Tij[ULU_3D][T12]);
    Bz[T13] = calc_sym_grad_1D(dz, Tij[ULL_3D][T13], Tij[ULU_3D][T13]);
    Bz[T22] = calc_sym_grad_1D(dz, Tij[ULL_3D][T22], Tij[ULU_3D][T22]);
    Bz[T23] = calc_sym_grad_1D(dz, Tij[ULL_3D][T23], Tij[ULU_3D][T23]);
    Bz[T33] = calc_sym_grad_1D(dz, Tij[ULL_3D][T33], Tij[ULU_3D][T33]);

    Cz[T11] = calc_sym_grad_1D(dz, Tij[LUL_3D][T11], Tij[LUU_3D][T11]);
    Cz[T12] = calc_sym_grad_1D(dz, Tij[LUL_3D][T12], Tij[LUU_3D][T12]);
    Cz[T13] = calc_sym_grad_1D(dz, Tij[LUL_3D][T13], Tij[LUU_3D][T13]);
    Cz[T22] = calc_sym_grad_1D(dz, Tij[LUL_3D][T22], Tij[LUU_3D][T22]);
    Cz[T23] = calc_sym_grad_1D(dz, Tij[LUL_3D][T23], Tij[LUU_3D][T23]);
    Cz[T33] = calc_sym_grad_1D(dz, Tij[LUL_3D][T33], Tij[LUU_3D][T33]);

    Dz[T11] = calc_sym_grad_1D(dz, Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
    Dz[T12] = calc_sym_grad_1D(dz, Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
    Dz[T13] = calc_sym_grad_1D(dz, Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
    Dz[T22] = calc_sym_grad_1D(dz, Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
    Dz[T23] = calc_sym_grad_1D(dz, Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
    Dz[T33] = calc_sym_grad_1D(dz, Tij[UUL_3D][T33], Tij[UUU_3D][T33]);

    dTdz[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, Az[T11], Bz[T11], Cz[T11], Dz[T11]);
    dTdz[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, Az[T12], Bz[T12], Cz[T12], Dz[T12]);
    dTdz[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, Az[T13], Bz[T13], Cz[T13], Dz[T13]); 
    dTdz[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, Az[T22], Bz[T22], Cz[T22], Dz[T22]);
    dTdz[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, Az[T23], Bz[T23], Cz[T23], Dz[T23]);
    dTdz[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, Az[T33], Bz[T33], Cz[T33], Dz[T33]);

    dTdz[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, Bz[T11], Az[T11], Cz[T11], Dz[T11]);
    dTdz[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, Bz[T12], Az[T12], Cz[T12], Dz[T12]);
    dTdz[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, Bz[T13], Az[T13], Cz[T13], Dz[T13]); 
    dTdz[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, Bz[T22], Az[T22], Cz[T22], Dz[T22]);
    dTdz[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, Bz[T23], Az[T23], Cz[T23], Dz[T23]);
    dTdz[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, Bz[T33], Az[T33], Cz[T33], Dz[T33]);

    dTdz[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, Cz[T11], Az[T11], Bz[T11], Dz[T11]);
    dTdz[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, Cz[T12], Az[T12], Bz[T12], Dz[T12]);
    dTdz[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, Cz[T13], Az[T13], Bz[T13], Dz[T13]); 
    dTdz[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, Cz[T22], Az[T22], Bz[T22], Dz[T22]);
    dTdz[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, Cz[T23], Az[T23], Bz[T23], Dz[T23]);
    dTdz[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, Cz[T33], Az[T33], Bz[T33], Dz[T33]);

    dTdz[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, Dz[T11], Az[T11], Bz[T11], Cz[T11]);
    dTdz[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, Dz[T12], Az[T12], Bz[T12], Cz[T12]);
    dTdz[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, Dz[T13], Az[T13], Bz[T13], Cz[T13]); 
    dTdz[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, Dz[T22], Az[T22], Bz[T22], Cz[T22]);
    dTdz[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, Dz[T23], Az[T23], Bz[T23], Cz[T23]);
    dTdz[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, Dz[T33], Az[T33], Bz[T33], Cz[T33]);

    for (int j = LL_2D; j <= UU_2D; ++j) {
      for (int k = 0; k < 6; ++k) {
        temp_grad_d[6*j*ndim + k] = alpha*n_avg*vth_avg*dTdx[j][k];
        temp_grad_d[6*(j*ndim + 1) + k] = alpha*n_avg*vth_avg*dTdy[j][k];
        temp_grad_d[6*(j*ndim + 2) + k] = alpha*n_avg*vth_avg*dTdz[j][k];
      }
    }

    double dmin = fmin(fmin(dx, dy), dz);
    cfla = dt/(dmin*dmin);
  }

  return fmax(alpha*vth_avg*cfla, cfl);
}

static void
calc_grad_closure_update(const gkyl_ten_moment_grad_closure *gces,
  const double *temp_grad_d[], double *rhs)
{
  const int ndim = gces->ndim;
  int ncomp = 1;
  double divQx[6] = {0.0};
  double divQy[6] = {0.0};
  double divQz[6] = {0.0};

  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double dTdx[2][6] = {0.0};
    double dTdy[2][6] = {0.0};
    double dTdz[2][6] = {0.0};
    double q[2][10] = {0.0};
    
    for (int j = L_1D; j <= U_1D; ++j) {
      for (int k = 0; k < 6; ++k) {
        dTdx[j][k] = temp_grad_d[j][k];
      }
      q[j][Q111] = (dTdx[j][T11] + dTdx[j][T11] + dTdx[j][T11])/3.0;
      q[j][Q112] = (dTdx[j][T12] + dTdx[j][T12] + dTdy[j][T11])/3.0;
      q[j][Q113] = (dTdx[j][T13] + dTdx[j][T13] + dTdz[j][T11])/3.0;
      q[j][Q122] = (dTdx[j][T22] + dTdy[j][T12] + dTdy[j][T12])/3.0;
      q[j][Q123] = (dTdx[j][T23] + dTdy[j][T13] + dTdz[j][T12])/3.0;
      q[j][Q133] = (dTdx[j][T33] + dTdz[j][T13] + dTdz[j][T13])/3.0;
      q[j][Q222] = (dTdy[j][T22] + dTdy[j][T22] + dTdy[j][T22])/3.0;
      q[j][Q223] = (dTdy[j][T23] + dTdy[j][T23] + dTdz[j][T22])/3.0;
      q[j][Q233] = (dTdy[j][T33] + dTdz[j][T23] + dTdz[j][T23])/3.0;
      q[j][Q333] = (dTdz[j][T33] + dTdz[j][T33] + dTdz[j][T33])/3.0;
    }

    divQx[0] = calc_sym_grad_1D(dx, q[L_1D][Q111], q[U_1D][Q111]);
    divQx[1] = calc_sym_grad_1D(dx, q[L_1D][Q112], q[U_1D][Q112]);
    divQx[2] = calc_sym_grad_1D(dx, q[L_1D][Q113], q[U_1D][Q113]);
    divQx[3] = calc_sym_grad_1D(dx, q[L_1D][Q122], q[U_1D][Q122]);
    divQx[4] = calc_sym_grad_1D(dx, q[L_1D][Q123], q[U_1D][Q123]);
    divQx[5] = calc_sym_grad_1D(dx, q[L_1D][Q133], q[U_1D][Q133]);
  }
  else if (ndim == 2) {
    ncomp = 2;
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double dTdx[4][6] = {0.0};
    double dTdy[4][6] = {0.0};
    double dTdz[4][6] = {0.0};
    double q[4][10] = {0.0};

    int compx[4] = { U_1D, L_1D, U_1D, L_1D };
    int compy[4] = { U_1D, U_1D, L_1D, L_1D };
    for (int j = LL_2D; j <= UU_2D; ++j) {
      for (int k = 0; k < 6; ++k) {
        dTdx[j][k] = temp_grad_d[j][6*compx[j]*ndim + k];
        dTdy[j][k] = temp_grad_d[j][6*(compy[j]*ndim + 1) + k];
      }
      q[j][Q111] = (dTdx[j][T11] + dTdx[j][T11] + dTdx[j][T11])/3.0;
      q[j][Q112] = (dTdx[j][T12] + dTdx[j][T12] + dTdy[j][T11])/3.0;
      q[j][Q113] = (dTdx[j][T13] + dTdx[j][T13] + dTdz[j][T11])/3.0;
      q[j][Q122] = (dTdx[j][T22] + dTdy[j][T12] + dTdy[j][T12])/3.0;
      q[j][Q123] = (dTdx[j][T23] + dTdy[j][T13] + dTdz[j][T12])/3.0;
      q[j][Q133] = (dTdx[j][T33] + dTdz[j][T13] + dTdz[j][T13])/3.0;
      q[j][Q222] = (dTdy[j][T22] + dTdy[j][T22] + dTdy[j][T22])/3.0;
      q[j][Q223] = (dTdy[j][T23] + dTdy[j][T23] + dTdz[j][T22])/3.0;
      q[j][Q233] = (dTdy[j][T33] + dTdz[j][T23] + dTdz[j][T23])/3.0;
      q[j][Q333] = (dTdz[j][T33] + dTdz[j][T33] + dTdz[j][T33])/3.0;
    }

    divQx[0] = calc_sym_gradx_2D(dx, q[LL_2D][Q111], q[LU_2D][Q111], q[UL_2D][Q111], q[UU_2D][Q111]);
    divQx[1] = calc_sym_gradx_2D(dx, q[LL_2D][Q112], q[LU_2D][Q112], q[UL_2D][Q112], q[UU_2D][Q112]);
    divQx[2] = calc_sym_gradx_2D(dx, q[LL_2D][Q113], q[LU_2D][Q113], q[UL_2D][Q113], q[UU_2D][Q113]);
    divQx[3] = calc_sym_gradx_2D(dx, q[LL_2D][Q122], q[LU_2D][Q122], q[UL_2D][Q122], q[UU_2D][Q122]);
    divQx[4] = calc_sym_gradx_2D(dx, q[LL_2D][Q123], q[LU_2D][Q123], q[UL_2D][Q123], q[UU_2D][Q123]);
    divQx[5] = calc_sym_gradx_2D(dx, q[LL_2D][Q133], q[LU_2D][Q133], q[UL_2D][Q133], q[UU_2D][Q133]);

    divQy[0] = calc_sym_grady_2D(dy, q[LL_2D][Q112], q[LU_2D][Q112], q[UL_2D][Q112], q[UU_2D][Q112]);
    divQy[1] = calc_sym_grady_2D(dy, q[LL_2D][Q122], q[LU_2D][Q122], q[UL_2D][Q122], q[UU_2D][Q122]);
    divQy[2] = calc_sym_grady_2D(dy, q[LL_2D][Q123], q[LU_2D][Q123], q[UL_2D][Q123], q[UU_2D][Q123]);
    divQy[3] = calc_sym_grady_2D(dy, q[LL_2D][Q222], q[LU_2D][Q222], q[UL_2D][Q222], q[UU_2D][Q222]);
    divQy[4] = calc_sym_grady_2D(dy, q[LL_2D][Q223], q[LU_2D][Q223], q[UL_2D][Q223], q[UU_2D][Q223]);
    divQy[5] = calc_sym_grady_2D(dy, q[LL_2D][Q233], q[LU_2D][Q233], q[UL_2D][Q233], q[UU_2D][Q233]);
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];

    double dTdx[8][6] = {0.0};
    double dTdy[8][6] = {0.0};
    double dTdz[8][6] = {0.0};
    double q[8][10] = {0.0};

    int compx[8] = { UU_2D, LU_2D, UL_2D, LL_2D, UU_2D, LU_2D, UL_2D, LL_2D };
    int compy[8] = { UU_2D, UL_2D, UU_2D, UL_2D, LU_2D, LL_2D, LU_2D, LL_2D };
    int compz[8] = { UU_2D, UU_2D, UL_2D, UL_2D, LU_2D, LU_2D, LL_2D, LL_2D };
    for (int j = LL_2D; j <= UU_2D; ++j) {
      for (int k = 0; k < 6; ++k) {
        dTdx[j][k] = temp_grad_d[j][6*compx[j]*ndim + k];
        dTdy[j][k] = temp_grad_d[j][6*(compy[j]*ndim + 1) + k];
        dTdz[j][k] = temp_grad_d[j][6*(compz[j]*ndim + 2) + k];
      }
      q[j][Q111] = (dTdx[j][T11] + dTdx[j][T11] + dTdx[j][T11])/3.0;
      q[j][Q112] = (dTdx[j][T12] + dTdx[j][T12] + dTdy[j][T11])/3.0;
      q[j][Q113] = (dTdx[j][T13] + dTdx[j][T13] + dTdz[j][T11])/3.0;
      q[j][Q122] = (dTdx[j][T22] + dTdy[j][T12] + dTdy[j][T12])/3.0;
      q[j][Q123] = (dTdx[j][T23] + dTdy[j][T13] + dTdz[j][T12])/3.0;
      q[j][Q133] = (dTdx[j][T33] + dTdz[j][T13] + dTdz[j][T13])/3.0;
      q[j][Q222] = (dTdy[j][T22] + dTdy[j][T22] + dTdy[j][T22])/3.0;
      q[j][Q223] = (dTdy[j][T23] + dTdy[j][T23] + dTdz[j][T22])/3.0;
      q[j][Q233] = (dTdy[j][T33] + dTdz[j][T23] + dTdz[j][T23])/3.0;
      q[j][Q333] = (dTdz[j][T33] + dTdz[j][T33] + dTdz[j][T33])/3.0;
    }

    divQx[0] = calc_sym_gradx_3D(dx, q[LLL_3D][Q111], q[LLU_3D][Q111], q[LUL_3D][Q111], q[LUU_3D][Q111],
                                     q[ULL_3D][Q111], q[ULU_3D][Q111], q[UUL_3D][Q111], q[UUU_3D][Q111]);
    divQx[1] = calc_sym_gradx_3D(dx, q[LLL_3D][Q112], q[LLU_3D][Q112], q[LUL_3D][Q112], q[LUU_3D][Q112],
                                     q[ULL_3D][Q112], q[ULU_3D][Q112], q[UUL_3D][Q112], q[UUU_3D][Q112]);
    divQx[2] = calc_sym_gradx_3D(dx, q[LLL_3D][Q113], q[LLU_3D][Q113], q[LUL_3D][Q113], q[LUU_3D][Q113],
                                     q[ULL_3D][Q113], q[ULU_3D][Q113], q[UUL_3D][Q113], q[UUU_3D][Q113]);
    divQx[3] = calc_sym_gradx_3D(dx, q[LLL_3D][Q122], q[LLU_3D][Q122], q[LUL_3D][Q122], q[LUU_3D][Q122],
                                     q[ULL_3D][Q122], q[ULU_3D][Q122], q[UUL_3D][Q122], q[UUU_3D][Q122]);
    divQx[4] = calc_sym_gradx_3D(dx, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQx[5] = calc_sym_gradx_3D(dx, q[LLL_3D][Q133], q[LLU_3D][Q133], q[LUL_3D][Q133], q[LUU_3D][Q133],
                                     q[ULL_3D][Q133], q[ULU_3D][Q133], q[UUL_3D][Q133], q[UUU_3D][Q133]);

    divQy[0] = calc_sym_grady_3D(dy, q[LLL_3D][Q112], q[LLU_3D][Q112], q[LUL_3D][Q112], q[LUU_3D][Q112],
                                     q[ULL_3D][Q112], q[ULU_3D][Q112], q[UUL_3D][Q112], q[UUU_3D][Q112]);
    divQy[1] = calc_sym_grady_3D(dy, q[LLL_3D][Q122], q[LLU_3D][Q122], q[LUL_3D][Q122], q[LUU_3D][Q122],
                                     q[ULL_3D][Q122], q[ULU_3D][Q122], q[UUL_3D][Q122], q[UUU_3D][Q122]);
    divQy[2] = calc_sym_grady_3D(dy, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQy[3] = calc_sym_grady_3D(dy, q[LLL_3D][Q222], q[LLU_3D][Q222], q[LUL_3D][Q222], q[LUU_3D][Q222],
                                     q[ULL_3D][Q222], q[ULU_3D][Q222], q[UUL_3D][Q222], q[UUU_3D][Q222]);
    divQy[4] = calc_sym_grady_3D(dy, q[LLL_3D][Q223], q[LLU_3D][Q223], q[LUL_3D][Q223], q[LUU_3D][Q223],
                                     q[ULL_3D][Q223], q[ULU_3D][Q223], q[UUL_3D][Q223], q[UUU_3D][Q223]);
    divQy[5] = calc_sym_grady_3D(dy, q[LLL_3D][Q233], q[LLU_3D][Q233], q[LUL_3D][Q233], q[LUU_3D][Q233],
                                     q[ULL_3D][Q233], q[ULU_3D][Q233], q[UUL_3D][Q233], q[UUU_3D][Q233]);

    divQz[0] = calc_sym_gradz_3D(dz, q[LLL_3D][Q113], q[LLU_3D][Q113], q[LUL_3D][Q113], q[LUU_3D][Q113],
                                     q[ULL_3D][Q113], q[ULU_3D][Q113], q[UUL_3D][Q113], q[UUU_3D][Q113]);
    divQz[1] = calc_sym_gradz_3D(dz, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQz[2] = calc_sym_gradz_3D(dz, q[LLL_3D][Q133], q[LLU_3D][Q133], q[LUL_3D][Q133], q[LUU_3D][Q133],
                                     q[ULL_3D][Q133], q[ULU_3D][Q133], q[UUL_3D][Q133], q[UUU_3D][Q133]);
    divQz[3] = calc_sym_gradz_3D(dz, q[LLL_3D][Q223], q[LLU_3D][Q223], q[LUL_3D][Q223], q[LUU_3D][Q223],
                                     q[ULL_3D][Q223], q[ULU_3D][Q223], q[UUL_3D][Q223], q[UUU_3D][Q223]);
    divQz[4] = calc_sym_gradz_3D(dz, q[LLL_3D][Q233], q[LLU_3D][Q233], q[LUL_3D][Q233], q[LUU_3D][Q233],
                                     q[ULL_3D][Q233], q[ULU_3D][Q233], q[UUL_3D][Q233], q[UUU_3D][Q233]);
    divQz[5] = calc_sym_gradz_3D(dz, q[LLL_3D][Q333], q[LLU_3D][Q333], q[LUL_3D][Q333], q[LUU_3D][Q333],
                                     q[ULL_3D][Q333], q[ULU_3D][Q333], q[UUL_3D][Q333], q[UUU_3D][Q333]);
  }
  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = divQx[0] + divQy[0] + divQz[0];
  rhs[P12] = divQx[1] + divQy[1] + divQz[1];
  rhs[P13] = divQx[2] + divQy[2] + divQz[2];
  rhs[P22] = divQx[3] + divQy[3] + divQz[3];
  rhs[P23] = divQx[4] + divQy[4] + divQz[4];
  rhs[P33] = divQx[5] + divQy[5] + divQz[5];
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
  struct gkyl_array *cflrate, double dt, struct gkyl_array *temp_grad,
  struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };

  double cfla = 0.0, cfl = gces->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // delibrately a double

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(update_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(heat_flux_range, offsets_centers);

  const double* fluid_d[sz[ndim-1]];
  const double* em_tot_d[sz[ndim-1]];
  double *temp_grad_d;
  const double* temp_grad_up[sz[ndim-1]];
  double *rhs_d;

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_range);
  int counter = 0;
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_range, iter_vertex.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
    }

    temp_grad_d = gkyl_array_fetch(temp_grad, linc_vertex);
    
    cfla = calc_unmag_heat_flux(gces, fluid_d, gkyl_array_fetch(cflrate, linc_center),
      temp_grad_d, cfla, dt);
    counter = counter + 1;
  }

  if (cfla > cflm)
    is_cfl_violated = 1.0;

  counter = 0;
  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, update_range);
  while (gkyl_range_iter_next(&iter_center)) {
    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_center.idx);
    long linc_center = gkyl_range_idx(update_range, iter_center.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      temp_grad_up[i] = gkyl_array_fetch(temp_grad, linc_vertex + offsets_centers[i]);
    }

    rhs_d = gkyl_array_fetch(rhs, linc_center);

    calc_grad_closure_update(gces, temp_grad_up, rhs_d);
    counter = counter + 1;
  }

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
