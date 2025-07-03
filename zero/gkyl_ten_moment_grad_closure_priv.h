#include <gkyl_moment_non_ideal_priv.h>

typedef double (*heat_flux_calc_t)(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate, double cfl,
  double dt, double *rhs_d[]);

struct gkyl_ten_moment_grad_closure {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  double k0; // damping coefficient
  double omega;
  double cfl; // CFL number to use
  bool mag;
  struct gkyl_comm *comm;

  heat_flux_calc_t calc_q;
};

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

static void
var_setup(const gkyl_ten_moment_grad_closure *gces,
  int start, int end,
  const double *fluid_d[],
  double rho[], double p[], double Tij[][6])
{
  for (int j = start; j <= end; ++j) {
    rho[j] = fluid_d[j][RHO];
    p[j] = (fluid_d[j][P11] - fluid_d[j][MX]*fluid_d[j][MX]/fluid_d[j][RHO]
          + fluid_d[j][P22] - fluid_d[j][MY]*fluid_d[j][MY]/fluid_d[j][RHO]
          + fluid_d[j][P33] - fluid_d[j][MZ]*fluid_d[j][MZ]/fluid_d[j][RHO])/3.0;
    Tij[j][T11] = (fluid_d[j][P11] - fluid_d[j][MX]*fluid_d[j][MX]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T12] = (fluid_d[j][P12] - fluid_d[j][MX]*fluid_d[j][MY]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T13] = (fluid_d[j][P13] - fluid_d[j][MX]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T22] = (fluid_d[j][P22] - fluid_d[j][MY]*fluid_d[j][MY]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T23] = (fluid_d[j][P23] - fluid_d[j][MY]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T33] = (fluid_d[j][P33] - fluid_d[j][MZ]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
  }
}

static inline double
calc_sym_grad_limiter_2D(double alpha, double a, double b)
{
  double avg = (a + b)/2;
  double min = fmin(alpha*a, a/alpha);
  double max = fmax(alpha*a, a/alpha);
  if (avg <= min) {
    return min;
  } else if (avg >= max) {
    return max;
  } else {
    return avg;
  }
}

static inline double
calc_sym_grad_limiter_3D(double alpha, double a, double b, double c, double d)
{
  double avg = (a + b + c + d)/4;
  double min = fmin(alpha*a, a/alpha);
  double max = fmax(alpha*a, a/alpha);
  if (avg <= min) {
    return min;
  } else if (avg >= max) {
    return max;
  } else {
    return avg;
  }
}

static void
calc_temp_grad_1D(const double Tij[][6], const double dx, double dTdx[6])
{
  dTdx[T11] = calc_sym_grad_1D(dx, Tij[L_1D][T11], Tij[U_1D][T11]);
  dTdx[T12] = calc_sym_grad_1D(dx, Tij[L_1D][T12], Tij[U_1D][T12]);
  dTdx[T13] = calc_sym_grad_1D(dx, Tij[L_1D][T13], Tij[U_1D][T13]);
  dTdx[T22] = calc_sym_grad_1D(dx, Tij[L_1D][T22], Tij[U_1D][T22]);
  dTdx[T23] = calc_sym_grad_1D(dx, Tij[L_1D][T23], Tij[U_1D][T23]);
  dTdx[T33] = calc_sym_grad_1D(dx, Tij[L_1D][T33], Tij[U_1D][T33]);
}

static void
calc_temp_grad_2D(const double Tij[][6], const double dx, const double dy,
  double dTdx[][6], double dTdy[][6])
{
  double dTx[2][6] = {0.0};
  double dTy[2][6] = {0.0};

  double limit = 0.75;

  dTx[L_1D][T11] = calc_sym_grad_1D(dx, Tij[LL_2D][T11], Tij[UL_2D][T11]);
  dTx[L_1D][T12] = calc_sym_grad_1D(dx, Tij[LL_2D][T12], Tij[UL_2D][T12]);
  dTx[L_1D][T13] = calc_sym_grad_1D(dx, Tij[LL_2D][T13], Tij[UL_2D][T13]);
  dTx[L_1D][T22] = calc_sym_grad_1D(dx, Tij[LL_2D][T22], Tij[UL_2D][T22]);
  dTx[L_1D][T23] = calc_sym_grad_1D(dx, Tij[LL_2D][T23], Tij[UL_2D][T23]);
  dTx[L_1D][T33] = calc_sym_grad_1D(dx, Tij[LL_2D][T33], Tij[UL_2D][T33]);

  dTx[U_1D][T11] = calc_sym_grad_1D(dx, Tij[LU_2D][T11], Tij[UU_2D][T11]);
  dTx[U_1D][T12] = calc_sym_grad_1D(dx, Tij[LU_2D][T12], Tij[UU_2D][T12]);
  dTx[U_1D][T13] = calc_sym_grad_1D(dx, Tij[LU_2D][T13], Tij[UU_2D][T13]);
  dTx[U_1D][T22] = calc_sym_grad_1D(dx, Tij[LU_2D][T22], Tij[UU_2D][T22]);
  dTx[U_1D][T23] = calc_sym_grad_1D(dx, Tij[LU_2D][T23], Tij[UU_2D][T23]);
  dTx[U_1D][T33] = calc_sym_grad_1D(dx, Tij[LU_2D][T33], Tij[UU_2D][T33]);

  dTdx[L_1D][T11] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T11], dTx[U_1D][T11]);
  dTdx[L_1D][T12] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T12], dTx[U_1D][T12]);
  dTdx[L_1D][T13] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T13], dTx[U_1D][T13]);
  dTdx[L_1D][T22] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T22], dTx[U_1D][T22]);
  dTdx[L_1D][T23] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T23], dTx[U_1D][T23]);
  dTdx[L_1D][T33] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T33], dTx[U_1D][T33]);

  dTdx[U_1D][T11] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T11], dTx[L_1D][T11]);
  dTdx[U_1D][T12] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T12], dTx[L_1D][T12]);
  dTdx[U_1D][T13] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T13], dTx[L_1D][T13]);
  dTdx[U_1D][T22] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T22], dTx[L_1D][T22]);
  dTdx[U_1D][T23] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T23], dTx[L_1D][T23]);
  dTdx[U_1D][T33] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T33], dTx[L_1D][T33]);

  dTy[L_1D][T11] = calc_sym_grad_1D(dy, Tij[LL_2D][T11], Tij[LU_2D][T11]);
  dTy[L_1D][T12] = calc_sym_grad_1D(dy, Tij[LL_2D][T12], Tij[LU_2D][T12]);
  dTy[L_1D][T13] = calc_sym_grad_1D(dy, Tij[LL_2D][T13], Tij[LU_2D][T13]);
  dTy[L_1D][T22] = calc_sym_grad_1D(dy, Tij[LL_2D][T22], Tij[LU_2D][T22]);
  dTy[L_1D][T23] = calc_sym_grad_1D(dy, Tij[LL_2D][T23], Tij[LU_2D][T23]);
  dTy[L_1D][T33] = calc_sym_grad_1D(dy, Tij[LL_2D][T33], Tij[LU_2D][T33]);

  dTy[U_1D][T11] = calc_sym_grad_1D(dy, Tij[UL_2D][T11], Tij[UU_2D][T11]);
  dTy[U_1D][T12] = calc_sym_grad_1D(dy, Tij[UL_2D][T12], Tij[UU_2D][T12]);
  dTy[U_1D][T13] = calc_sym_grad_1D(dy, Tij[UL_2D][T13], Tij[UU_2D][T13]);
  dTy[U_1D][T22] = calc_sym_grad_1D(dy, Tij[UL_2D][T22], Tij[UU_2D][T22]);
  dTy[U_1D][T23] = calc_sym_grad_1D(dy, Tij[UL_2D][T23], Tij[UU_2D][T23]);
  dTy[U_1D][T33] = calc_sym_grad_1D(dy, Tij[UL_2D][T33], Tij[UU_2D][T33]);

  dTdy[L_1D][T11] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T11], dTy[U_1D][T11]);
  dTdy[L_1D][T12] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T12], dTy[U_1D][T12]);
  dTdy[L_1D][T13] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T13], dTy[U_1D][T13]);
  dTdy[L_1D][T22] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T22], dTy[U_1D][T22]);
  dTdy[L_1D][T23] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T23], dTy[U_1D][T23]);
  dTdy[L_1D][T33] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T33], dTy[U_1D][T33]);

  dTdy[U_1D][T11] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T11], dTy[L_1D][T11]);
  dTdy[U_1D][T12] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T12], dTy[L_1D][T12]);
  dTdy[U_1D][T13] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T13], dTy[L_1D][T13]);
  dTdy[U_1D][T22] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T22], dTy[L_1D][T22]);
  dTdy[U_1D][T23] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T23], dTy[L_1D][T23]);
  dTdy[U_1D][T33] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T33], dTy[L_1D][T33]);
}

static void
calc_temp_grad_3D(const double Tij[][6], const double dx, const double dy,
  const double dz, double dTdx[][6], double dTdy[][6], double dTdz[][6])
{
  double dTx[4][6] = {0.0};
  double dTy[4][6] = {0.0};
  double dTz[4][6] = {0.0};

  double limit = 0.75;

  dTx[LL_2D][T11] = calc_sym_grad_1D(dx, Tij[LLL_3D][T11], Tij[ULL_3D][T11]);
  dTx[LL_2D][T12] = calc_sym_grad_1D(dx, Tij[LLL_3D][T12], Tij[ULL_3D][T12]);
  dTx[LL_2D][T13] = calc_sym_grad_1D(dx, Tij[LLL_3D][T13], Tij[ULL_3D][T13]);
  dTx[LL_2D][T22] = calc_sym_grad_1D(dx, Tij[LLL_3D][T22], Tij[ULL_3D][T22]);
  dTx[LL_2D][T23] = calc_sym_grad_1D(dx, Tij[LLL_3D][T23], Tij[ULL_3D][T23]);
  dTx[LL_2D][T33] = calc_sym_grad_1D(dx, Tij[LLL_3D][T33], Tij[ULL_3D][T33]);

  dTx[LU_2D][T11] = calc_sym_grad_1D(dx, Tij[LUL_3D][T11], Tij[UUL_3D][T11]);
  dTx[LU_2D][T12] = calc_sym_grad_1D(dx, Tij[LUL_3D][T12], Tij[UUL_3D][T12]);
  dTx[LU_2D][T13] = calc_sym_grad_1D(dx, Tij[LUL_3D][T13], Tij[UUL_3D][T13]);
  dTx[LU_2D][T22] = calc_sym_grad_1D(dx, Tij[LUL_3D][T22], Tij[UUL_3D][T22]);
  dTx[LU_2D][T23] = calc_sym_grad_1D(dx, Tij[LUL_3D][T23], Tij[UUL_3D][T23]);
  dTx[LU_2D][T33] = calc_sym_grad_1D(dx, Tij[LUL_3D][T33], Tij[UUL_3D][T33]);
  
  dTx[UL_2D][T11] = calc_sym_grad_1D(dx, Tij[LLU_3D][T11], Tij[ULU_3D][T11]);
  dTx[UL_2D][T12] = calc_sym_grad_1D(dx, Tij[LLU_3D][T12], Tij[ULU_3D][T12]);
  dTx[UL_2D][T13] = calc_sym_grad_1D(dx, Tij[LLU_3D][T13], Tij[ULU_3D][T13]);
  dTx[UL_2D][T22] = calc_sym_grad_1D(dx, Tij[LLU_3D][T22], Tij[ULU_3D][T22]);
  dTx[UL_2D][T23] = calc_sym_grad_1D(dx, Tij[LLU_3D][T23], Tij[ULU_3D][T23]);
  dTx[UL_2D][T33] = calc_sym_grad_1D(dx, Tij[LLU_3D][T33], Tij[ULU_3D][T33]);
  
  dTx[UU_2D][T11] = calc_sym_grad_1D(dx, Tij[LUU_3D][T11], Tij[UUU_3D][T11]);
  dTx[UU_2D][T12] = calc_sym_grad_1D(dx, Tij[LUU_3D][T12], Tij[UUU_3D][T12]);
  dTx[UU_2D][T13] = calc_sym_grad_1D(dx, Tij[LUU_3D][T13], Tij[UUU_3D][T13]);
  dTx[UU_2D][T22] = calc_sym_grad_1D(dx, Tij[LUU_3D][T22], Tij[UUU_3D][T22]);
  dTx[UU_2D][T23] = calc_sym_grad_1D(dx, Tij[LUU_3D][T23], Tij[UUU_3D][T23]);  
  dTx[UU_2D][T33] = calc_sym_grad_1D(dx, Tij[LUU_3D][T33], Tij[UUU_3D][T33]);

  dTdx[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T11], dTx[LU_2D][T11],
    dTx[UL_2D][T11], dTx[UU_2D][T11]);
  dTdx[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T12], dTx[LU_2D][T12],
    dTx[UL_2D][T12], dTx[UU_2D][T12]);
  dTdx[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T13], dTx[LU_2D][T13],
    dTx[UL_2D][T13], dTx[UU_2D][T13]);
  dTdx[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T22], dTx[LU_2D][T22],
    dTx[UL_2D][T22], dTx[UU_2D][T22]);
  dTdx[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T23], dTx[LU_2D][T23],
    dTx[UL_2D][T23], dTx[UU_2D][T23]);
  dTdx[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T33], dTx[LU_2D][T33],
    dTx[UL_2D][T33], dTx[UU_2D][T33]);

  dTdx[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T11], dTx[LL_2D][T11],
    dTx[UL_2D][T11], dTx[UU_2D][T11]);
  dTdx[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T12], dTx[LL_2D][T12],
    dTx[UL_2D][T12], dTx[UU_2D][T12]);
  dTdx[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T13], dTx[LL_2D][T13],
    dTx[UL_2D][T13], dTx[UU_2D][T13]);
  dTdx[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T22], dTx[LL_2D][T22],
    dTx[UL_2D][T22], dTx[UU_2D][T22]);
  dTdx[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T23], dTx[LL_2D][T23],
    dTx[UL_2D][T23], dTx[UU_2D][T23]);
  dTdx[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T33], dTx[LL_2D][T33],
    dTx[UL_2D][T33], dTx[UU_2D][T33]);

  dTdx[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T11], dTx[LL_2D][T11],
    dTx[LU_2D][T11], dTx[UU_2D][T11]);
  dTdx[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T12], dTx[LL_2D][T12],
    dTx[LU_2D][T12], dTx[UU_2D][T12]);
  dTdx[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T13], dTx[LL_2D][T13],
    dTx[LU_2D][T13], dTx[UU_2D][T13]);
  dTdx[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T22], dTx[LL_2D][T22],
    dTx[LU_2D][T22], dTx[UU_2D][T22]);
  dTdx[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T23], dTx[LL_2D][T23],
    dTx[LU_2D][T23], dTx[UU_2D][T23]);
  dTdx[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T33], dTx[LL_2D][T33],
    dTx[LU_2D][T33], dTx[UU_2D][T33]);

  dTdx[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T11], dTx[LL_2D][T11],
    dTx[LU_2D][T11], dTx[UU_2D][T11]);
  dTdx[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T12], dTx[LL_2D][T12],
    dTx[LU_2D][T12], dTx[UU_2D][T12]);
  dTdx[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T13], dTx[LL_2D][T13],
    dTx[LU_2D][T13], dTx[UU_2D][T13]);
  dTdx[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T22], dTx[LL_2D][T22],
    dTx[LU_2D][T22], dTx[UU_2D][T22]);
  dTdx[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T23], dTx[LL_2D][T23],
    dTx[LU_2D][T23], dTx[UU_2D][T23]);
  dTdx[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T33], dTx[LL_2D][T33],
    dTx[LU_2D][T33], dTx[UU_2D][T33]);

  dTy[LL_2D][T11] = calc_sym_grad_1D(dy, Tij[LLL_3D][T11], Tij[LUL_3D][T11]);
  dTy[LL_2D][T12] = calc_sym_grad_1D(dy, Tij[LLL_3D][T12], Tij[LUL_3D][T12]);
  dTy[LL_2D][T13] = calc_sym_grad_1D(dy, Tij[LLL_3D][T13], Tij[LUL_3D][T13]);
  dTy[LL_2D][T22] = calc_sym_grad_1D(dy, Tij[LLL_3D][T22], Tij[LUL_3D][T22]);
  dTy[LL_2D][T23] = calc_sym_grad_1D(dy, Tij[LLL_3D][T23], Tij[LUL_3D][T23]);
  dTy[LL_2D][T33] = calc_sym_grad_1D(dy, Tij[LLL_3D][T33], Tij[LUL_3D][T33]);

  dTy[LU_2D][T11] = calc_sym_grad_1D(dy, Tij[ULL_3D][T11], Tij[UUL_3D][T11]);
  dTy[LU_2D][T12] = calc_sym_grad_1D(dy, Tij[ULL_3D][T12], Tij[UUL_3D][T12]);
  dTy[LU_2D][T13] = calc_sym_grad_1D(dy, Tij[ULL_3D][T13], Tij[UUL_3D][T13]);
  dTy[LU_2D][T22] = calc_sym_grad_1D(dy, Tij[ULL_3D][T22], Tij[UUL_3D][T22]);
  dTy[LU_2D][T23] = calc_sym_grad_1D(dy, Tij[ULL_3D][T23], Tij[UUL_3D][T23]);
  dTy[LU_2D][T33] = calc_sym_grad_1D(dy, Tij[ULL_3D][T33], Tij[UUL_3D][T33]);

  dTy[UL_2D][T11] = calc_sym_grad_1D(dy, Tij[LLU_3D][T11], Tij[LUU_3D][T11]);
  dTy[UL_2D][T12] = calc_sym_grad_1D(dy, Tij[LLU_3D][T12], Tij[LUU_3D][T12]);
  dTy[UL_2D][T13] = calc_sym_grad_1D(dy, Tij[LLU_3D][T13], Tij[LUU_3D][T13]);
  dTy[UL_2D][T22] = calc_sym_grad_1D(dy, Tij[LLU_3D][T22], Tij[LUU_3D][T22]);
  dTy[UL_2D][T23] = calc_sym_grad_1D(dy, Tij[LLU_3D][T23], Tij[LUU_3D][T23]);
  dTy[UL_2D][T33] = calc_sym_grad_1D(dy, Tij[LLU_3D][T33], Tij[LUU_3D][T33]);

  dTy[UU_2D][T11] = calc_sym_grad_1D(dy, Tij[ULU_3D][T11], Tij[UUU_3D][T11]);
  dTy[UU_2D][T12] = calc_sym_grad_1D(dy, Tij[ULU_3D][T12], Tij[UUU_3D][T12]);
  dTy[UU_2D][T13] = calc_sym_grad_1D(dy, Tij[ULU_3D][T13], Tij[UUU_3D][T13]);
  dTy[UU_2D][T22] = calc_sym_grad_1D(dy, Tij[ULU_3D][T22], Tij[UUU_3D][T22]);
  dTy[UU_2D][T23] = calc_sym_grad_1D(dy, Tij[ULU_3D][T23], Tij[UUU_3D][T23]);
  dTy[UU_2D][T33] = calc_sym_grad_1D(dy, Tij[ULU_3D][T33], Tij[UUU_3D][T33]);

  dTdy[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T11], dTy[LU_2D][T11],
    dTy[UL_2D][T11], dTy[UU_2D][T11]);
  dTdy[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T12], dTy[LU_2D][T12],
    dTy[UL_2D][T12], dTy[UU_2D][T12]);
  dTdy[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T13], dTy[LU_2D][T13],
    dTy[UL_2D][T13], dTy[UU_2D][T13]);
  dTdy[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T22], dTy[LU_2D][T22],
    dTy[UL_2D][T22], dTy[UU_2D][T22]);
  dTdy[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T23], dTy[LU_2D][T23],
    dTy[UL_2D][T23], dTy[UU_2D][T23]);
  dTdy[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T33], dTy[LU_2D][T33],
    dTy[UL_2D][T33], dTy[UU_2D][T33]);

  dTdy[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T11], dTy[LL_2D][T11],
    dTy[UL_2D][T11], dTy[UU_2D][T11]);
  dTdy[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T12], dTy[LL_2D][T12],
    dTy[UL_2D][T12], dTy[UU_2D][T12]);
  dTdy[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T13], dTy[LL_2D][T13],
    dTy[UL_2D][T13], dTy[UU_2D][T13]);
  dTdy[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T22], dTy[LL_2D][T22],
    dTy[UL_2D][T22], dTy[UU_2D][T22]);
  dTdy[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T23], dTy[LL_2D][T23],
    dTy[UL_2D][T23], dTy[UU_2D][T23]);
  dTdy[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T33], dTy[LL_2D][T33],
    dTy[UL_2D][T33], dTy[UU_2D][T33]);

  dTdy[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T11], dTy[LL_2D][T11],
    dTy[LU_2D][T11], dTy[UU_2D][T11]);
  dTdy[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T12], dTy[LL_2D][T12],
    dTy[LU_2D][T12], dTy[UU_2D][T12]);
  dTdy[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T13], dTy[LL_2D][T13],
    dTy[LU_2D][T13], dTy[UU_2D][T13]);
  dTdy[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T22], dTy[LL_2D][T22],
    dTy[LU_2D][T22], dTy[UU_2D][T22]);
  dTdy[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T23], dTy[LL_2D][T23],
    dTy[LU_2D][T23], dTy[UU_2D][T23]);
  dTdy[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T33], dTy[LL_2D][T33],
    dTy[LU_2D][T33], dTy[UU_2D][T33]);

  dTdy[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T11], dTy[LL_2D][T11],
    dTy[LU_2D][T11], dTy[UU_2D][T11]);
  dTdy[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T12], dTy[LL_2D][T12],
    dTy[LU_2D][T12], dTy[UU_2D][T12]);
  dTdy[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T13], dTy[LL_2D][T13],
    dTy[LU_2D][T13], dTy[UU_2D][T13]);
  dTdy[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T22], dTy[LL_2D][T22],
    dTy[LU_2D][T22], dTy[UU_2D][T22]);
  dTdy[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T23], dTy[LL_2D][T23],
    dTy[LU_2D][T23], dTy[UU_2D][T23]);
  dTdy[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T33], dTy[LL_2D][T33],
    dTy[LU_2D][T33], dTy[UU_2D][T33]);

  dTz[LL_2D][T11] = calc_sym_grad_1D(dz, Tij[LLL_3D][T11], Tij[LLU_3D][T11]);
  dTz[LL_2D][T12] = calc_sym_grad_1D(dz, Tij[LLL_3D][T12], Tij[LLU_3D][T12]);
  dTz[LL_2D][T13] = calc_sym_grad_1D(dz, Tij[LLL_3D][T13], Tij[LLU_3D][T13]);
  dTz[LL_2D][T22] = calc_sym_grad_1D(dz, Tij[LLL_3D][T22], Tij[LLU_3D][T22]);
  dTz[LL_2D][T23] = calc_sym_grad_1D(dz, Tij[LLL_3D][T23], Tij[LLU_3D][T23]);
  dTz[LL_2D][T33] = calc_sym_grad_1D(dz, Tij[LLL_3D][T33], Tij[LLU_3D][T33]);

  dTz[LU_2D][T11] = calc_sym_grad_1D(dz, Tij[ULL_3D][T11], Tij[ULU_3D][T11]);
  dTz[LU_2D][T12] = calc_sym_grad_1D(dz, Tij[ULL_3D][T12], Tij[ULU_3D][T12]);
  dTz[LU_2D][T13] = calc_sym_grad_1D(dz, Tij[ULL_3D][T13], Tij[ULU_3D][T13]);
  dTz[LU_2D][T22] = calc_sym_grad_1D(dz, Tij[ULL_3D][T22], Tij[ULU_3D][T22]);
  dTz[LU_2D][T23] = calc_sym_grad_1D(dz, Tij[ULL_3D][T23], Tij[ULU_3D][T23]);
  dTz[LU_2D][T33] = calc_sym_grad_1D(dz, Tij[ULL_3D][T33], Tij[ULU_3D][T33]);

  dTz[UL_2D][T11] = calc_sym_grad_1D(dz, Tij[LUL_3D][T11], Tij[LUU_3D][T11]);
  dTz[UL_2D][T12] = calc_sym_grad_1D(dz, Tij[LUL_3D][T12], Tij[LUU_3D][T12]);
  dTz[UL_2D][T13] = calc_sym_grad_1D(dz, Tij[LUL_3D][T13], Tij[LUU_3D][T13]);
  dTz[UL_2D][T22] = calc_sym_grad_1D(dz, Tij[LUL_3D][T22], Tij[LUU_3D][T22]);
  dTz[UL_2D][T23] = calc_sym_grad_1D(dz, Tij[LUL_3D][T23], Tij[LUU_3D][T23]);
  dTz[UL_2D][T33] = calc_sym_grad_1D(dz, Tij[LUL_3D][T33], Tij[LUU_3D][T33]);

  dTz[UU_2D][T11] = calc_sym_grad_1D(dz, Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
  dTz[UU_2D][T12] = calc_sym_grad_1D(dz, Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
  dTz[UU_2D][T13] = calc_sym_grad_1D(dz, Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
  dTz[UU_2D][T22] = calc_sym_grad_1D(dz, Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
  dTz[UU_2D][T23] = calc_sym_grad_1D(dz, Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
  dTz[UU_2D][T33] = calc_sym_grad_1D(dz, Tij[UUL_3D][T33], Tij[UUU_3D][T33]);

  dTdz[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T11], dTz[LU_2D][T11],
    dTz[UL_2D][T11], dTz[UU_2D][T11]);
  dTdz[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T12], dTz[LU_2D][T12],
    dTz[UL_2D][T12], dTz[UU_2D][T12]);
  dTdz[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T13], dTz[LU_2D][T13],
    dTz[UL_2D][T13], dTz[UU_2D][T13]);
  dTdz[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T22], dTz[LU_2D][T22],
    dTz[UL_2D][T22], dTz[UU_2D][T22]);
  dTdz[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T23], dTz[LU_2D][T23],
    dTz[UL_2D][T23], dTz[UU_2D][T23]);
  dTdz[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T33], dTz[LU_2D][T33],
    dTz[UL_2D][T33], dTz[UU_2D][T33]);

  dTdz[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T11], dTz[LL_2D][T11],
    dTz[UL_2D][T11], dTz[UU_2D][T11]);
  dTdz[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T12], dTz[LL_2D][T12],
    dTz[UL_2D][T12], dTz[UU_2D][T12]);
  dTdz[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T13], dTz[LL_2D][T13],
    dTz[UL_2D][T13], dTz[UU_2D][T13]);
  dTdz[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T22], dTz[LL_2D][T22],
    dTz[UL_2D][T22], dTz[UU_2D][T22]);
  dTdz[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T23], dTz[LL_2D][T23],
    dTz[UL_2D][T23], dTz[UU_2D][T23]);
  dTdz[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T33], dTz[LL_2D][T33],
    dTz[UL_2D][T33], dTz[UU_2D][T33]);

  dTdz[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T11], dTz[LL_2D][T11],
    dTz[LU_2D][T11], dTz[UU_2D][T11]);
  dTdz[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T12], dTz[LL_2D][T12],
    dTz[LU_2D][T12], dTz[UU_2D][T12]);
  dTdz[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T13], dTz[LL_2D][T13],
    dTz[LU_2D][T13], dTz[UU_2D][T13]);
  dTdz[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T22], dTz[LL_2D][T22],
    dTz[LU_2D][T22], dTz[UU_2D][T22]);
  dTdz[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T23], dTz[LL_2D][T23],
    dTz[LU_2D][T23], dTz[UU_2D][T23]);
  dTdz[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T33], dTz[LL_2D][T33],
    dTz[LU_2D][T33], dTz[UU_2D][T33]);

  dTdz[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T11], dTz[LL_2D][T11],
    dTz[LU_2D][T11], dTz[UU_2D][T11]);
  dTdz[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T12], dTz[LL_2D][T12],
    dTz[LU_2D][T12], dTz[UU_2D][T12]);
  dTdz[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T13], dTz[LL_2D][T13],
    dTz[LU_2D][T13], dTz[UU_2D][T13]);
  dTdz[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T22], dTz[LL_2D][T22],
    dTz[LU_2D][T22], dTz[UU_2D][T22]);
  dTdz[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T23], dTz[LL_2D][T23],
    dTz[LU_2D][T23], dTz[UU_2D][T23]);
  dTdz[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T33], dTz[LL_2D][T33],
    dTz[LU_2D][T33], dTz[UU_2D][T33]);
}

static void
calc_heat_grad_1D(const double q[10], const double dx, double *rhs_d[])
{
  double signx[2] = { 1.0, -1.0 };
  
  rhs_d[L_1D][P11] += signx[L_1D]*q[Q111]/dx;
  rhs_d[L_1D][P12] += signx[L_1D]*q[Q112]/dx;
  rhs_d[L_1D][P13] += signx[L_1D]*q[Q113]/dx;
  rhs_d[L_1D][P22] += signx[L_1D]*q[Q122]/dx;
  rhs_d[L_1D][P23] += signx[L_1D]*q[Q123]/dx;
  rhs_d[L_1D][P33] += signx[L_1D]*q[Q133]/dx;

  rhs_d[U_1D][P11] += signx[U_1D]*q[Q111]/dx;
  rhs_d[U_1D][P12] += signx[U_1D]*q[Q112]/dx;
  rhs_d[U_1D][P13] += signx[U_1D]*q[Q113]/dx;
  rhs_d[U_1D][P22] += signx[U_1D]*q[Q122]/dx;
  rhs_d[U_1D][P23] += signx[U_1D]*q[Q123]/dx;
  rhs_d[U_1D][P33] += signx[U_1D]*q[Q133]/dx;
}

static void
calc_heat_grad_2D(const double q[][10], const double dx, const double dy,
  double *rhs_d[])
{
  double signx[4] = { 1.0, 1.0, -1.0, -1.0 };
  double signy[4] = { 1.0, -1.0, 1.0, -1.0 };

  rhs_d[LL_2D][P11] += signx[LL_2D]*q[LL_2D][Q111]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q112]/(2.0*dy);
  rhs_d[LL_2D][P12] += signx[LL_2D]*q[LL_2D][Q112]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q122]/(2.0*dy);
  rhs_d[LL_2D][P13] += signx[LL_2D]*q[LL_2D][Q113]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q123]/(2.0*dy);
  rhs_d[LL_2D][P22] += signx[LL_2D]*q[LL_2D][Q122]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q222]/(2.0*dy);
  rhs_d[LL_2D][P23] += signx[LL_2D]*q[LL_2D][Q123]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q223]/(2.0*dy);
  rhs_d[LL_2D][P33] += signx[LL_2D]*q[LL_2D][Q133]/(2.0*dx)
    + signy[LL_2D]*q[LL_2D][Q233]/(2.0*dy);

  rhs_d[LU_2D][P11] += signx[LU_2D]*q[LU_2D][Q111]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q112]/(2.0*dy);
  rhs_d[LU_2D][P12] += signx[LU_2D]*q[LU_2D][Q112]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q122]/(2.0*dy);
  rhs_d[LU_2D][P13] += signx[LU_2D]*q[LU_2D][Q113]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q123]/(2.0*dy);
  rhs_d[LU_2D][P22] += signx[LU_2D]*q[LU_2D][Q122]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q222]/(2.0*dy);
  rhs_d[LU_2D][P23] += signx[LU_2D]*q[LU_2D][Q123]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q223]/(2.0*dy);
  rhs_d[LU_2D][P33] += signx[LU_2D]*q[LU_2D][Q133]/(2.0*dx)
    + signy[LU_2D]*q[LU_2D][Q233]/(2.0*dy);

  rhs_d[UL_2D][P11] += signx[UL_2D]*q[UL_2D][Q111]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q112]/(2.0*dy);
  rhs_d[UL_2D][P12] += signx[UL_2D]*q[UL_2D][Q112]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q122]/(2.0*dy);
  rhs_d[UL_2D][P13] += signx[UL_2D]*q[UL_2D][Q113]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q123]/(2.0*dy);
  rhs_d[UL_2D][P22] += signx[UL_2D]*q[UL_2D][Q122]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q222]/(2.0*dy);
  rhs_d[UL_2D][P23] += signx[UL_2D]*q[UL_2D][Q123]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q223]/(2.0*dy);
  rhs_d[UL_2D][P33] += signx[UL_2D]*q[UL_2D][Q133]/(2.0*dx)
    + signy[UL_2D]*q[UL_2D][Q233]/(2.0*dy);

  rhs_d[UU_2D][P11] += signx[UU_2D]*q[UU_2D][Q111]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q112]/(2.0*dy);
  rhs_d[UU_2D][P12] += signx[UU_2D]*q[UU_2D][Q112]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q122]/(2.0*dy);
  rhs_d[UU_2D][P13] += signx[UU_2D]*q[UU_2D][Q113]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q123]/(2.0*dy);
  rhs_d[UU_2D][P22] += signx[UU_2D]*q[UU_2D][Q122]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q222]/(2.0*dy);
  rhs_d[UU_2D][P23] += signx[UU_2D]*q[UU_2D][Q123]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q223]/(2.0*dy);
  rhs_d[UU_2D][P33] += signx[UU_2D]*q[UU_2D][Q133]/(2.0*dx)
    + signy[UU_2D]*q[UU_2D][Q233]/(2.0*dy);
}

static void
calc_heat_grad_3D(const double q[][10], const double dx, const double dy,
  const double dz, double *rhs_d[])
{
  double signx[8] = { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 };
  double signy[8] = { 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0 };
  double signz[8] = { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };

  rhs_d[LLL_3D][P11] += signx[LLL_3D]*q[LLL_3D][Q111]/(4.0*dx)
    + signy[LLL_3D]*q[LLL_3D][Q112]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q113]/(4.0*dz);
  rhs_d[LLL_3D][P12] += signx[LLL_3D]*q[LLL_3D][Q112]/(4.0*dx)
    + signy[LLL_3D]*q[LLL_3D][Q122]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q123]/(4.0*dz);
  rhs_d[LLL_3D][P13] += signx[LLL_3D]*q[LLL_3D][Q113]/(4.0*dx)
    + signy[LLL_3D]*q[LLL_3D][Q123]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q133]/(4.0*dz);
  rhs_d[LLL_3D][P22] += signx[LLL_3D]*q[LLL_3D][Q122]/(4.0*dx)
    + signy[LLL_3D]*q[LLL_3D][Q222]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q223]/(4.0*dz);
  rhs_d[LLL_3D][P23] += signx[LLL_3D]*q[LLL_3D][Q123]/(4.0*dx) 
    + signy[LLL_3D]*q[LLL_3D][Q223]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q233]/(4.0*dz);
  rhs_d[LLL_3D][P33] += signx[LLL_3D]*q[LLL_3D][Q133]/(4.0*dx)
    + signy[LLL_3D]*q[LLL_3D][Q233]/(4.0*dy)
    + signz[LLL_3D]*q[LLL_3D][Q333]/(4.0*dz);

  rhs_d[LLU_3D][P11] += signx[LLU_3D]*q[LLU_3D][Q111]/(4.0*dx)
    + signy[LLU_3D]*q[LLU_3D][Q112]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q113]/(4.0*dz);
  rhs_d[LLU_3D][P12] += signx[LLU_3D]*q[LLU_3D][Q112]/(4.0*dx)
    + signy[LLU_3D]*q[LLU_3D][Q122]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q123]/(4.0*dz);
  rhs_d[LLU_3D][P13] += signx[LLU_3D]*q[LLU_3D][Q113]/(4.0*dx)
    + signy[LLU_3D]*q[LLU_3D][Q123]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q133]/(4.0*dz);
  rhs_d[LLU_3D][P22] += signx[LLU_3D]*q[LLU_3D][Q122]/(4.0*dx)
    + signy[LLU_3D]*q[LLU_3D][Q222]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q223]/(4.0*dz);
  rhs_d[LLU_3D][P23] += signx[LLU_3D]*q[LLU_3D][Q123]/(4.0*dx) 
    + signy[LLU_3D]*q[LLU_3D][Q223]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q233]/(4.0*dz);
  rhs_d[LLU_3D][P33] += signx[LLU_3D]*q[LLU_3D][Q133]/(4.0*dx)
    + signy[LLU_3D]*q[LLU_3D][Q233]/(4.0*dy)
    + signz[LLU_3D]*q[LLU_3D][Q333]/(4.0*dz);

  rhs_d[LUL_3D][P11] += signx[LUL_3D]*q[LUL_3D][Q111]/(4.0*dx)
    + signy[LUL_3D]*q[LUL_3D][Q112]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q113]/(4.0*dz);
  rhs_d[LUL_3D][P12] += signx[LUL_3D]*q[LUL_3D][Q112]/(4.0*dx)
    + signy[LUL_3D]*q[LUL_3D][Q122]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q123]/(4.0*dz);
  rhs_d[LUL_3D][P13] += signx[LUL_3D]*q[LUL_3D][Q113]/(4.0*dx)
    + signy[LUL_3D]*q[LUL_3D][Q123]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q133]/(4.0*dz);
  rhs_d[LUL_3D][P22] += signx[LUL_3D]*q[LUL_3D][Q122]/(4.0*dx)
    + signy[LUL_3D]*q[LUL_3D][Q222]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q223]/(4.0*dz);
  rhs_d[LUL_3D][P23] += signx[LUL_3D]*q[LUL_3D][Q123]/(4.0*dx) 
    + signy[LUL_3D]*q[LUL_3D][Q223]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q233]/(4.0*dz);
  rhs_d[LUL_3D][P33] += signx[LUL_3D]*q[LUL_3D][Q133]/(4.0*dx)
    + signy[LUL_3D]*q[LUL_3D][Q233]/(4.0*dy)
    + signz[LUL_3D]*q[LUL_3D][Q333]/(4.0*dz);

  rhs_d[LUU_3D][P11] += signx[LUU_3D]*q[LUU_3D][Q111]/(4.0*dx)
    + signy[LUU_3D]*q[LUU_3D][Q112]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q113]/(4.0*dz);
  rhs_d[LUU_3D][P12] += signx[LUU_3D]*q[LUU_3D][Q112]/(4.0*dx)
    + signy[LUU_3D]*q[LUU_3D][Q122]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q123]/(4.0*dz);
  rhs_d[LUU_3D][P13] += signx[LUU_3D]*q[LUU_3D][Q113]/(4.0*dx)
    + signy[LUU_3D]*q[LUU_3D][Q123]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q133]/(4.0*dz);
  rhs_d[LUU_3D][P22] += signx[LUU_3D]*q[LUU_3D][Q122]/(4.0*dx)
    + signy[LUU_3D]*q[LUU_3D][Q222]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q223]/(4.0*dz);
  rhs_d[LUU_3D][P23] += signx[LUU_3D]*q[LUU_3D][Q123]/(4.0*dx) 
    + signy[LUU_3D]*q[LUU_3D][Q223]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q233]/(4.0*dz);
  rhs_d[LUU_3D][P33] += signx[LUU_3D]*q[LUU_3D][Q133]/(4.0*dx)
    + signy[LUU_3D]*q[LUU_3D][Q233]/(4.0*dy)
    + signz[LUU_3D]*q[LUU_3D][Q333]/(4.0*dz);

  rhs_d[ULL_3D][P11] += signx[ULL_3D]*q[ULL_3D][Q111]/(4.0*dx)
    + signy[ULL_3D]*q[ULL_3D][Q112]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q113]/(4.0*dz);
  rhs_d[ULL_3D][P12] += signx[ULL_3D]*q[ULL_3D][Q112]/(4.0*dx)
    + signy[ULL_3D]*q[ULL_3D][Q122]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q123]/(4.0*dz);
  rhs_d[ULL_3D][P13] += signx[ULL_3D]*q[ULL_3D][Q113]/(4.0*dx)
    + signy[ULL_3D]*q[ULL_3D][Q123]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q133]/(4.0*dz);
  rhs_d[ULL_3D][P22] += signx[ULL_3D]*q[ULL_3D][Q122]/(4.0*dx)
    + signy[ULL_3D]*q[ULL_3D][Q222]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q223]/(4.0*dz);
  rhs_d[ULL_3D][P23] += signx[ULL_3D]*q[ULL_3D][Q123]/(4.0*dx) 
    + signy[ULL_3D]*q[ULL_3D][Q223]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q233]/(4.0*dz);
  rhs_d[ULL_3D][P33] += signx[ULL_3D]*q[ULL_3D][Q133]/(4.0*dx)
    + signy[ULL_3D]*q[ULL_3D][Q233]/(4.0*dy)
    + signz[ULL_3D]*q[ULL_3D][Q333]/(4.0*dz);

  rhs_d[ULU_3D][P11] += signx[ULU_3D]*q[ULU_3D][Q111]/(4.0*dx)
    + signy[ULU_3D]*q[ULU_3D][Q112]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q113]/(4.0*dz);
  rhs_d[ULU_3D][P12] += signx[ULU_3D]*q[ULU_3D][Q112]/(4.0*dx)
    + signy[ULU_3D]*q[ULU_3D][Q122]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q123]/(4.0*dz);
  rhs_d[ULU_3D][P13] += signx[ULU_3D]*q[ULU_3D][Q113]/(4.0*dx)
    + signy[ULU_3D]*q[ULU_3D][Q123]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q133]/(4.0*dz);
  rhs_d[ULU_3D][P22] += signx[ULU_3D]*q[ULU_3D][Q122]/(4.0*dx)
    + signy[ULU_3D]*q[ULU_3D][Q222]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q223]/(4.0*dz);
  rhs_d[ULU_3D][P23] += signx[ULU_3D]*q[ULU_3D][Q123]/(4.0*dx) 
    + signy[ULU_3D]*q[ULU_3D][Q223]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q233]/(4.0*dz);
  rhs_d[ULU_3D][P33] += signx[ULU_3D]*q[ULU_3D][Q133]/(4.0*dx)
    + signy[ULU_3D]*q[ULU_3D][Q233]/(4.0*dy)
    + signz[ULU_3D]*q[ULU_3D][Q333]/(4.0*dz);

  rhs_d[UUL_3D][P11] += signx[UUL_3D]*q[UUL_3D][Q111]/(4.0*dx)
    + signy[UUL_3D]*q[UUL_3D][Q112]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q113]/(4.0*dz);
  rhs_d[UUL_3D][P12] += signx[UUL_3D]*q[UUL_3D][Q112]/(4.0*dx)
    + signy[UUL_3D]*q[UUL_3D][Q122]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q123]/(4.0*dz);
  rhs_d[UUL_3D][P13] += signx[UUL_3D]*q[UUL_3D][Q113]/(4.0*dx)
    + signy[UUL_3D]*q[UUL_3D][Q123]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q133]/(4.0*dz);
  rhs_d[UUL_3D][P22] += signx[UUL_3D]*q[UUL_3D][Q122]/(4.0*dx)
    + signy[UUL_3D]*q[UUL_3D][Q222]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q223]/(4.0*dz);
  rhs_d[UUL_3D][P23] += signx[UUL_3D]*q[UUL_3D][Q123]/(4.0*dx) 
    + signy[UUL_3D]*q[UUL_3D][Q223]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q233]/(4.0*dz);
  rhs_d[UUL_3D][P33] += signx[UUL_3D]*q[UUL_3D][Q133]/(4.0*dx)
    + signy[UUL_3D]*q[UUL_3D][Q233]/(4.0*dy)
    + signz[UUL_3D]*q[UUL_3D][Q333]/(4.0*dz);

  rhs_d[UUU_3D][P11] += signx[UUU_3D]*q[UUU_3D][Q111]/(4.0*dx)
    + signy[UUU_3D]*q[UUU_3D][Q112]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q113]/(4.0*dz);
  rhs_d[UUU_3D][P12] += signx[UUU_3D]*q[UUU_3D][Q112]/(4.0*dx)
    + signy[UUU_3D]*q[UUU_3D][Q122]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q123]/(4.0*dz);
  rhs_d[UUU_3D][P13] += signx[UUU_3D]*q[UUU_3D][Q113]/(4.0*dx)
    + signy[UUU_3D]*q[UUU_3D][Q123]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q133]/(4.0*dz);
  rhs_d[UUU_3D][P22] += signx[UUU_3D]*q[UUU_3D][Q122]/(4.0*dx)
    + signy[UUU_3D]*q[UUU_3D][Q222]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q223]/(4.0*dz);
  rhs_d[UUU_3D][P23] += signx[UUU_3D]*q[UUU_3D][Q123]/(4.0*dx) 
    + signy[UUU_3D]*q[UUU_3D][Q223]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q233]/(4.0*dz);
  rhs_d[UUU_3D][P33] += signx[UUU_3D]*q[UUU_3D][Q133]/(4.0*dx)
    + signy[UUU_3D]*q[UUU_3D][Q233]/(4.0*dy)
    + signz[UUU_3D]*q[UUU_3D][Q333]/(4.0*dz);
}

static void
calc_source_1D(double dTdx[6], double sx, double sy, double sz)
{
  sx = (3.0*dTdx[T11] + dTdx[T22] + dTdx[T33])/6.0;
  sy = (dTdx[T12])/3.0;
  sz = (dTdx[T13])/3.0;
}

static void
calc_source_2D(double dTdx[][6], double dTdy[][6], double sx[], double sy[], double sz[])
{
  int compx[4] = { L_1D, U_1D, L_1D, U_1D };
  int compy[4] = { L_1D, L_1D, U_1D, U_1D };

  sx[LL_2D] = (3.0*dTdx[compx[LL_2D]][T11] + dTdx[compx[LL_2D]][T22]
    + dTdx[compx[LL_2D]][T33] + 2.0*dTdy[compy[LL_2D]][T12])/6.0;
  sy[LL_2D] = (3.0*dTdy[compy[LL_2D]][T22] + dTdy[compx[LL_2D]][T11]
    + dTdy[compy[LL_2D]][T33] + 2.0*dTdx[compx[LL_2D]][T12])/6.0;
  sz[LL_2D] = (2.0*dTdx[compx[LL_2D]][T13] + 2.0*dTdy[compy[LL_2D]][T23])/2.0;

  sx[LU_2D] = (3.0*dTdx[compx[LU_2D]][T11] + dTdx[compx[LU_2D]][T22]
    + dTdx[compx[LU_2D]][T33] + 2.0*dTdy[compy[LU_2D]][T12])/6.0;
  sy[LU_2D] = (3.0*dTdy[compy[LU_2D]][T22] + dTdy[compx[LU_2D]][T11]
    + dTdy[compy[LU_2D]][T33] + 2.0*dTdx[compx[LU_2D]][T12])/6.0;
  sz[LU_2D] = (2.0*dTdx[compx[LU_2D]][T13] + 2.0*dTdy[compy[LU_2D]][T23])/2.0;

  sx[LL_2D] = (3.0*dTdx[compx[LL_2D]][T11] + dTdx[compx[LL_2D]][T22]
    + dTdx[compx[LL_2D]][T33] + 2.0*dTdy[compy[LL_2D]][T12])/6.0;
  sy[LL_2D] = (3.0*dTdy[compy[LL_2D]][T22] + dTdy[compx[LL_2D]][T11]
    + dTdy[compy[LL_2D]][T33] + 2.0*dTdx[compx[LL_2D]][T12])/6.0;
  sz[LL_2D] = (2.0*dTdx[compx[LL_2D]][T13] + 2.0*dTdy[compy[LL_2D]][T23])/2.0;

  sx[LL_2D] = (3.0*dTdx[compx[LL_2D]][T11] + dTdx[compx[LL_2D]][T22]
    + dTdx[compx[LL_2D]][T33] + 2.0*dTdy[compy[LL_2D]][T12])/6.0;
  sy[LL_2D] = (3.0*dTdy[compy[LL_2D]][T22] + dTdy[compx[LL_2D]][T11]
    + dTdy[compy[LL_2D]][T33] + 2.0*dTdx[compx[LL_2D]][T12])/6.0;
  sz[LL_2D] = (2.0*dTdx[compx[LL_2D]][T13] + 2.0*dTdy[compy[LL_2D]][T23])/2.0;
}

static void
calc_source_3D(double dTdx[][6], double dTdy[][6], double dTdz[][6], double sx[],
  double sy[], double sz[])
{
  int compx[8] = { LL_2D, UL_2D, LU_2D, UU_2D, LL_2D, UL_2D, LU_2D, UU_2D };
  int compy[8] = { LL_2D, UL_2D, LL_2D, UL_2D, LU_2D, UU_2D, LU_2D, UU_2D };
  int compz[8] = { LL_2D, LL_2D, UL_2D, UL_2D, LU_2D, LU_2D, UU_2D, UU_2D };

  sx[LLL_3D] = (3.0*dTdx[compx[LLL_3D]][T11] + dTdx[compx[LLL_3D]][T22]
    + dTdx[compx[LLL_3D]][T33] + 2.0*dTdy[compy[LLL_3D]][T12]
    + 2.0*dTdz[compz[LLL_3D]][T13])/6.0;
  sy[LLL_3D] = (dTdy[compy[LLL_3D]][T11] + 3.0*dTdy[compy[LLL_3D]][T22]
    + dTdy[compy[LLL_3D]][T33] + 2.0*dTdx[compx[LLL_3D]][T12]
    + 2.0*dTdz[compz[LLL_3D]][T23])/6.0;
  sz[LLL_3D] = (dTdz[compz[LLL_3D]][T11] + dTdz[compz[LLL_3D]][T22]
    + 3.0*dTdz[compz[LLL_3D]][T33] + 2.0*dTdx[compx[LLL_3D]][T13]
    + 2.0*dTdy[compy[LLL_3D]][T23])/6.0;

  sx[LLU_3D] = (3.0*dTdx[compx[LLU_3D]][T11] + dTdx[compx[LLU_3D]][T22]
    + dTdx[compx[LLU_3D]][T33] + 2.0*dTdy[compy[LLU_3D]][T12]
    + 2.0*dTdz[compz[LLU_3D]][T13])/6.0;
  sy[LLU_3D] = (dTdy[compy[LLU_3D]][T11] + 3.0*dTdy[compy[LLU_3D]][T22]
    + dTdy[compy[LLU_3D]][T33] + 2.0*dTdx[compx[LLU_3D]][T12]
    + 2.0*dTdz[compz[LLU_3D]][T23])/6.0;
  sz[LLU_3D] = (dTdz[compz[LLU_3D]][T11] + dTdz[compz[LLU_3D]][T22]
    + 3.0*dTdz[compz[LLU_3D]][T33] + 2.0*dTdx[compx[LLU_3D]][T13]
    + 2.0*dTdy[compy[LLU_3D]][T23])/6.0;

  sx[LUL_3D] = (3.0*dTdx[compx[LUL_3D]][T11] + dTdx[compx[LUL_3D]][T22]
    + dTdx[compx[LUL_3D]][T33] + 2.0*dTdy[compy[LUL_3D]][T12]
    + 2.0*dTdz[compz[LUL_3D]][T13])/6.0;
  sy[LUL_3D] = (dTdy[compy[LUL_3D]][T11] + 3.0*dTdy[compy[LUL_3D]][T22]
    + dTdy[compy[LUL_3D]][T33] + 2.0*dTdx[compx[LUL_3D]][T12]
    + 2.0*dTdz[compz[LUL_3D]][T23])/6.0;
  sz[LUL_3D] = (dTdz[compz[LUL_3D]][T11] + dTdz[compz[LUL_3D]][T22]
    + 3.0*dTdz[compz[LUL_3D]][T33] + 2.0*dTdx[compx[LUL_3D]][T13]
    + 2.0*dTdy[compy[LUL_3D]][T23])/6.0;

  sx[LUU_3D] = (3.0*dTdx[compx[LUU_3D]][T11] + dTdx[compx[LUU_3D]][T22]
    + dTdx[compx[LUU_3D]][T33] + 2.0*dTdy[compy[LUU_3D]][T12]
    + 2.0*dTdz[compz[LUU_3D]][T13])/6.0;
  sy[LUU_3D] = (dTdy[compy[LUU_3D]][T11] + 3.0*dTdy[compy[LUU_3D]][T22]
    + dTdy[compy[LUU_3D]][T33] + 2.0*dTdx[compx[LUU_3D]][T12]
    + 2.0*dTdz[compz[LUU_3D]][T23])/6.0;
  sz[LUU_3D] = (dTdz[compz[LUU_3D]][T11] + dTdz[compz[LUU_3D]][T22]
    + 3.0*dTdz[compz[LUU_3D]][T33] + 2.0*dTdx[compx[LUU_3D]][T13]
    + 2.0*dTdy[compy[LUU_3D]][T23])/6.0;

  sx[ULL_3D] = (3.0*dTdx[compx[ULL_3D]][T11] + dTdx[compx[ULL_3D]][T22]
    + dTdx[compx[ULL_3D]][T33] + 2.0*dTdy[compy[ULL_3D]][T12]
    + 2.0*dTdz[compz[ULL_3D]][T13])/6.0;
  sy[ULL_3D] = (dTdy[compy[ULL_3D]][T11] + 3.0*dTdy[compy[ULL_3D]][T22]
    + dTdy[compy[ULL_3D]][T33] + 2.0*dTdx[compx[ULL_3D]][T12]
    + 2.0*dTdz[compz[ULL_3D]][T23])/6.0;
  sz[ULL_3D] = (dTdz[compz[ULL_3D]][T11] + dTdz[compz[ULL_3D]][T22]
    + 3.0*dTdz[compz[ULL_3D]][T33] + 2.0*dTdx[compx[ULL_3D]][T13]
    + 2.0*dTdy[compy[ULL_3D]][T23])/6.0;

  sx[ULU_3D] = (3.0*dTdx[compx[ULU_3D]][T11] + dTdx[compx[ULU_3D]][T22]
    + dTdx[compx[ULU_3D]][T33] + 2.0*dTdy[compy[ULU_3D]][T12]
    + 2.0*dTdz[compz[ULU_3D]][T13])/6.0;
  sy[ULU_3D] = (dTdy[compy[ULU_3D]][T11] + 3.0*dTdy[compy[ULU_3D]][T22]
    + dTdy[compy[ULU_3D]][T33] + 2.0*dTdx[compx[ULU_3D]][T12]
    + 2.0*dTdz[compz[ULU_3D]][T23])/6.0;
  sz[ULU_3D] = (dTdz[compz[ULU_3D]][T11] + dTdz[compz[ULU_3D]][T22]
    + 3.0*dTdz[compz[ULU_3D]][T33] + 2.0*dTdx[compx[ULU_3D]][T13]
    + 2.0*dTdy[compy[ULU_3D]][T23])/6.0;

  sx[UUL_3D] = (3.0*dTdx[compx[UUL_3D]][T11] + dTdx[compx[UUL_3D]][T22]
    + dTdx[compx[UUL_3D]][T33] + 2.0*dTdy[compy[UUL_3D]][T12]
    + 2.0*dTdz[compz[UUL_3D]][T13])/6.0;
  sy[UUL_3D] = (dTdy[compy[UUL_3D]][T11] + 3.0*dTdy[compy[UUL_3D]][T22]
    + dTdy[compy[UUL_3D]][T33] + 2.0*dTdx[compx[UUL_3D]][T12]
    + 2.0*dTdz[compz[UUL_3D]][T23])/6.0;
  sz[UUL_3D] = (dTdz[compz[UUL_3D]][T11] + dTdz[compz[UUL_3D]][T22]
    + 3.0*dTdz[compz[UUL_3D]][T33] + 2.0*dTdx[compx[UUL_3D]][T13]
    + 2.0*dTdy[compy[UUL_3D]][T23])/6.0;

  sx[UUU_3D] = (3.0*dTdx[compx[UUU_3D]][T11] + dTdx[compx[UUU_3D]][T22]
    + dTdx[compx[UUU_3D]][T33] + 2.0*dTdy[compy[UUU_3D]][T12]
    + 2.0*dTdz[compz[UUU_3D]][T13])/6.0;
  sy[UUU_3D] = (dTdy[compy[UUU_3D]][T11] + 3.0*dTdy[compy[UUU_3D]][T22]
    + dTdy[compy[UUU_3D]][T33] + 2.0*dTdx[compx[UUU_3D]][T12]
    + 2.0*dTdz[compz[UUU_3D]][T23])/6.0;
  sz[UUU_3D] = (dTdz[compz[UUU_3D]][T11] + dTdz[compz[UUU_3D]][T22]
    + 3.0*dTdz[compz[UUU_3D]][T33] + 2.0*dTdx[compx[UUU_3D]][T13]
    + 2.0*dTdy[compy[UUU_3D]][T23])/6.0;
}

static void
calc_heat_flux_1D(double qx, double qy, double qz, double q[10])
{
  q[Q111] = 1.2*qx;
  q[Q112] = 0.4*qy;
  q[Q113] = 0.4*qz;
  q[Q122] = 0.4*qx;
  q[Q123] = 0.0;
  q[Q133] = 0.4*qx;
}

static void
calc_heat_flux_2D(double qx[], double qy[], double qz[], double q[][10])
{
  q[LL_2D][Q111] = 1.2*qx[LL_2D];
  q[LL_2D][Q112] = 0.4*qy[LL_2D];
  q[LL_2D][Q113] = 0.4*qz[LL_2D];
  q[LL_2D][Q122] = 0.4*qx[LL_2D];
  q[LL_2D][Q123] = 0.0;
  q[LL_2D][Q133] = 0.4*qx[LL_2D];
  q[LL_2D][Q222] = 1.2*qy[LL_2D];
  q[LL_2D][Q223] = 0.4*qz[LL_2D];
  q[LL_2D][Q233] = 0.4*qy[LL_2D];

  q[LU_2D][Q111] = 1.2*qx[LU_2D];
  q[LU_2D][Q112] = 0.4*qy[LU_2D];
  q[LU_2D][Q113] = 0.4*qz[LU_2D];
  q[LU_2D][Q122] = 0.4*qx[LU_2D];
  q[LU_2D][Q123] = 0.0;
  q[LU_2D][Q133] = 0.4*qx[LU_2D];
  q[LU_2D][Q222] = 1.2*qy[LU_2D];
  q[LU_2D][Q223] = 0.4*qz[LU_2D];
  q[LU_2D][Q233] = 0.4*qy[LU_2D];

  q[UL_2D][Q111] = 1.2*qx[UL_2D];
  q[UL_2D][Q112] = 0.4*qy[UL_2D];
  q[UL_2D][Q113] = 0.4*qz[UL_2D];
  q[UL_2D][Q122] = 0.4*qx[UL_2D];
  q[UL_2D][Q123] = 0.0;
  q[UL_2D][Q133] = 0.4*qx[UL_2D];
  q[UL_2D][Q222] = 1.2*qy[UL_2D];
  q[UL_2D][Q223] = 0.4*qz[UL_2D];
  q[UL_2D][Q233] = 0.4*qy[UL_2D];

  q[UU_2D][Q111] = 1.2*qx[UU_2D];
  q[UU_2D][Q112] = 0.4*qy[UU_2D];
  q[UU_2D][Q113] = 0.4*qz[UU_2D];
  q[UU_2D][Q122] = 0.4*qx[UU_2D];
  q[UU_2D][Q123] = 0.0;
  q[UU_2D][Q133] = 0.4*qx[UU_2D];
  q[UU_2D][Q222] = 1.2*qy[UU_2D];
  q[UU_2D][Q223] = 0.4*qz[UU_2D];
  q[UU_2D][Q233] = 0.4*qy[UU_2D];
}

static void
calc_heat_flux_3D(double qx[], double qy[], double qz[], double q[][10])
{
  q[LLL_3D][Q111] = 1.2*qx[LLL_3D];
  q[LLL_3D][Q112] = 0.4*qy[LLL_3D];
  q[LLL_3D][Q113] = 0.4*qz[LLL_3D];
  q[LLL_3D][Q122] = 0.4*qx[LLL_3D];
  q[LLL_3D][Q123] = 0.0;
  q[LLL_3D][Q133] = 0.4*qx[LLL_3D];
  q[LLL_3D][Q222] = 1.2*qy[LLL_3D];
  q[LLL_3D][Q223] = 0.4*qz[LLL_3D];
  q[LLL_3D][Q233] = 0.4*qy[LLL_3D];
  q[LLL_3D][Q333] = 1.2*qx[LLL_3D];

  q[LLU_3D][Q111] = 1.2*qx[LLU_3D];
  q[LLU_3D][Q112] = 0.4*qy[LLU_3D];
  q[LLU_3D][Q113] = 0.4*qz[LLU_3D];
  q[LLU_3D][Q122] = 0.4*qx[LLU_3D];
  q[LLU_3D][Q123] = 0.0;
  q[LLU_3D][Q133] = 0.4*qx[LLU_3D];
  q[LLU_3D][Q222] = 1.2*qy[LLU_3D];
  q[LLU_3D][Q223] = 0.4*qz[LLU_3D];
  q[LLU_3D][Q233] = 0.4*qy[LLU_3D];
  q[LLU_3D][Q333] = 1.2*qx[LLU_3D];

  q[LUL_3D][Q111] = 1.2*qx[LUL_3D];
  q[LUL_3D][Q112] = 0.4*qy[LUL_3D];
  q[LUL_3D][Q113] = 0.4*qz[LUL_3D];
  q[LUL_3D][Q122] = 0.4*qx[LUL_3D];
  q[LUL_3D][Q123] = 0.0;
  q[LUL_3D][Q133] = 0.4*qx[LUL_3D];
  q[LUL_3D][Q222] = 1.2*qy[LUL_3D];
  q[LUL_3D][Q223] = 0.4*qz[LUL_3D];
  q[LUL_3D][Q233] = 0.4*qy[LUL_3D];
  q[LUL_3D][Q333] = 1.2*qx[LUL_3D];

  q[LUU_3D][Q111] = 1.2*qx[LUU_3D];
  q[LUU_3D][Q112] = 0.4*qy[LUU_3D];
  q[LUU_3D][Q113] = 0.4*qz[LUU_3D];
  q[LUU_3D][Q122] = 0.4*qx[LUU_3D];
  q[LUU_3D][Q123] = 0.0;
  q[LUU_3D][Q133] = 0.4*qx[LUU_3D];
  q[LUU_3D][Q222] = 1.2*qy[LUU_3D];
  q[LUU_3D][Q223] = 0.4*qz[LUU_3D];
  q[LUU_3D][Q233] = 0.4*qy[LUU_3D];
  q[LUU_3D][Q333] = 1.2*qx[LUU_3D];

  q[ULL_3D][Q111] = 1.2*qx[ULL_3D];
  q[ULL_3D][Q112] = 0.4*qy[ULL_3D];
  q[ULL_3D][Q113] = 0.4*qz[ULL_3D];
  q[ULL_3D][Q122] = 0.4*qx[ULL_3D];
  q[ULL_3D][Q123] = 0.0;
  q[ULL_3D][Q133] = 0.4*qx[ULL_3D];
  q[ULL_3D][Q222] = 1.2*qy[ULL_3D];
  q[ULL_3D][Q223] = 0.4*qz[ULL_3D];
  q[ULL_3D][Q233] = 0.4*qy[ULL_3D];
  q[ULL_3D][Q333] = 1.2*qx[ULL_3D];

  q[ULU_3D][Q111] = 1.2*qx[ULU_3D];
  q[ULU_3D][Q112] = 0.4*qy[ULU_3D];
  q[ULU_3D][Q113] = 0.4*qz[ULU_3D];
  q[ULU_3D][Q122] = 0.4*qx[ULU_3D];
  q[ULU_3D][Q123] = 0.0;
  q[ULU_3D][Q133] = 0.4*qx[ULU_3D];
  q[ULU_3D][Q222] = 1.2*qy[ULU_3D];
  q[ULU_3D][Q223] = 0.4*qz[ULU_3D];
  q[ULU_3D][Q233] = 0.4*qy[ULU_3D];
  q[ULU_3D][Q333] = 1.2*qx[ULU_3D];

  q[UUL_3D][Q111] = 1.2*qx[UUL_3D];
  q[UUL_3D][Q112] = 0.4*qy[UUL_3D];
  q[UUL_3D][Q113] = 0.4*qz[UUL_3D];
  q[UUL_3D][Q122] = 0.4*qx[UUL_3D];
  q[UUL_3D][Q123] = 0.0;
  q[UUL_3D][Q133] = 0.4*qx[UUL_3D];
  q[UUL_3D][Q222] = 1.2*qy[UUL_3D];
  q[UUL_3D][Q223] = 0.4*qz[UUL_3D];
  q[UUL_3D][Q233] = 0.4*qy[UUL_3D];
  q[UUL_3D][Q333] = 1.2*qx[UUL_3D];

  q[UUU_3D][Q111] = 1.2*qx[UUU_3D];
  q[UUU_3D][Q112] = 0.4*qy[UUU_3D];
  q[UUU_3D][Q113] = 0.4*qz[UUU_3D];
  q[UUU_3D][Q122] = 0.4*qx[UUU_3D];
  q[UUU_3D][Q123] = 0.0;
  q[UUU_3D][Q133] = 0.4*qx[UUU_3D];
  q[UUU_3D][Q222] = 1.2*qy[UUU_3D];
  q[UUU_3D][Q223] = 0.4*qz[UUU_3D];
  q[UUU_3D][Q233] = 0.4*qy[UUU_3D];
  q[UUU_3D][Q333] = 1.2*qx[UUU_3D];
}

double
calc_mag_heat_flux_1d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;

  const double dx = gces->grid.dx[0];
  double dTdx[6] = {0.0};
  double Tij[2][6] = {0.0};
  double rho[2] = {0.0};
  double p[2] = {0.0};
  double q[10] = {0.0};
  var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
  p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);

  calc_temp_grad_1D(Tij, dx, dTdx);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  double omega = gces->omega;

  double Bx = (em_tot_d[L_1D][BX] + em_tot_d[U_1D][BX])/2.0;
  double By = (em_tot_d[L_1D][BY] + em_tot_d[U_1D][BY])/2.0;
  double Bz = (em_tot_d[L_1D][BZ] + em_tot_d[U_1D][BZ])/2.0;

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;
  
  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi_par = alpha*vth_avg*rho_avg;
  double chi_perp = alpha*vth_avg*rho_avg/(1 + omega*omega);
  double chi_carat = alpha*vth_avg*rho_avg*omega/(1 + omega*omega);

  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  calc_source_1D(dTdx, sx, sy, sz);

  double qx = chi_perp*sx + (chi_par - chi_perp)*bx*(bx*sx + by*sy + bz*sz)
    - chi_carat*(by*sz - bz*sy);
  double qy = chi_perp*sy + (chi_par - chi_perp)*by*(bx*sx + by*sy + bz*sz)
    - chi_carat*(bx*sz - bz*sx);
  double qz = chi_perp*sz + (chi_par - chi_perp)*bz*(bx*sx + by*sy + bz*sz)
    - chi_carat*(bx*sy - by*sx);

  calc_heat_flux_1D(qx, qy, qz, q);

  calc_heat_grad_1D(q, dx, rhs_d);

  double cfla = dt/(dx*dx);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_mag_heat_flux_2d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  double dTdx[2][6] = {0.0};
  double dTdy[2][6] = {0.0};
  double Tij[4][6] = {0.0};
  double rho[4] = {0.0};
  double p[4] = {0.0};
  double sx[4] = {0.0};
  double sy[4] = {0.0};
  double sz[4] = {0.0};
  double qx[4] = {0.0};
  double qy[4] = {0.0};
  double qz[4] = {0.0};
  double q[4][10] = {0.0};
  var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
  p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);

  calc_temp_grad_2D(Tij, dx, dy, dTdx, dTdy);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  double omega = gces->omega;

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi_par = alpha*vth_avg*rho_avg;
  double chi_perp = alpha*vth_avg*rho_avg/(1 + omega*omega);
  double chi_carat = alpha*vth_avg*rho_avg*omega/(1 + omega*omega);

  double Bx = calc_harmonic_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX],
    em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
  double By = calc_harmonic_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY],
    em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
  double Bz = calc_harmonic_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ],
    em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;

  calc_source_2D(dTdx, dTdy, sx, sy, sz);

  qx[LL_2D] = chi_perp*sx[LL_2D] + (chi_par - chi_perp)*bx*(bx*sx[LL_2D]
    + by*sy[LL_2D] + bz*sz[LL_2D]) - chi_carat*(by*sz[LL_2D] - bz*sy[LL_2D]);
  qy[LL_2D] = chi_perp*sy[LL_2D] + (chi_par - chi_perp)*by*(bx*sx[LL_2D]
    + by*sy[LL_2D] + bz*sz[LL_2D]) - chi_carat*(bx*sz[LL_2D] - bz*sx[LL_2D]);
  qz[LL_2D] = chi_perp*sz[LL_2D] + (chi_par - chi_perp)*bz*(bx*sx[LL_2D]
    + by*sy[LL_2D] + bz*sz[LL_2D]) - chi_carat*(bx*sy[LL_2D] - by*sx[LL_2D]);

  qx[LU_2D] = chi_perp*sx[LU_2D] + (chi_par - chi_perp)*bx*(bx*sx[LU_2D]
    + by*sy[LU_2D] + bz*sz[LU_2D]) - chi_carat*(by*sz[LU_2D] - bz*sy[LU_2D]);
  qy[LU_2D] = chi_perp*sy[LU_2D] + (chi_par - chi_perp)*by*(bx*sx[LU_2D]
    + by*sy[LU_2D] + bz*sz[LU_2D]) - chi_carat*(bx*sz[LU_2D] - bz*sx[LU_2D]);
  qz[LU_2D] = chi_perp*sz[LU_2D] + (chi_par - chi_perp)*bz*(bx*sx[LU_2D]
    + by*sy[LU_2D] + bz*sz[LU_2D]) - chi_carat*(bx*sy[LU_2D] - by*sx[LU_2D]);

  qx[UL_2D] = chi_perp*sx[UL_2D] + (chi_par - chi_perp)*bx*(bx*sx[UL_2D]
    + by*sy[UL_2D] + bz*sz[UL_2D]) - chi_carat*(by*sz[UL_2D] - bz*sy[UL_2D]);
  qy[UL_2D] = chi_perp*sy[UL_2D] + (chi_par - chi_perp)*by*(bx*sx[UL_2D]
    + by*sy[UL_2D] + bz*sz[UL_2D]) - chi_carat*(bx*sz[UL_2D] - bz*sx[UL_2D]);
  qz[UL_2D] = chi_perp*sz[UL_2D] + (chi_par - chi_perp)*bz*(bx*sx[UL_2D]
    + by*sy[UL_2D] + bz*sz[UL_2D]) - chi_carat*(bx*sy[UL_2D] - by*sx[UL_2D]);

  qx[UU_2D] = chi_perp*sx[UU_2D] + (chi_par - chi_perp)*bx*(bx*sx[UU_2D]
    + by*sy[UU_2D] + bz*sz[UU_2D]) - chi_carat*(by*sz[UU_2D] - bz*sy[UU_2D]);
  qy[UU_2D] = chi_perp*sy[UU_2D] + (chi_par - chi_perp)*by*(bx*sx[UU_2D]
    + by*sy[UU_2D] + bz*sz[UU_2D]) - chi_carat*(bx*sz[UU_2D] - bz*sx[UU_2D]);
  qz[UU_2D] = chi_perp*sz[UU_2D] + (chi_par - chi_perp)*bz*(bx*sx[UU_2D]
    + by*sy[UU_2D] + bz*sz[UU_2D]) - chi_carat*(bx*sy[UU_2D] - by*sx[UU_2D]);

  calc_heat_flux_2D(qx, qy, qz, q);

  calc_heat_grad_2D(q, dx, dy, rhs_d);

  double da = fmin(dx, dy);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_mag_heat_flux_3d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  const double dz = gces->grid.dx[2];

  double dTx[4][6] = {0.0};
  double dTdx[4][6] = {0.0};
  double dTdy[4][6] = {0.0};
  double dTdz[4][6] = {0.0};
  double Tij[8][6] = {0.0};
  double rho[8] = {0.0};
  double p[8] = {0.0};
  double sx[8] = {0.0};
  double sy[8] = {0.0};
  double sz[8] = {0.0};
  double qx[8] = {0.0};
  double qy[8] = {0.0};
  double qz[8] = {0.0};
  double q[8][10] = {0.0};
  var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                 rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
  p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                               p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);

  calc_temp_grad_3D(Tij, dx, dy, dz, dTdx, dTdy, dTdz);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  double omega = gces->omega;

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi_par = alpha*vth_avg*rho_avg;
  double chi_perp = alpha*vth_avg*rho_avg/(1 + omega*omega);
  double chi_carat = alpha*vth_avg*rho_avg*omega/(1 + omega*omega);

  double Bx = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BX], em_tot_d[LLU_3D][BX],
    em_tot_d[LUL_3D][BX], em_tot_d[LUU_3D][BX], em_tot_d[ULL_3D][BX],
    em_tot_d[ULU_3D][BX], em_tot_d[UUL_3D][BX], em_tot_d[UUU_3D][BX]);
  double By = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BY], em_tot_d[LLU_3D][BY],
    em_tot_d[LUL_3D][BY], em_tot_d[LUU_3D][BY], em_tot_d[ULL_3D][BY],
    em_tot_d[ULU_3D][BY], em_tot_d[UUL_3D][BY], em_tot_d[UUU_3D][BY]);
  double Bz = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BZ], em_tot_d[LLU_3D][BZ],
    em_tot_d[LUL_3D][BZ], em_tot_d[LUU_3D][BZ], em_tot_d[ULL_3D][BZ],
    em_tot_d[ULU_3D][BZ], em_tot_d[UUL_3D][BZ], em_tot_d[UUU_3D][BZ]);

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;

  calc_source_3D(dTdx, dTdx, dTdz, sx, sy, sz);

  qx[LLL_3D] = chi_perp*sx[LLL_3D] + (chi_par - chi_perp)*bx*(bx*sx[LLL_3D]
    + by*sy[LLL_3D] + bz*sz[LLL_3D]) - chi_carat*(by*sz[LLL_3D] - bz*sy[LLL_3D]);
  qy[LLL_3D] = chi_perp*sy[LLL_3D] + (chi_par - chi_perp)*by*(bx*sx[LLL_3D]
    + by*sy[LLL_3D] + bz*sz[LLL_3D]) - chi_carat*(bx*sz[LLL_3D] - bz*sx[LLL_3D]);
  qz[LLL_3D] = chi_perp*sz[LLL_3D] + (chi_par - chi_perp)*bz*(bx*sx[LLL_3D]
    + by*sy[LLL_3D] + bz*sz[LLL_3D]) - chi_carat*(bx*sy[LLL_3D] - by*sx[LLL_3D]);

  qx[LLU_3D] = chi_perp*sx[LLU_3D] + (chi_par - chi_perp)*bx*(bx*sx[LLU_3D]
    + by*sy[LLU_3D] + bz*sz[LLU_3D]) - chi_carat*(by*sz[LLU_3D] - bz*sy[LLU_3D]);
  qy[LLU_3D] = chi_perp*sy[LLU_3D] + (chi_par - chi_perp)*by*(bx*sx[LLU_3D]
    + by*sy[LLU_3D] + bz*sz[LLU_3D]) - chi_carat*(bx*sz[LLU_3D] - bz*sx[LLU_3D]);
  qz[LLU_3D] = chi_perp*sz[LLU_3D] + (chi_par - chi_perp)*bz*(bx*sx[LLU_3D]
    + by*sy[LLU_3D] + bz*sz[LLU_3D]) - chi_carat*(bx*sy[LLU_3D] - by*sx[LLU_3D]);

  qx[LUL_3D] = chi_perp*sx[LUL_3D] + (chi_par - chi_perp)*bx*(bx*sx[LUL_3D]
    + by*sy[LUL_3D] + bz*sz[LUL_3D]) - chi_carat*(by*sz[LUL_3D] - bz*sy[LUL_3D]);
  qy[LUL_3D] = chi_perp*sy[LUL_3D] + (chi_par - chi_perp)*by*(bx*sx[LUL_3D]
    + by*sy[LUL_3D] + bz*sz[LUL_3D]) - chi_carat*(bx*sz[LUL_3D] - bz*sx[LUL_3D]);
  qz[LUL_3D] = chi_perp*sz[LUL_3D] + (chi_par - chi_perp)*bz*(bx*sx[LUL_3D]
    + by*sy[LUL_3D] + bz*sz[LUL_3D]) - chi_carat*(bx*sy[LUL_3D] - by*sx[LUL_3D]);
 
  qx[LUU_3D] = chi_perp*sx[LUU_3D] + (chi_par - chi_perp)*bx*(bx*sx[LUU_3D]
    + by*sy[LUU_3D] + bz*sz[LUU_3D]) - chi_carat*(by*sz[LUU_3D] - bz*sy[LUU_3D]);
  qy[LUU_3D] = chi_perp*sy[LUU_3D] + (chi_par - chi_perp)*by*(bx*sx[LUU_3D]
    + by*sy[LUU_3D] + bz*sz[LUU_3D]) - chi_carat*(bx*sz[LUU_3D] - bz*sx[LUU_3D]);
  qz[LUU_3D] = chi_perp*sz[LUU_3D] + (chi_par - chi_perp)*bz*(bx*sx[LUU_3D]
    + by*sy[LUU_3D] + bz*sz[LUU_3D]) - chi_carat*(bx*sy[LUU_3D] - by*sx[LUU_3D]);

  qx[ULL_3D] = chi_perp*sx[ULL_3D] + (chi_par - chi_perp)*bx*(bx*sx[ULL_3D]
    + by*sy[ULL_3D] + bz*sz[ULL_3D]) - chi_carat*(by*sz[ULL_3D] - bz*sy[ULL_3D]);
  qy[ULL_3D] = chi_perp*sy[ULL_3D] + (chi_par - chi_perp)*by*(bx*sx[ULL_3D]
    + by*sy[ULL_3D] + bz*sz[ULL_3D]) - chi_carat*(bx*sz[ULL_3D] - bz*sx[ULL_3D]);
  qz[ULL_3D] = chi_perp*sz[ULL_3D] + (chi_par - chi_perp)*bz*(bx*sx[ULL_3D]
    + by*sy[ULL_3D] + bz*sz[ULL_3D]) - chi_carat*(bx*sy[ULL_3D] - by*sx[ULL_3D]);

  qx[ULU_3D] = chi_perp*sx[ULU_3D] + (chi_par - chi_perp)*bx*(bx*sx[ULU_3D]
    + by*sy[ULU_3D] + bz*sz[ULU_3D]) - chi_carat*(by*sz[ULU_3D] - bz*sy[ULU_3D]);
  qy[ULU_3D] = chi_perp*sy[ULU_3D] + (chi_par - chi_perp)*by*(bx*sx[ULU_3D]
    + by*sy[ULU_3D] + bz*sz[ULU_3D]) - chi_carat*(bx*sz[ULU_3D] - bz*sx[ULU_3D]);
  qz[ULU_3D] = chi_perp*sz[ULU_3D] + (chi_par - chi_perp)*bz*(bx*sx[ULU_3D]
    + by*sy[ULU_3D] + bz*sz[ULU_3D]) - chi_carat*(bx*sy[ULU_3D] - by*sx[ULU_3D]);

  qx[UUL_3D] = chi_perp*sx[UUL_3D] + (chi_par - chi_perp)*bx*(bx*sx[UUL_3D]
    + by*sy[UUL_3D] + bz*sz[UUL_3D]) - chi_carat*(by*sz[UUL_3D] - bz*sy[UUL_3D]);
  qy[UUL_3D] = chi_perp*sy[UUL_3D] + (chi_par - chi_perp)*by*(bx*sx[UUL_3D]
    + by*sy[UUL_3D] + bz*sz[UUL_3D]) - chi_carat*(bx*sz[UUL_3D] - bz*sx[UUL_3D]);
  qz[UUL_3D] = chi_perp*sz[UUL_3D] + (chi_par - chi_perp)*bz*(bx*sx[UUL_3D]
    + by*sy[UUL_3D] + bz*sz[UUL_3D]) - chi_carat*(bx*sy[UUL_3D] - by*sx[UUL_3D]);

  qx[UUU_3D] = chi_perp*sx[UUU_3D] + (chi_par - chi_perp)*bx*(bx*sx[UUU_3D]
    + by*sy[UUU_3D] + bz*sz[UUU_3D]) - chi_carat*(by*sz[UUU_3D] - bz*sy[UUU_3D]);
  qy[UUU_3D] = chi_perp*sy[UUU_3D] + (chi_par - chi_perp)*by*(bx*sx[UUU_3D]
    + by*sy[UUU_3D] + bz*sz[UUU_3D]) - chi_carat*(bx*sz[UUU_3D] - bz*sx[UUU_3D]);
  qz[UUU_3D] = chi_perp*sz[UUU_3D] + (chi_par - chi_perp)*bz*(bx*sx[UUU_3D]
    + by*sy[UUU_3D] + bz*sz[UUU_3D]) - chi_carat*(bx*sy[UUU_3D] - by*sx[UUU_3D]);
  
  calc_heat_flux_3D(qx, qy, qz, q);

  calc_heat_grad_3D(q, dx, dy, dz, rhs_d);

  double da = fmin(fmin(dx, dy), dz);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_par_heat_flux_1d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;

  const double dx = gces->grid.dx[0];
  double dTdx[6] = {0.0};
  double Tij[2][6] = {0.0};
  double rho[2] = {0.0};
  double p[2] = {0.0};
  double q[10] = {0.0};
  var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
  p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);

  calc_temp_grad_1D(Tij, dx, dTdx);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  double Bx = (em_tot_d[L_1D][BX] + em_tot_d[U_1D][BX])/2.0;
  double By = (em_tot_d[L_1D][BY] + em_tot_d[U_1D][BY])/2.0;
  double Bz = (em_tot_d[L_1D][BZ] + em_tot_d[U_1D][BZ])/2.0;

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;
  
  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  calc_source_1D(dTdx, sx, sy, sz);

  double qx = chi*bx*(bx*sx + by*sy + bz*sz);
  double qy = chi*by*(bx*sx + by*sy + bz*sz);
  double qz = chi*bz*(bx*sx + by*sy + bz*sz);

  calc_heat_flux_1D(qx, qy, qz, q);

  calc_heat_grad_1D(q, dx, rhs_d);

  double cfla = dt/(dx*dx);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_par_heat_flux_2d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  double dTdx[2][6] = {0.0};
  double dTdy[2][6] = {0.0};
  double Tij[4][6] = {0.0};
  double rho[4] = {0.0};
  double p[4] = {0.0};
  double sx[4] = {0.0};
  double sy[4] = {0.0};
  double sz[4] = {0.0};
  double qx[4] = {0.0};
  double qy[4] = {0.0};
  double qz[4] = {0.0};
  double q[4][10] = {0.0};
  var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
  p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);

  calc_temp_grad_2D(Tij, dx, dy, dTdx, dTdy);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  double Bx = calc_harmonic_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX],
    em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
  double By = calc_harmonic_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY],
    em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
  double Bz = calc_harmonic_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ],
    em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;

  calc_source_2D(dTdx, dTdy, sx, sy, sz);

  qx[LL_2D] = chi*bx*(bx*sx[LL_2D] + by*sy[LL_2D] + bz*sz[LL_2D]);
  qy[LL_2D] = chi*by*(bx*sx[LL_2D] + by*sy[LL_2D] + bz*sz[LL_2D]);
  qz[LL_2D] = chi*bz*(bx*sx[LL_2D] + by*sy[LL_2D] + bz*sz[LL_2D]);

  qx[LU_2D] = chi*bx*(bx*sx[LU_2D] + by*sy[LU_2D] + bz*sz[LU_2D]);
  qy[LU_2D] = chi*by*(bx*sx[LU_2D] + by*sy[LU_2D] + bz*sz[LU_2D]);
  qz[LU_2D] = chi*bz*(bx*sx[LU_2D] + by*sy[LU_2D] + bz*sz[LU_2D]);

  qx[UL_2D] = chi*bx*(bx*sx[UL_2D] + by*sy[UL_2D] + bz*sz[UL_2D]);
  qy[UL_2D] = chi*by*(bx*sx[UL_2D] + by*sy[UL_2D] + bz*sz[UL_2D]);
  qz[UL_2D] = chi*bz*(bx*sx[UL_2D] + by*sy[UL_2D] + bz*sz[UL_2D]);

  qx[UU_2D] = chi*bx*(bx*sx[UU_2D] + by*sy[UU_2D] + bz*sz[UU_2D]);
  qy[UU_2D] = chi*by*(bx*sx[UU_2D] + by*sy[UU_2D] + bz*sz[UU_2D]);
  qz[UU_2D] = chi*bz*(bx*sx[UU_2D] + by*sy[UU_2D] + bz*sz[UU_2D]);

  calc_heat_flux_2D(qx, qy, qz, q);

  calc_heat_grad_2D(q, dx, dy, rhs_d);

  double da = fmin(dx, dy);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_par_heat_flux_3d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  const double dz = gces->grid.dx[2];

  double dTx[4][6] = {0.0};
  double dTdx[4][6] = {0.0};
  double dTdy[4][6] = {0.0};
  double dTdz[4][6] = {0.0};
  double Tij[8][6] = {0.0};
  double rho[8] = {0.0};
  double p[8] = {0.0};
  double sx[8] = {0.0};
  double sy[8] = {0.0};
  double sz[8] = {0.0};
  double qx[8] = {0.0};
  double qy[8] = {0.0};
  double qz[8] = {0.0};
  double q[8][10] = {0.0};
  var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                 rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
  p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                               p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);

  calc_temp_grad_3D(Tij, dx, dy, dz, dTdx, dTdy, dTdz);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  double Bx = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BX], em_tot_d[LLU_3D][BX],
    em_tot_d[LUL_3D][BX], em_tot_d[LUU_3D][BX], em_tot_d[ULL_3D][BX],
    em_tot_d[ULU_3D][BX], em_tot_d[UUL_3D][BX], em_tot_d[UUU_3D][BX]);
  double By = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BY], em_tot_d[LLU_3D][BY],
    em_tot_d[LUL_3D][BY], em_tot_d[LUU_3D][BY], em_tot_d[ULL_3D][BY],
    em_tot_d[ULU_3D][BY], em_tot_d[UUL_3D][BY], em_tot_d[UUU_3D][BY]);
  double Bz = calc_harmonic_avg_3D(em_tot_d[LLL_3D][BZ], em_tot_d[LLU_3D][BZ],
    em_tot_d[LUL_3D][BZ], em_tot_d[LUU_3D][BZ], em_tot_d[ULL_3D][BZ],
    em_tot_d[ULU_3D][BZ], em_tot_d[UUL_3D][BZ], em_tot_d[UUU_3D][BZ]);

  double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double bx = Bx/B_mag;
  double by = By/B_mag;
  double bz = Bz/B_mag;

  calc_source_3D(dTdx, dTdx, dTdz, sx, sy, sz);

  qx[LLL_3D] = chi*bx*(bx*sx[LLL_3D] + by*sy[LLL_3D] + bz*sz[LLL_3D]);
  qy[LLL_3D] = chi*by*(bx*sx[LLL_3D] + by*sy[LLL_3D] + bz*sz[LLL_3D]);
  qz[LLL_3D] = chi*bz*(bx*sx[LLL_3D] + by*sy[LLL_3D] + bz*sz[LLL_3D]);

  qx[LLU_3D] = chi*bx*(bx*sx[LLU_3D] + by*sy[LLU_3D] + bz*sz[LLU_3D]);
  qy[LLU_3D] = chi*by*(bx*sx[LLU_3D] + by*sy[LLU_3D] + bz*sz[LLU_3D]);
  qz[LLU_3D] = chi*bz*(bx*sx[LLU_3D] + by*sy[LLU_3D] + bz*sz[LLU_3D]);

  qx[LUL_3D] = chi*bx*(bx*sx[LUL_3D] + by*sy[LUL_3D] + bz*sz[LUL_3D]);
  qy[LUL_3D] = chi*by*(bx*sx[LUL_3D] + by*sy[LUL_3D] + bz*sz[LUL_3D]);
  qz[LUL_3D] = chi*bz*(bx*sx[LUL_3D] + by*sy[LUL_3D] + bz*sz[LUL_3D]);

  qx[LUU_3D] = chi*bx*(bx*sx[LUU_3D] + by*sy[LUU_3D] + bz*sz[LUU_3D]);
  qy[LUU_3D] = chi*by*(bx*sx[LUU_3D] + by*sy[LUU_3D] + bz*sz[LUU_3D]);
  qz[LUU_3D] = chi*bz*(bx*sx[LUU_3D] + by*sy[LUU_3D] + bz*sz[LUU_3D]);

  qx[ULL_3D] = chi*bx*(bx*sx[ULL_3D] + by*sy[ULL_3D] + bz*sz[ULL_3D]);
  qy[ULL_3D] = chi*by*(bx*sx[ULL_3D] + by*sy[ULL_3D] + bz*sz[ULL_3D]);
  qz[ULL_3D] = chi*bz*(bx*sx[ULL_3D] + by*sy[ULL_3D] + bz*sz[ULL_3D]);

  qx[ULU_3D] = chi*bx*(bx*sx[ULU_3D] + by*sy[ULU_3D] + bz*sz[ULU_3D]);
  qy[ULU_3D] = chi*by*(bx*sx[ULU_3D] + by*sy[ULU_3D] + bz*sz[ULU_3D]);
  qz[ULU_3D] = chi*bz*(bx*sx[ULU_3D] + by*sy[ULU_3D] + bz*sz[ULU_3D]);

  qx[UUL_3D] = chi*bx*(bx*sx[UUL_3D] + by*sy[UUL_3D] + bz*sz[UUL_3D]);
  qy[UUL_3D] = chi*by*(bx*sx[UUL_3D] + by*sy[UUL_3D] + bz*sz[UUL_3D]);
  qz[UUL_3D] = chi*bz*(bx*sx[UUL_3D] + by*sy[UUL_3D] + bz*sz[UUL_3D]);

  qx[UUU_3D] = chi*bx*(bx*sx[UUU_3D] + by*sy[UUU_3D] + bz*sz[UUU_3D]);
  qy[UUU_3D] = chi*by*(bx*sx[UUU_3D] + by*sy[UUU_3D] + bz*sz[UUU_3D]);
  qz[UUU_3D] = chi*bz*(bx*sx[UUU_3D] + by*sy[UUU_3D] + bz*sz[UUU_3D]);
  
  calc_heat_flux_3D(qx, qy, qz, q);

  calc_heat_grad_3D(q, dx, dy, dz, rhs_d);

  double da = fmin(fmin(dx, dy), dz);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_unmag_heat_flux_1d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;

  const double dx = gces->grid.dx[0];
  double dTdx[6] = {0.0};
  double Tij[2][6] = {0.0};
  double rho[2] = {0.0};
  double p[2] = {0.0};
  double q[10] = {0.0};
  var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
  p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);

  calc_temp_grad_1D(Tij, dx, dTdx);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  
  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  q[Q111] = chi*dTdx[T11];
  q[Q112] = chi*2.0*dTdx[T12]/3.0;
  q[Q113] = chi*2.0*dTdx[T13]/3.0;
  q[Q122] = chi*dTdx[T22]/3.0;
  q[Q123] = chi*dTdx[T23]/3.0;
  q[Q133] = chi*dTdx[T33]/3.0;

  calc_heat_grad_1D(q, dx, rhs_d);

  double cfla = dt/(dx*dx);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_unmag_heat_flux_2d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  double dTdx[2][6] = {0.0};
  double dTdy[2][6] = {0.0};
  double Tij[4][6] = {0.0};
  double rho[4] = {0.0};
  double p[4] = {0.0};
  double q[4][10] = {0.0};
  var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
  p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);

  calc_temp_grad_2D(Tij, dx, dy, dTdx, dTdy);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  int compx[4] = { L_1D, U_1D, L_1D, U_1D };
  int compy[4] = { L_1D, L_1D, U_1D, U_1D };

  q[LL_2D][Q111] = chi*dTdx[compx[LL_2D]][T11];
  q[LL_2D][Q112] = chi*(2.0*dTdx[compx[LL_2D]][T12] + dTdy[compy[LL_2D]][T11])/3.0;
  q[LL_2D][Q113] = chi*2.0*dTdx[compx[LL_2D]][T13]/3.0;
  q[LL_2D][Q122] = chi*(dTdx[compx[LL_2D]][T22] + 2.0*dTdy[compy[LL_2D]][T12])/3.0;
  q[LL_2D][Q123] = chi*(dTdx[compx[LL_2D]][T23] + dTdy[compy[LL_2D]][T13])/3.0;
  q[LL_2D][Q133] = chi*dTdx[compx[LL_2D]][T33]/3.0;
  q[LL_2D][Q222] = chi*dTdy[compy[LL_2D]][T22];
  q[LL_2D][Q223] = chi*2.0*dTdy[compy[LL_2D]][T23]/3.0;
  q[LL_2D][Q233] = chi*dTdy[compy[LL_2D]][T33]/3.0;

  q[LU_2D][Q111] = chi*dTdx[compx[LU_2D]][T11];
  q[LU_2D][Q112] = chi*(2.0*dTdx[compx[LU_2D]][T12] + dTdy[compy[LU_2D]][T11])/3.0;
  q[LU_2D][Q113] = chi*2.0*dTdx[compx[LU_2D]][T13]/3.0;
  q[LU_2D][Q122] = chi*(dTdx[compx[LU_2D]][T22] + 2.0*dTdy[compy[LU_2D]][T12])/3.0;
  q[LU_2D][Q123] = chi*(dTdx[compx[LU_2D]][T23] + dTdy[compy[LU_2D]][T13])/3.0;
  q[LU_2D][Q133] = chi*dTdx[compx[LU_2D]][T33]/3.0;
  q[LU_2D][Q222] = chi*dTdy[compy[LU_2D]][T22];
  q[LU_2D][Q223] = chi*2.0*dTdy[compy[LU_2D]][T23]/3.0;
  q[LU_2D][Q233] = chi*dTdy[compy[LU_2D]][T33]/3.0;

  q[UL_2D][Q111] = chi*dTdx[compx[UL_2D]][T11];
  q[UL_2D][Q112] = chi*(2.0*dTdx[compx[UL_2D]][T12] + dTdy[compy[UL_2D]][T11])/3.0;
  q[UL_2D][Q113] = chi*2.0*dTdx[compx[UL_2D]][T13]/3.0;
  q[UL_2D][Q122] = chi*(dTdx[compx[UL_2D]][T22] + 2.0*dTdy[compy[UL_2D]][T12])/3.0;
  q[UL_2D][Q123] = chi*(dTdx[compx[UL_2D]][T23] + dTdy[compy[UL_2D]][T13])/3.0;
  q[UL_2D][Q133] = chi*dTdx[compx[UL_2D]][T33]/3.0;
  q[UL_2D][Q222] = chi*dTdy[compy[UL_2D]][T22];
  q[UL_2D][Q223] = chi*2.0*dTdy[compy[UL_2D]][T23]/3.0;
  q[UL_2D][Q233] = chi*dTdy[compy[UL_2D]][T33]/3.0;

  q[UU_2D][Q111] = chi*dTdx[compx[UU_2D]][T11];
  q[UU_2D][Q112] = chi*(2.0*dTdx[compx[UU_2D]][T12] + dTdy[compy[UU_2D]][T11])/3.0;
  q[UU_2D][Q113] = chi*2.0*dTdx[compx[UU_2D]][T13]/3.0;
  q[UU_2D][Q122] = chi*(dTdx[compx[UU_2D]][T22] + 2.0*dTdy[compy[UU_2D]][T12])/3.0;
  q[UU_2D][Q123] = chi*(dTdx[compx[UU_2D]][T23] + dTdy[compy[UU_2D]][T13])/3.0;
  q[UU_2D][Q133] = chi*dTdx[compx[UU_2D]][T33]/3.0;
  q[UU_2D][Q222] = chi*dTdy[compy[UU_2D]][T22];
  q[UU_2D][Q223] = chi*2.0*dTdy[compy[UU_2D]][T23]/3.0;
  q[UU_2D][Q233] = chi*dTdy[compy[UU_2D]][T33]/3.0;

  calc_heat_grad_2D(q, dx, dy, rhs_d);
  
  double da = fmin(dx, dy);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

double
calc_unmag_heat_flux_3d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], const double *em_tot_d[], double *cflrate,
  double cfl, double dt, double *rhs_d[])
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  
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

  calc_temp_grad_3D(Tij, dx, dy, dz, dTdx, dTdy, dTdz);
  
  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  int compx[8] = { LL_2D, UL_2D, LU_2D, UU_2D, LL_2D, UL_2D, LU_2D, UU_2D };
  int compy[8] = { LL_2D, UL_2D, LL_2D, UL_2D, LU_2D, UU_2D, LU_2D, UU_2D };
  int compz[8] = { LL_2D, LL_2D, UL_2D, UL_2D, LU_2D, LU_2D, UU_2D, UU_2D };

  q[LLL_3D][Q111] = chi*dTdx[compx[LLL_3D]][T11];
  q[LLL_3D][Q112] = chi*(2.0*dTdx[compx[LLL_3D]][T12]
    + dTdy[compy[LLL_3D]][T11])/3.0;
  q[LLL_3D][Q113] = chi*(2.0*dTdx[compx[LLL_3D]][T13]
    + dTdz[compz[LLL_3D]][T11])/3.0;
  q[LLL_3D][Q122] = chi*(dTdx[compx[LLL_3D]][T22]
    + 2.0*dTdy[compy[LLL_3D]][T12])/3.0;
  q[LLL_3D][Q123] = chi*(dTdx[compx[LLL_3D]][T23]
    + dTdy[compy[LLL_3D]][T13] + dTdz[compz[LLL_3D]][T12])/3.0;
  q[LLL_3D][Q133] = chi*(dTdx[compx[LLL_3D]][T33]
    + 2.0*dTdz[compz[LLL_3D]][T13])/3.0;
  q[LLL_3D][Q222] = chi*dTdy[compy[LLL_3D]][T22];
  q[LLL_3D][Q223] = chi*(2.0*dTdy[compy[LLL_3D]][T23]
    + dTdz[compz[LLL_3D]][T22])/3.0;
  q[LLL_3D][Q233] = chi*(dTdy[compy[LLL_3D]][T33]
    + 2.0*dTdz[compz[LLL_3D]][T23])/3.0;
  q[LLL_3D][Q333] = chi*dTdz[compz[LLL_3D]][T33];

  q[LLU_3D][Q111] = chi*dTdx[compx[LLU_3D]][T11];
  q[LLU_3D][Q112] = chi*(2.0*dTdx[compx[LLU_3D]][T12]
    + dTdy[compy[LLU_3D]][T11])/3.0;
  q[LLU_3D][Q113] = chi*(2.0*dTdx[compx[LLU_3D]][T13]
    + dTdz[compz[LLU_3D]][T11])/3.0;
  q[LLU_3D][Q122] = chi*(dTdx[compx[LLU_3D]][T22]
    + 2.0*dTdy[compy[LLU_3D]][T12])/3.0;
  q[LLU_3D][Q123] = chi*(dTdx[compx[LLU_3D]][T23]
    + dTdy[compy[LLU_3D]][T13] + dTdz[compz[LLU_3D]][T12])/3.0;
  q[LLU_3D][Q133] = chi*(dTdx[compx[LLU_3D]][T33]
    + 2.0*dTdz[compz[LLU_3D]][T13])/3.0;
  q[LLU_3D][Q222] = chi*dTdy[compy[LLU_3D]][T22];
  q[LLU_3D][Q223] = chi*(2.0*dTdy[compy[LLU_3D]][T23]
    + dTdz[compz[LLU_3D]][T22])/3.0;
  q[LLU_3D][Q233] = chi*(dTdy[compy[LLU_3D]][T33]
    + 2.0*dTdz[compz[LLU_3D]][T23])/3.0;
  q[LLU_3D][Q333] = chi*dTdz[compz[LLU_3D]][T33];

 q[LUL_3D][Q111] = chi*dTdx[compx[LUL_3D]][T11];
  q[LUL_3D][Q112] = chi*(2.0*dTdx[compx[LUL_3D]][T12]
    + dTdy[compy[LUL_3D]][T11])/3.0;
  q[LUL_3D][Q113] = chi*(2.0*dTdx[compx[LUL_3D]][T13]
    + dTdz[compz[LUL_3D]][T11])/3.0;
  q[LUL_3D][Q122] = chi*(dTdx[compx[LUL_3D]][T22]
    + 2.0*dTdy[compy[LUL_3D]][T12])/3.0;
  q[LUL_3D][Q123] = chi*(dTdx[compx[LUL_3D]][T23]
    + dTdy[compy[LUL_3D]][T13] + dTdz[compz[LUL_3D]][T12])/3.0;
  q[LUL_3D][Q133] = chi*(dTdx[compx[LUL_3D]][T33]
    + 2.0*dTdz[compz[LUL_3D]][T13])/3.0;
  q[LUL_3D][Q222] = chi*dTdy[compy[LUL_3D]][T22];
  q[LUL_3D][Q223] = chi*(2.0*dTdy[compy[LUL_3D]][T23]
    + dTdz[compz[LUL_3D]][T22])/3.0;
  q[LUL_3D][Q233] = chi*(dTdy[compy[LUL_3D]][T33]
    + 2.0*dTdz[compz[LUL_3D]][T23])/3.0;
  q[LUL_3D][Q333] = chi*dTdz[compz[LUL_3D]][T33];

 q[LUU_3D][Q111] = chi*dTdx[compx[LUU_3D]][T11];
  q[LUU_3D][Q112] = chi*(2.0*dTdx[compx[LUU_3D]][T12]
    + dTdy[compy[LUU_3D]][T11])/3.0;
  q[LUU_3D][Q113] = chi*(2.0*dTdx[compx[LUU_3D]][T13]
    + dTdz[compz[LUU_3D]][T11])/3.0;
  q[LUU_3D][Q122] = chi*(dTdx[compx[LUU_3D]][T22]
    + 2.0*dTdy[compy[LUU_3D]][T12])/3.0;
  q[LUU_3D][Q123] = chi*(dTdx[compx[LUU_3D]][T23]
    + dTdy[compy[LUU_3D]][T13] + dTdz[compz[LUU_3D]][T12])/3.0;
  q[LUU_3D][Q133] = chi*(dTdx[compx[LUU_3D]][T33]
    + 2.0*dTdz[compz[LUU_3D]][T13])/3.0;
  q[LUU_3D][Q222] = chi*dTdy[compy[LUU_3D]][T22];
  q[LUU_3D][Q223] = chi*(2.0*dTdy[compy[LUU_3D]][T23]
    + dTdz[compz[LUU_3D]][T22])/3.0;
  q[LUU_3D][Q233] = chi*(dTdy[compy[LUU_3D]][T33]
    + 2.0*dTdz[compz[LUU_3D]][T23])/3.0;
  q[LUU_3D][Q333] = chi*dTdz[compz[LUU_3D]][T33];

 q[ULL_3D][Q111] = chi*dTdx[compx[ULL_3D]][T11];
  q[ULL_3D][Q112] = chi*(2.0*dTdx[compx[ULL_3D]][T12]
    + dTdy[compy[ULL_3D]][T11])/3.0;
  q[ULL_3D][Q113] = chi*(2.0*dTdx[compx[ULL_3D]][T13]
    + dTdz[compz[ULL_3D]][T11])/3.0;
  q[ULL_3D][Q122] = chi*(dTdx[compx[ULL_3D]][T22]
    + 2.0*dTdy[compy[ULL_3D]][T12])/3.0;
  q[ULL_3D][Q123] = chi*(dTdx[compx[ULL_3D]][T23]
    + dTdy[compy[ULL_3D]][T13] + dTdz[compz[ULL_3D]][T12])/3.0;
  q[ULL_3D][Q133] = chi*(dTdx[compx[ULL_3D]][T33]
    + 2.0*dTdz[compz[ULL_3D]][T13])/3.0;
  q[ULL_3D][Q222] = chi*dTdy[compy[ULL_3D]][T22];
  q[ULL_3D][Q223] = chi*(2.0*dTdy[compy[ULL_3D]][T23]
    + dTdz[compz[ULL_3D]][T22])/3.0;
  q[ULL_3D][Q233] = chi*(dTdy[compy[ULL_3D]][T33]
    + 2.0*dTdz[compz[ULL_3D]][T23])/3.0;
  q[ULL_3D][Q333] = chi*dTdz[compz[ULL_3D]][T33];

  q[ULU_3D][Q111] = chi*dTdx[compx[ULU_3D]][T11];
  q[ULU_3D][Q112] = chi*(2.0*dTdx[compx[ULU_3D]][T12]
    + dTdy[compy[ULU_3D]][T11])/3.0;
  q[ULU_3D][Q113] = chi*(2.0*dTdx[compx[ULU_3D]][T13]
    + dTdz[compz[ULU_3D]][T11])/3.0;
  q[ULU_3D][Q122] = chi*(dTdx[compx[ULU_3D]][T22]
    + 2.0*dTdy[compy[ULU_3D]][T12])/3.0;
  q[ULU_3D][Q123] = chi*(dTdx[compx[ULU_3D]][T23]
    + dTdy[compy[ULU_3D]][T13] + dTdz[compz[ULU_3D]][T12])/3.0;
  q[ULU_3D][Q133] = chi*(dTdx[compx[ULU_3D]][T33]
    + 2.0*dTdz[compz[ULU_3D]][T13])/3.0;
  q[ULU_3D][Q222] = chi*dTdy[compy[ULU_3D]][T22];
  q[ULU_3D][Q223] = chi*(2.0*dTdy[compy[ULU_3D]][T23]
    + dTdz[compz[ULU_3D]][T22])/3.0;
  q[ULU_3D][Q233] = chi*(dTdy[compy[ULU_3D]][T33]
    + 2.0*dTdz[compz[ULU_3D]][T23])/3.0;
  q[ULU_3D][Q333] = chi*dTdz[compz[ULU_3D]][T33];

  q[UUL_3D][Q111] = chi*dTdx[compx[UUL_3D]][T11];
  q[UUL_3D][Q112] = chi*(2.0*dTdx[compx[UUL_3D]][T12]
    + dTdy[compy[UUL_3D]][T11])/3.0;
  q[UUL_3D][Q113] = chi*(2.0*dTdx[compx[UUL_3D]][T13]
    + dTdz[compz[UUL_3D]][T11])/3.0;
  q[UUL_3D][Q122] = chi*(dTdx[compx[UUL_3D]][T22]
    + 2.0*dTdy[compy[UUL_3D]][T12])/3.0;
  q[UUL_3D][Q123] = chi*(dTdx[compx[UUL_3D]][T23]
    + dTdy[compy[UUL_3D]][T13] + dTdz[compz[UUL_3D]][T12])/3.0;
  q[UUL_3D][Q133] = chi*(dTdx[compx[UUL_3D]][T33]
    + 2.0*dTdz[compz[UUL_3D]][T13])/3.0;
  q[UUL_3D][Q222] = chi*dTdy[compy[UUL_3D]][T22];
  q[UUL_3D][Q223] = chi*(2.0*dTdy[compy[UUL_3D]][T23]
    + dTdz[compz[UUL_3D]][T22])/3.0;
  q[UUL_3D][Q233] = chi*(dTdy[compy[UUL_3D]][T33]
    + 2.0*dTdz[compz[UUL_3D]][T23])/3.0;
  q[UUL_3D][Q333] = chi*dTdz[compz[UUL_3D]][T33];

  q[UUU_3D][Q111] = chi*dTdx[compx[UUU_3D]][T11];
  q[UUU_3D][Q112] = chi*(2.0*dTdx[compx[UUU_3D]][T12]
    + dTdy[compy[UUU_3D]][T11])/3.0;
  q[UUU_3D][Q113] = chi*(2.0*dTdx[compx[UUU_3D]][T13]
    + dTdz[compz[UUU_3D]][T11])/3.0;
  q[UUU_3D][Q122] = chi*(dTdx[compx[UUU_3D]][T22]
    + 2.0*dTdy[compy[UUU_3D]][T12])/3.0;
  q[UUU_3D][Q123] = chi*(dTdx[compx[UUU_3D]][T23]
    + dTdy[compy[UUU_3D]][T13] + dTdz[compz[UUU_3D]][T12])/3.0;
  q[UUU_3D][Q133] = chi*(dTdx[compx[UUU_3D]][T33]
    + 2.0*dTdz[compz[UUU_3D]][T13])/3.0;
  q[UUU_3D][Q222] = chi*dTdy[compy[UUU_3D]][T22];
  q[UUU_3D][Q223] = chi*(2.0*dTdy[compy[UUU_3D]][T23]
    + dTdz[compz[UUU_3D]][T22])/3.0;
  q[UUU_3D][Q233] = chi*(dTdy[compy[UUU_3D]][T33]
    + 2.0*dTdz[compz[UUU_3D]][T23])/3.0;
  q[UUU_3D][Q333] = chi*dTdz[compz[UUU_3D]][T33];

  calc_heat_grad_3D(q, dx, dy, dz, rhs_d);

  double da = fmin(fmin(dx, dy), dz);
  double cfla = dt/(da*da);
  return fmax(alpha*vth_avg*cfla, cfl);
}

static const heat_flux_calc_t grad_closure_mag_funcs[3] = { calc_mag_heat_flux_1d,
  calc_mag_heat_flux_2d, calc_mag_heat_flux_3d };

static const heat_flux_calc_t grad_closure_par_funcs[3] = { calc_par_heat_flux_1d,
  calc_par_heat_flux_2d, calc_par_heat_flux_3d };

static const heat_flux_calc_t grad_closure_unmag_funcs[3] = { calc_unmag_heat_flux_1d,
  calc_unmag_heat_flux_2d, calc_unmag_heat_flux_3d };

void
grad_closure_calc_q_choose(struct gkyl_ten_moment_grad_closure *gces)
{
  if (gces->mag && gces->omega) {
    gces->calc_q = grad_closure_mag_funcs[gces->ndim - 1];
  }
  else if (gces->mag) {
    gces->calc_q = grad_closure_par_funcs[gces->ndim - 1];
  }
  else {
    gces->calc_q = grad_closure_unmag_funcs[gces->ndim - 1];
  }
}
