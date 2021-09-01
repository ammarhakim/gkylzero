#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_moment_braginskii_priv.h>

// Makes indexing cleaner
static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned P11 = 4;
static const unsigned P12 = 5;
static const unsigned P13 = 6;
static const unsigned P22 = 7;
static const unsigned P23 = 8;
static const unsigned P33 = 9;

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

// 1D stencil locations (L: lower, C: center, U: upper)
enum loc_1d {
  L_1D, C_1D, U_1D
};

// 2D stencil locations (L: lower, C: center, U: upper)
enum loc_2d {
  LL_2D, LC_2D, LU_2D,
  CL_2D, CC_2D, CU_2D,
  UL_2D, UC_2D, UU_2D
};

// 3D stencil locations (L: lower, C: center, U: upper)
enum loc_3d {
  LLL_3D, LLC_3D, LLU_3D,
  LCL_3D, LCC_3D, LCU_3D,
  LUL_3D, LUC_3D, LUU_3D,
  
  CLL_3D, CLC_3D, CLU_3D,
  CCL_3D, CCC_3D, CCU_3D,
  CUL_3D, CUC_3D, CUU_3D,

  ULL_3D, ULC_3D, ULU_3D,
  UCL_3D, UCC_3D, UCU_3D,
  UUL_3D, UUC_3D, UUU_3D
};

struct gkyl_ten_moment_grad_closure {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  double k0; // damping coefficient
};

static void
create_offsets(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
unmag_grad_closure_update(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double *rhs)
{
  const int ndim = gces->ndim;
  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double Tij[3][6] = {0.0};
    double rho[3] = {0.0};
    double p[3] = {0.0}; 
    for (int j = L_1D; j <= U_1D; ++j) {
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
    double rhoL = calc_harmonic_avg_1D(rho[L_1D], rho[C_1D]);
    double pL = calc_harmonic_avg_1D(p[L_1D], p[C_1D]);
    double rhoU = calc_harmonic_avg_1D(rho[C_1D], rho[U_1D]);
    double pU = calc_harmonic_avg_1D(p[C_1D], p[U_1D]);
    double vthL = sqrt(pL/rhoL);
    double vthU = sqrt(pU/rhoU);

    double dTdxL[6] = { calc_sym_grad_1D(dx, Tij[L_1D][T11], Tij[C_1D][T11]),
                        calc_sym_grad_1D(dx, Tij[L_1D][T12], Tij[C_1D][T12]),
                        calc_sym_grad_1D(dx, Tij[L_1D][T13], Tij[C_1D][T13]),
                        calc_sym_grad_1D(dx, Tij[L_1D][T22], Tij[C_1D][T22]),
                        calc_sym_grad_1D(dx, Tij[L_1D][T23], Tij[C_1D][T23]),
                        calc_sym_grad_1D(dx, Tij[L_1D][T33], Tij[C_1D][T33]) };

    double dTdxU[6] = { calc_sym_grad_1D(dx, Tij[C_1D][T11], Tij[U_1D][T11]),
                        calc_sym_grad_1D(dx, Tij[C_1D][T12], Tij[U_1D][T12]),
                        calc_sym_grad_1D(dx, Tij[C_1D][T13], Tij[U_1D][T13]),
                        calc_sym_grad_1D(dx, Tij[C_1D][T22], Tij[U_1D][T22]),
                        calc_sym_grad_1D(dx, Tij[C_1D][T23], Tij[U_1D][T23]),
                        calc_sym_grad_1D(dx, Tij[C_1D][T33], Tij[U_1D][T33]) };
                        
    double qL[9] = {0.0};
    double qU[9] = {0.0};
    double alpha = 1.0/gces->k0;
    
    qL[Q111] = alpha*vthL*rhoL*(dTdxL[T11] + dTdxL[T11] + dTdxL[T11])/3.0;
    qL[Q112] = alpha*vthL*rhoL*(dTdxL[T12] + dTdxL[T12])/3.0;
    qL[Q113] = alpha*vthL*rhoL*(dTdxL[T13] + dTdxL[T13])/3.0;
    qL[Q122] = alpha*vthL*rhoL*(dTdxL[T22])/3.0;
    qL[Q123] = alpha*vthL*rhoL*(dTdxL[T23])/3.0;
    qL[Q133] = alpha*vthL*rhoL*(dTdxL[T33])/3.0; 

    qU[Q111] = alpha*vthU*rhoU*(dTdxU[T11] + dTdxU[T11] + dTdxU[T11])/3.0;
    qU[Q112] = alpha*vthU*rhoU*(dTdxU[T12] + dTdxU[T12])/3.0;
    qU[Q113] = alpha*vthU*rhoU*(dTdxU[T13] + dTdxU[T13])/3.0;
    qU[Q122] = alpha*vthU*rhoU*(dTdxU[T22])/3.0;
    qU[Q123] = alpha*vthU*rhoU*(dTdxU[T23])/3.0;
    qU[Q133] = alpha*vthU*rhoU*(dTdxU[T33])/3.0; 

    rhs[0] = (qU[Q111] - qL[Q111])/dx;
    rhs[1] = (qU[Q112] - qL[Q112])/dx;
    rhs[2] = (qU[Q113] - qL[Q113])/dx;
    rhs[3] = (qU[Q122] - qL[Q122])/dx;
    rhs[4] = (qU[Q123] - qL[Q123])/dx;
    rhs[5] = (qU[Q133] - qL[Q133])/dx;
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[0];
    double Tij[9][6] = {0.0};
    double rho[9] = {0.0};
    double p[9] = {0.0}; 
    for (int j = LL_2D; j <= UU_2D; ++j) {
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
    double rhoLL = calc_harmonic_avg_2D(rho[LL_2D], rho[LC_2D], rho[CL_2D], rho[CC_2D]);
    double pLL = calc_harmonic_avg_2D(p[LL_2D], p[LC_2D], p[CL_2D], p[CC_2D]);
    double rhoLU = calc_harmonic_avg_2D(rho[LC_2D], rho[LU_2D], rho[CC_2D], rho[CU_2D]);
    double pLU = calc_harmonic_avg_2D(p[LC_2D], p[LU_2D], p[CC_2D], p[CU_2D]);
    double rhoUL = calc_harmonic_avg_2D(rho[CL_2D], rho[CC_2D], rho[UL_2D], rho[UC_2D]);
    double pUL = calc_harmonic_avg_2D(p[CL_2D], p[CC_2D], p[UL_2D], p[UC_2D]);
    double rhoUU = calc_harmonic_avg_2D(rho[CC_2D], rho[CU_2D], rho[UC_2D], rho[UU_2D]);
    double pUU = calc_harmonic_avg_2D(p[CC_2D], p[CU_2D], p[UC_2D], p[UU_2D]);

    double vthLL = sqrt(pLL/rhoLL);
    double vthLU = sqrt(pLU/rhoLU);
    double vthUL = sqrt(pUL/rhoUL);
    double vthUU = sqrt(pUU/rhoUU);

    double dTdxLL[6] = { calc_sym_gradx_2D(dx, Tij[LL_2D][T11], Tij[LC_2D][T11], Tij[CL_2D][T11], Tij[CC_2D][T11]),
                         calc_sym_gradx_2D(dx, Tij[LL_2D][T12], Tij[LC_2D][T12], Tij[CL_2D][T12], Tij[CC_2D][T12]),
                         calc_sym_gradx_2D(dx, Tij[LL_2D][T13], Tij[LC_2D][T13], Tij[CL_2D][T13], Tij[CC_2D][T13]),
                         calc_sym_gradx_2D(dx, Tij[LL_2D][T22], Tij[LC_2D][T22], Tij[CL_2D][T22], Tij[CC_2D][T22]), 
                         calc_sym_gradx_2D(dx, Tij[LL_2D][T23], Tij[LC_2D][T23], Tij[CL_2D][T23], Tij[CC_2D][T23]),
                         calc_sym_gradx_2D(dx, Tij[LL_2D][T33], Tij[LC_2D][T33], Tij[CL_2D][T33], Tij[CC_2D][T33]) };

    double dTdyLL[6] = { calc_sym_grady_2D(dy, Tij[LL_2D][T11], Tij[LC_2D][T11], Tij[CL_2D][T11], Tij[CC_2D][T11]),
                         calc_sym_grady_2D(dy, Tij[LL_2D][T12], Tij[LC_2D][T12], Tij[CL_2D][T12], Tij[CC_2D][T12]),
                         calc_sym_grady_2D(dy, Tij[LL_2D][T13], Tij[LC_2D][T13], Tij[CL_2D][T13], Tij[CC_2D][T13]),
                         calc_sym_grady_2D(dy, Tij[LL_2D][T22], Tij[LC_2D][T22], Tij[CL_2D][T22], Tij[CC_2D][T22]), 
                         calc_sym_grady_2D(dy, Tij[LL_2D][T23], Tij[LC_2D][T23], Tij[CL_2D][T23], Tij[CC_2D][T23]),
                         calc_sym_grady_2D(dy, Tij[LL_2D][T33], Tij[LC_2D][T33], Tij[CL_2D][T33], Tij[CC_2D][T33]) };

    double dTdxLU[6] = { calc_sym_gradx_2D(dx, Tij[LC_2D][T11], Tij[LU_2D][T11], Tij[CC_2D][T11], Tij[CU_2D][T11]),
                         calc_sym_gradx_2D(dx, Tij[LC_2D][T12], Tij[LU_2D][T12], Tij[CC_2D][T12], Tij[CU_2D][T12]),
                         calc_sym_gradx_2D(dx, Tij[LC_2D][T13], Tij[LU_2D][T13], Tij[CC_2D][T13], Tij[CU_2D][T13]),
                         calc_sym_gradx_2D(dx, Tij[LC_2D][T22], Tij[LU_2D][T22], Tij[CC_2D][T22], Tij[CU_2D][T22]),
                         calc_sym_gradx_2D(dx, Tij[LC_2D][T23], Tij[LU_2D][T23], Tij[CC_2D][T23], Tij[CU_2D][T23]),
                         calc_sym_gradx_2D(dx, Tij[LC_2D][T33], Tij[LU_2D][T33], Tij[CC_2D][T33], Tij[CU_2D][T33]) };

    double dTdyLU[6] = { calc_sym_grady_2D(dy, Tij[LC_2D][T11], Tij[LU_2D][T11], Tij[CC_2D][T11], Tij[CU_2D][T11]),
                         calc_sym_grady_2D(dy, Tij[LC_2D][T12], Tij[LU_2D][T12], Tij[CC_2D][T12], Tij[CU_2D][T12]),
                         calc_sym_grady_2D(dy, Tij[LC_2D][T13], Tij[LU_2D][T13], Tij[CC_2D][T13], Tij[CU_2D][T13]),
                         calc_sym_grady_2D(dy, Tij[LC_2D][T22], Tij[LU_2D][T22], Tij[CC_2D][T22], Tij[CU_2D][T22]),
                         calc_sym_grady_2D(dy, Tij[LC_2D][T23], Tij[LU_2D][T23], Tij[CC_2D][T23], Tij[CU_2D][T23]),
                         calc_sym_grady_2D(dy, Tij[LC_2D][T33], Tij[LU_2D][T33], Tij[CC_2D][T33], Tij[CU_2D][T33]) };

    double dTdxUL[6] = { calc_sym_gradx_2D(dx, Tij[CL_2D][T11], Tij[CC_2D][T11], Tij[UL_2D][T11], Tij[UC_2D][T11]),
                         calc_sym_gradx_2D(dx, Tij[CL_2D][T12], Tij[CC_2D][T12], Tij[UL_2D][T12], Tij[UC_2D][T12]),
                         calc_sym_gradx_2D(dx, Tij[CL_2D][T13], Tij[CC_2D][T13], Tij[UL_2D][T13], Tij[UC_2D][T13]),
                         calc_sym_gradx_2D(dx, Tij[CL_2D][T22], Tij[CC_2D][T22], Tij[UL_2D][T22], Tij[UC_2D][T22]), 
                         calc_sym_gradx_2D(dx, Tij[CL_2D][T23], Tij[CC_2D][T23], Tij[UL_2D][T23], Tij[UC_2D][T23]),
                         calc_sym_gradx_2D(dx, Tij[CL_2D][T33], Tij[CC_2D][T33], Tij[UL_2D][T33], Tij[UC_2D][T33]) };

    double dTdyUL[6] = { calc_sym_grady_2D(dy, Tij[CL_2D][T11], Tij[CC_2D][T11], Tij[UL_2D][T11], Tij[UC_2D][T11]),
                         calc_sym_grady_2D(dy, Tij[CL_2D][T12], Tij[CC_2D][T12], Tij[UL_2D][T12], Tij[UC_2D][T12]),
                         calc_sym_grady_2D(dy, Tij[CL_2D][T13], Tij[CC_2D][T13], Tij[UL_2D][T13], Tij[UC_2D][T13]),
                         calc_sym_grady_2D(dy, Tij[CL_2D][T22], Tij[CC_2D][T22], Tij[UL_2D][T22], Tij[UC_2D][T22]), 
                         calc_sym_grady_2D(dy, Tij[CL_2D][T23], Tij[CC_2D][T23], Tij[UL_2D][T23], Tij[UC_2D][T23]),
                         calc_sym_grady_2D(dy, Tij[CL_2D][T33], Tij[CC_2D][T33], Tij[UL_2D][T33], Tij[UC_2D][T33]) };

    double dTdxUU[6] = { calc_sym_gradx_2D(dx, Tij[CC_2D][T11], Tij[CU_2D][T11], Tij[UC_2D][T11], Tij[UU_2D][T11]),
                         calc_sym_gradx_2D(dx, Tij[CC_2D][T12], Tij[CU_2D][T12], Tij[UC_2D][T12], Tij[UU_2D][T12]),
                         calc_sym_gradx_2D(dx, Tij[CC_2D][T13], Tij[CU_2D][T13], Tij[UC_2D][T13], Tij[UU_2D][T13]),
                         calc_sym_gradx_2D(dx, Tij[CC_2D][T22], Tij[CU_2D][T22], Tij[UC_2D][T22], Tij[UU_2D][T22]),
                         calc_sym_gradx_2D(dx, Tij[CC_2D][T23], Tij[CU_2D][T23], Tij[UC_2D][T23], Tij[UU_2D][T23]),
                         calc_sym_gradx_2D(dx, Tij[CC_2D][T33], Tij[CU_2D][T33], Tij[UC_2D][T33], Tij[UU_2D][T33]) };

    double dTdyUU[6] = { calc_sym_grady_2D(dy, Tij[CC_2D][T11], Tij[CU_2D][T11], Tij[UC_2D][T11], Tij[UU_2D][T11]),
                         calc_sym_grady_2D(dy, Tij[CC_2D][T12], Tij[CU_2D][T12], Tij[UC_2D][T12], Tij[UU_2D][T12]),
                         calc_sym_grady_2D(dy, Tij[CC_2D][T13], Tij[CU_2D][T13], Tij[UC_2D][T13], Tij[UU_2D][T13]),
                         calc_sym_grady_2D(dy, Tij[CC_2D][T22], Tij[CU_2D][T22], Tij[UC_2D][T22], Tij[UU_2D][T22]),
                         calc_sym_grady_2D(dy, Tij[CC_2D][T23], Tij[CU_2D][T23], Tij[UC_2D][T23], Tij[UU_2D][T23]),
                         calc_sym_grady_2D(dy, Tij[CC_2D][T33], Tij[CU_2D][T33], Tij[UC_2D][T33], Tij[UU_2D][T33]) };

    double qLL[9] = {0.0};
    double qLU[9] = {0.0};
    double qUL[9] = {0.0};
    double qUU[9] = {0.0};
    double alpha = 1.0/gces->k0;
    
    qLL[Q111] = alpha*vthLL*rhoLL*(dTdxLL[T11] + dTdxLL[T11] + dTdxLL[T11])/3.0;
    qLL[Q112] = alpha*vthLL*rhoLL*(dTdxLL[T12] + dTdxLL[T12] + dTdyLL[T11])/3.0;
    qLL[Q113] = alpha*vthLL*rhoLL*(dTdxLL[T13] + dTdxLL[T13])/3.0;
    qLL[Q122] = alpha*vthLL*rhoLL*(dTdxLL[T22] + dTdyLL[T12] + dTdyLL[T12])/3.0;
    qLL[Q123] = alpha*vthLL*rhoLL*(dTdxLL[T23] + dTdyLL[T13])/3.0;
    qLL[Q133] = alpha*vthLL*rhoLL*(dTdxLL[T33])/3.0; 
    qLL[Q222] = alpha*vthLL*rhoLL*(dTdyLL[T22] + dTdyLL[T22] + dTdyLL[T22])/3.0;
    qLL[Q223] = alpha*vthLL*rhoLL*(dTdyLL[T23] + dTdyLL[T23])/3.0;
    qLL[Q233] = alpha*vthLL*rhoLL*(dTdyLL[T33])/3.0;

    qLU[Q111] = alpha*vthLU*rhoLU*(dTdxLU[T11] + dTdxLU[T11] + dTdxLU[T11])/3.0;
    qLU[Q112] = alpha*vthLU*rhoLU*(dTdxLU[T12] + dTdxLU[T12] + dTdyLU[T11])/3.0;
    qLU[Q113] = alpha*vthLU*rhoLU*(dTdxLU[T13] + dTdxLU[T13])/3.0;
    qLU[Q122] = alpha*vthLU*rhoLU*(dTdxLU[T22] + dTdyLU[T12] + dTdyLU[T12])/3.0;
    qLU[Q123] = alpha*vthLU*rhoLU*(dTdxLU[T23] + dTdyLU[T13])/3.0;
    qLU[Q133] = alpha*vthLU*rhoLU*(dTdxLU[T33])/3.0; 
    qLU[Q222] = alpha*vthLU*rhoLU*(dTdyLU[T22] + dTdyLU[T22] + dTdyLU[T22])/3.0;
    qLU[Q223] = alpha*vthLU*rhoLU*(dTdyLU[T23] + dTdyLU[T23])/3.0;
    qLU[Q233] = alpha*vthLU*rhoLU*(dTdyLU[T33])/3.0;

    qUL[Q111] = alpha*vthUL*rhoUL*(dTdxUL[T11] + dTdxUL[T11] + dTdxUL[T11])/3.0;
    qUL[Q112] = alpha*vthUL*rhoUL*(dTdxUL[T12] + dTdxUL[T12] + dTdyUL[T11])/3.0;
    qUL[Q113] = alpha*vthUL*rhoUL*(dTdxUL[T13] + dTdxUL[T13])/3.0;
    qUL[Q122] = alpha*vthUL*rhoUL*(dTdxUL[T22] + dTdyUL[T12] + dTdyUL[T12])/3.0;
    qUL[Q123] = alpha*vthUL*rhoUL*(dTdxUL[T23] + dTdyUL[T13])/3.0;
    qUL[Q133] = alpha*vthUL*rhoUL*(dTdxUL[T33])/3.0; 
    qUL[Q222] = alpha*vthUL*rhoUL*(dTdyUL[T22] + dTdyUL[T22] + dTdyUL[T22])/3.0;
    qUL[Q223] = alpha*vthUL*rhoUL*(dTdyUL[T23] + dTdyUL[T23])/3.0;
    qUL[Q233] = alpha*vthUL*rhoUL*(dTdyUL[T33])/3.0;

    qUU[Q111] = alpha*vthUU*rhoUU*(dTdxUU[T11] + dTdxUU[T11] + dTdxUU[T11])/3.0;
    qUU[Q112] = alpha*vthUU*rhoUU*(dTdxUU[T12] + dTdxUU[T12] + dTdyUU[T11])/3.0;
    qUU[Q113] = alpha*vthUU*rhoUU*(dTdxUU[T13] + dTdxUU[T13])/3.0;
    qUU[Q122] = alpha*vthUU*rhoUU*(dTdxUU[T22] + dTdyUU[T12] + dTdyUU[T12])/3.0;
    qUU[Q123] = alpha*vthUU*rhoUU*(dTdxUU[T23] + dTdyUU[T13])/3.0;
    qUU[Q133] = alpha*vthUU*rhoUU*(dTdxUU[T33])/3.0; 
    qUU[Q222] = alpha*vthUU*rhoUU*(dTdyUU[T22] + dTdyUU[T22] + dTdyUU[T22])/3.0;
    qUU[Q223] = alpha*vthUU*rhoUU*(dTdyUU[T23] + dTdyUU[T23])/3.0;
    qUU[Q233] = alpha*vthUU*rhoUU*(dTdyUU[T33])/3.0;

    double divQx[6] = { calc_sym_gradx_2D(dx, qLL[Q111], qLU[Q111], qUL[Q111], qUU[Q111]), 
                        calc_sym_gradx_2D(dx, qLL[Q112], qLU[Q112], qUL[Q112], qUU[Q112]), 
                        calc_sym_gradx_2D(dx, qLL[Q113], qLU[Q113], qUL[Q113], qUU[Q113]), 
                        calc_sym_gradx_2D(dx, qLL[Q122], qLU[Q122], qUL[Q122], qUU[Q122]), 
                        calc_sym_gradx_2D(dx, qLL[Q123], qLU[Q123], qUL[Q123], qUU[Q123]), 
                        calc_sym_gradx_2D(dx, qLL[Q133], qLU[Q133], qUL[Q133], qUU[Q133]) };

    double divQy[6] = { calc_sym_grady_2D(dy, qLL[Q112], qLU[Q112], qUL[Q112], qUU[Q112]), 
                        calc_sym_grady_2D(dy, qLL[Q122], qLU[Q122], qUL[Q122], qUU[Q122]), 
                        calc_sym_grady_2D(dy, qLL[Q123], qLU[Q123], qUL[Q123], qUU[Q123]), 
                        calc_sym_grady_2D(dy, qLL[Q222], qLU[Q222], qUL[Q222], qUU[Q222]), 
                        calc_sym_grady_2D(dy, qLL[Q223], qLU[Q223], qUL[Q223], qUU[Q223]), 
                        calc_sym_grady_2D(dy, qLL[Q233], qLU[Q233], qUL[Q233], qUU[Q233]) };

    rhs[0] = divQx[0] + divQy[0];
    rhs[1] = divQx[1] + divQy[1];
    rhs[2] = divQx[2] + divQy[2];
    rhs[3] = divQx[3] + divQy[3];
    rhs[4] = divQx[4] + divQy[4];
    rhs[5] = divQx[5] + divQy[5];
  }
}

gkyl_ten_moment_grad_closure*
gkyl_ten_moment_grad_closure_new(struct gkyl_ten_moment_grad_closure_inp inp)
{
  gkyl_ten_moment_grad_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_grad_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;

  return up;
}

void
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces, struct gkyl_range update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_range.ndim;
  long sz[] = { 3, 9, 27 };

  long offsets[sz[ndim-1]];
  create_offsets(&update_range, offsets);

  const double* fluid_d[sz[ndim-1]];
  double *rhs_d;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  while (gkyl_range_iter_next(&iter)) {
    
    long linc = gkyl_range_idx(&update_range, iter.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i)
      fluid_d[i] =  gkyl_array_cfetch(fluid, linc + offsets[i]); 

    rhs_d = gkyl_array_fetch(rhs, linc);
    unmag_grad_closure_update(gces, fluid_d, gkyl_array_fetch(cflrate, linc), rhs_d);
  }
}

void
gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure* up)
{
  free(up);
}
