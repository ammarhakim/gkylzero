#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
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
  const double *fluid_d[], double *cflrate, double *rhs_d)
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
                        
    double qL[10] = {0.0};
    double qU[10] = {0.0};
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

    rhs_d[RHO] = 0.0;
    rhs_d[MX] = 0.0;
    rhs_d[MY] = 0.0;
    rhs_d[MZ] = 0.0;
    rhs_d[P11] = (qU[Q111] - qL[Q111])/dx;
    rhs_d[P12] = (qU[Q112] - qL[Q112])/dx;
    rhs_d[P13] = (qU[Q113] - qL[Q113])/dx;
    rhs_d[P22] = (qU[Q122] - qL[Q122])/dx;
    rhs_d[P23] = (qU[Q123] - qL[Q123])/dx;
    rhs_d[P33] = (qU[Q133] - qL[Q133])/dx;
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
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

    double qLL[10] = {0.0};
    double qLU[10] = {0.0};
    double qUL[10] = {0.0};
    double qUU[10] = {0.0};
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

    rhs_d[RHO] = 0.0;
    rhs_d[MX] = 0.0;
    rhs_d[MY] = 0.0;
    rhs_d[MZ] = 0.0;
    rhs_d[P11] = divQx[0] + divQy[0];
    rhs_d[P12] = divQx[1] + divQy[1];
    rhs_d[P13] = divQx[2] + divQy[2];
    rhs_d[P22] = divQx[3] + divQy[3];
    rhs_d[P23] = divQx[4] + divQy[4];
    rhs_d[P33] = divQx[5] + divQy[5];
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];
    double Tij[27][6] = {0.0};
    double rho[27] = {0.0};
    double p[27] = {0.0}; 
    for (int j = LLL_3D; j <= UUU_3D; ++j) {
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

    double rhoLLL = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLC_3D], rho[LCL_3D], rho[LCC_3D], rho[CLL_3D], rho[CLC_3D], rho[CCL_3D], rho[CCC_3D]);
    double pLLL = calc_harmonic_avg_3D(p[LLL_3D], p[LLC_3D], p[LCL_3D], p[LCC_3D], p[CLL_3D], p[CLC_3D], p[CCL_3D], p[CCC_3D]);
    double rhoLLU = calc_harmonic_avg_3D(rho[LLC_3D], rho[LLU_3D], rho[LCC_3D], rho[LCU_3D], rho[CLC_3D], rho[CLU_3D], rho[CCC_3D], rho[CCU_3D]);
    double pLLU = calc_harmonic_avg_3D(p[LLC_3D], p[LLU_3D], p[LCC_3D], p[LCU_3D], p[CLC_3D], p[CLU_3D], p[CCC_3D], p[CCU_3D]);
    double rhoLUL = calc_harmonic_avg_3D(rho[LCL_3D], rho[LCC_3D], rho[LUL_3D], rho[LUC_3D], rho[CCL_3D], rho[CCC_3D], rho[CUL_3D], rho[CUC_3D]);
    double pLUL = calc_harmonic_avg_3D(p[LCL_3D], p[LCC_3D], p[LUL_3D], p[LUC_3D], p[CCL_3D], p[CCC_3D], p[CUL_3D], p[CUC_3D]);
    double rhoLUU = calc_harmonic_avg_3D(rho[LCC_3D], rho[LCU_3D], rho[LUC_3D], rho[LUU_3D], rho[CCC_3D], rho[CCU_3D], rho[CUC_3D], rho[CUU_3D]);
    double pLUU = calc_harmonic_avg_3D(p[LCC_3D], p[LCU_3D], p[LUC_3D], p[LUU_3D], p[CCC_3D], p[CCU_3D], p[CUC_3D], p[CUU_3D]);

    double rhoULL = calc_harmonic_avg_3D(rho[CLL_3D], rho[CLC_3D], rho[CCL_3D], rho[CCC_3D], rho[ULL_3D], rho[ULC_3D], rho[UCL_3D], rho[UCC_3D]);
    double pULL = calc_harmonic_avg_3D(p[CLL_3D], p[CLC_3D], p[CCL_3D], p[CCC_3D], p[ULL_3D], p[ULC_3D], p[UCL_3D], p[UCC_3D]);
    double rhoULU = calc_harmonic_avg_3D(rho[CLC_3D], rho[CLU_3D], rho[CCC_3D], rho[CCU_3D], rho[ULC_3D], rho[ULU_3D], rho[UCC_3D], rho[UCU_3D]);
    double pULU = calc_harmonic_avg_3D(p[CLC_3D], p[CLU_3D], p[CCC_3D], p[CCU_3D], p[ULC_3D], p[ULU_3D], p[UCC_3D], p[UCU_3D]);
    double rhoUUL = calc_harmonic_avg_3D(rho[CCL_3D], rho[CCC_3D], rho[CUL_3D], rho[CUC_3D], rho[UCL_3D], rho[UCC_3D], rho[UUL_3D], rho[UUC_3D]);
    double pUUL = calc_harmonic_avg_3D(p[CCL_3D], p[CCC_3D], p[CUL_3D], p[CUC_3D], p[UCL_3D], p[UCC_3D], p[UUL_3D], p[UUC_3D]);
    double rhoUUU = calc_harmonic_avg_3D(rho[CCC_3D], rho[CCU_3D], rho[CUC_3D], rho[CUU_3D], rho[UCC_3D], rho[UCU_3D], rho[UUC_3D], rho[UUU_3D]);
    double pUUU = calc_harmonic_avg_3D(p[CCC_3D], p[CCU_3D], p[CUC_3D], p[CUU_3D], p[UCC_3D], p[UCU_3D], p[UUC_3D], p[UUU_3D]);

    double vthLLL = sqrt(pLLL/rhoLLL);
    double vthLLU = sqrt(pLLU/rhoLLU);
    double vthLUL = sqrt(pLUL/rhoLUL);
    double vthLUU = sqrt(pLUU/rhoLUU);

    double vthULL = sqrt(pULL/rhoULL);
    double vthULU = sqrt(pULU/rhoULU);
    double vthUUL = sqrt(pUUL/rhoUUL);
    double vthUUU = sqrt(pUUU/rhoUUU);

    double dTdxLLL[6] = { calc_sym_gradx_3D(dx, Tij[LLL_3D][T11], Tij[LLC_3D][T11], Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[LLL_3D][T12], Tij[LLC_3D][T12], Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[LLL_3D][T13], Tij[LLC_3D][T13], Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[LLL_3D][T22], Tij[LLC_3D][T22], Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[LLL_3D][T23], Tij[LLC_3D][T23], Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[LLL_3D][T33], Tij[LLC_3D][T33], Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33]) };

    double dTdyLLL[6] = { calc_sym_grady_3D(dy, Tij[LLL_3D][T11], Tij[LLC_3D][T11], Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[LLL_3D][T12], Tij[LLC_3D][T12], Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[LLL_3D][T13], Tij[LLC_3D][T13], Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[LLL_3D][T22], Tij[LLC_3D][T22], Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[LLL_3D][T23], Tij[LLC_3D][T23], Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[LLL_3D][T33], Tij[LLC_3D][T33], Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33]) };

    double dTdzLLL[6] = { calc_sym_gradz_3D(dz, Tij[LLL_3D][T11], Tij[LLC_3D][T11], Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[LLL_3D][T12], Tij[LLC_3D][T12], Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[LLL_3D][T13], Tij[LLC_3D][T13], Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[LLL_3D][T22], Tij[LLC_3D][T22], Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[LLL_3D][T23], Tij[LLC_3D][T23], Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[LLL_3D][T33], Tij[LLC_3D][T33], Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33]) };

    double dTdxLLU[6] = { calc_sym_gradx_3D(dx, Tij[LLC_3D][T11], Tij[LLU_3D][T11], Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[LLC_3D][T12], Tij[LLU_3D][T12], Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[LLC_3D][T13], Tij[LLU_3D][T13], Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[LLC_3D][T22], Tij[LLU_3D][T22], Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[LLC_3D][T23], Tij[LLU_3D][T23], Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[LLC_3D][T33], Tij[LLU_3D][T33], Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33]) };

    double dTdyLLU[6] = { calc_sym_grady_3D(dy, Tij[LLC_3D][T11], Tij[LLU_3D][T11], Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[LLC_3D][T12], Tij[LLU_3D][T12], Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[LLC_3D][T13], Tij[LLU_3D][T13], Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[LLC_3D][T22], Tij[LLU_3D][T22], Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[LLC_3D][T23], Tij[LLU_3D][T23], Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[LLC_3D][T33], Tij[LLU_3D][T33], Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33]) };

    double dTdzLLU[6] = { calc_sym_gradz_3D(dz, Tij[LLC_3D][T11], Tij[LLU_3D][T11], Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[LLC_3D][T12], Tij[LLU_3D][T12], Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[LLC_3D][T13], Tij[LLU_3D][T13], Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[LLC_3D][T22], Tij[LLU_3D][T22], Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[LLC_3D][T23], Tij[LLU_3D][T23], Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[LLC_3D][T33], Tij[LLU_3D][T33], Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33]) };

    double dTdxLUL[6] = { calc_sym_gradx_3D(dx, Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[LUL_3D][T11], Tij[LUC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[LUL_3D][T12], Tij[LUC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[LUL_3D][T13], Tij[LUC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[LUL_3D][T22], Tij[LUC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[LUL_3D][T23], Tij[LUC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[LUL_3D][T33], Tij[LUC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33]) };

    double dTdyLUL[6] = { calc_sym_grady_3D(dy, Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[LUL_3D][T11], Tij[LUC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[LUL_3D][T12], Tij[LUC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[LUL_3D][T13], Tij[LUC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[LUL_3D][T22], Tij[LUC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[LUL_3D][T23], Tij[LUC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[LUL_3D][T33], Tij[LUC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33]) };

    double dTdzLUL[6] = { calc_sym_gradz_3D(dz, Tij[LCL_3D][T11], Tij[LCC_3D][T11], Tij[LUL_3D][T11], Tij[LUC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[LCL_3D][T12], Tij[LCC_3D][T12], Tij[LUL_3D][T12], Tij[LUC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[LCL_3D][T13], Tij[LCC_3D][T13], Tij[LUL_3D][T13], Tij[LUC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[LCL_3D][T22], Tij[LCC_3D][T22], Tij[LUL_3D][T22], Tij[LUC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[LCL_3D][T23], Tij[LCC_3D][T23], Tij[LUL_3D][T23], Tij[LUC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[LCL_3D][T33], Tij[LCC_3D][T33], Tij[LUL_3D][T33], Tij[LUC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33]) };

    double dTdxLUU[6] = { calc_sym_gradx_3D(dx, Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[LUC_3D][T11], Tij[LUU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[LUC_3D][T12], Tij[LUU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[LUC_3D][T13], Tij[LUU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[LUC_3D][T22], Tij[LUU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[LUC_3D][T23], Tij[LUU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[LUC_3D][T33], Tij[LUU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33]) };

    double dTdyLUU[6] = { calc_sym_grady_3D(dy, Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[LUC_3D][T11], Tij[LUU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[LUC_3D][T12], Tij[LUU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[LUC_3D][T13], Tij[LUU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[LUC_3D][T22], Tij[LUU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[LUC_3D][T23], Tij[LUU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[LUC_3D][T33], Tij[LUU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33]) };

    double dTdzLUU[6] = { calc_sym_gradz_3D(dz, Tij[LCC_3D][T11], Tij[LCU_3D][T11], Tij[LUC_3D][T11], Tij[LUU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[LCC_3D][T12], Tij[LCU_3D][T12], Tij[LUC_3D][T12], Tij[LUU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[LCC_3D][T13], Tij[LCU_3D][T13], Tij[LUC_3D][T13], Tij[LUU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[LCC_3D][T22], Tij[LCU_3D][T22], Tij[LUC_3D][T22], Tij[LUU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[LCC_3D][T23], Tij[LCU_3D][T23], Tij[LUC_3D][T23], Tij[LUU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[LCC_3D][T33], Tij[LCU_3D][T33], Tij[LUC_3D][T33], Tij[LUU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33]) };

    double dTdxULL[6] = { calc_sym_gradx_3D(dx, Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[ULL_3D][T11], Tij[ULC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[ULL_3D][T12], Tij[ULC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[ULL_3D][T13], Tij[ULC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[ULL_3D][T22], Tij[ULC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[ULL_3D][T23], Tij[ULC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[ULL_3D][T33], Tij[ULC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33]) };

    double dTdyULL[6] = { calc_sym_grady_3D(dy, Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[ULL_3D][T11], Tij[ULC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[ULL_3D][T12], Tij[ULC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[ULL_3D][T13], Tij[ULC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[ULL_3D][T22], Tij[ULC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[ULL_3D][T23], Tij[ULC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[ULL_3D][T33], Tij[ULC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33]) };

    double dTdzULL[6] = { calc_sym_gradz_3D(dz, Tij[CLL_3D][T11], Tij[CLC_3D][T11], Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[ULL_3D][T11], Tij[ULC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[CLL_3D][T12], Tij[CLC_3D][T12], Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[ULL_3D][T12], Tij[ULC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[CLL_3D][T13], Tij[CLC_3D][T13], Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[ULL_3D][T13], Tij[ULC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[CLL_3D][T22], Tij[CLC_3D][T22], Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[ULL_3D][T22], Tij[ULC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[CLL_3D][T23], Tij[CLC_3D][T23], Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[ULL_3D][T23], Tij[ULC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[CLL_3D][T33], Tij[CLC_3D][T33], Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[ULL_3D][T33], Tij[ULC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33]) };

    double dTdxULU[6] = { calc_sym_gradx_3D(dx, Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[ULC_3D][T11], Tij[ULU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[ULC_3D][T12], Tij[ULU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[ULC_3D][T13], Tij[ULU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[ULC_3D][T22], Tij[ULU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[ULC_3D][T23], Tij[ULU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[ULC_3D][T33], Tij[ULU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33]) };

    double dTdyULU[6] = { calc_sym_grady_3D(dy, Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[ULC_3D][T11], Tij[ULU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[ULC_3D][T12], Tij[ULU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[ULC_3D][T13], Tij[ULU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[ULC_3D][T22], Tij[ULU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[ULC_3D][T23], Tij[ULU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[ULC_3D][T33], Tij[ULU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33]) };

    double dTdzULU[6] = { calc_sym_gradz_3D(dz, Tij[CLC_3D][T11], Tij[CLU_3D][T11], Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[ULC_3D][T11], Tij[ULU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[CLC_3D][T12], Tij[CLU_3D][T12], Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[ULC_3D][T12], Tij[ULU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[CLC_3D][T13], Tij[CLU_3D][T13], Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[ULC_3D][T13], Tij[ULU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[CLC_3D][T22], Tij[CLU_3D][T22], Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[ULC_3D][T22], Tij[ULU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[CLC_3D][T23], Tij[CLU_3D][T23], Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[ULC_3D][T23], Tij[ULU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[CLC_3D][T33], Tij[CLU_3D][T33], Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[ULC_3D][T33], Tij[ULU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33]) };

    double dTdxUUL[6] = { calc_sym_gradx_3D(dx, Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11], Tij[UUL_3D][T11], Tij[UUC_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12], Tij[UUL_3D][T12], Tij[UUC_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13], Tij[UUL_3D][T13], Tij[UUC_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22], Tij[UUL_3D][T22], Tij[UUC_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23], Tij[UUL_3D][T23], Tij[UUC_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33], Tij[UUL_3D][T33], Tij[UUC_3D][T33]) };

    double dTdyUUL[6] = { calc_sym_grady_3D(dy, Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11], Tij[UUL_3D][T11], Tij[UUC_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12], Tij[UUL_3D][T12], Tij[UUC_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13], Tij[UUL_3D][T13], Tij[UUC_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22], Tij[UUL_3D][T22], Tij[UUC_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23], Tij[UUL_3D][T23], Tij[UUC_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33], Tij[UUL_3D][T33], Tij[UUC_3D][T33]) };

    double dTdzUUL[6] = { calc_sym_gradz_3D(dz, Tij[CCL_3D][T11], Tij[CCC_3D][T11], Tij[CUL_3D][T11], Tij[CUC_3D][T11], Tij[UCL_3D][T11], Tij[UCC_3D][T11], Tij[UUL_3D][T11], Tij[UUC_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[CCL_3D][T12], Tij[CCC_3D][T12], Tij[CUL_3D][T12], Tij[CUC_3D][T12], Tij[UCL_3D][T12], Tij[UCC_3D][T12], Tij[UUL_3D][T12], Tij[UUC_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[CCL_3D][T13], Tij[CCC_3D][T13], Tij[CUL_3D][T13], Tij[CUC_3D][T13], Tij[UCL_3D][T13], Tij[UCC_3D][T13], Tij[UUL_3D][T13], Tij[UUC_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[CCL_3D][T22], Tij[CCC_3D][T22], Tij[CUL_3D][T22], Tij[CUC_3D][T22], Tij[UCL_3D][T22], Tij[UCC_3D][T22], Tij[UUL_3D][T22], Tij[UUC_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[CCL_3D][T23], Tij[CCC_3D][T23], Tij[CUL_3D][T23], Tij[CUC_3D][T23], Tij[UCL_3D][T23], Tij[UCC_3D][T23], Tij[UUL_3D][T23], Tij[UUC_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[CCL_3D][T33], Tij[CCC_3D][T33], Tij[CUL_3D][T33], Tij[CUC_3D][T33], Tij[UCL_3D][T33], Tij[UCC_3D][T33], Tij[UUL_3D][T33], Tij[UUC_3D][T33]) };

    double dTdxUUU[6] = { calc_sym_gradx_3D(dx, Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11], Tij[UUC_3D][T11], Tij[UUU_3D][T11]),
                          calc_sym_gradx_3D(dx, Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12], Tij[UUC_3D][T12], Tij[UUU_3D][T12]),
                          calc_sym_gradx_3D(dx, Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13], Tij[UUC_3D][T13], Tij[UUU_3D][T13]),
                          calc_sym_gradx_3D(dx, Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22], Tij[UUC_3D][T22], Tij[UUU_3D][T22]),
                          calc_sym_gradx_3D(dx, Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23], Tij[UUC_3D][T23], Tij[UUU_3D][T23]),
                          calc_sym_gradx_3D(dx, Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33], Tij[UUC_3D][T33], Tij[UUU_3D][T33]) };

    double dTdyUUU[6] = { calc_sym_grady_3D(dy, Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11], Tij[UUC_3D][T11], Tij[UUU_3D][T11]),
                          calc_sym_grady_3D(dy, Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12], Tij[UUC_3D][T12], Tij[UUU_3D][T12]),
                          calc_sym_grady_3D(dy, Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13], Tij[UUC_3D][T13], Tij[UUU_3D][T13]),
                          calc_sym_grady_3D(dy, Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22], Tij[UUC_3D][T22], Tij[UUU_3D][T22]),
                          calc_sym_grady_3D(dy, Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23], Tij[UUC_3D][T23], Tij[UUU_3D][T23]),
                          calc_sym_grady_3D(dy, Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33], Tij[UUC_3D][T33], Tij[UUU_3D][T33]) };

    double dTdzUUU[6] = { calc_sym_gradz_3D(dz, Tij[CCC_3D][T11], Tij[CCU_3D][T11], Tij[CUC_3D][T11], Tij[CUU_3D][T11], Tij[UCC_3D][T11], Tij[UCU_3D][T11], Tij[UUC_3D][T11], Tij[UUU_3D][T11]),
                          calc_sym_gradz_3D(dz, Tij[CCC_3D][T12], Tij[CCU_3D][T12], Tij[CUC_3D][T12], Tij[CUU_3D][T12], Tij[UCC_3D][T12], Tij[UCU_3D][T12], Tij[UUC_3D][T12], Tij[UUU_3D][T12]),
                          calc_sym_gradz_3D(dz, Tij[CCC_3D][T13], Tij[CCU_3D][T13], Tij[CUC_3D][T13], Tij[CUU_3D][T13], Tij[UCC_3D][T13], Tij[UCU_3D][T13], Tij[UUC_3D][T13], Tij[UUU_3D][T13]),
                          calc_sym_gradz_3D(dz, Tij[CCC_3D][T22], Tij[CCU_3D][T22], Tij[CUC_3D][T22], Tij[CUU_3D][T22], Tij[UCC_3D][T22], Tij[UCU_3D][T22], Tij[UUC_3D][T22], Tij[UUU_3D][T22]),
                          calc_sym_gradz_3D(dz, Tij[CCC_3D][T23], Tij[CCU_3D][T23], Tij[CUC_3D][T23], Tij[CUU_3D][T23], Tij[UCC_3D][T23], Tij[UCU_3D][T23], Tij[UUC_3D][T23], Tij[UUU_3D][T23]),
                          calc_sym_gradz_3D(dz, Tij[CCC_3D][T33], Tij[CCU_3D][T33], Tij[CUC_3D][T33], Tij[CUU_3D][T33], Tij[UCC_3D][T33], Tij[UCU_3D][T33], Tij[UUC_3D][T33], Tij[UUU_3D][T33]) };

    double qLLL[10] = {0.0};
    double qLLU[10] = {0.0};
    double qLUL[10] = {0.0};
    double qLUU[10] = {0.0};
    double qULL[10] = {0.0};
    double qULU[10] = {0.0};
    double qUUL[10] = {0.0};
    double qUUU[10] = {0.0};
    double alpha = 1.0/gces->k0;

    qLLL[Q111] = alpha*vthLLL*rhoLLL*(dTdxLLL[T11] + dTdxLLL[T11] + dTdxLLL[T11])/3.0;
    qLLL[Q112] = alpha*vthLLL*rhoLLL*(dTdxLLL[T12] + dTdxLLL[T12] + dTdyLLL[T11])/3.0;
    qLLL[Q113] = alpha*vthLLL*rhoLLL*(dTdxLLL[T13] + dTdxLLL[T13] + dTdzLLL[T11])/3.0;
    qLLL[Q122] = alpha*vthLLL*rhoLLL*(dTdxLLL[T22] + dTdyLLL[T12] + dTdyLLL[T12])/3.0;
    qLLL[Q123] = alpha*vthLLL*rhoLLL*(dTdxLLL[T23] + dTdyLLL[T13] + dTdzLLL[T12])/3.0;
    qLLL[Q133] = alpha*vthLLL*rhoLLL*(dTdxLLL[T33] + dTdzLLL[T13] + dTdzLLL[T13])/3.0;
    qLLL[Q222] = alpha*vthLLL*rhoLLL*(dTdyLLL[T22] + dTdyLLL[T22] + dTdyLLL[T22])/3.0;
    qLLL[Q223] = alpha*vthLLL*rhoLLL*(dTdyLLL[T23] + dTdyLLL[T23] + dTdzLLL[T22])/3.0;
    qLLL[Q233] = alpha*vthLLL*rhoLLL*(dTdyLLL[T33] + dTdzLLL[T23] + dTdzLLL[T23])/3.0;
    qLLL[Q333] = alpha*vthLLL*rhoLLL*(dTdzLLL[T33] + dTdzLLL[T33] + dTdzLLL[T33])/3.0;

    qLLU[Q111] = alpha*vthLLU*rhoLLU*(dTdxLLU[T11] + dTdxLLU[T11] + dTdxLLU[T11])/3.0;
    qLLU[Q112] = alpha*vthLLU*rhoLLU*(dTdxLLU[T12] + dTdxLLU[T12] + dTdyLLU[T11])/3.0;
    qLLU[Q113] = alpha*vthLLU*rhoLLU*(dTdxLLU[T13] + dTdxLLU[T13] + dTdzLLU[T11])/3.0;
    qLLU[Q122] = alpha*vthLLU*rhoLLU*(dTdxLLU[T22] + dTdyLLU[T12] + dTdyLLU[T12])/3.0;
    qLLU[Q123] = alpha*vthLLU*rhoLLU*(dTdxLLU[T23] + dTdyLLU[T13] + dTdzLLU[T12])/3.0;
    qLLU[Q133] = alpha*vthLLU*rhoLLU*(dTdxLLU[T33] + dTdzLLU[T13] + dTdzLLU[T13])/3.0;
    qLLU[Q222] = alpha*vthLLU*rhoLLU*(dTdyLLU[T22] + dTdyLLU[T22] + dTdyLLU[T22])/3.0;
    qLLU[Q223] = alpha*vthLLU*rhoLLU*(dTdyLLU[T23] + dTdyLLU[T23] + dTdzLLU[T22])/3.0;
    qLLU[Q233] = alpha*vthLLU*rhoLLU*(dTdyLLU[T33] + dTdzLLU[T23] + dTdzLLU[T23])/3.0;
    qLLU[Q333] = alpha*vthLLU*rhoLLU*(dTdzLLU[T33] + dTdzLLU[T33] + dTdzLLU[T33])/3.0;

    qLUL[Q111] = alpha*vthLUL*rhoLUL*(dTdxLUL[T11] + dTdxLUL[T11] + dTdxLUL[T11])/3.0;
    qLUL[Q112] = alpha*vthLUL*rhoLUL*(dTdxLUL[T12] + dTdxLUL[T12] + dTdyLUL[T11])/3.0;
    qLUL[Q113] = alpha*vthLUL*rhoLUL*(dTdxLUL[T13] + dTdxLUL[T13] + dTdzLUL[T11])/3.0;
    qLUL[Q122] = alpha*vthLUL*rhoLUL*(dTdxLUL[T22] + dTdyLUL[T12] + dTdyLUL[T12])/3.0;
    qLUL[Q123] = alpha*vthLUL*rhoLUL*(dTdxLUL[T23] + dTdyLUL[T13] + dTdzLUL[T12])/3.0;
    qLUL[Q133] = alpha*vthLUL*rhoLUL*(dTdxLUL[T33] + dTdzLUL[T13] + dTdzLUL[T13])/3.0;
    qLUL[Q222] = alpha*vthLUL*rhoLUL*(dTdyLUL[T22] + dTdyLUL[T22] + dTdyLUL[T22])/3.0;
    qLUL[Q223] = alpha*vthLUL*rhoLUL*(dTdyLUL[T23] + dTdyLUL[T23] + dTdzLUL[T22])/3.0;
    qLUL[Q233] = alpha*vthLUL*rhoLUL*(dTdyLUL[T33] + dTdzLUL[T23] + dTdzLUL[T23])/3.0;
    qLUL[Q333] = alpha*vthLUL*rhoLUL*(dTdzLUL[T33] + dTdzLUL[T33] + dTdzLUL[T33])/3.0;

    qLUU[Q111] = alpha*vthLUU*rhoLUU*(dTdxLUU[T11] + dTdxLUU[T11] + dTdxLUU[T11])/3.0;
    qLUU[Q112] = alpha*vthLUU*rhoLUU*(dTdxLUU[T12] + dTdxLUU[T12] + dTdyLUU[T11])/3.0;
    qLUU[Q113] = alpha*vthLUU*rhoLUU*(dTdxLUU[T13] + dTdxLUU[T13] + dTdzLUU[T11])/3.0;
    qLUU[Q122] = alpha*vthLUU*rhoLUU*(dTdxLUU[T22] + dTdyLUU[T12] + dTdyLUU[T12])/3.0;
    qLUU[Q123] = alpha*vthLUU*rhoLUU*(dTdxLUU[T23] + dTdyLUU[T13] + dTdzLUU[T12])/3.0;
    qLUU[Q133] = alpha*vthLUU*rhoLUU*(dTdxLUU[T33] + dTdzLUU[T13] + dTdzLUU[T13])/3.0;
    qLUU[Q222] = alpha*vthLUU*rhoLUU*(dTdyLUU[T22] + dTdyLUU[T22] + dTdyLUU[T22])/3.0;
    qLUU[Q223] = alpha*vthLUU*rhoLUU*(dTdyLUU[T23] + dTdyLUU[T23] + dTdzLUU[T22])/3.0;
    qLUU[Q233] = alpha*vthLUU*rhoLUU*(dTdyLUU[T33] + dTdzLUU[T23] + dTdzLUU[T23])/3.0;
    qLUU[Q333] = alpha*vthLUU*rhoLUU*(dTdzLUU[T33] + dTdzLUU[T33] + dTdzLUU[T33])/3.0;

    qULL[Q111] = alpha*vthULL*rhoULL*(dTdxULL[T11] + dTdxULL[T11] + dTdxULL[T11])/3.0;
    qULL[Q112] = alpha*vthULL*rhoULL*(dTdxULL[T12] + dTdxULL[T12] + dTdyULL[T11])/3.0;
    qULL[Q113] = alpha*vthULL*rhoULL*(dTdxULL[T13] + dTdxULL[T13] + dTdzULL[T11])/3.0;
    qULL[Q122] = alpha*vthULL*rhoULL*(dTdxULL[T22] + dTdyULL[T12] + dTdyULL[T12])/3.0;
    qULL[Q123] = alpha*vthULL*rhoULL*(dTdxULL[T23] + dTdyULL[T13] + dTdzULL[T12])/3.0;
    qULL[Q133] = alpha*vthULL*rhoULL*(dTdxULL[T33] + dTdzULL[T13] + dTdzULL[T13])/3.0;
    qULL[Q222] = alpha*vthULL*rhoULL*(dTdyULL[T22] + dTdyULL[T22] + dTdyULL[T22])/3.0;
    qULL[Q223] = alpha*vthULL*rhoULL*(dTdyULL[T23] + dTdyULL[T23] + dTdzULL[T22])/3.0;
    qULL[Q233] = alpha*vthULL*rhoULL*(dTdyULL[T33] + dTdzULL[T23] + dTdzULL[T23])/3.0;
    qULL[Q333] = alpha*vthULL*rhoULL*(dTdzULL[T33] + dTdzULL[T33] + dTdzULL[T33])/3.0;

    qULU[Q111] = alpha*vthULU*rhoULU*(dTdxULU[T11] + dTdxULU[T11] + dTdxULU[T11])/3.0;
    qULU[Q112] = alpha*vthULU*rhoULU*(dTdxULU[T12] + dTdxULU[T12] + dTdyULU[T11])/3.0;
    qULU[Q113] = alpha*vthULU*rhoULU*(dTdxULU[T13] + dTdxULU[T13] + dTdzULU[T11])/3.0;
    qULU[Q122] = alpha*vthULU*rhoULU*(dTdxULU[T22] + dTdyULU[T12] + dTdyULU[T12])/3.0;
    qULU[Q123] = alpha*vthULU*rhoULU*(dTdxULU[T23] + dTdyULU[T13] + dTdzULU[T12])/3.0;
    qULU[Q133] = alpha*vthULU*rhoULU*(dTdxULU[T33] + dTdzULU[T13] + dTdzULU[T13])/3.0;
    qULU[Q222] = alpha*vthULU*rhoULU*(dTdyULU[T22] + dTdyULU[T22] + dTdyULU[T22])/3.0;
    qULU[Q223] = alpha*vthULU*rhoULU*(dTdyULU[T23] + dTdyULU[T23] + dTdzULU[T22])/3.0;
    qULU[Q233] = alpha*vthULU*rhoULU*(dTdyULU[T33] + dTdzULU[T23] + dTdzULU[T23])/3.0;
    qULU[Q333] = alpha*vthULU*rhoULU*(dTdzULU[T33] + dTdzULU[T33] + dTdzULU[T33])/3.0;

    qUUL[Q111] = alpha*vthUUL*rhoUUL*(dTdxUUL[T11] + dTdxUUL[T11] + dTdxUUL[T11])/3.0;
    qUUL[Q112] = alpha*vthUUL*rhoUUL*(dTdxUUL[T12] + dTdxUUL[T12] + dTdyUUL[T11])/3.0;
    qUUL[Q113] = alpha*vthUUL*rhoUUL*(dTdxUUL[T13] + dTdxUUL[T13] + dTdzUUL[T11])/3.0;
    qUUL[Q122] = alpha*vthUUL*rhoUUL*(dTdxUUL[T22] + dTdyUUL[T12] + dTdyUUL[T12])/3.0;
    qUUL[Q123] = alpha*vthUUL*rhoUUL*(dTdxUUL[T23] + dTdyUUL[T13] + dTdzUUL[T12])/3.0;
    qUUL[Q133] = alpha*vthUUL*rhoUUL*(dTdxUUL[T33] + dTdzUUL[T13] + dTdzUUL[T13])/3.0;
    qUUL[Q222] = alpha*vthUUL*rhoUUL*(dTdyUUL[T22] + dTdyUUL[T22] + dTdyUUL[T22])/3.0;
    qUUL[Q223] = alpha*vthUUL*rhoUUL*(dTdyUUL[T23] + dTdyUUL[T23] + dTdzUUL[T22])/3.0;
    qUUL[Q233] = alpha*vthUUL*rhoUUL*(dTdyUUL[T33] + dTdzUUL[T23] + dTdzUUL[T23])/3.0;
    qUUL[Q333] = alpha*vthUUL*rhoUUL*(dTdzUUL[T33] + dTdzUUL[T33] + dTdzUUL[T33])/3.0;

    qUUU[Q111] = alpha*vthUUU*rhoUUU*(dTdxUUU[T11] + dTdxUUU[T11] + dTdxUUU[T11])/3.0;
    qUUU[Q112] = alpha*vthUUU*rhoUUU*(dTdxUUU[T12] + dTdxUUU[T12] + dTdyUUU[T11])/3.0;
    qUUU[Q113] = alpha*vthUUU*rhoUUU*(dTdxUUU[T13] + dTdxUUU[T13] + dTdzUUU[T11])/3.0;
    qUUU[Q122] = alpha*vthUUU*rhoUUU*(dTdxUUU[T22] + dTdyUUU[T12] + dTdyUUU[T12])/3.0;
    qUUU[Q123] = alpha*vthUUU*rhoUUU*(dTdxUUU[T23] + dTdyUUU[T13] + dTdzUUU[T12])/3.0;
    qUUU[Q133] = alpha*vthUUU*rhoUUU*(dTdxUUU[T33] + dTdzUUU[T13] + dTdzUUU[T13])/3.0;
    qUUU[Q222] = alpha*vthUUU*rhoUUU*(dTdyUUU[T22] + dTdyUUU[T22] + dTdyUUU[T22])/3.0;
    qUUU[Q223] = alpha*vthUUU*rhoUUU*(dTdyUUU[T23] + dTdyUUU[T23] + dTdzUUU[T22])/3.0;
    qUUU[Q233] = alpha*vthUUU*rhoUUU*(dTdyUUU[T33] + dTdzUUU[T23] + dTdzUUU[T23])/3.0;
    qUUU[Q333] = alpha*vthUUU*rhoUUU*(dTdzUUU[T33] + dTdzUUU[T33] + dTdzUUU[T33])/3.0;

    double divQx[6] = { calc_sym_gradx_3D(dx, qLLL[Q111], qLLU[Q111], qLUL[Q111], qLUU[Q111], qULL[Q111], qULU[Q111], qUUL[Q111], qUUU[Q111]), 
                        calc_sym_gradx_3D(dx, qLLL[Q112], qLLU[Q112], qLUL[Q112], qLUU[Q112], qULL[Q112], qULU[Q112], qUUL[Q112], qUUU[Q112]), 
                        calc_sym_gradx_3D(dx, qLLL[Q113], qLLU[Q113], qLUL[Q113], qLUU[Q113], qULL[Q113], qULU[Q113], qUUL[Q113], qUUU[Q113]), 
                        calc_sym_gradx_3D(dx, qLLL[Q122], qLLU[Q122], qLUL[Q122], qLUU[Q122], qULL[Q122], qULU[Q122], qUUL[Q122], qUUU[Q122]), 
                        calc_sym_gradx_3D(dx, qLLL[Q123], qLLU[Q123], qLUL[Q123], qLUU[Q123], qULL[Q123], qULU[Q123], qUUL[Q123], qUUU[Q123]), 
                        calc_sym_gradx_3D(dx, qLLL[Q133], qLLU[Q133], qLUL[Q133], qLUU[Q133], qULL[Q133], qULU[Q133], qUUL[Q133], qUUU[Q133]) };

    double divQy[6] = { calc_sym_grady_3D(dy, qLLL[Q112], qLLU[Q112], qLUL[Q112], qLUU[Q112], qULL[Q112], qULU[Q112], qUUL[Q112], qUUU[Q112]), 
                        calc_sym_grady_3D(dy, qLLL[Q122], qLLU[Q122], qLUL[Q122], qLUU[Q122], qULL[Q122], qULU[Q122], qUUL[Q122], qUUU[Q122]), 
                        calc_sym_grady_3D(dy, qLLL[Q123], qLLU[Q123], qLUL[Q123], qLUU[Q123], qULL[Q123], qULU[Q123], qUUL[Q123], qUUU[Q123]), 
                        calc_sym_grady_3D(dy, qLLL[Q222], qLLU[Q222], qLUL[Q222], qLUU[Q222], qULL[Q222], qULU[Q222], qUUL[Q222], qUUU[Q222]), 
                        calc_sym_grady_3D(dy, qLLL[Q223], qLLU[Q223], qLUL[Q223], qLUU[Q223], qULL[Q223], qULU[Q223], qUUL[Q223], qUUU[Q223]), 
                        calc_sym_grady_3D(dy, qLLL[Q233], qLLU[Q233], qLUL[Q233], qLUU[Q233], qULL[Q233], qULU[Q233], qUUL[Q233], qUUU[Q233]) };

    double divQz[6] = { calc_sym_gradz_3D(dz, qLLL[Q113], qLLU[Q113], qLUL[Q113], qLUU[Q113], qULL[Q113], qULU[Q113], qUUL[Q113], qUUU[Q113]), 
                        calc_sym_gradz_3D(dz, qLLL[Q123], qLLU[Q123], qLUL[Q123], qLUU[Q123], qULL[Q123], qULU[Q123], qUUL[Q123], qUUU[Q123]), 
                        calc_sym_gradz_3D(dz, qLLL[Q133], qLLU[Q133], qLUL[Q133], qLUU[Q133], qULL[Q133], qULU[Q133], qUUL[Q133], qUUU[Q133]), 
                        calc_sym_gradz_3D(dz, qLLL[Q223], qLLU[Q223], qLUL[Q223], qLUU[Q223], qULL[Q223], qULU[Q223], qUUL[Q223], qUUU[Q223]), 
                        calc_sym_gradz_3D(dz, qLLL[Q233], qLLU[Q233], qLUL[Q233], qLUU[Q233], qULL[Q233], qULU[Q233], qUUL[Q233], qUUU[Q233]), 
                        calc_sym_gradz_3D(dz, qLLL[Q333], qLLU[Q333], qLUL[Q333], qLUU[Q333], qULL[Q333], qULU[Q333], qUUL[Q333], qUUU[Q333]) };

    rhs_d[RHO] = 0.0;
    rhs_d[MX] = 0.0;
    rhs_d[MY] = 0.0;
    rhs_d[MZ] = 0.0;
    rhs_d[P11] = divQx[0] + divQy[0] + divQz[0];
    rhs_d[P12] = divQx[1] + divQy[1] + divQz[1];
    rhs_d[P13] = divQx[2] + divQy[2] + divQz[2];
    rhs_d[P22] = divQx[3] + divQy[3] + divQz[3];
    rhs_d[P23] = divQx[4] + divQy[4] + divQz[4];
    rhs_d[P33] = divQx[5] + divQy[5] + divQz[5];
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
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces, 
  const struct gkyl_range *update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 3, 9, 27 };

  long offsets[sz[ndim-1]];
  create_offsets(update_range, offsets);

  const double* fluid_d[sz[ndim-1]];
  double *rhs_d;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    
    long linc = gkyl_range_idx(update_range, iter.idx);
    
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
