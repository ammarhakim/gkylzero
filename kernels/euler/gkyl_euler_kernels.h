#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void fluid_vars_pressure_1x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_limiter_1x_ser_p1(double gas_gamma, double *p_c, double *fluid_l, double *fluid_c, double *fluid_r);
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_1x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_1x_ser_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_limiter_1x_ser_p2(double gas_gamma, double *p_c, double *fluid_l, double *fluid_c, double *fluid_r);
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_1x_ser_p2(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_2x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH int fluid_vars_u_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_3x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH int fluid_vars_u_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfz_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_2x_tensor_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH int fluid_vars_u_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 

inline static void
roe_avg_euler(const double *ql, const double *qr, const double *pl, const double *pr, double *avg)
{
  double rhol = ql[0], rhor = qr[0];

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[1]*ravgl1 + qr[1]*ravgr1)*ravg2;
  double v = (ql[2]*ravgl1 + qr[2]*ravgr1)*ravg2;
  double w = (ql[3]*ravgl1 + qr[3]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl[0])*ravgl1 + (qr[4]+pr[0])*ravgr1)*ravg2;

  avg[0] = u; avg[1] = v; avg[2] = w; avg[3] = enth;
}

// evaluated at avg = {u, v, w, enth}
inline static void
waves_euler(double gas_gamma, const double *delta, const double *avg, double *waves, double *s)
{
  double g1 = gas_gamma - 1;

  double u = avg[0], v = avg[1], w = avg[2], enth = avg[3];

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;

  // Compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  double a2 = delta[2] - v*delta[0];
  double a3 = delta[3] - w*delta[0];
  double a5 = 0.5*(delta[1] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  double *wv;
  // Wave 1: eigenvalue is u-c
  wv = &waves[0];
  wv[0] = a1;
  wv[1] = a1*(u-a);
  wv[2] = a1*v;
  wv[3] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a4;
  wv[1] = a4*u;
  wv[2] = a4*v + a2;
  wv[3] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wv = &waves[10];
  wv[0] = a5;
  wv[1] = a5*(u+a);
  wv[2] = a5*v;
  wv[3] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
}

inline static void
qfluct_euler(const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5], *w2 = &waves[10];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
}

inline static double
minmod(double a, double b, double c)
{
  double sa = GKYL_SGN(a);
  double sb = GKYL_SGN(b);
  double sc = GKYL_SGN(c);
  if( (sa==sb) && (sb==sc) ) {
    if (sa<0)
      return GKYL_MAX(GKYL_MAX(a,b),c);
    else
      return GKYL_MIN(GKYL_MIN(a,b),c);
  }
  else {
     return 0;
  }
}

EXTERN_C_END 
