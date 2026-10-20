#include <gkyl_pkpm_vars_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH int pkpm_vars_closure_2x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *xi, const double *grad_xi,
  const double *theta, const double *dtB_over_magB, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl,
  double* GKYL_RESTRICT cond) 
{ 
  // w_lo, w_up:   Velocity domain edges; dv: velocity cell size. 
  // norm_moms:    Input [M0, M2, M0_G, M1, M3, M1_G, mu (x, y, z)] moments of the normalized distributions (stage state). 
  // rhs_moms:     Input moments of the stage-A right-hand sides, same layout. 
  // bvar, grad_b, u, grad_u, xi, grad_xi, theta, dtB_over_magB: fields of the forces. 
  // alpha_0_prev, c_h_prev: closures of the previous stage (may alias the outputs), for the edge speeds. 
  // f_lo, f_up:   [F, G, M0_perp, N_perp] of the lowest and uppermost velocity cells of this configuration cell. 
  // c_0, alpha_0, c_h, a_perp: Output closures at the 4 Gauss-Legendre nodes, interpolated to p=1. 
  // flux_lo, flux_up: Output edge fluxes [8 blocks] in the stored face format. 
  // cfl:          Output [lower, upper] edge-flux time-step rates. 
  // cond:         Output [1]: the smallest determinant of the 3x3 system over the nodes relative to 
  //               its edge-free value (1 when the edge terms vanish; small or negative when the 
  //               velocity-domain edges carry a non-negligible trace: too small velocity extents). 
  // Returns the number of nodes where the closures are undefined (everything is zero in the cell then). 
 
  const double dv1 = 2.0/dv; 
  const double *n0 = &norm_moms[0]; 
  const double *n2 = &norm_moms[4]; 
  const double *nG = &norm_moms[8]; 
  const double *n1 = &norm_moms[12]; 
  const double *r0 = &rhs_moms[0]; 
  const double *r2 = &rhs_moms[4]; 
  const double *rG = &rhs_moms[8]; 
  const double *r1 = &rhs_moms[12]; 
  const double *mu0 = &norm_moms[24]; 
  const double *mu1 = &norm_moms[28]; 
  const double *mu2 = &norm_moms[32]; 
  const double *rM0 = &rhs_moms[24]; 
  const double *rM1 = &rhs_moms[28]; 
  const double *rM2 = &rhs_moms[32]; 
  const double *b0 = &bvar[0]; 
  const double *b1 = &bvar[4]; 
  const double *b2 = &bvar[8]; 
  const double *u0 = &u[0]; 
  const double *u1 = &u[4]; 
  const double *u2 = &u[8]; 
  const double *dtb0 = &dtB_over_magB[0]; 
  const double *dtb1 = &dtB_over_magB[4]; 
  const double *dtb2 = &dtB_over_magB[8]; 
  const double *gux0 = &grad_u[0]; 
  const double *gux1 = &grad_u[4]; 
  const double *gux2 = &grad_u[8]; 
  const double *guy0 = &grad_u[12]; 
  const double *guy1 = &grad_u[16]; 
  const double *guy2 = &grad_u[20]; 
  const double *gbx0 = &grad_b[0]; 
  const double *gbx1 = &grad_b[4]; 
  const double *gbx2 = &grad_b[8]; 
  const double *gby0 = &grad_b[12]; 
  const double *gby1 = &grad_b[16]; 
  const double *gby2 = &grad_b[20]; 
  const double *gxix = &grad_xi[0]; 
  const double *gxiy = &grad_xi[12]; 
  const double *xi_1 = &xi[0]; 
  const double *lo0 = &f_lo[0], *up0 = &f_up[0]; 
  const double *lo1 = &f_lo[12], *up1 = &f_up[12]; 
  const double *lo2 = &f_lo[24], *up2 = &f_up[24]; 
  const double *lo3 = &f_lo[36], *up3 = &f_up[36]; 
  const double *lo4 = &f_lo[48], *up4 = &f_up[48]; 
  const double *lo5 = &f_lo[60], *up5 = &f_up[60]; 
  const double *lo6 = &f_lo[72], *up6 = &f_up[72]; 
  const double *lo7 = &f_lo[84], *up7 = &f_up[84]; 
 
  double c0_n[4] = {0.0}, a0_n[4] = {0.0}, c_n[4] = {0.0}; 
  double ap_n[3][4] = {{0.0}}, flo_n[8][4] = {{0.0}}, fup_n[8][4] = {{0.0}}; 
  double cfl_lo = 0.0, cfl_up = 0.0; 
  int num_null = 0; 
  double cond_min = 1.0; 
 
  { 
    const double n0_q = 0.5*n0[3]-0.5*n0[2]-0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = 0.5*n1[3]-0.5*n1[2]-0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = 0.5*n2[3]-0.5*n2[2]-0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = 0.5*nG[3]-0.5*nG[2]-0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = 0.5*r0[3]-0.5*r0[2]-0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = 0.5*r1[3]-0.5*r1[2]-0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = 0.5*r2[3]-0.5*r2[2]-0.5*r2[1]+0.5*r2[0] + 2.0*(0.5*rG[3]-0.5*rG[2]-0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {0.5*mu0[3]-0.5*mu0[2]-0.5*mu0[1]+0.5*mu0[0], 0.5*mu1[3]-0.5*mu1[2]-0.5*mu1[1]+0.5*mu1[0], 0.5*mu2[3]-0.5*mu2[2]-0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {0.5*rM0[3]-0.5*rM0[2]-0.5*rM0[1]+0.5*rM0[0], 0.5*rM1[3]-0.5*rM1[2]-0.5*rM1[1]+0.5*rM1[0], 0.5*rM2[3]-0.5*rM2[2]-0.5*rM2[1]+0.5*rM2[0]}; 
    const double a0p_q = 0.5*alpha_0_prev[3]-0.5*alpha_0_prev[2]-0.5*alpha_0_prev[1]+0.5*alpha_0_prev[0]; 
    const double cp_q = 0.5*c_h_prev[3]-0.5*c_h_prev[2]-0.5*c_h_prev[1]+0.5*c_h_prev[0]; 
    const double xi_q = 0.5*xi_1[3]-0.5*xi_1[2]-0.5*xi_1[1]+0.5*xi_1[0]; 
    const double th_q = 0.5*theta[3]-0.5*theta[2]-0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = 0.5*b0[3]-0.5*b0[2]-0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = 0.5*b1[3]-0.5*b1[2]-0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = 0.5*b2[3]-0.5*b2[2]-0.5*b2[1]+0.5*b2[0]; 
    const double u_q[3] = {0.5*u0[3]-0.5*u0[2]-0.5*u0[1]+0.5*u0[0], 0.5*u1[3]-0.5*u1[2]-0.5*u1[1]+0.5*u1[0], 0.5*u2[3]-0.5*u2[2]-0.5*u2[1]+0.5*u2[0]}; 
    const double d_q[3] = {0.5*dtb0[3]-0.5*dtb0[2]-0.5*dtb0[1]+0.5*dtb0[0], 0.5*dtb1[3]-0.5*dtb1[2]-0.5*dtb1[1]+0.5*dtb1[0], 0.5*dtb2[3]-0.5*dtb2[2]-0.5*dtb2[1]+0.5*dtb2[0]}; 
    const double gu00_q = 0.5*gux0[3]-0.5*gux0[2]-0.5*gux0[1]+0.5*gux0[0]; 
    const double gu01_q = 0.5*gux1[3]-0.5*gux1[2]-0.5*gux1[1]+0.5*gux1[0]; 
    const double gu02_q = 0.5*gux2[3]-0.5*gux2[2]-0.5*gux2[1]+0.5*gux2[0]; 
    const double gu10_q = 0.5*guy0[3]-0.5*guy0[2]-0.5*guy0[1]+0.5*guy0[0]; 
    const double gu11_q = 0.5*guy1[3]-0.5*guy1[2]-0.5*guy1[1]+0.5*guy1[0]; 
    const double gu12_q = 0.5*guy2[3]-0.5*guy2[2]-0.5*guy2[1]+0.5*guy2[0]; 
    const double gxi0_q = 0.5*gxix[3]-0.5*gxix[2]-0.5*gxix[1]+0.5*gxix[0]; 
    const double gxi1_q = 0.5*gxiy[3]-0.5*gxiy[2]-0.5*gxiy[1]+0.5*gxiy[0]; 
    const double gb_q[6] = {0.5*gbx0[3]-0.5*gbx0[2]-0.5*gbx0[1]+0.5*gbx0[0], 0.5*gbx1[3]-0.5*gbx1[2]-0.5*gbx1[1]+0.5*gbx1[0], 0.5*gbx2[3]-0.5*gbx2[2]-0.5*gbx2[1]+0.5*gbx2[0], 0.5*gby0[3]-0.5*gby0[2]-0.5*gby0[1]+0.5*gby0[0], 0.5*gby1[3]-0.5*gby1[2]-0.5*gby1[1]+0.5*gby1[0], 0.5*gby2[3]-0.5*gby2[2]-0.5*gby2[1]+0.5*gby2[0]}; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q + b1_q*b0_q*gu10_q + b1_q*b1_q*gu11_q + b1_q*b2_q*gu12_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q + b1_q*gxi1_q; 
    const double Ulo_q[8] = {0.7905694150420947*lo0[11]-0.7905694150420948*lo0[10]-0.7905694150420948*lo0[9]-0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]+0.6123724356957944*lo0[6]+0.6123724356957944*lo0[5]+0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]-0.3535533905932737*lo0[2]-0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], 0.7905694150420947*lo1[11]-0.7905694150420948*lo1[10]-0.7905694150420948*lo1[9]-0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]+0.6123724356957944*lo1[6]+0.6123724356957944*lo1[5]+0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]-0.3535533905932737*lo1[2]-0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], 0.7905694150420947*lo2[11]-0.7905694150420948*lo2[10]-0.7905694150420948*lo2[9]-0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]+0.6123724356957944*lo2[6]+0.6123724356957944*lo2[5]+0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]-0.3535533905932737*lo2[2]-0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], 0.7905694150420947*lo3[11]-0.7905694150420948*lo3[10]-0.7905694150420948*lo3[9]-0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]+0.6123724356957944*lo3[6]+0.6123724356957944*lo3[5]+0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]-0.3535533905932737*lo3[2]-0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], 0.7905694150420947*lo4[11]-0.7905694150420948*lo4[10]-0.7905694150420948*lo4[9]-0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]+0.6123724356957944*lo4[6]+0.6123724356957944*lo4[5]+0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]-0.3535533905932737*lo4[2]-0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], 0.7905694150420947*lo5[11]-0.7905694150420948*lo5[10]-0.7905694150420948*lo5[9]-0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]+0.6123724356957944*lo5[6]+0.6123724356957944*lo5[5]+0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]-0.3535533905932737*lo5[2]-0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], 0.7905694150420947*lo6[11]-0.7905694150420948*lo6[10]-0.7905694150420948*lo6[9]-0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]+0.6123724356957944*lo6[6]+0.6123724356957944*lo6[5]+0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]-0.3535533905932737*lo6[2]-0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], 0.7905694150420947*lo7[11]-0.7905694150420948*lo7[10]-0.7905694150420948*lo7[9]-0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]+0.6123724356957944*lo7[6]+0.6123724356957944*lo7[5]+0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]-0.3535533905932737*lo7[2]-0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {0.7905694150420947*up0[11]-0.7905694150420948*up0[10]-0.7905694150420948*up0[9]+0.6123724356957944*up0[8]+0.7905694150420947*up0[7]-0.6123724356957944*up0[6]-0.6123724356957944*up0[5]+0.3535533905932737*up0[4]+0.6123724356957944*up0[3]-0.3535533905932737*up0[2]-0.3535533905932737*up0[1]+0.3535533905932737*up0[0], 0.7905694150420947*up1[11]-0.7905694150420948*up1[10]-0.7905694150420948*up1[9]+0.6123724356957944*up1[8]+0.7905694150420947*up1[7]-0.6123724356957944*up1[6]-0.6123724356957944*up1[5]+0.3535533905932737*up1[4]+0.6123724356957944*up1[3]-0.3535533905932737*up1[2]-0.3535533905932737*up1[1]+0.3535533905932737*up1[0], 0.7905694150420947*up2[11]-0.7905694150420948*up2[10]-0.7905694150420948*up2[9]+0.6123724356957944*up2[8]+0.7905694150420947*up2[7]-0.6123724356957944*up2[6]-0.6123724356957944*up2[5]+0.3535533905932737*up2[4]+0.6123724356957944*up2[3]-0.3535533905932737*up2[2]-0.3535533905932737*up2[1]+0.3535533905932737*up2[0], 0.7905694150420947*up3[11]-0.7905694150420948*up3[10]-0.7905694150420948*up3[9]+0.6123724356957944*up3[8]+0.7905694150420947*up3[7]-0.6123724356957944*up3[6]-0.6123724356957944*up3[5]+0.3535533905932737*up3[4]+0.6123724356957944*up3[3]-0.3535533905932737*up3[2]-0.3535533905932737*up3[1]+0.3535533905932737*up3[0], 0.7905694150420947*up4[11]-0.7905694150420948*up4[10]-0.7905694150420948*up4[9]+0.6123724356957944*up4[8]+0.7905694150420947*up4[7]-0.6123724356957944*up4[6]-0.6123724356957944*up4[5]+0.3535533905932737*up4[4]+0.6123724356957944*up4[3]-0.3535533905932737*up4[2]-0.3535533905932737*up4[1]+0.3535533905932737*up4[0], 0.7905694150420947*up5[11]-0.7905694150420948*up5[10]-0.7905694150420948*up5[9]+0.6123724356957944*up5[8]+0.7905694150420947*up5[7]-0.6123724356957944*up5[6]-0.6123724356957944*up5[5]+0.3535533905932737*up5[4]+0.6123724356957944*up5[3]-0.3535533905932737*up5[2]-0.3535533905932737*up5[1]+0.3535533905932737*up5[0], 0.7905694150420947*up6[11]-0.7905694150420948*up6[10]-0.7905694150420948*up6[9]+0.6123724356957944*up6[8]+0.7905694150420947*up6[7]-0.6123724356957944*up6[6]-0.6123724356957944*up6[5]+0.3535533905932737*up6[4]+0.6123724356957944*up6[3]-0.3535533905932737*up6[2]-0.3535533905932737*up6[1]+0.3535533905932737*up6[0], 0.7905694150420947*up7[11]-0.7905694150420948*up7[10]-0.7905694150420948*up7[9]+0.6123724356957944*up7[8]+0.7905694150420947*up7[7]-0.6123724356957944*up7[6]-0.6123724356957944*up7[5]+0.3535533905932737*up7[4]+0.6123724356957944*up7[3]-0.3535533905932737*up7[2]-0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    // first-harmonic coefficients at the node: coupling flux T and its speed bounds at the edges 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    const double gu_q[6] = {gu00_q, gu01_q, gu02_q, gu10_q, gu11_q, gu12_q}; 
    const double gx_q[2] = {gxi0_q, gxi1_q}; 
    struct pkpm_harm_coeffs hc; 
    pkpm_harm_coeffs_set(&hc, 2, b_q, u_q, xi_q, th_q, d_q, gu_q, gb_q, gx_q, 0, 0.0, 0); 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_coupling_flux(&hc, w_lo, Ulo_q, Tlo); 
    pkpm_harm_coupling_flux(&hc, w_up, Uup_q, Tup); 
    pkpm_harm_speed_bounds(&hc, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_speed_bounds(&hc, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double ap_lo = a0p_q + cp_q*w_lo + k_lo, ap_up = a0p_q + cp_q*w_up + k_up; 
    double wf_lo = 0.0, wu_lo = 0.0, wf_up = 0.0, wu_up = 0.0; 
    { 
      const double sl = ap_lo + llo_lo, sr = ap_lo + lhi_lo; 
      if (sr < 0.0) { 
        wf_lo = 1.0; 
      } else if (sl < 0.0) { 
        wf_lo = -sl/(sr - sl); wu_lo = sl*sr/(sr - sl); 
      } 
    } 
    { 
      const double sl = ap_up + llo_up, sr = ap_up + lhi_up; 
      if (sl > 0.0) { 
        wf_up = 1.0; 
      } else if (sr > 0.0) { 
        wf_up = sr/(sr - sl); wu_up = -sl*sr/(sr - sl); 
      } 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) wf_e U_e, known part kU (k_e wf_e U_e) and 
    // x (coupling and HLL terms wf_e T_e + wu_e U_e) 
    const double eF0 = wf_lo*Flo_q - wf_up*Fup_q; 
    const double eF1 = wf_lo*w_lo*Flo_q - wf_up*w_up*Fup_q; 
    const double eF2 = wf_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*Fup_q; 
    const double eF3 = wf_lo*w_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = wf_lo*Glo_q - wf_up*Gup_q; 
    const double eG1 = wf_lo*w_lo*Glo_q - wf_up*w_up*Gup_q; 
    const double kF0 = wf_lo*k_lo*Flo_q - wf_up*k_up*Fup_q; 
    const double kF1 = wf_lo*w_lo*k_lo*Flo_q - wf_up*w_up*k_up*Fup_q; 
    const double kF2 = wf_lo*w_lo*w_lo*k_lo*Flo_q - wf_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = wf_lo*k_lo*Glo_q - wf_up*k_up*Gup_q; 
    double xlo[8], xup[8]; 
    for (int c = 0; c < 8; ++c) { 
      xlo[c] = wf_lo*Tlo[c] + wu_lo*Ulo_q[c]; 
      xup[c] = wf_up*Tup[c] + wu_up*Uup_q[c]; 
    } 
    const double xF0 = xlo[0] - xup[0]; 
    const double xF1 = w_lo*xlo[0] - w_up*xup[0]; 
    const double xF2 = w_lo*w_lo*xlo[0] - w_up*w_up*xup[0]; 
    const double xG0 = xlo[1] - xup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) 
      v[a] = rM_q[a] + (wf_lo*k_lo*Ulo_q[2 + a] + xlo[2 + a]) - (wf_up*k_up*Uup_q[2 + a] + xup[2 + a]); 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*Ulo_q[2 + a] - wf_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*w_lo*Ulo_q[2 + a] - wf_up*w_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - kF0 - xF0, b1 = -r1_q - kF1 - xF1, b2 = -r2_q - kF2 - 2.0*kG0 - xF2 - 2.0*xG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[0] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[0] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[0] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_n[a][0] = -(P[a] + (c0_n[0] + c_n[0])*Q[a] + Aa[a]*a0_n[0] + Ac[a]*c_n[0])/n0_q; 
      const double al_lo = a0_n[0] + c_n[0]*w_lo + k_lo, al_up = a0_n[0] + c_n[0]*w_up + k_up; 
      for (int c = 0; c < 8; ++c) { 
        flo_n[c][0] = wf_lo*al_lo*Ulo_q[c] + xlo[c]; 
        fup_n[c][0] = wf_up*al_up*Uup_q[c] + xup[c]; 
      } 
      cfl_lo = fmax(cfl_lo, fmax(fabs(al_lo + llo_lo), fabs(al_lo + lhi_lo))); 
      cfl_up = fmax(cfl_up, fmax(fabs(al_up + llo_up), fabs(al_up + lhi_up))); 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = -(0.5*n0[3])+0.5*n0[2]-0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = -(0.5*n1[3])+0.5*n1[2]-0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = -(0.5*n2[3])+0.5*n2[2]-0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = -(0.5*nG[3])+0.5*nG[2]-0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = -(0.5*r0[3])+0.5*r0[2]-0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = -(0.5*r1[3])+0.5*r1[2]-0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = -(0.5*r2[3])+0.5*r2[2]-0.5*r2[1]+0.5*r2[0] + 2.0*(-(0.5*rG[3])+0.5*rG[2]-0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {-(0.5*mu0[3])+0.5*mu0[2]-0.5*mu0[1]+0.5*mu0[0], -(0.5*mu1[3])+0.5*mu1[2]-0.5*mu1[1]+0.5*mu1[0], -(0.5*mu2[3])+0.5*mu2[2]-0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {-(0.5*rM0[3])+0.5*rM0[2]-0.5*rM0[1]+0.5*rM0[0], -(0.5*rM1[3])+0.5*rM1[2]-0.5*rM1[1]+0.5*rM1[0], -(0.5*rM2[3])+0.5*rM2[2]-0.5*rM2[1]+0.5*rM2[0]}; 
    const double a0p_q = -(0.5*alpha_0_prev[3])+0.5*alpha_0_prev[2]-0.5*alpha_0_prev[1]+0.5*alpha_0_prev[0]; 
    const double cp_q = -(0.5*c_h_prev[3])+0.5*c_h_prev[2]-0.5*c_h_prev[1]+0.5*c_h_prev[0]; 
    const double xi_q = -(0.5*xi_1[3])+0.5*xi_1[2]-0.5*xi_1[1]+0.5*xi_1[0]; 
    const double th_q = -(0.5*theta[3])+0.5*theta[2]-0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = -(0.5*b0[3])+0.5*b0[2]-0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = -(0.5*b1[3])+0.5*b1[2]-0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = -(0.5*b2[3])+0.5*b2[2]-0.5*b2[1]+0.5*b2[0]; 
    const double u_q[3] = {-(0.5*u0[3])+0.5*u0[2]-0.5*u0[1]+0.5*u0[0], -(0.5*u1[3])+0.5*u1[2]-0.5*u1[1]+0.5*u1[0], -(0.5*u2[3])+0.5*u2[2]-0.5*u2[1]+0.5*u2[0]}; 
    const double d_q[3] = {-(0.5*dtb0[3])+0.5*dtb0[2]-0.5*dtb0[1]+0.5*dtb0[0], -(0.5*dtb1[3])+0.5*dtb1[2]-0.5*dtb1[1]+0.5*dtb1[0], -(0.5*dtb2[3])+0.5*dtb2[2]-0.5*dtb2[1]+0.5*dtb2[0]}; 
    const double gu00_q = -(0.5*gux0[3])+0.5*gux0[2]-0.5*gux0[1]+0.5*gux0[0]; 
    const double gu01_q = -(0.5*gux1[3])+0.5*gux1[2]-0.5*gux1[1]+0.5*gux1[0]; 
    const double gu02_q = -(0.5*gux2[3])+0.5*gux2[2]-0.5*gux2[1]+0.5*gux2[0]; 
    const double gu10_q = -(0.5*guy0[3])+0.5*guy0[2]-0.5*guy0[1]+0.5*guy0[0]; 
    const double gu11_q = -(0.5*guy1[3])+0.5*guy1[2]-0.5*guy1[1]+0.5*guy1[0]; 
    const double gu12_q = -(0.5*guy2[3])+0.5*guy2[2]-0.5*guy2[1]+0.5*guy2[0]; 
    const double gxi0_q = -(0.5*gxix[3])+0.5*gxix[2]-0.5*gxix[1]+0.5*gxix[0]; 
    const double gxi1_q = -(0.5*gxiy[3])+0.5*gxiy[2]-0.5*gxiy[1]+0.5*gxiy[0]; 
    const double gb_q[6] = {-(0.5*gbx0[3])+0.5*gbx0[2]-0.5*gbx0[1]+0.5*gbx0[0], -(0.5*gbx1[3])+0.5*gbx1[2]-0.5*gbx1[1]+0.5*gbx1[0], -(0.5*gbx2[3])+0.5*gbx2[2]-0.5*gbx2[1]+0.5*gbx2[0], -(0.5*gby0[3])+0.5*gby0[2]-0.5*gby0[1]+0.5*gby0[0], -(0.5*gby1[3])+0.5*gby1[2]-0.5*gby1[1]+0.5*gby1[0], -(0.5*gby2[3])+0.5*gby2[2]-0.5*gby2[1]+0.5*gby2[0]}; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q + b1_q*b0_q*gu10_q + b1_q*b1_q*gu11_q + b1_q*b2_q*gu12_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q + b1_q*gxi1_q; 
    const double Ulo_q[8] = {-(0.7905694150420947*lo0[11])+0.7905694150420948*lo0[10]-0.7905694150420948*lo0[9]+0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]-0.6123724356957944*lo0[6]+0.6123724356957944*lo0[5]-0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]+0.3535533905932737*lo0[2]-0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], -(0.7905694150420947*lo1[11])+0.7905694150420948*lo1[10]-0.7905694150420948*lo1[9]+0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]-0.6123724356957944*lo1[6]+0.6123724356957944*lo1[5]-0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]+0.3535533905932737*lo1[2]-0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], -(0.7905694150420947*lo2[11])+0.7905694150420948*lo2[10]-0.7905694150420948*lo2[9]+0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]-0.6123724356957944*lo2[6]+0.6123724356957944*lo2[5]-0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]+0.3535533905932737*lo2[2]-0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], -(0.7905694150420947*lo3[11])+0.7905694150420948*lo3[10]-0.7905694150420948*lo3[9]+0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]-0.6123724356957944*lo3[6]+0.6123724356957944*lo3[5]-0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]+0.3535533905932737*lo3[2]-0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], -(0.7905694150420947*lo4[11])+0.7905694150420948*lo4[10]-0.7905694150420948*lo4[9]+0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]-0.6123724356957944*lo4[6]+0.6123724356957944*lo4[5]-0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]+0.3535533905932737*lo4[2]-0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], -(0.7905694150420947*lo5[11])+0.7905694150420948*lo5[10]-0.7905694150420948*lo5[9]+0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]-0.6123724356957944*lo5[6]+0.6123724356957944*lo5[5]-0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]+0.3535533905932737*lo5[2]-0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], -(0.7905694150420947*lo6[11])+0.7905694150420948*lo6[10]-0.7905694150420948*lo6[9]+0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]-0.6123724356957944*lo6[6]+0.6123724356957944*lo6[5]-0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]+0.3535533905932737*lo6[2]-0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], -(0.7905694150420947*lo7[11])+0.7905694150420948*lo7[10]-0.7905694150420948*lo7[9]+0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]-0.6123724356957944*lo7[6]+0.6123724356957944*lo7[5]-0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]+0.3535533905932737*lo7[2]-0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {-(0.7905694150420947*up0[11])+0.7905694150420948*up0[10]-0.7905694150420948*up0[9]-0.6123724356957944*up0[8]+0.7905694150420947*up0[7]+0.6123724356957944*up0[6]-0.6123724356957944*up0[5]-0.3535533905932737*up0[4]+0.6123724356957944*up0[3]+0.3535533905932737*up0[2]-0.3535533905932737*up0[1]+0.3535533905932737*up0[0], -(0.7905694150420947*up1[11])+0.7905694150420948*up1[10]-0.7905694150420948*up1[9]-0.6123724356957944*up1[8]+0.7905694150420947*up1[7]+0.6123724356957944*up1[6]-0.6123724356957944*up1[5]-0.3535533905932737*up1[4]+0.6123724356957944*up1[3]+0.3535533905932737*up1[2]-0.3535533905932737*up1[1]+0.3535533905932737*up1[0], -(0.7905694150420947*up2[11])+0.7905694150420948*up2[10]-0.7905694150420948*up2[9]-0.6123724356957944*up2[8]+0.7905694150420947*up2[7]+0.6123724356957944*up2[6]-0.6123724356957944*up2[5]-0.3535533905932737*up2[4]+0.6123724356957944*up2[3]+0.3535533905932737*up2[2]-0.3535533905932737*up2[1]+0.3535533905932737*up2[0], -(0.7905694150420947*up3[11])+0.7905694150420948*up3[10]-0.7905694150420948*up3[9]-0.6123724356957944*up3[8]+0.7905694150420947*up3[7]+0.6123724356957944*up3[6]-0.6123724356957944*up3[5]-0.3535533905932737*up3[4]+0.6123724356957944*up3[3]+0.3535533905932737*up3[2]-0.3535533905932737*up3[1]+0.3535533905932737*up3[0], -(0.7905694150420947*up4[11])+0.7905694150420948*up4[10]-0.7905694150420948*up4[9]-0.6123724356957944*up4[8]+0.7905694150420947*up4[7]+0.6123724356957944*up4[6]-0.6123724356957944*up4[5]-0.3535533905932737*up4[4]+0.6123724356957944*up4[3]+0.3535533905932737*up4[2]-0.3535533905932737*up4[1]+0.3535533905932737*up4[0], -(0.7905694150420947*up5[11])+0.7905694150420948*up5[10]-0.7905694150420948*up5[9]-0.6123724356957944*up5[8]+0.7905694150420947*up5[7]+0.6123724356957944*up5[6]-0.6123724356957944*up5[5]-0.3535533905932737*up5[4]+0.6123724356957944*up5[3]+0.3535533905932737*up5[2]-0.3535533905932737*up5[1]+0.3535533905932737*up5[0], -(0.7905694150420947*up6[11])+0.7905694150420948*up6[10]-0.7905694150420948*up6[9]-0.6123724356957944*up6[8]+0.7905694150420947*up6[7]+0.6123724356957944*up6[6]-0.6123724356957944*up6[5]-0.3535533905932737*up6[4]+0.6123724356957944*up6[3]+0.3535533905932737*up6[2]-0.3535533905932737*up6[1]+0.3535533905932737*up6[0], -(0.7905694150420947*up7[11])+0.7905694150420948*up7[10]-0.7905694150420948*up7[9]-0.6123724356957944*up7[8]+0.7905694150420947*up7[7]+0.6123724356957944*up7[6]-0.6123724356957944*up7[5]-0.3535533905932737*up7[4]+0.6123724356957944*up7[3]+0.3535533905932737*up7[2]-0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    // first-harmonic coefficients at the node: coupling flux T and its speed bounds at the edges 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    const double gu_q[6] = {gu00_q, gu01_q, gu02_q, gu10_q, gu11_q, gu12_q}; 
    const double gx_q[2] = {gxi0_q, gxi1_q}; 
    struct pkpm_harm_coeffs hc; 
    pkpm_harm_coeffs_set(&hc, 2, b_q, u_q, xi_q, th_q, d_q, gu_q, gb_q, gx_q, 0, 0.0, 0); 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_coupling_flux(&hc, w_lo, Ulo_q, Tlo); 
    pkpm_harm_coupling_flux(&hc, w_up, Uup_q, Tup); 
    pkpm_harm_speed_bounds(&hc, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_speed_bounds(&hc, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double ap_lo = a0p_q + cp_q*w_lo + k_lo, ap_up = a0p_q + cp_q*w_up + k_up; 
    double wf_lo = 0.0, wu_lo = 0.0, wf_up = 0.0, wu_up = 0.0; 
    { 
      const double sl = ap_lo + llo_lo, sr = ap_lo + lhi_lo; 
      if (sr < 0.0) { 
        wf_lo = 1.0; 
      } else if (sl < 0.0) { 
        wf_lo = -sl/(sr - sl); wu_lo = sl*sr/(sr - sl); 
      } 
    } 
    { 
      const double sl = ap_up + llo_up, sr = ap_up + lhi_up; 
      if (sl > 0.0) { 
        wf_up = 1.0; 
      } else if (sr > 0.0) { 
        wf_up = sr/(sr - sl); wu_up = -sl*sr/(sr - sl); 
      } 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) wf_e U_e, known part kU (k_e wf_e U_e) and 
    // x (coupling and HLL terms wf_e T_e + wu_e U_e) 
    const double eF0 = wf_lo*Flo_q - wf_up*Fup_q; 
    const double eF1 = wf_lo*w_lo*Flo_q - wf_up*w_up*Fup_q; 
    const double eF2 = wf_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*Fup_q; 
    const double eF3 = wf_lo*w_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = wf_lo*Glo_q - wf_up*Gup_q; 
    const double eG1 = wf_lo*w_lo*Glo_q - wf_up*w_up*Gup_q; 
    const double kF0 = wf_lo*k_lo*Flo_q - wf_up*k_up*Fup_q; 
    const double kF1 = wf_lo*w_lo*k_lo*Flo_q - wf_up*w_up*k_up*Fup_q; 
    const double kF2 = wf_lo*w_lo*w_lo*k_lo*Flo_q - wf_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = wf_lo*k_lo*Glo_q - wf_up*k_up*Gup_q; 
    double xlo[8], xup[8]; 
    for (int c = 0; c < 8; ++c) { 
      xlo[c] = wf_lo*Tlo[c] + wu_lo*Ulo_q[c]; 
      xup[c] = wf_up*Tup[c] + wu_up*Uup_q[c]; 
    } 
    const double xF0 = xlo[0] - xup[0]; 
    const double xF1 = w_lo*xlo[0] - w_up*xup[0]; 
    const double xF2 = w_lo*w_lo*xlo[0] - w_up*w_up*xup[0]; 
    const double xG0 = xlo[1] - xup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) 
      v[a] = rM_q[a] + (wf_lo*k_lo*Ulo_q[2 + a] + xlo[2 + a]) - (wf_up*k_up*Uup_q[2 + a] + xup[2 + a]); 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*Ulo_q[2 + a] - wf_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*w_lo*Ulo_q[2 + a] - wf_up*w_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - kF0 - xF0, b1 = -r1_q - kF1 - xF1, b2 = -r2_q - kF2 - 2.0*kG0 - xF2 - 2.0*xG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[1] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[1] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[1] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_n[a][1] = -(P[a] + (c0_n[1] + c_n[1])*Q[a] + Aa[a]*a0_n[1] + Ac[a]*c_n[1])/n0_q; 
      const double al_lo = a0_n[1] + c_n[1]*w_lo + k_lo, al_up = a0_n[1] + c_n[1]*w_up + k_up; 
      for (int c = 0; c < 8; ++c) { 
        flo_n[c][1] = wf_lo*al_lo*Ulo_q[c] + xlo[c]; 
        fup_n[c][1] = wf_up*al_up*Uup_q[c] + xup[c]; 
      } 
      cfl_lo = fmax(cfl_lo, fmax(fabs(al_lo + llo_lo), fabs(al_lo + lhi_lo))); 
      cfl_up = fmax(cfl_up, fmax(fabs(al_up + llo_up), fabs(al_up + lhi_up))); 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = -(0.5*n0[3])-0.5*n0[2]+0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = -(0.5*n1[3])-0.5*n1[2]+0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = -(0.5*n2[3])-0.5*n2[2]+0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = -(0.5*nG[3])-0.5*nG[2]+0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = -(0.5*r0[3])-0.5*r0[2]+0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = -(0.5*r1[3])-0.5*r1[2]+0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = -(0.5*r2[3])-0.5*r2[2]+0.5*r2[1]+0.5*r2[0] + 2.0*(-(0.5*rG[3])-0.5*rG[2]+0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {-(0.5*mu0[3])-0.5*mu0[2]+0.5*mu0[1]+0.5*mu0[0], -(0.5*mu1[3])-0.5*mu1[2]+0.5*mu1[1]+0.5*mu1[0], -(0.5*mu2[3])-0.5*mu2[2]+0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {-(0.5*rM0[3])-0.5*rM0[2]+0.5*rM0[1]+0.5*rM0[0], -(0.5*rM1[3])-0.5*rM1[2]+0.5*rM1[1]+0.5*rM1[0], -(0.5*rM2[3])-0.5*rM2[2]+0.5*rM2[1]+0.5*rM2[0]}; 
    const double a0p_q = -(0.5*alpha_0_prev[3])-0.5*alpha_0_prev[2]+0.5*alpha_0_prev[1]+0.5*alpha_0_prev[0]; 
    const double cp_q = -(0.5*c_h_prev[3])-0.5*c_h_prev[2]+0.5*c_h_prev[1]+0.5*c_h_prev[0]; 
    const double xi_q = -(0.5*xi_1[3])-0.5*xi_1[2]+0.5*xi_1[1]+0.5*xi_1[0]; 
    const double th_q = -(0.5*theta[3])-0.5*theta[2]+0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = -(0.5*b0[3])-0.5*b0[2]+0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = -(0.5*b1[3])-0.5*b1[2]+0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = -(0.5*b2[3])-0.5*b2[2]+0.5*b2[1]+0.5*b2[0]; 
    const double u_q[3] = {-(0.5*u0[3])-0.5*u0[2]+0.5*u0[1]+0.5*u0[0], -(0.5*u1[3])-0.5*u1[2]+0.5*u1[1]+0.5*u1[0], -(0.5*u2[3])-0.5*u2[2]+0.5*u2[1]+0.5*u2[0]}; 
    const double d_q[3] = {-(0.5*dtb0[3])-0.5*dtb0[2]+0.5*dtb0[1]+0.5*dtb0[0], -(0.5*dtb1[3])-0.5*dtb1[2]+0.5*dtb1[1]+0.5*dtb1[0], -(0.5*dtb2[3])-0.5*dtb2[2]+0.5*dtb2[1]+0.5*dtb2[0]}; 
    const double gu00_q = -(0.5*gux0[3])-0.5*gux0[2]+0.5*gux0[1]+0.5*gux0[0]; 
    const double gu01_q = -(0.5*gux1[3])-0.5*gux1[2]+0.5*gux1[1]+0.5*gux1[0]; 
    const double gu02_q = -(0.5*gux2[3])-0.5*gux2[2]+0.5*gux2[1]+0.5*gux2[0]; 
    const double gu10_q = -(0.5*guy0[3])-0.5*guy0[2]+0.5*guy0[1]+0.5*guy0[0]; 
    const double gu11_q = -(0.5*guy1[3])-0.5*guy1[2]+0.5*guy1[1]+0.5*guy1[0]; 
    const double gu12_q = -(0.5*guy2[3])-0.5*guy2[2]+0.5*guy2[1]+0.5*guy2[0]; 
    const double gxi0_q = -(0.5*gxix[3])-0.5*gxix[2]+0.5*gxix[1]+0.5*gxix[0]; 
    const double gxi1_q = -(0.5*gxiy[3])-0.5*gxiy[2]+0.5*gxiy[1]+0.5*gxiy[0]; 
    const double gb_q[6] = {-(0.5*gbx0[3])-0.5*gbx0[2]+0.5*gbx0[1]+0.5*gbx0[0], -(0.5*gbx1[3])-0.5*gbx1[2]+0.5*gbx1[1]+0.5*gbx1[0], -(0.5*gbx2[3])-0.5*gbx2[2]+0.5*gbx2[1]+0.5*gbx2[0], -(0.5*gby0[3])-0.5*gby0[2]+0.5*gby0[1]+0.5*gby0[0], -(0.5*gby1[3])-0.5*gby1[2]+0.5*gby1[1]+0.5*gby1[0], -(0.5*gby2[3])-0.5*gby2[2]+0.5*gby2[1]+0.5*gby2[0]}; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q + b1_q*b0_q*gu10_q + b1_q*b1_q*gu11_q + b1_q*b2_q*gu12_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q + b1_q*gxi1_q; 
    const double Ulo_q[8] = {-(0.7905694150420947*lo0[11])-0.7905694150420948*lo0[10]+0.7905694150420948*lo0[9]+0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]+0.6123724356957944*lo0[6]-0.6123724356957944*lo0[5]-0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]-0.3535533905932737*lo0[2]+0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], -(0.7905694150420947*lo1[11])-0.7905694150420948*lo1[10]+0.7905694150420948*lo1[9]+0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]+0.6123724356957944*lo1[6]-0.6123724356957944*lo1[5]-0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]-0.3535533905932737*lo1[2]+0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], -(0.7905694150420947*lo2[11])-0.7905694150420948*lo2[10]+0.7905694150420948*lo2[9]+0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]+0.6123724356957944*lo2[6]-0.6123724356957944*lo2[5]-0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]-0.3535533905932737*lo2[2]+0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], -(0.7905694150420947*lo3[11])-0.7905694150420948*lo3[10]+0.7905694150420948*lo3[9]+0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]+0.6123724356957944*lo3[6]-0.6123724356957944*lo3[5]-0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]-0.3535533905932737*lo3[2]+0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], -(0.7905694150420947*lo4[11])-0.7905694150420948*lo4[10]+0.7905694150420948*lo4[9]+0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]+0.6123724356957944*lo4[6]-0.6123724356957944*lo4[5]-0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]-0.3535533905932737*lo4[2]+0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], -(0.7905694150420947*lo5[11])-0.7905694150420948*lo5[10]+0.7905694150420948*lo5[9]+0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]+0.6123724356957944*lo5[6]-0.6123724356957944*lo5[5]-0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]-0.3535533905932737*lo5[2]+0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], -(0.7905694150420947*lo6[11])-0.7905694150420948*lo6[10]+0.7905694150420948*lo6[9]+0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]+0.6123724356957944*lo6[6]-0.6123724356957944*lo6[5]-0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]-0.3535533905932737*lo6[2]+0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], -(0.7905694150420947*lo7[11])-0.7905694150420948*lo7[10]+0.7905694150420948*lo7[9]+0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]+0.6123724356957944*lo7[6]-0.6123724356957944*lo7[5]-0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]-0.3535533905932737*lo7[2]+0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {-(0.7905694150420947*up0[11])-0.7905694150420948*up0[10]+0.7905694150420948*up0[9]-0.6123724356957944*up0[8]+0.7905694150420947*up0[7]-0.6123724356957944*up0[6]+0.6123724356957944*up0[5]-0.3535533905932737*up0[4]+0.6123724356957944*up0[3]-0.3535533905932737*up0[2]+0.3535533905932737*up0[1]+0.3535533905932737*up0[0], -(0.7905694150420947*up1[11])-0.7905694150420948*up1[10]+0.7905694150420948*up1[9]-0.6123724356957944*up1[8]+0.7905694150420947*up1[7]-0.6123724356957944*up1[6]+0.6123724356957944*up1[5]-0.3535533905932737*up1[4]+0.6123724356957944*up1[3]-0.3535533905932737*up1[2]+0.3535533905932737*up1[1]+0.3535533905932737*up1[0], -(0.7905694150420947*up2[11])-0.7905694150420948*up2[10]+0.7905694150420948*up2[9]-0.6123724356957944*up2[8]+0.7905694150420947*up2[7]-0.6123724356957944*up2[6]+0.6123724356957944*up2[5]-0.3535533905932737*up2[4]+0.6123724356957944*up2[3]-0.3535533905932737*up2[2]+0.3535533905932737*up2[1]+0.3535533905932737*up2[0], -(0.7905694150420947*up3[11])-0.7905694150420948*up3[10]+0.7905694150420948*up3[9]-0.6123724356957944*up3[8]+0.7905694150420947*up3[7]-0.6123724356957944*up3[6]+0.6123724356957944*up3[5]-0.3535533905932737*up3[4]+0.6123724356957944*up3[3]-0.3535533905932737*up3[2]+0.3535533905932737*up3[1]+0.3535533905932737*up3[0], -(0.7905694150420947*up4[11])-0.7905694150420948*up4[10]+0.7905694150420948*up4[9]-0.6123724356957944*up4[8]+0.7905694150420947*up4[7]-0.6123724356957944*up4[6]+0.6123724356957944*up4[5]-0.3535533905932737*up4[4]+0.6123724356957944*up4[3]-0.3535533905932737*up4[2]+0.3535533905932737*up4[1]+0.3535533905932737*up4[0], -(0.7905694150420947*up5[11])-0.7905694150420948*up5[10]+0.7905694150420948*up5[9]-0.6123724356957944*up5[8]+0.7905694150420947*up5[7]-0.6123724356957944*up5[6]+0.6123724356957944*up5[5]-0.3535533905932737*up5[4]+0.6123724356957944*up5[3]-0.3535533905932737*up5[2]+0.3535533905932737*up5[1]+0.3535533905932737*up5[0], -(0.7905694150420947*up6[11])-0.7905694150420948*up6[10]+0.7905694150420948*up6[9]-0.6123724356957944*up6[8]+0.7905694150420947*up6[7]-0.6123724356957944*up6[6]+0.6123724356957944*up6[5]-0.3535533905932737*up6[4]+0.6123724356957944*up6[3]-0.3535533905932737*up6[2]+0.3535533905932737*up6[1]+0.3535533905932737*up6[0], -(0.7905694150420947*up7[11])-0.7905694150420948*up7[10]+0.7905694150420948*up7[9]-0.6123724356957944*up7[8]+0.7905694150420947*up7[7]-0.6123724356957944*up7[6]+0.6123724356957944*up7[5]-0.3535533905932737*up7[4]+0.6123724356957944*up7[3]-0.3535533905932737*up7[2]+0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    // first-harmonic coefficients at the node: coupling flux T and its speed bounds at the edges 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    const double gu_q[6] = {gu00_q, gu01_q, gu02_q, gu10_q, gu11_q, gu12_q}; 
    const double gx_q[2] = {gxi0_q, gxi1_q}; 
    struct pkpm_harm_coeffs hc; 
    pkpm_harm_coeffs_set(&hc, 2, b_q, u_q, xi_q, th_q, d_q, gu_q, gb_q, gx_q, 0, 0.0, 0); 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_coupling_flux(&hc, w_lo, Ulo_q, Tlo); 
    pkpm_harm_coupling_flux(&hc, w_up, Uup_q, Tup); 
    pkpm_harm_speed_bounds(&hc, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_speed_bounds(&hc, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double ap_lo = a0p_q + cp_q*w_lo + k_lo, ap_up = a0p_q + cp_q*w_up + k_up; 
    double wf_lo = 0.0, wu_lo = 0.0, wf_up = 0.0, wu_up = 0.0; 
    { 
      const double sl = ap_lo + llo_lo, sr = ap_lo + lhi_lo; 
      if (sr < 0.0) { 
        wf_lo = 1.0; 
      } else if (sl < 0.0) { 
        wf_lo = -sl/(sr - sl); wu_lo = sl*sr/(sr - sl); 
      } 
    } 
    { 
      const double sl = ap_up + llo_up, sr = ap_up + lhi_up; 
      if (sl > 0.0) { 
        wf_up = 1.0; 
      } else if (sr > 0.0) { 
        wf_up = sr/(sr - sl); wu_up = -sl*sr/(sr - sl); 
      } 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) wf_e U_e, known part kU (k_e wf_e U_e) and 
    // x (coupling and HLL terms wf_e T_e + wu_e U_e) 
    const double eF0 = wf_lo*Flo_q - wf_up*Fup_q; 
    const double eF1 = wf_lo*w_lo*Flo_q - wf_up*w_up*Fup_q; 
    const double eF2 = wf_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*Fup_q; 
    const double eF3 = wf_lo*w_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = wf_lo*Glo_q - wf_up*Gup_q; 
    const double eG1 = wf_lo*w_lo*Glo_q - wf_up*w_up*Gup_q; 
    const double kF0 = wf_lo*k_lo*Flo_q - wf_up*k_up*Fup_q; 
    const double kF1 = wf_lo*w_lo*k_lo*Flo_q - wf_up*w_up*k_up*Fup_q; 
    const double kF2 = wf_lo*w_lo*w_lo*k_lo*Flo_q - wf_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = wf_lo*k_lo*Glo_q - wf_up*k_up*Gup_q; 
    double xlo[8], xup[8]; 
    for (int c = 0; c < 8; ++c) { 
      xlo[c] = wf_lo*Tlo[c] + wu_lo*Ulo_q[c]; 
      xup[c] = wf_up*Tup[c] + wu_up*Uup_q[c]; 
    } 
    const double xF0 = xlo[0] - xup[0]; 
    const double xF1 = w_lo*xlo[0] - w_up*xup[0]; 
    const double xF2 = w_lo*w_lo*xlo[0] - w_up*w_up*xup[0]; 
    const double xG0 = xlo[1] - xup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) 
      v[a] = rM_q[a] + (wf_lo*k_lo*Ulo_q[2 + a] + xlo[2 + a]) - (wf_up*k_up*Uup_q[2 + a] + xup[2 + a]); 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*Ulo_q[2 + a] - wf_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*w_lo*Ulo_q[2 + a] - wf_up*w_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - kF0 - xF0, b1 = -r1_q - kF1 - xF1, b2 = -r2_q - kF2 - 2.0*kG0 - xF2 - 2.0*xG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[2] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[2] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[2] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_n[a][2] = -(P[a] + (c0_n[2] + c_n[2])*Q[a] + Aa[a]*a0_n[2] + Ac[a]*c_n[2])/n0_q; 
      const double al_lo = a0_n[2] + c_n[2]*w_lo + k_lo, al_up = a0_n[2] + c_n[2]*w_up + k_up; 
      for (int c = 0; c < 8; ++c) { 
        flo_n[c][2] = wf_lo*al_lo*Ulo_q[c] + xlo[c]; 
        fup_n[c][2] = wf_up*al_up*Uup_q[c] + xup[c]; 
      } 
      cfl_lo = fmax(cfl_lo, fmax(fabs(al_lo + llo_lo), fabs(al_lo + lhi_lo))); 
      cfl_up = fmax(cfl_up, fmax(fabs(al_up + llo_up), fabs(al_up + lhi_up))); 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = 0.5*n0[3]+0.5*n0[2]+0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = 0.5*n1[3]+0.5*n1[2]+0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = 0.5*n2[3]+0.5*n2[2]+0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = 0.5*nG[3]+0.5*nG[2]+0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = 0.5*r0[3]+0.5*r0[2]+0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = 0.5*r1[3]+0.5*r1[2]+0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = 0.5*r2[3]+0.5*r2[2]+0.5*r2[1]+0.5*r2[0] + 2.0*(0.5*rG[3]+0.5*rG[2]+0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {0.5*mu0[3]+0.5*mu0[2]+0.5*mu0[1]+0.5*mu0[0], 0.5*mu1[3]+0.5*mu1[2]+0.5*mu1[1]+0.5*mu1[0], 0.5*mu2[3]+0.5*mu2[2]+0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {0.5*rM0[3]+0.5*rM0[2]+0.5*rM0[1]+0.5*rM0[0], 0.5*rM1[3]+0.5*rM1[2]+0.5*rM1[1]+0.5*rM1[0], 0.5*rM2[3]+0.5*rM2[2]+0.5*rM2[1]+0.5*rM2[0]}; 
    const double a0p_q = 0.5*alpha_0_prev[3]+0.5*alpha_0_prev[2]+0.5*alpha_0_prev[1]+0.5*alpha_0_prev[0]; 
    const double cp_q = 0.5*c_h_prev[3]+0.5*c_h_prev[2]+0.5*c_h_prev[1]+0.5*c_h_prev[0]; 
    const double xi_q = 0.5*xi_1[3]+0.5*xi_1[2]+0.5*xi_1[1]+0.5*xi_1[0]; 
    const double th_q = 0.5*theta[3]+0.5*theta[2]+0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = 0.5*b0[3]+0.5*b0[2]+0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = 0.5*b1[3]+0.5*b1[2]+0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = 0.5*b2[3]+0.5*b2[2]+0.5*b2[1]+0.5*b2[0]; 
    const double u_q[3] = {0.5*u0[3]+0.5*u0[2]+0.5*u0[1]+0.5*u0[0], 0.5*u1[3]+0.5*u1[2]+0.5*u1[1]+0.5*u1[0], 0.5*u2[3]+0.5*u2[2]+0.5*u2[1]+0.5*u2[0]}; 
    const double d_q[3] = {0.5*dtb0[3]+0.5*dtb0[2]+0.5*dtb0[1]+0.5*dtb0[0], 0.5*dtb1[3]+0.5*dtb1[2]+0.5*dtb1[1]+0.5*dtb1[0], 0.5*dtb2[3]+0.5*dtb2[2]+0.5*dtb2[1]+0.5*dtb2[0]}; 
    const double gu00_q = 0.5*gux0[3]+0.5*gux0[2]+0.5*gux0[1]+0.5*gux0[0]; 
    const double gu01_q = 0.5*gux1[3]+0.5*gux1[2]+0.5*gux1[1]+0.5*gux1[0]; 
    const double gu02_q = 0.5*gux2[3]+0.5*gux2[2]+0.5*gux2[1]+0.5*gux2[0]; 
    const double gu10_q = 0.5*guy0[3]+0.5*guy0[2]+0.5*guy0[1]+0.5*guy0[0]; 
    const double gu11_q = 0.5*guy1[3]+0.5*guy1[2]+0.5*guy1[1]+0.5*guy1[0]; 
    const double gu12_q = 0.5*guy2[3]+0.5*guy2[2]+0.5*guy2[1]+0.5*guy2[0]; 
    const double gxi0_q = 0.5*gxix[3]+0.5*gxix[2]+0.5*gxix[1]+0.5*gxix[0]; 
    const double gxi1_q = 0.5*gxiy[3]+0.5*gxiy[2]+0.5*gxiy[1]+0.5*gxiy[0]; 
    const double gb_q[6] = {0.5*gbx0[3]+0.5*gbx0[2]+0.5*gbx0[1]+0.5*gbx0[0], 0.5*gbx1[3]+0.5*gbx1[2]+0.5*gbx1[1]+0.5*gbx1[0], 0.5*gbx2[3]+0.5*gbx2[2]+0.5*gbx2[1]+0.5*gbx2[0], 0.5*gby0[3]+0.5*gby0[2]+0.5*gby0[1]+0.5*gby0[0], 0.5*gby1[3]+0.5*gby1[2]+0.5*gby1[1]+0.5*gby1[0], 0.5*gby2[3]+0.5*gby2[2]+0.5*gby2[1]+0.5*gby2[0]}; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q + b1_q*b0_q*gu10_q + b1_q*b1_q*gu11_q + b1_q*b2_q*gu12_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q + b1_q*gxi1_q; 
    const double Ulo_q[8] = {0.7905694150420947*lo0[11]+0.7905694150420948*lo0[10]+0.7905694150420948*lo0[9]-0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]-0.6123724356957944*lo0[6]-0.6123724356957944*lo0[5]+0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]+0.3535533905932737*lo0[2]+0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], 0.7905694150420947*lo1[11]+0.7905694150420948*lo1[10]+0.7905694150420948*lo1[9]-0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]-0.6123724356957944*lo1[6]-0.6123724356957944*lo1[5]+0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]+0.3535533905932737*lo1[2]+0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], 0.7905694150420947*lo2[11]+0.7905694150420948*lo2[10]+0.7905694150420948*lo2[9]-0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]-0.6123724356957944*lo2[6]-0.6123724356957944*lo2[5]+0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]+0.3535533905932737*lo2[2]+0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], 0.7905694150420947*lo3[11]+0.7905694150420948*lo3[10]+0.7905694150420948*lo3[9]-0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]-0.6123724356957944*lo3[6]-0.6123724356957944*lo3[5]+0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]+0.3535533905932737*lo3[2]+0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], 0.7905694150420947*lo4[11]+0.7905694150420948*lo4[10]+0.7905694150420948*lo4[9]-0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]-0.6123724356957944*lo4[6]-0.6123724356957944*lo4[5]+0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]+0.3535533905932737*lo4[2]+0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], 0.7905694150420947*lo5[11]+0.7905694150420948*lo5[10]+0.7905694150420948*lo5[9]-0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]-0.6123724356957944*lo5[6]-0.6123724356957944*lo5[5]+0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]+0.3535533905932737*lo5[2]+0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], 0.7905694150420947*lo6[11]+0.7905694150420948*lo6[10]+0.7905694150420948*lo6[9]-0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]-0.6123724356957944*lo6[6]-0.6123724356957944*lo6[5]+0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]+0.3535533905932737*lo6[2]+0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], 0.7905694150420947*lo7[11]+0.7905694150420948*lo7[10]+0.7905694150420948*lo7[9]-0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]-0.6123724356957944*lo7[6]-0.6123724356957944*lo7[5]+0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]+0.3535533905932737*lo7[2]+0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {0.7905694150420947*up0[11]+0.7905694150420948*up0[10]+0.7905694150420948*up0[9]+0.6123724356957944*up0[8]+0.7905694150420947*up0[7]+0.6123724356957944*up0[6]+0.6123724356957944*up0[5]+0.3535533905932737*up0[4]+0.6123724356957944*up0[3]+0.3535533905932737*up0[2]+0.3535533905932737*up0[1]+0.3535533905932737*up0[0], 0.7905694150420947*up1[11]+0.7905694150420948*up1[10]+0.7905694150420948*up1[9]+0.6123724356957944*up1[8]+0.7905694150420947*up1[7]+0.6123724356957944*up1[6]+0.6123724356957944*up1[5]+0.3535533905932737*up1[4]+0.6123724356957944*up1[3]+0.3535533905932737*up1[2]+0.3535533905932737*up1[1]+0.3535533905932737*up1[0], 0.7905694150420947*up2[11]+0.7905694150420948*up2[10]+0.7905694150420948*up2[9]+0.6123724356957944*up2[8]+0.7905694150420947*up2[7]+0.6123724356957944*up2[6]+0.6123724356957944*up2[5]+0.3535533905932737*up2[4]+0.6123724356957944*up2[3]+0.3535533905932737*up2[2]+0.3535533905932737*up2[1]+0.3535533905932737*up2[0], 0.7905694150420947*up3[11]+0.7905694150420948*up3[10]+0.7905694150420948*up3[9]+0.6123724356957944*up3[8]+0.7905694150420947*up3[7]+0.6123724356957944*up3[6]+0.6123724356957944*up3[5]+0.3535533905932737*up3[4]+0.6123724356957944*up3[3]+0.3535533905932737*up3[2]+0.3535533905932737*up3[1]+0.3535533905932737*up3[0], 0.7905694150420947*up4[11]+0.7905694150420948*up4[10]+0.7905694150420948*up4[9]+0.6123724356957944*up4[8]+0.7905694150420947*up4[7]+0.6123724356957944*up4[6]+0.6123724356957944*up4[5]+0.3535533905932737*up4[4]+0.6123724356957944*up4[3]+0.3535533905932737*up4[2]+0.3535533905932737*up4[1]+0.3535533905932737*up4[0], 0.7905694150420947*up5[11]+0.7905694150420948*up5[10]+0.7905694150420948*up5[9]+0.6123724356957944*up5[8]+0.7905694150420947*up5[7]+0.6123724356957944*up5[6]+0.6123724356957944*up5[5]+0.3535533905932737*up5[4]+0.6123724356957944*up5[3]+0.3535533905932737*up5[2]+0.3535533905932737*up5[1]+0.3535533905932737*up5[0], 0.7905694150420947*up6[11]+0.7905694150420948*up6[10]+0.7905694150420948*up6[9]+0.6123724356957944*up6[8]+0.7905694150420947*up6[7]+0.6123724356957944*up6[6]+0.6123724356957944*up6[5]+0.3535533905932737*up6[4]+0.6123724356957944*up6[3]+0.3535533905932737*up6[2]+0.3535533905932737*up6[1]+0.3535533905932737*up6[0], 0.7905694150420947*up7[11]+0.7905694150420948*up7[10]+0.7905694150420948*up7[9]+0.6123724356957944*up7[8]+0.7905694150420947*up7[7]+0.6123724356957944*up7[6]+0.6123724356957944*up7[5]+0.3535533905932737*up7[4]+0.6123724356957944*up7[3]+0.3535533905932737*up7[2]+0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    // first-harmonic coefficients at the node: coupling flux T and its speed bounds at the edges 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    const double gu_q[6] = {gu00_q, gu01_q, gu02_q, gu10_q, gu11_q, gu12_q}; 
    const double gx_q[2] = {gxi0_q, gxi1_q}; 
    struct pkpm_harm_coeffs hc; 
    pkpm_harm_coeffs_set(&hc, 2, b_q, u_q, xi_q, th_q, d_q, gu_q, gb_q, gx_q, 0, 0.0, 0); 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_coupling_flux(&hc, w_lo, Ulo_q, Tlo); 
    pkpm_harm_coupling_flux(&hc, w_up, Uup_q, Tup); 
    pkpm_harm_speed_bounds(&hc, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_speed_bounds(&hc, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double ap_lo = a0p_q + cp_q*w_lo + k_lo, ap_up = a0p_q + cp_q*w_up + k_up; 
    double wf_lo = 0.0, wu_lo = 0.0, wf_up = 0.0, wu_up = 0.0; 
    { 
      const double sl = ap_lo + llo_lo, sr = ap_lo + lhi_lo; 
      if (sr < 0.0) { 
        wf_lo = 1.0; 
      } else if (sl < 0.0) { 
        wf_lo = -sl/(sr - sl); wu_lo = sl*sr/(sr - sl); 
      } 
    } 
    { 
      const double sl = ap_up + llo_up, sr = ap_up + lhi_up; 
      if (sl > 0.0) { 
        wf_up = 1.0; 
      } else if (sr > 0.0) { 
        wf_up = sr/(sr - sl); wu_up = -sl*sr/(sr - sl); 
      } 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) wf_e U_e, known part kU (k_e wf_e U_e) and 
    // x (coupling and HLL terms wf_e T_e + wu_e U_e) 
    const double eF0 = wf_lo*Flo_q - wf_up*Fup_q; 
    const double eF1 = wf_lo*w_lo*Flo_q - wf_up*w_up*Fup_q; 
    const double eF2 = wf_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*Fup_q; 
    const double eF3 = wf_lo*w_lo*w_lo*w_lo*Flo_q - wf_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = wf_lo*Glo_q - wf_up*Gup_q; 
    const double eG1 = wf_lo*w_lo*Glo_q - wf_up*w_up*Gup_q; 
    const double kF0 = wf_lo*k_lo*Flo_q - wf_up*k_up*Fup_q; 
    const double kF1 = wf_lo*w_lo*k_lo*Flo_q - wf_up*w_up*k_up*Fup_q; 
    const double kF2 = wf_lo*w_lo*w_lo*k_lo*Flo_q - wf_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = wf_lo*k_lo*Glo_q - wf_up*k_up*Gup_q; 
    double xlo[8], xup[8]; 
    for (int c = 0; c < 8; ++c) { 
      xlo[c] = wf_lo*Tlo[c] + wu_lo*Ulo_q[c]; 
      xup[c] = wf_up*Tup[c] + wu_up*Uup_q[c]; 
    } 
    const double xF0 = xlo[0] - xup[0]; 
    const double xF1 = w_lo*xlo[0] - w_up*xup[0]; 
    const double xF2 = w_lo*w_lo*xlo[0] - w_up*w_up*xup[0]; 
    const double xG0 = xlo[1] - xup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) 
      v[a] = rM_q[a] + (wf_lo*k_lo*Ulo_q[2 + a] + xlo[2 + a]) - (wf_up*k_up*Uup_q[2 + a] + xup[2 + a]); 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*Ulo_q[2 + a] - wf_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = wf_lo*w_lo*Ulo_q[2 + a] - wf_up*w_up*Uup_q[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - kF0 - xF0, b1 = -r1_q - kF1 - xF1, b2 = -r2_q - kF2 - 2.0*kG0 - xF2 - 2.0*xG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[3] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[3] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[3] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_n[a][3] = -(P[a] + (c0_n[3] + c_n[3])*Q[a] + Aa[a]*a0_n[3] + Ac[a]*c_n[3])/n0_q; 
      const double al_lo = a0_n[3] + c_n[3]*w_lo + k_lo, al_up = a0_n[3] + c_n[3]*w_up + k_up; 
      for (int c = 0; c < 8; ++c) { 
        flo_n[c][3] = wf_lo*al_lo*Ulo_q[c] + xlo[c]; 
        fup_n[c][3] = wf_up*al_up*Uup_q[c] + xup[c]; 
      } 
      cfl_lo = fmax(cfl_lo, fmax(fabs(al_lo + llo_lo), fabs(al_lo + lhi_lo))); 
      cfl_up = fmax(cfl_up, fmax(fabs(al_up + llo_up), fabs(al_up + lhi_up))); 
    } else { 
      num_null += 1; 
    } 
  } 
 
  cfl[0] = 2.5*dv1*cfl_lo; 
  cfl[1] = 2.5*dv1*cfl_up; 
  cond[0] = cond_min; 
  if (num_null > 0) { 
    // a node with M0 <= 0 or a non-positive determinant zeroes everything in the cell, its 
    // CFL rate included (the closures it would have needed are not applied) 
    cfl[0] = 0.0; cfl[1] = 0.0; 
    for (int k = 0; k < 4; ++k) { 
      c_0[k] = 0.0; alpha_0[k] = 0.0; c_h[k] = 0.0; 
    } 
    for (int k = 0; k < 12; ++k) { 
      a_perp[k] = 0.0; 
    } 
    for (int k = 0; k < 32; ++k) { 
      flux_lo[k] = 0.0; flux_up[k] = 0.0; 
    } 
    return num_null; 
  } 
 
  c_0[0] = 0.5*c0_n[3]+0.5*c0_n[2]+0.5*c0_n[1]+0.5*c0_n[0]; 
  c_0[1] = 0.5*c0_n[3]+0.5*c0_n[2]-0.5*c0_n[1]-0.5*c0_n[0]; 
  c_0[2] = 0.5*c0_n[3]-0.5*c0_n[2]+0.5*c0_n[1]-0.5*c0_n[0]; 
  c_0[3] = 0.5*c0_n[3]-0.5*c0_n[2]-0.5*c0_n[1]+0.5*c0_n[0]; 
  alpha_0[0] = 0.5*a0_n[3]+0.5*a0_n[2]+0.5*a0_n[1]+0.5*a0_n[0]; 
  alpha_0[1] = 0.5*a0_n[3]+0.5*a0_n[2]-0.5*a0_n[1]-0.5*a0_n[0]; 
  alpha_0[2] = 0.5*a0_n[3]-0.5*a0_n[2]+0.5*a0_n[1]-0.5*a0_n[0]; 
  alpha_0[3] = 0.5*a0_n[3]-0.5*a0_n[2]-0.5*a0_n[1]+0.5*a0_n[0]; 
  c_h[0] = 0.5*c_n[3]+0.5*c_n[2]+0.5*c_n[1]+0.5*c_n[0]; 
  c_h[1] = 0.5*c_n[3]+0.5*c_n[2]-0.5*c_n[1]-0.5*c_n[0]; 
  c_h[2] = 0.5*c_n[3]-0.5*c_n[2]+0.5*c_n[1]-0.5*c_n[0]; 
  c_h[3] = 0.5*c_n[3]-0.5*c_n[2]-0.5*c_n[1]+0.5*c_n[0]; 
 
  for (int a = 0; a < 3; ++a) { 
    a_perp[a*4 + 0] = (0.5)*ap_n[a][0] + (0.5)*ap_n[a][1] + (0.5)*ap_n[a][2] + (0.5)*ap_n[a][3]; 
    a_perp[a*4 + 1] = (-0.5)*ap_n[a][0] + (-0.5)*ap_n[a][1] + (0.5)*ap_n[a][2] + (0.5)*ap_n[a][3]; 
    a_perp[a*4 + 2] = (-0.5)*ap_n[a][0] + (0.5)*ap_n[a][1] + (-0.5)*ap_n[a][2] + (0.5)*ap_n[a][3]; 
    a_perp[a*4 + 3] = (0.5)*ap_n[a][0] + (-0.5)*ap_n[a][1] + (-0.5)*ap_n[a][2] + (0.5)*ap_n[a][3]; 
  } 
  for (int c = 0; c < 8; ++c) { 
    flux_lo[c*4 + 0] = (1.0)*flo_n[c][0] + (1.0)*flo_n[c][1] + (1.0)*flo_n[c][2] + (1.0)*flo_n[c][3]; 
    flux_up[c*4 + 0] = (1.0)*fup_n[c][0] + (1.0)*fup_n[c][1] + (1.0)*fup_n[c][2] + (1.0)*fup_n[c][3]; 
    flux_lo[c*4 + 1] = (-1.4433756729740645)*flo_n[c][0] + (-1.4433756729740645)*flo_n[c][1] + (1.4433756729740645)*flo_n[c][2] + (1.4433756729740645)*flo_n[c][3]; 
    flux_up[c*4 + 1] = (-1.4433756729740645)*fup_n[c][0] + (-1.4433756729740645)*fup_n[c][1] + (1.4433756729740645)*fup_n[c][2] + (1.4433756729740645)*fup_n[c][3]; 
    flux_lo[c*4 + 2] = (-1.4433756729740645)*flo_n[c][0] + (1.4433756729740645)*flo_n[c][1] + (-1.4433756729740645)*flo_n[c][2] + (1.4433756729740645)*flo_n[c][3]; 
    flux_up[c*4 + 2] = (-1.4433756729740645)*fup_n[c][0] + (1.4433756729740645)*fup_n[c][1] + (-1.4433756729740645)*fup_n[c][2] + (1.4433756729740645)*fup_n[c][3]; 
    flux_lo[c*4 + 3] = (2.0833333333333335)*flo_n[c][0] + (-2.0833333333333335)*flo_n[c][1] + (-2.0833333333333335)*flo_n[c][2] + (2.0833333333333335)*flo_n[c][3]; 
    flux_up[c*4 + 3] = (2.0833333333333335)*fup_n[c][0] + (-2.0833333333333335)*fup_n[c][1] + (-2.0833333333333335)*fup_n[c][2] + (2.0833333333333335)*fup_n[c][3]; 
  } 
 
  return 0; 
} 
