#include <gkyl_pkpm_vars_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH int pkpm_vars_closure_1x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *coef, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl,
  double* GKYL_RESTRICT cond) 
{ 
  // w_lo, w_up:   Velocity domain edges; dv: velocity cell size. 
  // norm_moms:    Input [M0, M2, M0_G, M1, M3, M1_G, mu (x, y, z)] moments of the normalized distributions (stage state). 
  // rhs_moms:     Input moments of the stage-A right-hand sides, same layout. 
  // bvar:         Magnetic field unit vector [bx, by, bz]. 
  // coef:         Stored operator of the cell (gkyl_pkpm_harm_ops.h). 
  // theta:        T_perp/T. 
  // alpha_0_prev, c_h_prev: closures of the previous stage (may alias the outputs), for the edge speeds. 
  // f_lo, f_up:   [F, G, M0_perp, N_perp] of the lowest and uppermost velocity cells of this configuration cell. 
  // c_0, alpha_0, c_h, a_perp: Output closures at the 2 Gauss-Legendre nodes, interpolated to p=1. 
  // flux_lo, flux_up: Output edge fluxes [8 blocks] in the stored face format. 
  // cfl:          Output [lower, upper] edge-flux time-step rates. 
  // cond:         Output [1]: the smallest determinant of the 3x3 system over the nodes relative to 
  //               its edge-free value (1 when the edge terms vanish; small or negative when the 
  //               velocity-domain edges carry a non-negligible trace: too small velocity extents). 
  // Returns the number of nodes where the closures are undefined (everything is zero in the cell then). 
 
  const double dv1 = 2.0/dv; 
  const double *n0 = &norm_moms[0]; 
  const double *n2 = &norm_moms[2]; 
  const double *nG = &norm_moms[4]; 
  const double *n1 = &norm_moms[6]; 
  const double *r0 = &rhs_moms[0]; 
  const double *r2 = &rhs_moms[2]; 
  const double *rG = &rhs_moms[4]; 
  const double *r1 = &rhs_moms[6]; 
  const double *mu0 = &norm_moms[12]; 
  const double *mu1 = &norm_moms[14]; 
  const double *mu2 = &norm_moms[16]; 
  const double *rM0 = &rhs_moms[12]; 
  const double *rM1 = &rhs_moms[14]; 
  const double *rM2 = &rhs_moms[16]; 
  const double *b0 = &bvar[0]; 
  const double *b1 = &bvar[2]; 
  const double *b2 = &bvar[4]; 
  const double *lo0 = &f_lo[0], *up0 = &f_up[0]; 
  const double *lo1 = &f_lo[6], *up1 = &f_up[6]; 
  const double *lo2 = &f_lo[12], *up2 = &f_up[12]; 
  const double *lo3 = &f_lo[18], *up3 = &f_up[18]; 
  const double *lo4 = &f_lo[24], *up4 = &f_up[24]; 
  const double *lo5 = &f_lo[30], *up5 = &f_up[30]; 
  const double *lo6 = &f_lo[36], *up6 = &f_up[36]; 
  const double *lo7 = &f_lo[42], *up7 = &f_up[42]; 
 
  double c0_n[2] = {0.0}, a0_n[2] = {0.0}, c_n[2] = {0.0}; 
  double ap_n[3][2] = {{0.0}}, flo_n[8][2] = {{0.0}}, fup_n[8][2] = {{0.0}}; 
  double cfl_lo = 0.0, cfl_up = 0.0; 
  int num_null = 0; 
  double cond_min = 1.0; 
 
  { 
    const double n0_q = 0.7071067811865475*n0[0]-0.7071067811865475*n0[1]; 
    const double n1_q = 0.7071067811865475*n1[0]-0.7071067811865475*n1[1]; 
    const double n2_q = 0.7071067811865475*n2[0]-0.7071067811865475*n2[1]; 
    const double nG_q = 0.7071067811865475*nG[0]-0.7071067811865475*nG[1]; 
    const double r0_q = 0.7071067811865475*r0[0]-0.7071067811865475*r0[1]; 
    const double r1_q = 0.7071067811865475*r1[0]-0.7071067811865475*r1[1]; 
    const double r2_q = 0.7071067811865475*r2[0]-0.7071067811865475*r2[1] + 2.0*(0.7071067811865475*rG[0]-0.7071067811865475*rG[1]); 
    const double mu_q[3] = {0.7071067811865475*mu0[0]-0.7071067811865475*mu0[1], 0.7071067811865475*mu1[0]-0.7071067811865475*mu1[1], 0.7071067811865475*mu2[0]-0.7071067811865475*mu2[1]}; 
    const double rM_q[3] = {0.7071067811865475*rM0[0]-0.7071067811865475*rM0[1], 0.7071067811865475*rM1[0]-0.7071067811865475*rM1[1], 0.7071067811865475*rM2[0]-0.7071067811865475*rM2[1]}; 
    const double a0p_q = 0.7071067811865475*alpha_0_prev[0]-0.7071067811865475*alpha_0_prev[1]; 
    const double cp_q = 0.7071067811865475*c_h_prev[0]-0.7071067811865475*c_h_prev[1]; 
    const double th_q = 0.7071067811865475*theta[0]-0.7071067811865475*theta[1]; 
    const double b0_q = 0.7071067811865475*b0[0]-0.7071067811865475*b0[1]; 
    const double b1_q = 0.7071067811865475*b1[0]-0.7071067811865475*b1[1]; 
    const double b2_q = 0.7071067811865475*b2[0]-0.7071067811865475*b2[1]; 
    const double Ulo_q[8] = {-(1.1180339887498951*lo0[5])+1.118033988749895*lo0[4]+0.8660254037844386*lo0[3]-0.8660254037844386*lo0[2]-0.5*lo0[1]+0.5*lo0[0], -(1.1180339887498951*lo1[5])+1.118033988749895*lo1[4]+0.8660254037844386*lo1[3]-0.8660254037844386*lo1[2]-0.5*lo1[1]+0.5*lo1[0], -(1.1180339887498951*lo2[5])+1.118033988749895*lo2[4]+0.8660254037844386*lo2[3]-0.8660254037844386*lo2[2]-0.5*lo2[1]+0.5*lo2[0], -(1.1180339887498951*lo3[5])+1.118033988749895*lo3[4]+0.8660254037844386*lo3[3]-0.8660254037844386*lo3[2]-0.5*lo3[1]+0.5*lo3[0], -(1.1180339887498951*lo4[5])+1.118033988749895*lo4[4]+0.8660254037844386*lo4[3]-0.8660254037844386*lo4[2]-0.5*lo4[1]+0.5*lo4[0], -(1.1180339887498951*lo5[5])+1.118033988749895*lo5[4]+0.8660254037844386*lo5[3]-0.8660254037844386*lo5[2]-0.5*lo5[1]+0.5*lo5[0], -(1.1180339887498951*lo6[5])+1.118033988749895*lo6[4]+0.8660254037844386*lo6[3]-0.8660254037844386*lo6[2]-0.5*lo6[1]+0.5*lo6[0], -(1.1180339887498951*lo7[5])+1.118033988749895*lo7[4]+0.8660254037844386*lo7[3]-0.8660254037844386*lo7[2]-0.5*lo7[1]+0.5*lo7[0]}; 
    const double Uup_q[8] = {-(1.1180339887498951*up0[5])+1.118033988749895*up0[4]-0.8660254037844386*up0[3]+0.8660254037844386*up0[2]-0.5*up0[1]+0.5*up0[0], -(1.1180339887498951*up1[5])+1.118033988749895*up1[4]-0.8660254037844386*up1[3]+0.8660254037844386*up1[2]-0.5*up1[1]+0.5*up1[0], -(1.1180339887498951*up2[5])+1.118033988749895*up2[4]-0.8660254037844386*up2[3]+0.8660254037844386*up2[2]-0.5*up2[1]+0.5*up2[0], -(1.1180339887498951*up3[5])+1.118033988749895*up3[4]-0.8660254037844386*up3[3]+0.8660254037844386*up3[2]-0.5*up3[1]+0.5*up3[0], -(1.1180339887498951*up4[5])+1.118033988749895*up4[4]-0.8660254037844386*up4[3]+0.8660254037844386*up4[2]-0.5*up4[1]+0.5*up4[0], -(1.1180339887498951*up5[5])+1.118033988749895*up5[4]-0.8660254037844386*up5[3]+0.8660254037844386*up5[2]-0.5*up5[1]+0.5*up5[0], -(1.1180339887498951*up6[5])+1.118033988749895*up6[4]-0.8660254037844386*up6[3]+0.8660254037844386*up6[2]-0.5*up6[1]+0.5*up6[0], -(1.1180339887498951*up7[5])+1.118033988749895*up7[4]-0.8660254037844386*up7[3]+0.8660254037844386*up7[2]-0.5*up7[1]+0.5*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = 0.7071067811865475*coef[0]-0.7071067811865475*coef[1]; 
    nv[1] = 0.7071067811865475*coef[3]-0.7071067811865475*coef[4]; 
    nv[2] = 0.7071067811865475*coef[66]-0.7071067811865475*coef[67]; 
    nv[3] = 0.7071067811865475*coef[81]-0.7071067811865475*coef[82]; 
    nv[4] = 0.7071067811865475*coef[192]-0.7071067811865475*coef[193]; 
    nv[5] = 0.7071067811865475*coef[207]-0.7071067811865475*coef[208]; 
    nv[6] = 0.7071067811865475*coef[318]-0.7071067811865475*coef[319]; 
    nv[7] = 0.7071067811865475*coef[333]-0.7071067811865475*coef[334]; 
    nv[8] = 0.7071067811865475*coef[69]-0.7071067811865475*coef[70]; 
    nv[9] = 0.7071067811865475*coef[84]-0.7071067811865475*coef[85]; 
    nv[10] = 0.7071067811865475*coef[195]-0.7071067811865475*coef[196]; 
    nv[11] = 0.7071067811865475*coef[210]-0.7071067811865475*coef[211]; 
    nv[12] = 0.7071067811865475*coef[321]-0.7071067811865475*coef[322]; 
    nv[13] = 0.7071067811865475*coef[336]-0.7071067811865475*coef[337]; 
    nv[14] = 0.7071067811865475*coef[72]-0.7071067811865475*coef[73]; 
    nv[15] = 0.7071067811865475*coef[87]-0.7071067811865475*coef[88]; 
    nv[16] = 0.7071067811865475*coef[198]-0.7071067811865475*coef[199]; 
    nv[17] = 0.7071067811865475*coef[213]-0.7071067811865475*coef[214]; 
    nv[18] = 0.7071067811865475*coef[324]-0.7071067811865475*coef[325]; 
    nv[19] = 0.7071067811865475*coef[339]-0.7071067811865475*coef[340]; 
    nv[20] = 0.7071067811865475*coef[75]-0.7071067811865475*coef[76]; 
    nv[21] = 0.7071067811865475*coef[90]-0.7071067811865475*coef[91]; 
    nv[22] = 0.7071067811865475*coef[201]-0.7071067811865475*coef[202]; 
    nv[23] = 0.7071067811865475*coef[216]-0.7071067811865475*coef[217]; 
    nv[24] = 0.7071067811865475*coef[327]-0.7071067811865475*coef[328]; 
    nv[25] = 0.7071067811865475*coef[342]-0.7071067811865475*coef[343]; 
    nv[26] = 0.7071067811865475*coef[78]-0.7071067811865475*coef[79]; 
    nv[27] = 0.7071067811865475*coef[93]-0.7071067811865475*coef[94]; 
    nv[28] = 0.7071067811865475*coef[204]-0.7071067811865475*coef[205]; 
    nv[29] = 0.7071067811865475*coef[219]-0.7071067811865475*coef[220]; 
    nv[30] = 0.7071067811865475*coef[330]-0.7071067811865475*coef[331]; 
    nv[31] = 0.7071067811865475*coef[345]-0.7071067811865475*coef[346]; 
    nv[32] = 0.7071067811865475*coef[96]-0.7071067811865475*coef[97]; 
    nv[33] = 0.7071067811865475*coef[105]-0.7071067811865475*coef[106]; 
    nv[34] = 0.7071067811865475*coef[114]-0.7071067811865475*coef[115]; 
    nv[35] = 0.7071067811865475*coef[222]-0.7071067811865475*coef[223]; 
    nv[36] = 0.7071067811865475*coef[231]-0.7071067811865475*coef[232]; 
    nv[37] = 0.7071067811865475*coef[240]-0.7071067811865475*coef[241]; 
    nv[38] = 0.7071067811865475*coef[348]-0.7071067811865475*coef[349]; 
    nv[39] = 0.7071067811865475*coef[357]-0.7071067811865475*coef[358]; 
    nv[40] = 0.7071067811865475*coef[366]-0.7071067811865475*coef[367]; 
    nv[41] = 0.7071067811865475*coef[99]-0.7071067811865475*coef[100]; 
    nv[42] = 0.7071067811865475*coef[108]-0.7071067811865475*coef[109]; 
    nv[43] = 0.7071067811865475*coef[117]-0.7071067811865475*coef[118]; 
    nv[44] = 0.7071067811865475*coef[225]-0.7071067811865475*coef[226]; 
    nv[45] = 0.7071067811865475*coef[234]-0.7071067811865475*coef[235]; 
    nv[46] = 0.7071067811865475*coef[243]-0.7071067811865475*coef[244]; 
    nv[47] = 0.7071067811865475*coef[351]-0.7071067811865475*coef[352]; 
    nv[48] = 0.7071067811865475*coef[360]-0.7071067811865475*coef[361]; 
    nv[49] = 0.7071067811865475*coef[369]-0.7071067811865475*coef[370]; 
    nv[50] = 0.7071067811865475*coef[102]-0.7071067811865475*coef[103]; 
    nv[51] = 0.7071067811865475*coef[111]-0.7071067811865475*coef[112]; 
    nv[52] = 0.7071067811865475*coef[120]-0.7071067811865475*coef[121]; 
    nv[53] = 0.7071067811865475*coef[228]-0.7071067811865475*coef[229]; 
    nv[54] = 0.7071067811865475*coef[237]-0.7071067811865475*coef[238]; 
    nv[55] = 0.7071067811865475*coef[246]-0.7071067811865475*coef[247]; 
    nv[56] = 0.7071067811865475*coef[354]-0.7071067811865475*coef[355]; 
    nv[57] = 0.7071067811865475*coef[363]-0.7071067811865475*coef[364]; 
    nv[58] = 0.7071067811865475*coef[372]-0.7071067811865475*coef[373]; 
    nv[59] = 0.7071067811865475*coef[12]-0.7071067811865475*coef[13]; 
    nv[60] = 0.7071067811865475*coef[15]-0.7071067811865475*coef[16]; 
    nv[61] = 0.7071067811865475*coef[18]-0.7071067811865475*coef[19]; 
    nv[PKPM_HS_NODE_THETA] = th_q; 
    nv[PKPM_HS_NODE_B] = b0_q; nv[PKPM_HS_NODE_B + 1] = b1_q; nv[PKPM_HS_NODE_B + 2] = b2_q; 
    nv[PKPM_HS_NODE_A0] = a0p_q; nv[PKPM_HS_NODE_CH] = cp_q; 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_node_coupling(nv, w_lo, Ulo_q, Tlo); 
    pkpm_harm_node_coupling(nv, w_up, Uup_q, Tup); 
    pkpm_harm_node_bounds(nv, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_node_bounds(nv, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = pkpm_harm_node_alpha_k(nv, w_lo), k_up = pkpm_harm_node_alpha_k(nv, w_up); 
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
    const double n0_q = 0.7071067811865475*n0[1]+0.7071067811865475*n0[0]; 
    const double n1_q = 0.7071067811865475*n1[1]+0.7071067811865475*n1[0]; 
    const double n2_q = 0.7071067811865475*n2[1]+0.7071067811865475*n2[0]; 
    const double nG_q = 0.7071067811865475*nG[1]+0.7071067811865475*nG[0]; 
    const double r0_q = 0.7071067811865475*r0[1]+0.7071067811865475*r0[0]; 
    const double r1_q = 0.7071067811865475*r1[1]+0.7071067811865475*r1[0]; 
    const double r2_q = 0.7071067811865475*r2[1]+0.7071067811865475*r2[0] + 2.0*(0.7071067811865475*rG[1]+0.7071067811865475*rG[0]); 
    const double mu_q[3] = {0.7071067811865475*mu0[1]+0.7071067811865475*mu0[0], 0.7071067811865475*mu1[1]+0.7071067811865475*mu1[0], 0.7071067811865475*mu2[1]+0.7071067811865475*mu2[0]}; 
    const double rM_q[3] = {0.7071067811865475*rM0[1]+0.7071067811865475*rM0[0], 0.7071067811865475*rM1[1]+0.7071067811865475*rM1[0], 0.7071067811865475*rM2[1]+0.7071067811865475*rM2[0]}; 
    const double a0p_q = 0.7071067811865475*alpha_0_prev[1]+0.7071067811865475*alpha_0_prev[0]; 
    const double cp_q = 0.7071067811865475*c_h_prev[1]+0.7071067811865475*c_h_prev[0]; 
    const double th_q = 0.7071067811865475*theta[1]+0.7071067811865475*theta[0]; 
    const double b0_q = 0.7071067811865475*b0[1]+0.7071067811865475*b0[0]; 
    const double b1_q = 0.7071067811865475*b1[1]+0.7071067811865475*b1[0]; 
    const double b2_q = 0.7071067811865475*b2[1]+0.7071067811865475*b2[0]; 
    const double Ulo_q[8] = {1.1180339887498951*lo0[5]+1.118033988749895*lo0[4]-0.8660254037844386*lo0[3]-0.8660254037844386*lo0[2]+0.5*lo0[1]+0.5*lo0[0], 1.1180339887498951*lo1[5]+1.118033988749895*lo1[4]-0.8660254037844386*lo1[3]-0.8660254037844386*lo1[2]+0.5*lo1[1]+0.5*lo1[0], 1.1180339887498951*lo2[5]+1.118033988749895*lo2[4]-0.8660254037844386*lo2[3]-0.8660254037844386*lo2[2]+0.5*lo2[1]+0.5*lo2[0], 1.1180339887498951*lo3[5]+1.118033988749895*lo3[4]-0.8660254037844386*lo3[3]-0.8660254037844386*lo3[2]+0.5*lo3[1]+0.5*lo3[0], 1.1180339887498951*lo4[5]+1.118033988749895*lo4[4]-0.8660254037844386*lo4[3]-0.8660254037844386*lo4[2]+0.5*lo4[1]+0.5*lo4[0], 1.1180339887498951*lo5[5]+1.118033988749895*lo5[4]-0.8660254037844386*lo5[3]-0.8660254037844386*lo5[2]+0.5*lo5[1]+0.5*lo5[0], 1.1180339887498951*lo6[5]+1.118033988749895*lo6[4]-0.8660254037844386*lo6[3]-0.8660254037844386*lo6[2]+0.5*lo6[1]+0.5*lo6[0], 1.1180339887498951*lo7[5]+1.118033988749895*lo7[4]-0.8660254037844386*lo7[3]-0.8660254037844386*lo7[2]+0.5*lo7[1]+0.5*lo7[0]}; 
    const double Uup_q[8] = {1.1180339887498951*up0[5]+1.118033988749895*up0[4]+0.8660254037844386*up0[3]+0.8660254037844386*up0[2]+0.5*up0[1]+0.5*up0[0], 1.1180339887498951*up1[5]+1.118033988749895*up1[4]+0.8660254037844386*up1[3]+0.8660254037844386*up1[2]+0.5*up1[1]+0.5*up1[0], 1.1180339887498951*up2[5]+1.118033988749895*up2[4]+0.8660254037844386*up2[3]+0.8660254037844386*up2[2]+0.5*up2[1]+0.5*up2[0], 1.1180339887498951*up3[5]+1.118033988749895*up3[4]+0.8660254037844386*up3[3]+0.8660254037844386*up3[2]+0.5*up3[1]+0.5*up3[0], 1.1180339887498951*up4[5]+1.118033988749895*up4[4]+0.8660254037844386*up4[3]+0.8660254037844386*up4[2]+0.5*up4[1]+0.5*up4[0], 1.1180339887498951*up5[5]+1.118033988749895*up5[4]+0.8660254037844386*up5[3]+0.8660254037844386*up5[2]+0.5*up5[1]+0.5*up5[0], 1.1180339887498951*up6[5]+1.118033988749895*up6[4]+0.8660254037844386*up6[3]+0.8660254037844386*up6[2]+0.5*up6[1]+0.5*up6[0], 1.1180339887498951*up7[5]+1.118033988749895*up7[4]+0.8660254037844386*up7[3]+0.8660254037844386*up7[2]+0.5*up7[1]+0.5*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = 0.7071067811865475*coef[1]+0.7071067811865475*coef[0]; 
    nv[1] = 0.7071067811865475*coef[4]+0.7071067811865475*coef[3]; 
    nv[2] = 0.7071067811865475*coef[67]+0.7071067811865475*coef[66]; 
    nv[3] = 0.7071067811865475*coef[82]+0.7071067811865475*coef[81]; 
    nv[4] = 0.7071067811865475*coef[193]+0.7071067811865475*coef[192]; 
    nv[5] = 0.7071067811865475*coef[208]+0.7071067811865475*coef[207]; 
    nv[6] = 0.7071067811865475*coef[319]+0.7071067811865475*coef[318]; 
    nv[7] = 0.7071067811865475*coef[334]+0.7071067811865475*coef[333]; 
    nv[8] = 0.7071067811865475*coef[70]+0.7071067811865475*coef[69]; 
    nv[9] = 0.7071067811865475*coef[85]+0.7071067811865475*coef[84]; 
    nv[10] = 0.7071067811865475*coef[196]+0.7071067811865475*coef[195]; 
    nv[11] = 0.7071067811865475*coef[211]+0.7071067811865475*coef[210]; 
    nv[12] = 0.7071067811865475*coef[322]+0.7071067811865475*coef[321]; 
    nv[13] = 0.7071067811865475*coef[337]+0.7071067811865475*coef[336]; 
    nv[14] = 0.7071067811865475*coef[73]+0.7071067811865475*coef[72]; 
    nv[15] = 0.7071067811865475*coef[88]+0.7071067811865475*coef[87]; 
    nv[16] = 0.7071067811865475*coef[199]+0.7071067811865475*coef[198]; 
    nv[17] = 0.7071067811865475*coef[214]+0.7071067811865475*coef[213]; 
    nv[18] = 0.7071067811865475*coef[325]+0.7071067811865475*coef[324]; 
    nv[19] = 0.7071067811865475*coef[340]+0.7071067811865475*coef[339]; 
    nv[20] = 0.7071067811865475*coef[76]+0.7071067811865475*coef[75]; 
    nv[21] = 0.7071067811865475*coef[91]+0.7071067811865475*coef[90]; 
    nv[22] = 0.7071067811865475*coef[202]+0.7071067811865475*coef[201]; 
    nv[23] = 0.7071067811865475*coef[217]+0.7071067811865475*coef[216]; 
    nv[24] = 0.7071067811865475*coef[328]+0.7071067811865475*coef[327]; 
    nv[25] = 0.7071067811865475*coef[343]+0.7071067811865475*coef[342]; 
    nv[26] = 0.7071067811865475*coef[79]+0.7071067811865475*coef[78]; 
    nv[27] = 0.7071067811865475*coef[94]+0.7071067811865475*coef[93]; 
    nv[28] = 0.7071067811865475*coef[205]+0.7071067811865475*coef[204]; 
    nv[29] = 0.7071067811865475*coef[220]+0.7071067811865475*coef[219]; 
    nv[30] = 0.7071067811865475*coef[331]+0.7071067811865475*coef[330]; 
    nv[31] = 0.7071067811865475*coef[346]+0.7071067811865475*coef[345]; 
    nv[32] = 0.7071067811865475*coef[97]+0.7071067811865475*coef[96]; 
    nv[33] = 0.7071067811865475*coef[106]+0.7071067811865475*coef[105]; 
    nv[34] = 0.7071067811865475*coef[115]+0.7071067811865475*coef[114]; 
    nv[35] = 0.7071067811865475*coef[223]+0.7071067811865475*coef[222]; 
    nv[36] = 0.7071067811865475*coef[232]+0.7071067811865475*coef[231]; 
    nv[37] = 0.7071067811865475*coef[241]+0.7071067811865475*coef[240]; 
    nv[38] = 0.7071067811865475*coef[349]+0.7071067811865475*coef[348]; 
    nv[39] = 0.7071067811865475*coef[358]+0.7071067811865475*coef[357]; 
    nv[40] = 0.7071067811865475*coef[367]+0.7071067811865475*coef[366]; 
    nv[41] = 0.7071067811865475*coef[100]+0.7071067811865475*coef[99]; 
    nv[42] = 0.7071067811865475*coef[109]+0.7071067811865475*coef[108]; 
    nv[43] = 0.7071067811865475*coef[118]+0.7071067811865475*coef[117]; 
    nv[44] = 0.7071067811865475*coef[226]+0.7071067811865475*coef[225]; 
    nv[45] = 0.7071067811865475*coef[235]+0.7071067811865475*coef[234]; 
    nv[46] = 0.7071067811865475*coef[244]+0.7071067811865475*coef[243]; 
    nv[47] = 0.7071067811865475*coef[352]+0.7071067811865475*coef[351]; 
    nv[48] = 0.7071067811865475*coef[361]+0.7071067811865475*coef[360]; 
    nv[49] = 0.7071067811865475*coef[370]+0.7071067811865475*coef[369]; 
    nv[50] = 0.7071067811865475*coef[103]+0.7071067811865475*coef[102]; 
    nv[51] = 0.7071067811865475*coef[112]+0.7071067811865475*coef[111]; 
    nv[52] = 0.7071067811865475*coef[121]+0.7071067811865475*coef[120]; 
    nv[53] = 0.7071067811865475*coef[229]+0.7071067811865475*coef[228]; 
    nv[54] = 0.7071067811865475*coef[238]+0.7071067811865475*coef[237]; 
    nv[55] = 0.7071067811865475*coef[247]+0.7071067811865475*coef[246]; 
    nv[56] = 0.7071067811865475*coef[355]+0.7071067811865475*coef[354]; 
    nv[57] = 0.7071067811865475*coef[364]+0.7071067811865475*coef[363]; 
    nv[58] = 0.7071067811865475*coef[373]+0.7071067811865475*coef[372]; 
    nv[59] = 0.7071067811865475*coef[13]+0.7071067811865475*coef[12]; 
    nv[60] = 0.7071067811865475*coef[16]+0.7071067811865475*coef[15]; 
    nv[61] = 0.7071067811865475*coef[19]+0.7071067811865475*coef[18]; 
    nv[PKPM_HS_NODE_THETA] = th_q; 
    nv[PKPM_HS_NODE_B] = b0_q; nv[PKPM_HS_NODE_B + 1] = b1_q; nv[PKPM_HS_NODE_B + 2] = b2_q; 
    nv[PKPM_HS_NODE_A0] = a0p_q; nv[PKPM_HS_NODE_CH] = cp_q; 
    double Tlo[8], Tup[8], llo_lo, lhi_lo, llo_up, lhi_up; 
    pkpm_harm_node_coupling(nv, w_lo, Ulo_q, Tlo); 
    pkpm_harm_node_coupling(nv, w_up, Uup_q, Tup); 
    pkpm_harm_node_bounds(nv, w_lo, &llo_lo, &lhi_lo); 
    pkpm_harm_node_bounds(nv, w_up, &llo_up, &lhi_up); 
    // known part of alpha at the edges, and the HLL weights against the empty exterior from the 
    // previous closures: edge flux = wf (alpha U + T(U)) + wu U 
    const double k_lo = pkpm_harm_node_alpha_k(nv, w_lo), k_up = pkpm_harm_node_alpha_k(nv, w_up); 
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
 
  cfl[0] = 2.5*dv1*cfl_lo; 
  cfl[1] = 2.5*dv1*cfl_up; 
  cond[0] = cond_min; 
  if (num_null > 0) { 
    // a node with M0 <= 0 or a non-positive determinant zeroes everything in the cell, its 
    // CFL rate included (the closures it would have needed are not applied) 
    cfl[0] = 0.0; cfl[1] = 0.0; 
    for (int k = 0; k < 2; ++k) { 
      c_0[k] = 0.0; alpha_0[k] = 0.0; c_h[k] = 0.0; 
    } 
    for (int k = 0; k < 6; ++k) { 
      a_perp[k] = 0.0; 
    } 
    for (int k = 0; k < 16; ++k) { 
      flux_lo[k] = 0.0; flux_up[k] = 0.0; 
    } 
    return num_null; 
  } 
 
  c_0[0] = 0.7071067811865475*c0_n[1]+0.7071067811865475*c0_n[0]; 
  c_0[1] = 0.7071067811865475*c0_n[1]-0.7071067811865475*c0_n[0]; 
  alpha_0[0] = 0.7071067811865475*a0_n[1]+0.7071067811865475*a0_n[0]; 
  alpha_0[1] = 0.7071067811865475*a0_n[1]-0.7071067811865475*a0_n[0]; 
  c_h[0] = 0.7071067811865475*c_n[1]+0.7071067811865475*c_n[0]; 
  c_h[1] = 0.7071067811865475*c_n[1]-0.7071067811865475*c_n[0]; 
 
  for (int a = 0; a < 3; ++a) { 
    a_perp[a*2 + 0] = (0.7071067811865475)*ap_n[a][0] + (0.7071067811865475)*ap_n[a][1]; 
    a_perp[a*2 + 1] = (-0.7071067811865475)*ap_n[a][0] + (0.7071067811865475)*ap_n[a][1]; 
  } 
  for (int c = 0; c < 8; ++c) { 
    flux_lo[c*2 + 0] = (1.0)*flo_n[c][0] + (1.0)*flo_n[c][1]; 
    flux_up[c*2 + 0] = (1.0)*fup_n[c][0] + (1.0)*fup_n[c][1]; 
    flux_lo[c*2 + 1] = (-1.4433756729740645)*flo_n[c][0] + (1.4433756729740645)*flo_n[c][1]; 
    flux_up[c*2 + 1] = (-1.4433756729740645)*fup_n[c][0] + (1.4433756729740645)*fup_n[c][1]; 
  } 
 
  return 0; 
} 
