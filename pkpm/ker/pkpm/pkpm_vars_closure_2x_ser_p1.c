#include <gkyl_pkpm_vars_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH int pkpm_vars_closure_2x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
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
    const double th_q = 0.5*theta[3]-0.5*theta[2]-0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = 0.5*b0[3]-0.5*b0[2]-0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = 0.5*b1[3]-0.5*b1[2]-0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = 0.5*b2[3]-0.5*b2[2]-0.5*b2[1]+0.5*b2[0]; 
    const double Ulo_q[8] = {0.7905694150420947*lo0[11]-0.7905694150420948*lo0[10]-0.7905694150420948*lo0[9]-0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]+0.6123724356957944*lo0[6]+0.6123724356957944*lo0[5]+0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]-0.3535533905932737*lo0[2]-0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], 0.7905694150420947*lo1[11]-0.7905694150420948*lo1[10]-0.7905694150420948*lo1[9]-0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]+0.6123724356957944*lo1[6]+0.6123724356957944*lo1[5]+0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]-0.3535533905932737*lo1[2]-0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], 0.7905694150420947*lo2[11]-0.7905694150420948*lo2[10]-0.7905694150420948*lo2[9]-0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]+0.6123724356957944*lo2[6]+0.6123724356957944*lo2[5]+0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]-0.3535533905932737*lo2[2]-0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], 0.7905694150420947*lo3[11]-0.7905694150420948*lo3[10]-0.7905694150420948*lo3[9]-0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]+0.6123724356957944*lo3[6]+0.6123724356957944*lo3[5]+0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]-0.3535533905932737*lo3[2]-0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], 0.7905694150420947*lo4[11]-0.7905694150420948*lo4[10]-0.7905694150420948*lo4[9]-0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]+0.6123724356957944*lo4[6]+0.6123724356957944*lo4[5]+0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]-0.3535533905932737*lo4[2]-0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], 0.7905694150420947*lo5[11]-0.7905694150420948*lo5[10]-0.7905694150420948*lo5[9]-0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]+0.6123724356957944*lo5[6]+0.6123724356957944*lo5[5]+0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]-0.3535533905932737*lo5[2]-0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], 0.7905694150420947*lo6[11]-0.7905694150420948*lo6[10]-0.7905694150420948*lo6[9]-0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]+0.6123724356957944*lo6[6]+0.6123724356957944*lo6[5]+0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]-0.3535533905932737*lo6[2]-0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], 0.7905694150420947*lo7[11]-0.7905694150420948*lo7[10]-0.7905694150420948*lo7[9]-0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]+0.6123724356957944*lo7[6]+0.6123724356957944*lo7[5]+0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]-0.3535533905932737*lo7[2]-0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {0.7905694150420947*up0[11]-0.7905694150420948*up0[10]-0.7905694150420948*up0[9]+0.6123724356957944*up0[8]+0.7905694150420947*up0[7]-0.6123724356957944*up0[6]-0.6123724356957944*up0[5]+0.3535533905932737*up0[4]+0.6123724356957944*up0[3]-0.3535533905932737*up0[2]-0.3535533905932737*up0[1]+0.3535533905932737*up0[0], 0.7905694150420947*up1[11]-0.7905694150420948*up1[10]-0.7905694150420948*up1[9]+0.6123724356957944*up1[8]+0.7905694150420947*up1[7]-0.6123724356957944*up1[6]-0.6123724356957944*up1[5]+0.3535533905932737*up1[4]+0.6123724356957944*up1[3]-0.3535533905932737*up1[2]-0.3535533905932737*up1[1]+0.3535533905932737*up1[0], 0.7905694150420947*up2[11]-0.7905694150420948*up2[10]-0.7905694150420948*up2[9]+0.6123724356957944*up2[8]+0.7905694150420947*up2[7]-0.6123724356957944*up2[6]-0.6123724356957944*up2[5]+0.3535533905932737*up2[4]+0.6123724356957944*up2[3]-0.3535533905932737*up2[2]-0.3535533905932737*up2[1]+0.3535533905932737*up2[0], 0.7905694150420947*up3[11]-0.7905694150420948*up3[10]-0.7905694150420948*up3[9]+0.6123724356957944*up3[8]+0.7905694150420947*up3[7]-0.6123724356957944*up3[6]-0.6123724356957944*up3[5]+0.3535533905932737*up3[4]+0.6123724356957944*up3[3]-0.3535533905932737*up3[2]-0.3535533905932737*up3[1]+0.3535533905932737*up3[0], 0.7905694150420947*up4[11]-0.7905694150420948*up4[10]-0.7905694150420948*up4[9]+0.6123724356957944*up4[8]+0.7905694150420947*up4[7]-0.6123724356957944*up4[6]-0.6123724356957944*up4[5]+0.3535533905932737*up4[4]+0.6123724356957944*up4[3]-0.3535533905932737*up4[2]-0.3535533905932737*up4[1]+0.3535533905932737*up4[0], 0.7905694150420947*up5[11]-0.7905694150420948*up5[10]-0.7905694150420948*up5[9]+0.6123724356957944*up5[8]+0.7905694150420947*up5[7]-0.6123724356957944*up5[6]-0.6123724356957944*up5[5]+0.3535533905932737*up5[4]+0.6123724356957944*up5[3]-0.3535533905932737*up5[2]-0.3535533905932737*up5[1]+0.3535533905932737*up5[0], 0.7905694150420947*up6[11]-0.7905694150420948*up6[10]-0.7905694150420948*up6[9]+0.6123724356957944*up6[8]+0.7905694150420947*up6[7]-0.6123724356957944*up6[6]-0.6123724356957944*up6[5]+0.3535533905932737*up6[4]+0.6123724356957944*up6[3]-0.3535533905932737*up6[2]-0.3535533905932737*up6[1]+0.3535533905932737*up6[0], 0.7905694150420947*up7[11]-0.7905694150420948*up7[10]-0.7905694150420948*up7[9]+0.6123724356957944*up7[8]+0.7905694150420947*up7[7]-0.6123724356957944*up7[6]-0.6123724356957944*up7[5]+0.3535533905932737*up7[4]+0.6123724356957944*up7[3]-0.3535533905932737*up7[2]-0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = 0.5*coef[4]-0.5*coef[3]-0.5*coef[1]+0.5*coef[0]; 
    nv[1] = 0.5*coef[13]-0.5*coef[12]-0.5*coef[10]+0.5*coef[9]; 
    nv[2] = 0.5*coef[202]-0.5*coef[201]-0.5*coef[199]+0.5*coef[198]; 
    nv[3] = 0.5*coef[247]-0.5*coef[246]-0.5*coef[244]+0.5*coef[243]; 
    nv[4] = 0.5*coef[580]-0.5*coef[579]-0.5*coef[577]+0.5*coef[576]; 
    nv[5] = 0.5*coef[625]-0.5*coef[624]-0.5*coef[622]+0.5*coef[621]; 
    nv[6] = 0.5*coef[958]-0.5*coef[957]-0.5*coef[955]+0.5*coef[954]; 
    nv[7] = 0.5*coef[1003]-0.5*coef[1002]-0.5*coef[1000]+0.5*coef[999]; 
    nv[8] = 0.5*coef[211]-0.5*coef[210]-0.5*coef[208]+0.5*coef[207]; 
    nv[9] = 0.5*coef[256]-0.5*coef[255]-0.5*coef[253]+0.5*coef[252]; 
    nv[10] = 0.5*coef[589]-0.5*coef[588]-0.5*coef[586]+0.5*coef[585]; 
    nv[11] = 0.5*coef[634]-0.5*coef[633]-0.5*coef[631]+0.5*coef[630]; 
    nv[12] = 0.5*coef[967]-0.5*coef[966]-0.5*coef[964]+0.5*coef[963]; 
    nv[13] = 0.5*coef[1012]-0.5*coef[1011]-0.5*coef[1009]+0.5*coef[1008]; 
    nv[14] = 0.5*coef[220]-0.5*coef[219]-0.5*coef[217]+0.5*coef[216]; 
    nv[15] = 0.5*coef[265]-0.5*coef[264]-0.5*coef[262]+0.5*coef[261]; 
    nv[16] = 0.5*coef[598]-0.5*coef[597]-0.5*coef[595]+0.5*coef[594]; 
    nv[17] = 0.5*coef[643]-0.5*coef[642]-0.5*coef[640]+0.5*coef[639]; 
    nv[18] = 0.5*coef[976]-0.5*coef[975]-0.5*coef[973]+0.5*coef[972]; 
    nv[19] = 0.5*coef[1021]-0.5*coef[1020]-0.5*coef[1018]+0.5*coef[1017]; 
    nv[20] = 0.5*coef[229]-0.5*coef[228]-0.5*coef[226]+0.5*coef[225]; 
    nv[21] = 0.5*coef[274]-0.5*coef[273]-0.5*coef[271]+0.5*coef[270]; 
    nv[22] = 0.5*coef[607]-0.5*coef[606]-0.5*coef[604]+0.5*coef[603]; 
    nv[23] = 0.5*coef[652]-0.5*coef[651]-0.5*coef[649]+0.5*coef[648]; 
    nv[24] = 0.5*coef[985]-0.5*coef[984]-0.5*coef[982]+0.5*coef[981]; 
    nv[25] = 0.5*coef[1030]-0.5*coef[1029]-0.5*coef[1027]+0.5*coef[1026]; 
    nv[26] = 0.5*coef[238]-0.5*coef[237]-0.5*coef[235]+0.5*coef[234]; 
    nv[27] = 0.5*coef[283]-0.5*coef[282]-0.5*coef[280]+0.5*coef[279]; 
    nv[28] = 0.5*coef[616]-0.5*coef[615]-0.5*coef[613]+0.5*coef[612]; 
    nv[29] = 0.5*coef[661]-0.5*coef[660]-0.5*coef[658]+0.5*coef[657]; 
    nv[30] = 0.5*coef[994]-0.5*coef[993]-0.5*coef[991]+0.5*coef[990]; 
    nv[31] = 0.5*coef[1039]-0.5*coef[1038]-0.5*coef[1036]+0.5*coef[1035]; 
    nv[32] = 0.5*coef[292]-0.5*coef[291]-0.5*coef[289]+0.5*coef[288]; 
    nv[33] = 0.5*coef[319]-0.5*coef[318]-0.5*coef[316]+0.5*coef[315]; 
    nv[34] = 0.5*coef[346]-0.5*coef[345]-0.5*coef[343]+0.5*coef[342]; 
    nv[35] = 0.5*coef[670]-0.5*coef[669]-0.5*coef[667]+0.5*coef[666]; 
    nv[36] = 0.5*coef[697]-0.5*coef[696]-0.5*coef[694]+0.5*coef[693]; 
    nv[37] = 0.5*coef[724]-0.5*coef[723]-0.5*coef[721]+0.5*coef[720]; 
    nv[38] = 0.5*coef[1048]-0.5*coef[1047]-0.5*coef[1045]+0.5*coef[1044]; 
    nv[39] = 0.5*coef[1075]-0.5*coef[1074]-0.5*coef[1072]+0.5*coef[1071]; 
    nv[40] = 0.5*coef[1102]-0.5*coef[1101]-0.5*coef[1099]+0.5*coef[1098]; 
    nv[41] = 0.5*coef[301]-0.5*coef[300]-0.5*coef[298]+0.5*coef[297]; 
    nv[42] = 0.5*coef[328]-0.5*coef[327]-0.5*coef[325]+0.5*coef[324]; 
    nv[43] = 0.5*coef[355]-0.5*coef[354]-0.5*coef[352]+0.5*coef[351]; 
    nv[44] = 0.5*coef[679]-0.5*coef[678]-0.5*coef[676]+0.5*coef[675]; 
    nv[45] = 0.5*coef[706]-0.5*coef[705]-0.5*coef[703]+0.5*coef[702]; 
    nv[46] = 0.5*coef[733]-0.5*coef[732]-0.5*coef[730]+0.5*coef[729]; 
    nv[47] = 0.5*coef[1057]-0.5*coef[1056]-0.5*coef[1054]+0.5*coef[1053]; 
    nv[48] = 0.5*coef[1084]-0.5*coef[1083]-0.5*coef[1081]+0.5*coef[1080]; 
    nv[49] = 0.5*coef[1111]-0.5*coef[1110]-0.5*coef[1108]+0.5*coef[1107]; 
    nv[50] = 0.5*coef[310]-0.5*coef[309]-0.5*coef[307]+0.5*coef[306]; 
    nv[51] = 0.5*coef[337]-0.5*coef[336]-0.5*coef[334]+0.5*coef[333]; 
    nv[52] = 0.5*coef[364]-0.5*coef[363]-0.5*coef[361]+0.5*coef[360]; 
    nv[53] = 0.5*coef[688]-0.5*coef[687]-0.5*coef[685]+0.5*coef[684]; 
    nv[54] = 0.5*coef[715]-0.5*coef[714]-0.5*coef[712]+0.5*coef[711]; 
    nv[55] = 0.5*coef[742]-0.5*coef[741]-0.5*coef[739]+0.5*coef[738]; 
    nv[56] = 0.5*coef[1066]-0.5*coef[1065]-0.5*coef[1063]+0.5*coef[1062]; 
    nv[57] = 0.5*coef[1093]-0.5*coef[1092]-0.5*coef[1090]+0.5*coef[1089]; 
    nv[58] = 0.5*coef[1120]-0.5*coef[1119]-0.5*coef[1117]+0.5*coef[1116]; 
    nv[59] = 0.5*coef[40]-0.5*coef[39]-0.5*coef[37]+0.5*coef[36]; 
    nv[60] = 0.5*coef[49]-0.5*coef[48]-0.5*coef[46]+0.5*coef[45]; 
    nv[61] = 0.5*coef[58]-0.5*coef[57]-0.5*coef[55]+0.5*coef[54]; 
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
    const double th_q = -(0.5*theta[3])+0.5*theta[2]-0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = -(0.5*b0[3])+0.5*b0[2]-0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = -(0.5*b1[3])+0.5*b1[2]-0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = -(0.5*b2[3])+0.5*b2[2]-0.5*b2[1]+0.5*b2[0]; 
    const double Ulo_q[8] = {-(0.7905694150420947*lo0[11])+0.7905694150420948*lo0[10]-0.7905694150420948*lo0[9]+0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]-0.6123724356957944*lo0[6]+0.6123724356957944*lo0[5]-0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]+0.3535533905932737*lo0[2]-0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], -(0.7905694150420947*lo1[11])+0.7905694150420948*lo1[10]-0.7905694150420948*lo1[9]+0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]-0.6123724356957944*lo1[6]+0.6123724356957944*lo1[5]-0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]+0.3535533905932737*lo1[2]-0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], -(0.7905694150420947*lo2[11])+0.7905694150420948*lo2[10]-0.7905694150420948*lo2[9]+0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]-0.6123724356957944*lo2[6]+0.6123724356957944*lo2[5]-0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]+0.3535533905932737*lo2[2]-0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], -(0.7905694150420947*lo3[11])+0.7905694150420948*lo3[10]-0.7905694150420948*lo3[9]+0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]-0.6123724356957944*lo3[6]+0.6123724356957944*lo3[5]-0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]+0.3535533905932737*lo3[2]-0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], -(0.7905694150420947*lo4[11])+0.7905694150420948*lo4[10]-0.7905694150420948*lo4[9]+0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]-0.6123724356957944*lo4[6]+0.6123724356957944*lo4[5]-0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]+0.3535533905932737*lo4[2]-0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], -(0.7905694150420947*lo5[11])+0.7905694150420948*lo5[10]-0.7905694150420948*lo5[9]+0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]-0.6123724356957944*lo5[6]+0.6123724356957944*lo5[5]-0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]+0.3535533905932737*lo5[2]-0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], -(0.7905694150420947*lo6[11])+0.7905694150420948*lo6[10]-0.7905694150420948*lo6[9]+0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]-0.6123724356957944*lo6[6]+0.6123724356957944*lo6[5]-0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]+0.3535533905932737*lo6[2]-0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], -(0.7905694150420947*lo7[11])+0.7905694150420948*lo7[10]-0.7905694150420948*lo7[9]+0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]-0.6123724356957944*lo7[6]+0.6123724356957944*lo7[5]-0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]+0.3535533905932737*lo7[2]-0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {-(0.7905694150420947*up0[11])+0.7905694150420948*up0[10]-0.7905694150420948*up0[9]-0.6123724356957944*up0[8]+0.7905694150420947*up0[7]+0.6123724356957944*up0[6]-0.6123724356957944*up0[5]-0.3535533905932737*up0[4]+0.6123724356957944*up0[3]+0.3535533905932737*up0[2]-0.3535533905932737*up0[1]+0.3535533905932737*up0[0], -(0.7905694150420947*up1[11])+0.7905694150420948*up1[10]-0.7905694150420948*up1[9]-0.6123724356957944*up1[8]+0.7905694150420947*up1[7]+0.6123724356957944*up1[6]-0.6123724356957944*up1[5]-0.3535533905932737*up1[4]+0.6123724356957944*up1[3]+0.3535533905932737*up1[2]-0.3535533905932737*up1[1]+0.3535533905932737*up1[0], -(0.7905694150420947*up2[11])+0.7905694150420948*up2[10]-0.7905694150420948*up2[9]-0.6123724356957944*up2[8]+0.7905694150420947*up2[7]+0.6123724356957944*up2[6]-0.6123724356957944*up2[5]-0.3535533905932737*up2[4]+0.6123724356957944*up2[3]+0.3535533905932737*up2[2]-0.3535533905932737*up2[1]+0.3535533905932737*up2[0], -(0.7905694150420947*up3[11])+0.7905694150420948*up3[10]-0.7905694150420948*up3[9]-0.6123724356957944*up3[8]+0.7905694150420947*up3[7]+0.6123724356957944*up3[6]-0.6123724356957944*up3[5]-0.3535533905932737*up3[4]+0.6123724356957944*up3[3]+0.3535533905932737*up3[2]-0.3535533905932737*up3[1]+0.3535533905932737*up3[0], -(0.7905694150420947*up4[11])+0.7905694150420948*up4[10]-0.7905694150420948*up4[9]-0.6123724356957944*up4[8]+0.7905694150420947*up4[7]+0.6123724356957944*up4[6]-0.6123724356957944*up4[5]-0.3535533905932737*up4[4]+0.6123724356957944*up4[3]+0.3535533905932737*up4[2]-0.3535533905932737*up4[1]+0.3535533905932737*up4[0], -(0.7905694150420947*up5[11])+0.7905694150420948*up5[10]-0.7905694150420948*up5[9]-0.6123724356957944*up5[8]+0.7905694150420947*up5[7]+0.6123724356957944*up5[6]-0.6123724356957944*up5[5]-0.3535533905932737*up5[4]+0.6123724356957944*up5[3]+0.3535533905932737*up5[2]-0.3535533905932737*up5[1]+0.3535533905932737*up5[0], -(0.7905694150420947*up6[11])+0.7905694150420948*up6[10]-0.7905694150420948*up6[9]-0.6123724356957944*up6[8]+0.7905694150420947*up6[7]+0.6123724356957944*up6[6]-0.6123724356957944*up6[5]-0.3535533905932737*up6[4]+0.6123724356957944*up6[3]+0.3535533905932737*up6[2]-0.3535533905932737*up6[1]+0.3535533905932737*up6[0], -(0.7905694150420947*up7[11])+0.7905694150420948*up7[10]-0.7905694150420948*up7[9]-0.6123724356957944*up7[8]+0.7905694150420947*up7[7]+0.6123724356957944*up7[6]-0.6123724356957944*up7[5]-0.3535533905932737*up7[4]+0.6123724356957944*up7[3]+0.3535533905932737*up7[2]-0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = -(0.5*coef[4])+0.5*coef[3]-0.5*coef[1]+0.5*coef[0]; 
    nv[1] = -(0.5*coef[13])+0.5*coef[12]-0.5*coef[10]+0.5*coef[9]; 
    nv[2] = -(0.5*coef[202])+0.5*coef[201]-0.5*coef[199]+0.5*coef[198]; 
    nv[3] = -(0.5*coef[247])+0.5*coef[246]-0.5*coef[244]+0.5*coef[243]; 
    nv[4] = -(0.5*coef[580])+0.5*coef[579]-0.5*coef[577]+0.5*coef[576]; 
    nv[5] = -(0.5*coef[625])+0.5*coef[624]-0.5*coef[622]+0.5*coef[621]; 
    nv[6] = -(0.5*coef[958])+0.5*coef[957]-0.5*coef[955]+0.5*coef[954]; 
    nv[7] = -(0.5*coef[1003])+0.5*coef[1002]-0.5*coef[1000]+0.5*coef[999]; 
    nv[8] = -(0.5*coef[211])+0.5*coef[210]-0.5*coef[208]+0.5*coef[207]; 
    nv[9] = -(0.5*coef[256])+0.5*coef[255]-0.5*coef[253]+0.5*coef[252]; 
    nv[10] = -(0.5*coef[589])+0.5*coef[588]-0.5*coef[586]+0.5*coef[585]; 
    nv[11] = -(0.5*coef[634])+0.5*coef[633]-0.5*coef[631]+0.5*coef[630]; 
    nv[12] = -(0.5*coef[967])+0.5*coef[966]-0.5*coef[964]+0.5*coef[963]; 
    nv[13] = -(0.5*coef[1012])+0.5*coef[1011]-0.5*coef[1009]+0.5*coef[1008]; 
    nv[14] = -(0.5*coef[220])+0.5*coef[219]-0.5*coef[217]+0.5*coef[216]; 
    nv[15] = -(0.5*coef[265])+0.5*coef[264]-0.5*coef[262]+0.5*coef[261]; 
    nv[16] = -(0.5*coef[598])+0.5*coef[597]-0.5*coef[595]+0.5*coef[594]; 
    nv[17] = -(0.5*coef[643])+0.5*coef[642]-0.5*coef[640]+0.5*coef[639]; 
    nv[18] = -(0.5*coef[976])+0.5*coef[975]-0.5*coef[973]+0.5*coef[972]; 
    nv[19] = -(0.5*coef[1021])+0.5*coef[1020]-0.5*coef[1018]+0.5*coef[1017]; 
    nv[20] = -(0.5*coef[229])+0.5*coef[228]-0.5*coef[226]+0.5*coef[225]; 
    nv[21] = -(0.5*coef[274])+0.5*coef[273]-0.5*coef[271]+0.5*coef[270]; 
    nv[22] = -(0.5*coef[607])+0.5*coef[606]-0.5*coef[604]+0.5*coef[603]; 
    nv[23] = -(0.5*coef[652])+0.5*coef[651]-0.5*coef[649]+0.5*coef[648]; 
    nv[24] = -(0.5*coef[985])+0.5*coef[984]-0.5*coef[982]+0.5*coef[981]; 
    nv[25] = -(0.5*coef[1030])+0.5*coef[1029]-0.5*coef[1027]+0.5*coef[1026]; 
    nv[26] = -(0.5*coef[238])+0.5*coef[237]-0.5*coef[235]+0.5*coef[234]; 
    nv[27] = -(0.5*coef[283])+0.5*coef[282]-0.5*coef[280]+0.5*coef[279]; 
    nv[28] = -(0.5*coef[616])+0.5*coef[615]-0.5*coef[613]+0.5*coef[612]; 
    nv[29] = -(0.5*coef[661])+0.5*coef[660]-0.5*coef[658]+0.5*coef[657]; 
    nv[30] = -(0.5*coef[994])+0.5*coef[993]-0.5*coef[991]+0.5*coef[990]; 
    nv[31] = -(0.5*coef[1039])+0.5*coef[1038]-0.5*coef[1036]+0.5*coef[1035]; 
    nv[32] = -(0.5*coef[292])+0.5*coef[291]-0.5*coef[289]+0.5*coef[288]; 
    nv[33] = -(0.5*coef[319])+0.5*coef[318]-0.5*coef[316]+0.5*coef[315]; 
    nv[34] = -(0.5*coef[346])+0.5*coef[345]-0.5*coef[343]+0.5*coef[342]; 
    nv[35] = -(0.5*coef[670])+0.5*coef[669]-0.5*coef[667]+0.5*coef[666]; 
    nv[36] = -(0.5*coef[697])+0.5*coef[696]-0.5*coef[694]+0.5*coef[693]; 
    nv[37] = -(0.5*coef[724])+0.5*coef[723]-0.5*coef[721]+0.5*coef[720]; 
    nv[38] = -(0.5*coef[1048])+0.5*coef[1047]-0.5*coef[1045]+0.5*coef[1044]; 
    nv[39] = -(0.5*coef[1075])+0.5*coef[1074]-0.5*coef[1072]+0.5*coef[1071]; 
    nv[40] = -(0.5*coef[1102])+0.5*coef[1101]-0.5*coef[1099]+0.5*coef[1098]; 
    nv[41] = -(0.5*coef[301])+0.5*coef[300]-0.5*coef[298]+0.5*coef[297]; 
    nv[42] = -(0.5*coef[328])+0.5*coef[327]-0.5*coef[325]+0.5*coef[324]; 
    nv[43] = -(0.5*coef[355])+0.5*coef[354]-0.5*coef[352]+0.5*coef[351]; 
    nv[44] = -(0.5*coef[679])+0.5*coef[678]-0.5*coef[676]+0.5*coef[675]; 
    nv[45] = -(0.5*coef[706])+0.5*coef[705]-0.5*coef[703]+0.5*coef[702]; 
    nv[46] = -(0.5*coef[733])+0.5*coef[732]-0.5*coef[730]+0.5*coef[729]; 
    nv[47] = -(0.5*coef[1057])+0.5*coef[1056]-0.5*coef[1054]+0.5*coef[1053]; 
    nv[48] = -(0.5*coef[1084])+0.5*coef[1083]-0.5*coef[1081]+0.5*coef[1080]; 
    nv[49] = -(0.5*coef[1111])+0.5*coef[1110]-0.5*coef[1108]+0.5*coef[1107]; 
    nv[50] = -(0.5*coef[310])+0.5*coef[309]-0.5*coef[307]+0.5*coef[306]; 
    nv[51] = -(0.5*coef[337])+0.5*coef[336]-0.5*coef[334]+0.5*coef[333]; 
    nv[52] = -(0.5*coef[364])+0.5*coef[363]-0.5*coef[361]+0.5*coef[360]; 
    nv[53] = -(0.5*coef[688])+0.5*coef[687]-0.5*coef[685]+0.5*coef[684]; 
    nv[54] = -(0.5*coef[715])+0.5*coef[714]-0.5*coef[712]+0.5*coef[711]; 
    nv[55] = -(0.5*coef[742])+0.5*coef[741]-0.5*coef[739]+0.5*coef[738]; 
    nv[56] = -(0.5*coef[1066])+0.5*coef[1065]-0.5*coef[1063]+0.5*coef[1062]; 
    nv[57] = -(0.5*coef[1093])+0.5*coef[1092]-0.5*coef[1090]+0.5*coef[1089]; 
    nv[58] = -(0.5*coef[1120])+0.5*coef[1119]-0.5*coef[1117]+0.5*coef[1116]; 
    nv[59] = -(0.5*coef[40])+0.5*coef[39]-0.5*coef[37]+0.5*coef[36]; 
    nv[60] = -(0.5*coef[49])+0.5*coef[48]-0.5*coef[46]+0.5*coef[45]; 
    nv[61] = -(0.5*coef[58])+0.5*coef[57]-0.5*coef[55]+0.5*coef[54]; 
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
    const double th_q = -(0.5*theta[3])-0.5*theta[2]+0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = -(0.5*b0[3])-0.5*b0[2]+0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = -(0.5*b1[3])-0.5*b1[2]+0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = -(0.5*b2[3])-0.5*b2[2]+0.5*b2[1]+0.5*b2[0]; 
    const double Ulo_q[8] = {-(0.7905694150420947*lo0[11])-0.7905694150420948*lo0[10]+0.7905694150420948*lo0[9]+0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]+0.6123724356957944*lo0[6]-0.6123724356957944*lo0[5]-0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]-0.3535533905932737*lo0[2]+0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], -(0.7905694150420947*lo1[11])-0.7905694150420948*lo1[10]+0.7905694150420948*lo1[9]+0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]+0.6123724356957944*lo1[6]-0.6123724356957944*lo1[5]-0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]-0.3535533905932737*lo1[2]+0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], -(0.7905694150420947*lo2[11])-0.7905694150420948*lo2[10]+0.7905694150420948*lo2[9]+0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]+0.6123724356957944*lo2[6]-0.6123724356957944*lo2[5]-0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]-0.3535533905932737*lo2[2]+0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], -(0.7905694150420947*lo3[11])-0.7905694150420948*lo3[10]+0.7905694150420948*lo3[9]+0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]+0.6123724356957944*lo3[6]-0.6123724356957944*lo3[5]-0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]-0.3535533905932737*lo3[2]+0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], -(0.7905694150420947*lo4[11])-0.7905694150420948*lo4[10]+0.7905694150420948*lo4[9]+0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]+0.6123724356957944*lo4[6]-0.6123724356957944*lo4[5]-0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]-0.3535533905932737*lo4[2]+0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], -(0.7905694150420947*lo5[11])-0.7905694150420948*lo5[10]+0.7905694150420948*lo5[9]+0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]+0.6123724356957944*lo5[6]-0.6123724356957944*lo5[5]-0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]-0.3535533905932737*lo5[2]+0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], -(0.7905694150420947*lo6[11])-0.7905694150420948*lo6[10]+0.7905694150420948*lo6[9]+0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]+0.6123724356957944*lo6[6]-0.6123724356957944*lo6[5]-0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]-0.3535533905932737*lo6[2]+0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], -(0.7905694150420947*lo7[11])-0.7905694150420948*lo7[10]+0.7905694150420948*lo7[9]+0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]+0.6123724356957944*lo7[6]-0.6123724356957944*lo7[5]-0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]-0.3535533905932737*lo7[2]+0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {-(0.7905694150420947*up0[11])-0.7905694150420948*up0[10]+0.7905694150420948*up0[9]-0.6123724356957944*up0[8]+0.7905694150420947*up0[7]-0.6123724356957944*up0[6]+0.6123724356957944*up0[5]-0.3535533905932737*up0[4]+0.6123724356957944*up0[3]-0.3535533905932737*up0[2]+0.3535533905932737*up0[1]+0.3535533905932737*up0[0], -(0.7905694150420947*up1[11])-0.7905694150420948*up1[10]+0.7905694150420948*up1[9]-0.6123724356957944*up1[8]+0.7905694150420947*up1[7]-0.6123724356957944*up1[6]+0.6123724356957944*up1[5]-0.3535533905932737*up1[4]+0.6123724356957944*up1[3]-0.3535533905932737*up1[2]+0.3535533905932737*up1[1]+0.3535533905932737*up1[0], -(0.7905694150420947*up2[11])-0.7905694150420948*up2[10]+0.7905694150420948*up2[9]-0.6123724356957944*up2[8]+0.7905694150420947*up2[7]-0.6123724356957944*up2[6]+0.6123724356957944*up2[5]-0.3535533905932737*up2[4]+0.6123724356957944*up2[3]-0.3535533905932737*up2[2]+0.3535533905932737*up2[1]+0.3535533905932737*up2[0], -(0.7905694150420947*up3[11])-0.7905694150420948*up3[10]+0.7905694150420948*up3[9]-0.6123724356957944*up3[8]+0.7905694150420947*up3[7]-0.6123724356957944*up3[6]+0.6123724356957944*up3[5]-0.3535533905932737*up3[4]+0.6123724356957944*up3[3]-0.3535533905932737*up3[2]+0.3535533905932737*up3[1]+0.3535533905932737*up3[0], -(0.7905694150420947*up4[11])-0.7905694150420948*up4[10]+0.7905694150420948*up4[9]-0.6123724356957944*up4[8]+0.7905694150420947*up4[7]-0.6123724356957944*up4[6]+0.6123724356957944*up4[5]-0.3535533905932737*up4[4]+0.6123724356957944*up4[3]-0.3535533905932737*up4[2]+0.3535533905932737*up4[1]+0.3535533905932737*up4[0], -(0.7905694150420947*up5[11])-0.7905694150420948*up5[10]+0.7905694150420948*up5[9]-0.6123724356957944*up5[8]+0.7905694150420947*up5[7]-0.6123724356957944*up5[6]+0.6123724356957944*up5[5]-0.3535533905932737*up5[4]+0.6123724356957944*up5[3]-0.3535533905932737*up5[2]+0.3535533905932737*up5[1]+0.3535533905932737*up5[0], -(0.7905694150420947*up6[11])-0.7905694150420948*up6[10]+0.7905694150420948*up6[9]-0.6123724356957944*up6[8]+0.7905694150420947*up6[7]-0.6123724356957944*up6[6]+0.6123724356957944*up6[5]-0.3535533905932737*up6[4]+0.6123724356957944*up6[3]-0.3535533905932737*up6[2]+0.3535533905932737*up6[1]+0.3535533905932737*up6[0], -(0.7905694150420947*up7[11])-0.7905694150420948*up7[10]+0.7905694150420948*up7[9]-0.6123724356957944*up7[8]+0.7905694150420947*up7[7]-0.6123724356957944*up7[6]+0.6123724356957944*up7[5]-0.3535533905932737*up7[4]+0.6123724356957944*up7[3]-0.3535533905932737*up7[2]+0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = -(0.5*coef[4])-0.5*coef[3]+0.5*coef[1]+0.5*coef[0]; 
    nv[1] = -(0.5*coef[13])-0.5*coef[12]+0.5*coef[10]+0.5*coef[9]; 
    nv[2] = -(0.5*coef[202])-0.5*coef[201]+0.5*coef[199]+0.5*coef[198]; 
    nv[3] = -(0.5*coef[247])-0.5*coef[246]+0.5*coef[244]+0.5*coef[243]; 
    nv[4] = -(0.5*coef[580])-0.5*coef[579]+0.5*coef[577]+0.5*coef[576]; 
    nv[5] = -(0.5*coef[625])-0.5*coef[624]+0.5*coef[622]+0.5*coef[621]; 
    nv[6] = -(0.5*coef[958])-0.5*coef[957]+0.5*coef[955]+0.5*coef[954]; 
    nv[7] = -(0.5*coef[1003])-0.5*coef[1002]+0.5*coef[1000]+0.5*coef[999]; 
    nv[8] = -(0.5*coef[211])-0.5*coef[210]+0.5*coef[208]+0.5*coef[207]; 
    nv[9] = -(0.5*coef[256])-0.5*coef[255]+0.5*coef[253]+0.5*coef[252]; 
    nv[10] = -(0.5*coef[589])-0.5*coef[588]+0.5*coef[586]+0.5*coef[585]; 
    nv[11] = -(0.5*coef[634])-0.5*coef[633]+0.5*coef[631]+0.5*coef[630]; 
    nv[12] = -(0.5*coef[967])-0.5*coef[966]+0.5*coef[964]+0.5*coef[963]; 
    nv[13] = -(0.5*coef[1012])-0.5*coef[1011]+0.5*coef[1009]+0.5*coef[1008]; 
    nv[14] = -(0.5*coef[220])-0.5*coef[219]+0.5*coef[217]+0.5*coef[216]; 
    nv[15] = -(0.5*coef[265])-0.5*coef[264]+0.5*coef[262]+0.5*coef[261]; 
    nv[16] = -(0.5*coef[598])-0.5*coef[597]+0.5*coef[595]+0.5*coef[594]; 
    nv[17] = -(0.5*coef[643])-0.5*coef[642]+0.5*coef[640]+0.5*coef[639]; 
    nv[18] = -(0.5*coef[976])-0.5*coef[975]+0.5*coef[973]+0.5*coef[972]; 
    nv[19] = -(0.5*coef[1021])-0.5*coef[1020]+0.5*coef[1018]+0.5*coef[1017]; 
    nv[20] = -(0.5*coef[229])-0.5*coef[228]+0.5*coef[226]+0.5*coef[225]; 
    nv[21] = -(0.5*coef[274])-0.5*coef[273]+0.5*coef[271]+0.5*coef[270]; 
    nv[22] = -(0.5*coef[607])-0.5*coef[606]+0.5*coef[604]+0.5*coef[603]; 
    nv[23] = -(0.5*coef[652])-0.5*coef[651]+0.5*coef[649]+0.5*coef[648]; 
    nv[24] = -(0.5*coef[985])-0.5*coef[984]+0.5*coef[982]+0.5*coef[981]; 
    nv[25] = -(0.5*coef[1030])-0.5*coef[1029]+0.5*coef[1027]+0.5*coef[1026]; 
    nv[26] = -(0.5*coef[238])-0.5*coef[237]+0.5*coef[235]+0.5*coef[234]; 
    nv[27] = -(0.5*coef[283])-0.5*coef[282]+0.5*coef[280]+0.5*coef[279]; 
    nv[28] = -(0.5*coef[616])-0.5*coef[615]+0.5*coef[613]+0.5*coef[612]; 
    nv[29] = -(0.5*coef[661])-0.5*coef[660]+0.5*coef[658]+0.5*coef[657]; 
    nv[30] = -(0.5*coef[994])-0.5*coef[993]+0.5*coef[991]+0.5*coef[990]; 
    nv[31] = -(0.5*coef[1039])-0.5*coef[1038]+0.5*coef[1036]+0.5*coef[1035]; 
    nv[32] = -(0.5*coef[292])-0.5*coef[291]+0.5*coef[289]+0.5*coef[288]; 
    nv[33] = -(0.5*coef[319])-0.5*coef[318]+0.5*coef[316]+0.5*coef[315]; 
    nv[34] = -(0.5*coef[346])-0.5*coef[345]+0.5*coef[343]+0.5*coef[342]; 
    nv[35] = -(0.5*coef[670])-0.5*coef[669]+0.5*coef[667]+0.5*coef[666]; 
    nv[36] = -(0.5*coef[697])-0.5*coef[696]+0.5*coef[694]+0.5*coef[693]; 
    nv[37] = -(0.5*coef[724])-0.5*coef[723]+0.5*coef[721]+0.5*coef[720]; 
    nv[38] = -(0.5*coef[1048])-0.5*coef[1047]+0.5*coef[1045]+0.5*coef[1044]; 
    nv[39] = -(0.5*coef[1075])-0.5*coef[1074]+0.5*coef[1072]+0.5*coef[1071]; 
    nv[40] = -(0.5*coef[1102])-0.5*coef[1101]+0.5*coef[1099]+0.5*coef[1098]; 
    nv[41] = -(0.5*coef[301])-0.5*coef[300]+0.5*coef[298]+0.5*coef[297]; 
    nv[42] = -(0.5*coef[328])-0.5*coef[327]+0.5*coef[325]+0.5*coef[324]; 
    nv[43] = -(0.5*coef[355])-0.5*coef[354]+0.5*coef[352]+0.5*coef[351]; 
    nv[44] = -(0.5*coef[679])-0.5*coef[678]+0.5*coef[676]+0.5*coef[675]; 
    nv[45] = -(0.5*coef[706])-0.5*coef[705]+0.5*coef[703]+0.5*coef[702]; 
    nv[46] = -(0.5*coef[733])-0.5*coef[732]+0.5*coef[730]+0.5*coef[729]; 
    nv[47] = -(0.5*coef[1057])-0.5*coef[1056]+0.5*coef[1054]+0.5*coef[1053]; 
    nv[48] = -(0.5*coef[1084])-0.5*coef[1083]+0.5*coef[1081]+0.5*coef[1080]; 
    nv[49] = -(0.5*coef[1111])-0.5*coef[1110]+0.5*coef[1108]+0.5*coef[1107]; 
    nv[50] = -(0.5*coef[310])-0.5*coef[309]+0.5*coef[307]+0.5*coef[306]; 
    nv[51] = -(0.5*coef[337])-0.5*coef[336]+0.5*coef[334]+0.5*coef[333]; 
    nv[52] = -(0.5*coef[364])-0.5*coef[363]+0.5*coef[361]+0.5*coef[360]; 
    nv[53] = -(0.5*coef[688])-0.5*coef[687]+0.5*coef[685]+0.5*coef[684]; 
    nv[54] = -(0.5*coef[715])-0.5*coef[714]+0.5*coef[712]+0.5*coef[711]; 
    nv[55] = -(0.5*coef[742])-0.5*coef[741]+0.5*coef[739]+0.5*coef[738]; 
    nv[56] = -(0.5*coef[1066])-0.5*coef[1065]+0.5*coef[1063]+0.5*coef[1062]; 
    nv[57] = -(0.5*coef[1093])-0.5*coef[1092]+0.5*coef[1090]+0.5*coef[1089]; 
    nv[58] = -(0.5*coef[1120])-0.5*coef[1119]+0.5*coef[1117]+0.5*coef[1116]; 
    nv[59] = -(0.5*coef[40])-0.5*coef[39]+0.5*coef[37]+0.5*coef[36]; 
    nv[60] = -(0.5*coef[49])-0.5*coef[48]+0.5*coef[46]+0.5*coef[45]; 
    nv[61] = -(0.5*coef[58])-0.5*coef[57]+0.5*coef[55]+0.5*coef[54]; 
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
    const double th_q = 0.5*theta[3]+0.5*theta[2]+0.5*theta[1]+0.5*theta[0]; 
    const double b0_q = 0.5*b0[3]+0.5*b0[2]+0.5*b0[1]+0.5*b0[0]; 
    const double b1_q = 0.5*b1[3]+0.5*b1[2]+0.5*b1[1]+0.5*b1[0]; 
    const double b2_q = 0.5*b2[3]+0.5*b2[2]+0.5*b2[1]+0.5*b2[0]; 
    const double Ulo_q[8] = {0.7905694150420947*lo0[11]+0.7905694150420948*lo0[10]+0.7905694150420948*lo0[9]-0.6123724356957944*lo0[8]+0.7905694150420947*lo0[7]-0.6123724356957944*lo0[6]-0.6123724356957944*lo0[5]+0.3535533905932737*lo0[4]-0.6123724356957944*lo0[3]+0.3535533905932737*lo0[2]+0.3535533905932737*lo0[1]+0.3535533905932737*lo0[0], 0.7905694150420947*lo1[11]+0.7905694150420948*lo1[10]+0.7905694150420948*lo1[9]-0.6123724356957944*lo1[8]+0.7905694150420947*lo1[7]-0.6123724356957944*lo1[6]-0.6123724356957944*lo1[5]+0.3535533905932737*lo1[4]-0.6123724356957944*lo1[3]+0.3535533905932737*lo1[2]+0.3535533905932737*lo1[1]+0.3535533905932737*lo1[0], 0.7905694150420947*lo2[11]+0.7905694150420948*lo2[10]+0.7905694150420948*lo2[9]-0.6123724356957944*lo2[8]+0.7905694150420947*lo2[7]-0.6123724356957944*lo2[6]-0.6123724356957944*lo2[5]+0.3535533905932737*lo2[4]-0.6123724356957944*lo2[3]+0.3535533905932737*lo2[2]+0.3535533905932737*lo2[1]+0.3535533905932737*lo2[0], 0.7905694150420947*lo3[11]+0.7905694150420948*lo3[10]+0.7905694150420948*lo3[9]-0.6123724356957944*lo3[8]+0.7905694150420947*lo3[7]-0.6123724356957944*lo3[6]-0.6123724356957944*lo3[5]+0.3535533905932737*lo3[4]-0.6123724356957944*lo3[3]+0.3535533905932737*lo3[2]+0.3535533905932737*lo3[1]+0.3535533905932737*lo3[0], 0.7905694150420947*lo4[11]+0.7905694150420948*lo4[10]+0.7905694150420948*lo4[9]-0.6123724356957944*lo4[8]+0.7905694150420947*lo4[7]-0.6123724356957944*lo4[6]-0.6123724356957944*lo4[5]+0.3535533905932737*lo4[4]-0.6123724356957944*lo4[3]+0.3535533905932737*lo4[2]+0.3535533905932737*lo4[1]+0.3535533905932737*lo4[0], 0.7905694150420947*lo5[11]+0.7905694150420948*lo5[10]+0.7905694150420948*lo5[9]-0.6123724356957944*lo5[8]+0.7905694150420947*lo5[7]-0.6123724356957944*lo5[6]-0.6123724356957944*lo5[5]+0.3535533905932737*lo5[4]-0.6123724356957944*lo5[3]+0.3535533905932737*lo5[2]+0.3535533905932737*lo5[1]+0.3535533905932737*lo5[0], 0.7905694150420947*lo6[11]+0.7905694150420948*lo6[10]+0.7905694150420948*lo6[9]-0.6123724356957944*lo6[8]+0.7905694150420947*lo6[7]-0.6123724356957944*lo6[6]-0.6123724356957944*lo6[5]+0.3535533905932737*lo6[4]-0.6123724356957944*lo6[3]+0.3535533905932737*lo6[2]+0.3535533905932737*lo6[1]+0.3535533905932737*lo6[0], 0.7905694150420947*lo7[11]+0.7905694150420948*lo7[10]+0.7905694150420948*lo7[9]-0.6123724356957944*lo7[8]+0.7905694150420947*lo7[7]-0.6123724356957944*lo7[6]-0.6123724356957944*lo7[5]+0.3535533905932737*lo7[4]-0.6123724356957944*lo7[3]+0.3535533905932737*lo7[2]+0.3535533905932737*lo7[1]+0.3535533905932737*lo7[0]}; 
    const double Uup_q[8] = {0.7905694150420947*up0[11]+0.7905694150420948*up0[10]+0.7905694150420948*up0[9]+0.6123724356957944*up0[8]+0.7905694150420947*up0[7]+0.6123724356957944*up0[6]+0.6123724356957944*up0[5]+0.3535533905932737*up0[4]+0.6123724356957944*up0[3]+0.3535533905932737*up0[2]+0.3535533905932737*up0[1]+0.3535533905932737*up0[0], 0.7905694150420947*up1[11]+0.7905694150420948*up1[10]+0.7905694150420948*up1[9]+0.6123724356957944*up1[8]+0.7905694150420947*up1[7]+0.6123724356957944*up1[6]+0.6123724356957944*up1[5]+0.3535533905932737*up1[4]+0.6123724356957944*up1[3]+0.3535533905932737*up1[2]+0.3535533905932737*up1[1]+0.3535533905932737*up1[0], 0.7905694150420947*up2[11]+0.7905694150420948*up2[10]+0.7905694150420948*up2[9]+0.6123724356957944*up2[8]+0.7905694150420947*up2[7]+0.6123724356957944*up2[6]+0.6123724356957944*up2[5]+0.3535533905932737*up2[4]+0.6123724356957944*up2[3]+0.3535533905932737*up2[2]+0.3535533905932737*up2[1]+0.3535533905932737*up2[0], 0.7905694150420947*up3[11]+0.7905694150420948*up3[10]+0.7905694150420948*up3[9]+0.6123724356957944*up3[8]+0.7905694150420947*up3[7]+0.6123724356957944*up3[6]+0.6123724356957944*up3[5]+0.3535533905932737*up3[4]+0.6123724356957944*up3[3]+0.3535533905932737*up3[2]+0.3535533905932737*up3[1]+0.3535533905932737*up3[0], 0.7905694150420947*up4[11]+0.7905694150420948*up4[10]+0.7905694150420948*up4[9]+0.6123724356957944*up4[8]+0.7905694150420947*up4[7]+0.6123724356957944*up4[6]+0.6123724356957944*up4[5]+0.3535533905932737*up4[4]+0.6123724356957944*up4[3]+0.3535533905932737*up4[2]+0.3535533905932737*up4[1]+0.3535533905932737*up4[0], 0.7905694150420947*up5[11]+0.7905694150420948*up5[10]+0.7905694150420948*up5[9]+0.6123724356957944*up5[8]+0.7905694150420947*up5[7]+0.6123724356957944*up5[6]+0.6123724356957944*up5[5]+0.3535533905932737*up5[4]+0.6123724356957944*up5[3]+0.3535533905932737*up5[2]+0.3535533905932737*up5[1]+0.3535533905932737*up5[0], 0.7905694150420947*up6[11]+0.7905694150420948*up6[10]+0.7905694150420948*up6[9]+0.6123724356957944*up6[8]+0.7905694150420947*up6[7]+0.6123724356957944*up6[6]+0.6123724356957944*up6[5]+0.3535533905932737*up6[4]+0.6123724356957944*up6[3]+0.3535533905932737*up6[2]+0.3535533905932737*up6[1]+0.3535533905932737*up6[0], 0.7905694150420947*up7[11]+0.7905694150420948*up7[10]+0.7905694150420948*up7[9]+0.6123724356957944*up7[8]+0.7905694150420947*up7[7]+0.6123724356957944*up7[6]+0.6123724356957944*up7[5]+0.3535533905932737*up7[4]+0.6123724356957944*up7[3]+0.3535533905932737*up7[2]+0.3535533905932737*up7[1]+0.3535533905932737*up7[0]}; 
    const double Flo_q = Ulo_q[0], Glo_q = Ulo_q[1], Fup_q = Uup_q[0], Gup_q = Uup_q[1]; 
    const double b_q[3] = {b0_q, b1_q, b2_q}; 
    // node values of the w flux (the stored coefficients at the node, theta, b and the previous 
    // closures): coupling flux T, its speed bounds and the known part of alpha at the edges 
    double nv[PKPM_HS_NODE_NV]; 
    nv[0] = 0.5*coef[4]+0.5*coef[3]+0.5*coef[1]+0.5*coef[0]; 
    nv[1] = 0.5*coef[13]+0.5*coef[12]+0.5*coef[10]+0.5*coef[9]; 
    nv[2] = 0.5*coef[202]+0.5*coef[201]+0.5*coef[199]+0.5*coef[198]; 
    nv[3] = 0.5*coef[247]+0.5*coef[246]+0.5*coef[244]+0.5*coef[243]; 
    nv[4] = 0.5*coef[580]+0.5*coef[579]+0.5*coef[577]+0.5*coef[576]; 
    nv[5] = 0.5*coef[625]+0.5*coef[624]+0.5*coef[622]+0.5*coef[621]; 
    nv[6] = 0.5*coef[958]+0.5*coef[957]+0.5*coef[955]+0.5*coef[954]; 
    nv[7] = 0.5*coef[1003]+0.5*coef[1002]+0.5*coef[1000]+0.5*coef[999]; 
    nv[8] = 0.5*coef[211]+0.5*coef[210]+0.5*coef[208]+0.5*coef[207]; 
    nv[9] = 0.5*coef[256]+0.5*coef[255]+0.5*coef[253]+0.5*coef[252]; 
    nv[10] = 0.5*coef[589]+0.5*coef[588]+0.5*coef[586]+0.5*coef[585]; 
    nv[11] = 0.5*coef[634]+0.5*coef[633]+0.5*coef[631]+0.5*coef[630]; 
    nv[12] = 0.5*coef[967]+0.5*coef[966]+0.5*coef[964]+0.5*coef[963]; 
    nv[13] = 0.5*coef[1012]+0.5*coef[1011]+0.5*coef[1009]+0.5*coef[1008]; 
    nv[14] = 0.5*coef[220]+0.5*coef[219]+0.5*coef[217]+0.5*coef[216]; 
    nv[15] = 0.5*coef[265]+0.5*coef[264]+0.5*coef[262]+0.5*coef[261]; 
    nv[16] = 0.5*coef[598]+0.5*coef[597]+0.5*coef[595]+0.5*coef[594]; 
    nv[17] = 0.5*coef[643]+0.5*coef[642]+0.5*coef[640]+0.5*coef[639]; 
    nv[18] = 0.5*coef[976]+0.5*coef[975]+0.5*coef[973]+0.5*coef[972]; 
    nv[19] = 0.5*coef[1021]+0.5*coef[1020]+0.5*coef[1018]+0.5*coef[1017]; 
    nv[20] = 0.5*coef[229]+0.5*coef[228]+0.5*coef[226]+0.5*coef[225]; 
    nv[21] = 0.5*coef[274]+0.5*coef[273]+0.5*coef[271]+0.5*coef[270]; 
    nv[22] = 0.5*coef[607]+0.5*coef[606]+0.5*coef[604]+0.5*coef[603]; 
    nv[23] = 0.5*coef[652]+0.5*coef[651]+0.5*coef[649]+0.5*coef[648]; 
    nv[24] = 0.5*coef[985]+0.5*coef[984]+0.5*coef[982]+0.5*coef[981]; 
    nv[25] = 0.5*coef[1030]+0.5*coef[1029]+0.5*coef[1027]+0.5*coef[1026]; 
    nv[26] = 0.5*coef[238]+0.5*coef[237]+0.5*coef[235]+0.5*coef[234]; 
    nv[27] = 0.5*coef[283]+0.5*coef[282]+0.5*coef[280]+0.5*coef[279]; 
    nv[28] = 0.5*coef[616]+0.5*coef[615]+0.5*coef[613]+0.5*coef[612]; 
    nv[29] = 0.5*coef[661]+0.5*coef[660]+0.5*coef[658]+0.5*coef[657]; 
    nv[30] = 0.5*coef[994]+0.5*coef[993]+0.5*coef[991]+0.5*coef[990]; 
    nv[31] = 0.5*coef[1039]+0.5*coef[1038]+0.5*coef[1036]+0.5*coef[1035]; 
    nv[32] = 0.5*coef[292]+0.5*coef[291]+0.5*coef[289]+0.5*coef[288]; 
    nv[33] = 0.5*coef[319]+0.5*coef[318]+0.5*coef[316]+0.5*coef[315]; 
    nv[34] = 0.5*coef[346]+0.5*coef[345]+0.5*coef[343]+0.5*coef[342]; 
    nv[35] = 0.5*coef[670]+0.5*coef[669]+0.5*coef[667]+0.5*coef[666]; 
    nv[36] = 0.5*coef[697]+0.5*coef[696]+0.5*coef[694]+0.5*coef[693]; 
    nv[37] = 0.5*coef[724]+0.5*coef[723]+0.5*coef[721]+0.5*coef[720]; 
    nv[38] = 0.5*coef[1048]+0.5*coef[1047]+0.5*coef[1045]+0.5*coef[1044]; 
    nv[39] = 0.5*coef[1075]+0.5*coef[1074]+0.5*coef[1072]+0.5*coef[1071]; 
    nv[40] = 0.5*coef[1102]+0.5*coef[1101]+0.5*coef[1099]+0.5*coef[1098]; 
    nv[41] = 0.5*coef[301]+0.5*coef[300]+0.5*coef[298]+0.5*coef[297]; 
    nv[42] = 0.5*coef[328]+0.5*coef[327]+0.5*coef[325]+0.5*coef[324]; 
    nv[43] = 0.5*coef[355]+0.5*coef[354]+0.5*coef[352]+0.5*coef[351]; 
    nv[44] = 0.5*coef[679]+0.5*coef[678]+0.5*coef[676]+0.5*coef[675]; 
    nv[45] = 0.5*coef[706]+0.5*coef[705]+0.5*coef[703]+0.5*coef[702]; 
    nv[46] = 0.5*coef[733]+0.5*coef[732]+0.5*coef[730]+0.5*coef[729]; 
    nv[47] = 0.5*coef[1057]+0.5*coef[1056]+0.5*coef[1054]+0.5*coef[1053]; 
    nv[48] = 0.5*coef[1084]+0.5*coef[1083]+0.5*coef[1081]+0.5*coef[1080]; 
    nv[49] = 0.5*coef[1111]+0.5*coef[1110]+0.5*coef[1108]+0.5*coef[1107]; 
    nv[50] = 0.5*coef[310]+0.5*coef[309]+0.5*coef[307]+0.5*coef[306]; 
    nv[51] = 0.5*coef[337]+0.5*coef[336]+0.5*coef[334]+0.5*coef[333]; 
    nv[52] = 0.5*coef[364]+0.5*coef[363]+0.5*coef[361]+0.5*coef[360]; 
    nv[53] = 0.5*coef[688]+0.5*coef[687]+0.5*coef[685]+0.5*coef[684]; 
    nv[54] = 0.5*coef[715]+0.5*coef[714]+0.5*coef[712]+0.5*coef[711]; 
    nv[55] = 0.5*coef[742]+0.5*coef[741]+0.5*coef[739]+0.5*coef[738]; 
    nv[56] = 0.5*coef[1066]+0.5*coef[1065]+0.5*coef[1063]+0.5*coef[1062]; 
    nv[57] = 0.5*coef[1093]+0.5*coef[1092]+0.5*coef[1090]+0.5*coef[1089]; 
    nv[58] = 0.5*coef[1120]+0.5*coef[1119]+0.5*coef[1117]+0.5*coef[1116]; 
    nv[59] = 0.5*coef[40]+0.5*coef[39]+0.5*coef[37]+0.5*coef[36]; 
    nv[60] = 0.5*coef[49]+0.5*coef[48]+0.5*coef[46]+0.5*coef[45]; 
    nv[61] = 0.5*coef[58]+0.5*coef[57]+0.5*coef[55]+0.5*coef[54]; 
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
