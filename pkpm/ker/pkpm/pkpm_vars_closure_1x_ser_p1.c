#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_closure_1x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *grad_b, const double *grad_u, const double *xi, const double *grad_xi, const double *theta,
  const double *alpha_0_prev, const double *c_h_prev, const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl) 
{ 
  // w_lo, w_up:   Velocity domain edges; dv: velocity cell size. 
  // norm_moms:    Input [M0, M2, M0_G, M1, ...] moments of the normalized distributions (stage state). 
  // rhs_moms:     Input [M0, M2, M0_G, M1, ...] moments of the stage-A right-hand sides. 
  // bvar, grad_b, grad_u, xi, grad_xi, theta: fields of the forces (pkpm_vars_grad_u_* layouts). 
  // alpha_0_prev, c_h_prev: closures of the previous stage (may alias the outputs), for the upwind switch. 
  // f_lo, f_up:   [F, G] of the lowest and uppermost velocity cells of this configuration cell. 
  // c_0, alpha_0, c_h: Output closures at the 2 Gauss-Legendre nodes, interpolated to p=1. 
  // flux_lo, flux_up: Output edge fluxes [F flux, G flux] in the stored face format. 
  // cfl:          Output [lower, upper] edge-flux time-step rates. 
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
  const double *b0 = &bvar[0]; 
  const double *b1 = &bvar[2]; 
  const double *b2 = &bvar[4]; 
  const double *gux0 = &grad_u[0]; 
  const double *gux1 = &grad_u[2]; 
  const double *gux2 = &grad_u[4]; 
  const double *gbx = &grad_b[0]; 
  const double *gxix = &grad_xi[0]; 
  const double *xi_1 = &xi[0]; 
  const double *F_lo = &f_lo[0]; 
  const double *G_lo = &f_lo[6]; 
  const double *F_up = &f_up[0]; 
  const double *G_up = &f_up[6]; 
 
  double c0_n[2] = {0.0}, a0_n[2] = {0.0}, c_n[2] = {0.0}; 
  double flo_n[2] = {0.0}, glo_n[2] = {0.0}, fup_n[2] = {0.0}, gup_n[2] = {0.0}; 
  double cfl_lo = 0.0, cfl_up = 0.0; 
  int num_null = 0; 
 
  { 
    const double n0_q = 0.7071067811865475*n0[0]-0.7071067811865475*n0[1]; 
    const double n1_q = 0.7071067811865475*n1[0]-0.7071067811865475*n1[1]; 
    const double n2_q = 0.7071067811865475*n2[0]-0.7071067811865475*n2[1]; 
    const double nG_q = 0.7071067811865475*nG[0]-0.7071067811865475*nG[1]; 
    const double r0_q = 0.7071067811865475*r0[0]-0.7071067811865475*r0[1]; 
    const double r1_q = 0.7071067811865475*r1[0]-0.7071067811865475*r1[1]; 
    const double r2_q = 0.7071067811865475*r2[0]-0.7071067811865475*r2[1] + 2.0*(0.7071067811865475*rG[0]-0.7071067811865475*rG[1]); 
    const double a0p_q = 0.7071067811865475*alpha_0_prev[0]-0.7071067811865475*alpha_0_prev[1]; 
    const double cp_q = 0.7071067811865475*c_h_prev[0]-0.7071067811865475*c_h_prev[1]; 
    const double xi_q = 0.7071067811865475*xi_1[0]-0.7071067811865475*xi_1[1]; 
    const double th_q = 0.7071067811865475*theta[0]-0.7071067811865475*theta[1]; 
    const double b0_q = 0.7071067811865475*b0[0]-0.7071067811865475*b0[1]; 
    const double b1_q = 0.7071067811865475*b1[0]-0.7071067811865475*b1[1]; 
    const double b2_q = 0.7071067811865475*b2[0]-0.7071067811865475*b2[1]; 
    const double gu00_q = 0.7071067811865475*gux0[0]-0.7071067811865475*gux0[1]; 
    const double gu01_q = 0.7071067811865475*gux1[0]-0.7071067811865475*gux1[1]; 
    const double gu02_q = 0.7071067811865475*gux2[0]-0.7071067811865475*gux2[1]; 
    const double gxi0_q = 0.7071067811865475*gxix[0]-0.7071067811865475*gxix[1]; 
    const double gb0_q = 0.7071067811865475*gbx[0]-0.7071067811865475*gbx[1]; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q; 
    const double beta_q = xi_q*(0.0 + gb0_q); 
    const double Flo_q = -(1.1180339887498951*F_lo[5])+1.118033988749895*F_lo[4]+0.8660254037844386*F_lo[3]-0.8660254037844386*F_lo[2]-0.5*F_lo[1]+0.5*F_lo[0]; 
    const double Glo_q = -(1.1180339887498951*G_lo[5])+1.118033988749895*G_lo[4]+0.8660254037844386*G_lo[3]-0.8660254037844386*G_lo[2]-0.5*G_lo[1]+0.5*G_lo[0]; 
    const double Fup_q = -(1.1180339887498951*F_up[5])+1.118033988749895*F_up[4]-0.8660254037844386*F_up[3]+0.8660254037844386*F_up[2]-0.5*F_up[1]+0.5*F_up[0]; 
    const double Gup_q = -(1.1180339887498951*G_up[5])+1.118033988749895*G_up[4]-0.8660254037844386*G_up[3]+0.8660254037844386*G_up[2]-0.5*G_up[1]+0.5*G_up[0]; 
    const double Hlo_q = 4.0*th_q*Glo_q - 2.0*th_q*th_q*Flo_q, Hup_q = 4.0*th_q*Gup_q - 2.0*th_q*th_q*Fup_q; 
    // known part of the force at the edges and the upwind switches (previous closures) 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double s_lo = (a0p_q + cp_q*w_lo + k_lo) < 0.0 ? 1.0 : 0.0, s_up = (a0p_q + cp_q*w_up + k_up) > 0.0 ? 1.0 : 0.0; 
    const double sb_lo = beta_q < 0.0 ? 1.0 : 0.0, sb_up = beta_q > 0.0 ? 1.0 : 0.0; 
    // edge contributions to the moment rates: the k-th moment gains -n_e w_e^k Phi_e per edge, 
    // n_e = +1 up, -1 lo, Phi_e the outward-switched +w flux (closure part times these, known part kF/bH) 
    const double eF0 = s_lo*Flo_q - s_up*Fup_q; 
    const double eF1 = s_lo*w_lo*Flo_q - s_up*w_up*Fup_q; 
    const double eF2 = s_lo*w_lo*w_lo*Flo_q - s_up*w_up*w_up*Fup_q; 
    const double eF3 = s_lo*w_lo*w_lo*w_lo*Flo_q - s_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = s_lo*Glo_q - s_up*Gup_q; 
    const double eG1 = s_lo*w_lo*Glo_q - s_up*w_up*Gup_q; 
    const double kF0 = s_lo*k_lo*Flo_q - s_up*k_up*Fup_q; 
    const double kF1 = s_lo*w_lo*k_lo*Flo_q - s_up*w_up*k_up*Fup_q; 
    const double kF2 = s_lo*w_lo*w_lo*k_lo*Flo_q - s_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = s_lo*k_lo*Glo_q - s_up*k_up*Gup_q; 
    const double bH0 = sb_lo*beta_q*Hlo_q - sb_up*beta_q*Hup_q; 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q, A21 = 2.0*n1_q + eF2 + 2.0*eG0, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1; 
    const double b0 = -r0_q - kF0, b1 = -r1_q - kF1, b2 = -r2_q - kF2 - 2.0*kG0 - 2.0*bH0; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[0] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[0] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[0] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      const double al_lo = a0_n[0] + c_n[0]*w_lo + k_lo, al_up = a0_n[0] + c_n[0]*w_up + k_up; 
      flo_n[0] = s_lo*al_lo*Flo_q; 
      glo_n[0] = s_lo*al_lo*Glo_q + sb_lo*beta_q*Hlo_q; 
      fup_n[0] = s_up*al_up*Fup_q; 
      gup_n[0] = s_up*al_up*Gup_q + sb_up*beta_q*Hup_q; 
      cfl_lo = fmax(cfl_lo, fabs(al_lo) + 3.414213562373095*fabs(th_q*beta_q)); 
      cfl_up = fmax(cfl_up, fabs(al_up) + 3.414213562373095*fabs(th_q*beta_q)); 
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
    const double a0p_q = 0.7071067811865475*alpha_0_prev[1]+0.7071067811865475*alpha_0_prev[0]; 
    const double cp_q = 0.7071067811865475*c_h_prev[1]+0.7071067811865475*c_h_prev[0]; 
    const double xi_q = 0.7071067811865475*xi_1[1]+0.7071067811865475*xi_1[0]; 
    const double th_q = 0.7071067811865475*theta[1]+0.7071067811865475*theta[0]; 
    const double b0_q = 0.7071067811865475*b0[1]+0.7071067811865475*b0[0]; 
    const double b1_q = 0.7071067811865475*b1[1]+0.7071067811865475*b1[0]; 
    const double b2_q = 0.7071067811865475*b2[1]+0.7071067811865475*b2[0]; 
    const double gu00_q = 0.7071067811865475*gux0[1]+0.7071067811865475*gux0[0]; 
    const double gu01_q = 0.7071067811865475*gux1[1]+0.7071067811865475*gux1[0]; 
    const double gu02_q = 0.7071067811865475*gux2[1]+0.7071067811865475*gux2[0]; 
    const double gxi0_q = 0.7071067811865475*gxix[1]+0.7071067811865475*gxix[0]; 
    const double gb0_q = 0.7071067811865475*gbx[1]+0.7071067811865475*gbx[0]; 
    const double K_q = 0.0 + b0_q*b0_q*gu00_q + b0_q*b1_q*gu01_q + b0_q*b2_q*gu02_q; 
    const double bgx_q = 0.0 + b0_q*gxi0_q; 
    const double beta_q = xi_q*(0.0 + gb0_q); 
    const double Flo_q = 1.1180339887498951*F_lo[5]+1.118033988749895*F_lo[4]-0.8660254037844386*F_lo[3]-0.8660254037844386*F_lo[2]+0.5*F_lo[1]+0.5*F_lo[0]; 
    const double Glo_q = 1.1180339887498951*G_lo[5]+1.118033988749895*G_lo[4]-0.8660254037844386*G_lo[3]-0.8660254037844386*G_lo[2]+0.5*G_lo[1]+0.5*G_lo[0]; 
    const double Fup_q = 1.1180339887498951*F_up[5]+1.118033988749895*F_up[4]+0.8660254037844386*F_up[3]+0.8660254037844386*F_up[2]+0.5*F_up[1]+0.5*F_up[0]; 
    const double Gup_q = 1.1180339887498951*G_up[5]+1.118033988749895*G_up[4]+0.8660254037844386*G_up[3]+0.8660254037844386*G_up[2]+0.5*G_up[1]+0.5*G_up[0]; 
    const double Hlo_q = 4.0*th_q*Glo_q - 2.0*th_q*th_q*Flo_q, Hup_q = 4.0*th_q*Gup_q - 2.0*th_q*th_q*Fup_q; 
    // known part of the force at the edges and the upwind switches (previous closures) 
    const double k_lo = -K_q*w_lo - bgx_q*w_lo*w_lo, k_up = -K_q*w_up - bgx_q*w_up*w_up; 
    const double s_lo = (a0p_q + cp_q*w_lo + k_lo) < 0.0 ? 1.0 : 0.0, s_up = (a0p_q + cp_q*w_up + k_up) > 0.0 ? 1.0 : 0.0; 
    const double sb_lo = beta_q < 0.0 ? 1.0 : 0.0, sb_up = beta_q > 0.0 ? 1.0 : 0.0; 
    // edge contributions to the moment rates: the k-th moment gains -n_e w_e^k Phi_e per edge, 
    // n_e = +1 up, -1 lo, Phi_e the outward-switched +w flux (closure part times these, known part kF/bH) 
    const double eF0 = s_lo*Flo_q - s_up*Fup_q; 
    const double eF1 = s_lo*w_lo*Flo_q - s_up*w_up*Fup_q; 
    const double eF2 = s_lo*w_lo*w_lo*Flo_q - s_up*w_up*w_up*Fup_q; 
    const double eF3 = s_lo*w_lo*w_lo*w_lo*Flo_q - s_up*w_up*w_up*w_up*Fup_q; 
    const double eG0 = s_lo*Glo_q - s_up*Gup_q; 
    const double eG1 = s_lo*w_lo*Glo_q - s_up*w_up*Gup_q; 
    const double kF0 = s_lo*k_lo*Flo_q - s_up*k_up*Fup_q; 
    const double kF1 = s_lo*w_lo*k_lo*Flo_q - s_up*w_up*k_up*Fup_q; 
    const double kF2 = s_lo*w_lo*w_lo*k_lo*Flo_q - s_up*w_up*w_up*k_up*Fup_q; 
    const double kG0 = s_lo*k_lo*Glo_q - s_up*k_up*Gup_q; 
    const double bH0 = sb_lo*beta_q*Hlo_q - sb_up*beta_q*Hup_q; 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q, A21 = 2.0*n1_q + eF2 + 2.0*eG0, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1; 
    const double b0 = -r0_q - kF0, b1 = -r1_q - kF1, b2 = -r2_q - kF2 - 2.0*kG0 - 2.0*bH0; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_n[1] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_n[1] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_n[1] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      const double al_lo = a0_n[1] + c_n[1]*w_lo + k_lo, al_up = a0_n[1] + c_n[1]*w_up + k_up; 
      flo_n[1] = s_lo*al_lo*Flo_q; 
      glo_n[1] = s_lo*al_lo*Glo_q + sb_lo*beta_q*Hlo_q; 
      fup_n[1] = s_up*al_up*Fup_q; 
      gup_n[1] = s_up*al_up*Gup_q + sb_up*beta_q*Hup_q; 
      cfl_lo = fmax(cfl_lo, fabs(al_lo) + 3.414213562373095*fabs(th_q*beta_q)); 
      cfl_up = fmax(cfl_up, fabs(al_up) + 3.414213562373095*fabs(th_q*beta_q)); 
    } else { 
      num_null += 1; 
    } 
  } 
 
  cfl[0] = 2.5*dv1*cfl_lo; 
  cfl[1] = 2.5*dv1*cfl_up; 
  if (num_null > 0) { 
    c_0[0] = 0.0; alpha_0[0] = 0.0; c_h[0] = 0.0; 
    c_0[1] = 0.0; alpha_0[1] = 0.0; c_h[1] = 0.0; 
    flux_lo[0] = 0.0; flux_up[0] = 0.0; 
    flux_lo[1] = 0.0; flux_up[1] = 0.0; 
    flux_lo[2] = 0.0; flux_up[2] = 0.0; 
    flux_lo[3] = 0.0; flux_up[3] = 0.0; 
    return num_null; 
  } 
 
  c_0[0] = 0.7071067811865475*c0_n[1]+0.7071067811865475*c0_n[0]; 
  c_0[1] = 0.7071067811865475*c0_n[1]-0.7071067811865475*c0_n[0]; 
  alpha_0[0] = 0.7071067811865475*a0_n[1]+0.7071067811865475*a0_n[0]; 
  alpha_0[1] = 0.7071067811865475*a0_n[1]-0.7071067811865475*a0_n[0]; 
  c_h[0] = 0.7071067811865475*c_n[1]+0.7071067811865475*c_n[0]; 
  c_h[1] = 0.7071067811865475*c_n[1]-0.7071067811865475*c_n[0]; 
 
  flux_lo[0] = flo_n[1]+flo_n[0]; 
  flux_lo[1] = 1.4433756729740645*flo_n[1]-1.4433756729740645*flo_n[0]; 
  flux_lo[2] = glo_n[1]+glo_n[0]; 
  flux_lo[3] = 1.4433756729740645*glo_n[1]-1.4433756729740645*glo_n[0]; 
  flux_up[0] = fup_n[1]+fup_n[0]; 
  flux_up[1] = 1.4433756729740645*fup_n[1]-1.4433756729740645*fup_n[0]; 
  flux_up[2] = gup_n[1]+gup_n[0]; 
  flux_up[3] = 1.4433756729740645*gup_n[1]-1.4433756729740645*gup_n[0]; 
 
  return 0; 
} 
