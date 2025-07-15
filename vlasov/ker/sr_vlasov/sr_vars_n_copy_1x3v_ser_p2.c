#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_sqrt.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void sr_vars_n_copy_1x3v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // M0:    Lab frame density = Gamma*n.
  // n:     Rest-frame density computed as Gamma_inv*M0 where Gamma_inv = sqrt(1 - |V_drift|^2). 
 
  struct gkyl_mat x0 = gkyl_nmat_get(x, count+0); 
  double V_0[3] = {0.0}; 
  V_0[0] = gkyl_mat_get(&x0,0,0); 
  V_0[1] = gkyl_mat_get(&x0,1,0); 
  V_0[2] = gkyl_mat_get(&x0,2,0); 
  double V_0_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_0, V_0_sq); 
 
  struct gkyl_mat x1 = gkyl_nmat_get(x, count+1); 
  double V_1[3] = {0.0}; 
  V_1[0] = gkyl_mat_get(&x1,0,0); 
  V_1[1] = gkyl_mat_get(&x1,1,0); 
  V_1[2] = gkyl_mat_get(&x1,2,0); 
  double V_1_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_1, V_1_sq); 
 
  struct gkyl_mat x2 = gkyl_nmat_get(x, count+2); 
  double V_2[3] = {0.0}; 
  V_2[0] = gkyl_mat_get(&x2,0,0); 
  V_2[1] = gkyl_mat_get(&x2,1,0); 
  V_2[2] = gkyl_mat_get(&x2,2,0); 
  double V_2_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_2, V_2_sq); 
 
  double Gamma2_inv[3] = {0.0}; 
  Gamma2_inv[0] = (-1.0*V_2_sq[0])-1.0*V_1_sq[0]-1.0*V_0_sq[0]+1.414213562373095; 
  Gamma2_inv[1] = (-1.0*V_2_sq[1])-1.0*V_1_sq[1]-1.0*V_0_sq[1]; 
  Gamma2_inv[2] = (-1.0*V_2_sq[2])-1.0*V_1_sq[2]-1.0*V_0_sq[2]; 

  int cell_avg = 0;
  if (0.6324555320336759*Gamma2_inv[2]-0.9486832980505137*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*Gamma2_inv[0]-0.7905694150420947*Gamma2_inv[2] < 0.0) cell_avg = 1; 
  if (0.6324555320336759*Gamma2_inv[2]+0.9486832980505137*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    double Gamma2_inv_lobatto[2] = {0.0}; 
    double Gamma2_inv_p1[2] = {0.0}; 
    double Gamma_inv_p1[2] = {0.0}; 
    double V0_quad = 0.0; 
    double V1_quad = 0.0; 
    double V2_quad = 0.0; 

    V0_quad = 1.58113883008419*V_0[2]-1.224744871391589*V_0[1]+0.7071067811865475*V_0[0]; 
    V1_quad = 1.58113883008419*V_1[2]-1.224744871391589*V_1[1]+0.7071067811865475*V_1[0]; 
    V2_quad = 1.58113883008419*V_2[2]-1.224744871391589*V_2[1]+0.7071067811865475*V_2[0]; 
    Gamma2_inv_lobatto[0] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad - V2_quad*V2_quad; 
    if (Gamma2_inv_lobatto[0] < 0.0) Gamma2_inv_lobatto[0] = 1.0e-16; 

    V0_quad = 1.58113883008419*V_0[2]+1.224744871391589*V_0[1]+0.7071067811865475*V_0[0]; 
    V1_quad = 1.58113883008419*V_1[2]+1.224744871391589*V_1[1]+0.7071067811865475*V_1[0]; 
    V2_quad = 1.58113883008419*V_2[2]+1.224744871391589*V_2[1]+0.7071067811865475*V_2[0]; 
    Gamma2_inv_lobatto[1] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad - V2_quad*V2_quad; 
    if (Gamma2_inv_lobatto[1] < 0.0) Gamma2_inv_lobatto[1] = 1.0e-16; 

    Gamma2_inv_p1[0] = 0.7071067811865475*Gamma2_inv_lobatto[1]+0.7071067811865475*Gamma2_inv_lobatto[0]; 
    Gamma2_inv_p1[1] = 0.408248290463863*Gamma2_inv_lobatto[1]-0.408248290463863*Gamma2_inv_lobatto[0]; 
    ser_1x_p1_sqrt(Gamma2_inv_p1, Gamma_inv_p1); 
    n[0] = 0.7071067811865475*Gamma_inv_p1[1]*M0[1]+0.7071067811865475*Gamma_inv_p1[0]*M0[0]; 
    n[1] = 0.6324555320336759*Gamma_inv_p1[1]*M0[2]+0.7071067811865475*Gamma_inv_p1[0]*M0[1]+0.7071067811865475*M0[0]*Gamma_inv_p1[1]; 
    n[2] = 0.7071067811865475*Gamma_inv_p1[0]*M0[2]+0.6324555320336759*Gamma_inv_p1[1]*M0[1]; 
  } 
  else { 
    double Gamma_inv[3] = {0.0}; 
    ser_1x_p2_sqrt(Gamma2_inv, Gamma_inv); 
    n[0] = 0.7071067811865475*Gamma_inv[2]*M0[2]+0.7071067811865475*Gamma_inv[1]*M0[1]+0.7071067811865475*Gamma_inv[0]*M0[0]; 
    n[1] = 0.6324555320336759*Gamma_inv[1]*M0[2]+0.6324555320336759*M0[1]*Gamma_inv[2]+0.7071067811865475*Gamma_inv[0]*M0[1]+0.7071067811865475*M0[0]*Gamma_inv[1]; 
    n[2] = 0.4517539514526256*Gamma_inv[2]*M0[2]+0.7071067811865475*Gamma_inv[0]*M0[2]+0.7071067811865475*M0[0]*Gamma_inv[2]+0.6324555320336759*Gamma_inv[1]*M0[1]; 
  } 

} 
 
