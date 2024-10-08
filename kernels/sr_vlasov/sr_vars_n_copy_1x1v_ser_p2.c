#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_sqrt.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void sr_vars_n_copy_1x1v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n) 
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
 
  double Gamma2_inv[3] = {0.0}; 
  Gamma2_inv[0] = 1.414213562373095-1.0*V_0_sq[0]; 
  Gamma2_inv[1] = -1.0*V_0_sq[1]; 
  Gamma2_inv[2] = -1.0*V_0_sq[2]; 

  int cell_avg = 0;
  if (0.6324555320336759*Gamma2_inv[2]-0.9486832980505137*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*Gamma2_inv[0]-0.7905694150420947*Gamma2_inv[2] < 0.0) cell_avg = 1; 
  if (0.6324555320336759*Gamma2_inv[2]+0.9486832980505137*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    double Gamma2_inv_lobatto[2] = {0.0}; 
    double Gamma2_inv_p1[2] = {0.0}; 
    double Gamma_inv_p1[2] = {0.0}; 
    double V0_quad = 0.0; 

    V0_quad = 1.58113883008419*V_0[2]-1.224744871391589*V_0[1]+0.7071067811865475*V_0[0]; 
    Gamma2_inv_lobatto[0] = 1.0 - V0_quad*V0_quad; 
    if (Gamma2_inv_lobatto[0] < 0.0) Gamma2_inv_lobatto[0] = 1.0e-16; 

    V0_quad = 1.58113883008419*V_0[2]+1.224744871391589*V_0[1]+0.7071067811865475*V_0[0]; 
    Gamma2_inv_lobatto[1] = 1.0 - V0_quad*V0_quad; 
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
 
