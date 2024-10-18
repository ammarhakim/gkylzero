#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double *pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:       integer to indicate which matrix being fetched. 
  // x:           Input solution vector (nodal basis representation of solution). 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model (modal basis) at old time t^n.
  // pkpm_u:      [ux, uy, uz], Input flow velocity (modal basis) at old time t^n.
  // euler_pkpm:  [rho ux, rho uy, rho uz], Fluid output state vector (modal basis) at time t^{n+1}.
  // em:          [Ex, Ey, Ez, Bx, By, Bz], EM output state vector (modal basis) at time t^{n+1}.
  //              Source solve only updates Ex, Ey, Ez. 

  double ux_new[8] = {0.0}; 
  double uy_new[8] = {0.0}; 
  double uz_new[8] = {0.0}; 
  double Ex_new[8] = {0.0}; 
  double Ey_new[8] = {0.0}; 
  double Ez_new[8] = {0.0}; 

  double ux_new_nodal[GKYL_MAX_SPECIES][8]; 
  double uy_new_nodal[GKYL_MAX_SPECIES][8]; 
  double uz_new_nodal[GKYL_MAX_SPECIES][8]; 
  double Ex_new_nodal[8]; 
  double Ey_new_nodal[8]; 
  double Ez_new_nodal[8]; 

  // Factor of 2.0 is because solution at new time-step is A^{n+1} = 2.0*A_bar - A^{n}. 
  struct gkyl_mat sol0 = gkyl_nmat_get(x, count+0); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][0] = 2.0*gkyl_mat_get(&sol0, 0 + s*3, 0); 
    uy_new_nodal[s][0] = 2.0*gkyl_mat_get(&sol0, 1 + s*3, 0); 
    uz_new_nodal[s][0] = 2.0*gkyl_mat_get(&sol0, 2 + s*3, 0); 
  } 
  Ex_new_nodal[0] = 2.0*gkyl_mat_get(&sol0, 0 + num_species*3, 0); 
  Ey_new_nodal[0] = 2.0*gkyl_mat_get(&sol0, 1 + num_species*3, 0); 
  Ez_new_nodal[0] = 2.0*gkyl_mat_get(&sol0, 2 + num_species*3, 0); 

  struct gkyl_mat sol1 = gkyl_nmat_get(x, count+1); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][1] = 2.0*gkyl_mat_get(&sol1, 0 + s*3, 0); 
    uy_new_nodal[s][1] = 2.0*gkyl_mat_get(&sol1, 1 + s*3, 0); 
    uz_new_nodal[s][1] = 2.0*gkyl_mat_get(&sol1, 2 + s*3, 0); 
  } 
  Ex_new_nodal[1] = 2.0*gkyl_mat_get(&sol1, 0 + num_species*3, 0); 
  Ey_new_nodal[1] = 2.0*gkyl_mat_get(&sol1, 1 + num_species*3, 0); 
  Ez_new_nodal[1] = 2.0*gkyl_mat_get(&sol1, 2 + num_species*3, 0); 

  struct gkyl_mat sol2 = gkyl_nmat_get(x, count+2); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][2] = 2.0*gkyl_mat_get(&sol2, 0 + s*3, 0); 
    uy_new_nodal[s][2] = 2.0*gkyl_mat_get(&sol2, 1 + s*3, 0); 
    uz_new_nodal[s][2] = 2.0*gkyl_mat_get(&sol2, 2 + s*3, 0); 
  } 
  Ex_new_nodal[2] = 2.0*gkyl_mat_get(&sol2, 0 + num_species*3, 0); 
  Ey_new_nodal[2] = 2.0*gkyl_mat_get(&sol2, 1 + num_species*3, 0); 
  Ez_new_nodal[2] = 2.0*gkyl_mat_get(&sol2, 2 + num_species*3, 0); 

  struct gkyl_mat sol3 = gkyl_nmat_get(x, count+3); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][3] = 2.0*gkyl_mat_get(&sol3, 0 + s*3, 0); 
    uy_new_nodal[s][3] = 2.0*gkyl_mat_get(&sol3, 1 + s*3, 0); 
    uz_new_nodal[s][3] = 2.0*gkyl_mat_get(&sol3, 2 + s*3, 0); 
  } 
  Ex_new_nodal[3] = 2.0*gkyl_mat_get(&sol3, 0 + num_species*3, 0); 
  Ey_new_nodal[3] = 2.0*gkyl_mat_get(&sol3, 1 + num_species*3, 0); 
  Ez_new_nodal[3] = 2.0*gkyl_mat_get(&sol3, 2 + num_species*3, 0); 

  struct gkyl_mat sol4 = gkyl_nmat_get(x, count+4); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][4] = 2.0*gkyl_mat_get(&sol4, 0 + s*3, 0); 
    uy_new_nodal[s][4] = 2.0*gkyl_mat_get(&sol4, 1 + s*3, 0); 
    uz_new_nodal[s][4] = 2.0*gkyl_mat_get(&sol4, 2 + s*3, 0); 
  } 
  Ex_new_nodal[4] = 2.0*gkyl_mat_get(&sol4, 0 + num_species*3, 0); 
  Ey_new_nodal[4] = 2.0*gkyl_mat_get(&sol4, 1 + num_species*3, 0); 
  Ez_new_nodal[4] = 2.0*gkyl_mat_get(&sol4, 2 + num_species*3, 0); 

  struct gkyl_mat sol5 = gkyl_nmat_get(x, count+5); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][5] = 2.0*gkyl_mat_get(&sol5, 0 + s*3, 0); 
    uy_new_nodal[s][5] = 2.0*gkyl_mat_get(&sol5, 1 + s*3, 0); 
    uz_new_nodal[s][5] = 2.0*gkyl_mat_get(&sol5, 2 + s*3, 0); 
  } 
  Ex_new_nodal[5] = 2.0*gkyl_mat_get(&sol5, 0 + num_species*3, 0); 
  Ey_new_nodal[5] = 2.0*gkyl_mat_get(&sol5, 1 + num_species*3, 0); 
  Ez_new_nodal[5] = 2.0*gkyl_mat_get(&sol5, 2 + num_species*3, 0); 

  struct gkyl_mat sol6 = gkyl_nmat_get(x, count+6); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][6] = 2.0*gkyl_mat_get(&sol6, 0 + s*3, 0); 
    uy_new_nodal[s][6] = 2.0*gkyl_mat_get(&sol6, 1 + s*3, 0); 
    uz_new_nodal[s][6] = 2.0*gkyl_mat_get(&sol6, 2 + s*3, 0); 
  } 
  Ex_new_nodal[6] = 2.0*gkyl_mat_get(&sol6, 0 + num_species*3, 0); 
  Ey_new_nodal[6] = 2.0*gkyl_mat_get(&sol6, 1 + num_species*3, 0); 
  Ez_new_nodal[6] = 2.0*gkyl_mat_get(&sol6, 2 + num_species*3, 0); 

  struct gkyl_mat sol7 = gkyl_nmat_get(x, count+7); 
  for (int s = 0; s < num_species; ++s) { 
    ux_new_nodal[s][7] = 2.0*gkyl_mat_get(&sol7, 0 + s*3, 0); 
    uy_new_nodal[s][7] = 2.0*gkyl_mat_get(&sol7, 1 + s*3, 0); 
    uz_new_nodal[s][7] = 2.0*gkyl_mat_get(&sol7, 2 + s*3, 0); 
  } 
  Ex_new_nodal[7] = 2.0*gkyl_mat_get(&sol7, 0 + num_species*3, 0); 
  Ey_new_nodal[7] = 2.0*gkyl_mat_get(&sol7, 1 + num_species*3, 0); 
  Ez_new_nodal[7] = 2.0*gkyl_mat_get(&sol7, 2 + num_species*3, 0); 

  // Project 2.0*nodal solution onto modal bases and subtract solution at old time step. 
  // We then convert back to momentum with rho^n. 
  for (int s = 0; s < num_species; ++s) { 
    quad_nodal_to_modal_3d_ser_p1(ux_new_nodal[s], ux_new);
    quad_nodal_to_modal_3d_ser_p1(uy_new_nodal[s], uy_new);
    quad_nodal_to_modal_3d_ser_p1(uz_new_nodal[s], uz_new);
 
    const double *rho_old = &vlasov_pkpm_moms[s][0]; 
    const double *ux_old = &pkpm_u[s][0]; 
    const double *uy_old = &pkpm_u[s][8]; 
    const double *uz_old = &pkpm_u[s][16]; 
    ux_new[0] -= ux_old[0]; 
    uy_new[0] -= uy_old[0]; 
    uz_new[0] -= uz_old[0]; 

    ux_new[1] -= ux_old[1]; 
    uy_new[1] -= uy_old[1]; 
    uz_new[1] -= uz_old[1]; 

    ux_new[2] -= ux_old[2]; 
    uy_new[2] -= uy_old[2]; 
    uz_new[2] -= uz_old[2]; 

    ux_new[3] -= ux_old[3]; 
    uy_new[3] -= uy_old[3]; 
    uz_new[3] -= uz_old[3]; 

    ux_new[4] -= ux_old[4]; 
    uy_new[4] -= uy_old[4]; 
    uz_new[4] -= uz_old[4]; 

    ux_new[5] -= ux_old[5]; 
    uy_new[5] -= uy_old[5]; 
    uz_new[5] -= uz_old[5]; 

    ux_new[6] -= ux_old[6]; 
    uy_new[6] -= uy_old[6]; 
    uz_new[6] -= uz_old[6]; 

    ux_new[7] -= ux_old[7]; 
    uy_new[7] -= uy_old[7]; 
    uz_new[7] -= uz_old[7]; 

    double *out_rhoux = &euler_pkpm[s][0]; 
    double *out_rhouy = &euler_pkpm[s][8]; 
    double *out_rhouz = &euler_pkpm[s][16]; 

    binop_mul_3d_ser_p1(rho_old, ux_new, out_rhoux); 
    binop_mul_3d_ser_p1(rho_old, uy_new, out_rhouy); 
    binop_mul_3d_ser_p1(rho_old, uz_new, out_rhouz); 
  } 

  quad_nodal_to_modal_3d_ser_p1(Ex_new_nodal, Ex_new);
  quad_nodal_to_modal_3d_ser_p1(Ey_new_nodal, Ey_new);
  quad_nodal_to_modal_3d_ser_p1(Ez_new_nodal, Ez_new);

  double *out_Ex = &em[0]; 
  double *out_Ey = &em[8]; 
  double *out_Ez = &em[16]; 

  out_Ex[0] = 2.0*Ex_new[0]/epsilon0 - out_Ex[0]; 
  out_Ey[0] = 2.0*Ey_new[0]/epsilon0 - out_Ey[0]; 
  out_Ez[0] = 2.0*Ez_new[0]/epsilon0 - out_Ez[0]; 

  out_Ex[1] = 2.0*Ex_new[1]/epsilon0 - out_Ex[1]; 
  out_Ey[1] = 2.0*Ey_new[1]/epsilon0 - out_Ey[1]; 
  out_Ez[1] = 2.0*Ez_new[1]/epsilon0 - out_Ez[1]; 

  out_Ex[2] = 2.0*Ex_new[2]/epsilon0 - out_Ex[2]; 
  out_Ey[2] = 2.0*Ey_new[2]/epsilon0 - out_Ey[2]; 
  out_Ez[2] = 2.0*Ez_new[2]/epsilon0 - out_Ez[2]; 

  out_Ex[3] = 2.0*Ex_new[3]/epsilon0 - out_Ex[3]; 
  out_Ey[3] = 2.0*Ey_new[3]/epsilon0 - out_Ey[3]; 
  out_Ez[3] = 2.0*Ez_new[3]/epsilon0 - out_Ez[3]; 

  out_Ex[4] = 2.0*Ex_new[4]/epsilon0 - out_Ex[4]; 
  out_Ey[4] = 2.0*Ey_new[4]/epsilon0 - out_Ey[4]; 
  out_Ez[4] = 2.0*Ez_new[4]/epsilon0 - out_Ez[4]; 

  out_Ex[5] = 2.0*Ex_new[5]/epsilon0 - out_Ex[5]; 
  out_Ey[5] = 2.0*Ey_new[5]/epsilon0 - out_Ey[5]; 
  out_Ez[5] = 2.0*Ez_new[5]/epsilon0 - out_Ez[5]; 

  out_Ex[6] = 2.0*Ex_new[6]/epsilon0 - out_Ex[6]; 
  out_Ey[6] = 2.0*Ey_new[6]/epsilon0 - out_Ey[6]; 
  out_Ez[6] = 2.0*Ez_new[6]/epsilon0 - out_Ez[6]; 

  out_Ex[7] = 2.0*Ex_new[7]/epsilon0 - out_Ex[7]; 
  out_Ey[7] = 2.0*Ey_new[7]/epsilon0 - out_Ey[7]; 
  out_Ez[7] = 2.0*Ez_new[7]/epsilon0 - out_Ez[7]; 

} 
