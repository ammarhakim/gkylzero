#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:      integer to indicate which matrix being fetched. 
  // x:          Input solution vector. 
  // euler_pkpm: [rho ux, rho uy, rho uz], PKPM Fluid momentum output state vector.
  // em:         [Ex, Ey, Ez, Bx, By, Bz], EM output state vector.
  //             Source solve only updates Ex, Ey, Ez. 

  double rhoux_new[4] = {0.0}; 
  double rhouy_new[4] = {0.0}; 
  double rhouz_new[4] = {0.0}; 
  double Ex_new[4] = {0.0}; 
  double Ey_new[4] = {0.0}; 
  double Ez_new[4] = {0.0}; 

  double rhoux_new_nodal[GKYL_MAX_SPECIES][4]; 
  double rhouy_new_nodal[GKYL_MAX_SPECIES][4]; 
  double rhouz_new_nodal[GKYL_MAX_SPECIES][4]; 
  double Ex_new_nodal[4]; 
  double Ey_new_nodal[4]; 
  double Ez_new_nodal[4]; 

  // Factor of 2.0 is because solution at new time-step is A^{n+1} = 2.0*A_bar - A^{n}. 
  struct gkyl_mat sol0 = gkyl_nmat_get(x, count+0); 
  for (int s = 0; s < num_species; ++s) { 
    rhoux_new_nodal[s][0] = gkyl_mat_get(&sol0, 0 + s*3, 0); 
    rhouy_new_nodal[s][0] = gkyl_mat_get(&sol0, 1 + s*3, 0); 
    rhouz_new_nodal[s][0] = gkyl_mat_get(&sol0, 2 + s*3, 0); 
  } 
  Ex_new_nodal[0] = gkyl_mat_get(&sol0, 0 + num_species*3, 0); 
  Ey_new_nodal[0] = gkyl_mat_get(&sol0, 1 + num_species*3, 0); 
  Ez_new_nodal[0] = gkyl_mat_get(&sol0, 2 + num_species*3, 0); 

  struct gkyl_mat sol1 = gkyl_nmat_get(x, count+1); 
  for (int s = 0; s < num_species; ++s) { 
    rhoux_new_nodal[s][1] = gkyl_mat_get(&sol1, 0 + s*3, 0); 
    rhouy_new_nodal[s][1] = gkyl_mat_get(&sol1, 1 + s*3, 0); 
    rhouz_new_nodal[s][1] = gkyl_mat_get(&sol1, 2 + s*3, 0); 
  } 
  Ex_new_nodal[1] = gkyl_mat_get(&sol1, 0 + num_species*3, 0); 
  Ey_new_nodal[1] = gkyl_mat_get(&sol1, 1 + num_species*3, 0); 
  Ez_new_nodal[1] = gkyl_mat_get(&sol1, 2 + num_species*3, 0); 

  struct gkyl_mat sol2 = gkyl_nmat_get(x, count+2); 
  for (int s = 0; s < num_species; ++s) { 
    rhoux_new_nodal[s][2] = gkyl_mat_get(&sol2, 0 + s*3, 0); 
    rhouy_new_nodal[s][2] = gkyl_mat_get(&sol2, 1 + s*3, 0); 
    rhouz_new_nodal[s][2] = gkyl_mat_get(&sol2, 2 + s*3, 0); 
  } 
  Ex_new_nodal[2] = gkyl_mat_get(&sol2, 0 + num_species*3, 0); 
  Ey_new_nodal[2] = gkyl_mat_get(&sol2, 1 + num_species*3, 0); 
  Ez_new_nodal[2] = gkyl_mat_get(&sol2, 2 + num_species*3, 0); 

  struct gkyl_mat sol3 = gkyl_nmat_get(x, count+3); 
  for (int s = 0; s < num_species; ++s) { 
    rhoux_new_nodal[s][3] = gkyl_mat_get(&sol3, 0 + s*3, 0); 
    rhouy_new_nodal[s][3] = gkyl_mat_get(&sol3, 1 + s*3, 0); 
    rhouz_new_nodal[s][3] = gkyl_mat_get(&sol3, 2 + s*3, 0); 
  } 
  Ex_new_nodal[3] = gkyl_mat_get(&sol3, 0 + num_species*3, 0); 
  Ey_new_nodal[3] = gkyl_mat_get(&sol3, 1 + num_species*3, 0); 
  Ez_new_nodal[3] = gkyl_mat_get(&sol3, 2 + num_species*3, 0); 

  // Project nodal solution onto modal bases and subtract solution at old time step. 
  for (int s = 0; s < num_species; ++s) { 
    quad_nodal_to_modal_2d_ser_p1(rhoux_new_nodal[s], rhoux_new);
    quad_nodal_to_modal_2d_ser_p1(rhouy_new_nodal[s], rhouy_new);
    quad_nodal_to_modal_2d_ser_p1(rhouz_new_nodal[s], rhouz_new);
 
    double *out_rhoux = &euler_pkpm[s][0]; 
    double *out_rhouy = &euler_pkpm[s][4]; 
    double *out_rhouz = &euler_pkpm[s][8]; 

    out_rhoux[0] = 2.0*rhoux_new[0]/qbym[s] - out_rhoux[0]; 
    out_rhouy[0] = 2.0*rhouy_new[0]/qbym[s] - out_rhouy[0]; 
    out_rhouz[0] = 2.0*rhouz_new[0]/qbym[s] - out_rhouz[0]; 

    out_rhoux[1] = 2.0*rhoux_new[1]/qbym[s] - out_rhoux[1]; 
    out_rhouy[1] = 2.0*rhouy_new[1]/qbym[s] - out_rhouy[1]; 
    out_rhouz[1] = 2.0*rhouz_new[1]/qbym[s] - out_rhouz[1]; 

    out_rhoux[2] = 2.0*rhoux_new[2]/qbym[s] - out_rhoux[2]; 
    out_rhouy[2] = 2.0*rhouy_new[2]/qbym[s] - out_rhouy[2]; 
    out_rhouz[2] = 2.0*rhouz_new[2]/qbym[s] - out_rhouz[2]; 

    out_rhoux[3] = 2.0*rhoux_new[3]/qbym[s] - out_rhoux[3]; 
    out_rhouy[3] = 2.0*rhouy_new[3]/qbym[s] - out_rhouy[3]; 
    out_rhouz[3] = 2.0*rhouz_new[3]/qbym[s] - out_rhouz[3]; 

  } 

  quad_nodal_to_modal_2d_ser_p1(Ex_new_nodal, Ex_new);
  quad_nodal_to_modal_2d_ser_p1(Ey_new_nodal, Ey_new);
  quad_nodal_to_modal_2d_ser_p1(Ez_new_nodal, Ez_new);

  double *out_Ex = &em[0]; 
  double *out_Ey = &em[4]; 
  double *out_Ez = &em[8]; 

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

} 
