#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES],
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // app_accel:        Applied accelerations (external forces).
  // ext_em:           Externally applied EM fields.
  // app_current:      Applied external currents.
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.
  // em:               [Ex, Ey, Ez, Bx, By, Bz], EM input state vector.

  struct gkyl_mat lhs = gkyl_nmat_get(A_n, count); 
  struct gkyl_mat rhs = gkyl_nmat_get(rhs_n, count); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0); 

  double rho[GKYL_MAX_SPECIES][4]; 
  double rhoux[GKYL_MAX_SPECIES][4]; 
  double rhouy[GKYL_MAX_SPECIES][4]; 
  double rhouz[GKYL_MAX_SPECIES][4]; 

  double app_accel_x[GKYL_MAX_SPECIES][4]; 
  double app_accel_y[GKYL_MAX_SPECIES][4]; 
  double app_accel_z[GKYL_MAX_SPECIES][4]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_fluid = euler_pkpm[i]; 
    const double *inp_app_accel = app_accel[i]; 
    const double *inp_vlasov_pkpm_moms = vlasov_pkpm_moms[i]; 

    rho[i][0] = inp_vlasov_pkpm_moms[0]; 
    rhoux[i][0] = inp_fluid[0]; 
    rhouy[i][0] = inp_fluid[4]; 
    rhouz[i][0] = inp_fluid[8]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[4]; 
    app_accel_z[i][0] = inp_app_accel[8]; 

    rho[i][1] = inp_vlasov_pkpm_moms[1]; 
    rhoux[i][1] = inp_fluid[1]; 
    rhouy[i][1] = inp_fluid[5]; 
    rhouz[i][1] = inp_fluid[9]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[5]; 
    app_accel_z[i][1] = inp_app_accel[9]; 

    rho[i][2] = inp_vlasov_pkpm_moms[2]; 
    rhoux[i][2] = inp_fluid[2]; 
    rhouy[i][2] = inp_fluid[6]; 
    rhouz[i][2] = inp_fluid[10]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[6]; 
    app_accel_z[i][2] = inp_app_accel[10]; 

    rho[i][3] = inp_vlasov_pkpm_moms[3]; 
    rhoux[i][3] = inp_fluid[3]; 
    rhouy[i][3] = inp_fluid[7]; 
    rhouz[i][3] = inp_fluid[11]; 

    app_accel_x[i][3] = inp_app_accel[3]; 
    app_accel_y[i][3] = inp_app_accel[7]; 
    app_accel_z[i][3] = inp_app_accel[11]; 

  } 

  double *Ex = &em[0]; 
  double *Ey = &em[4]; 
  double *Ez = &em[8]; 
  double *Bx = &em[12]; 
  double *By = &em[16]; 
  double *Bz = &em[20]; 

  const double *ext_Ex = &ext_em[0]; 
  const double *ext_Ey = &ext_em[4]; 
  const double *ext_Ez = &ext_em[8]; 
  const double *ext_Bx = &ext_em[12]; 
  const double *ext_By = &ext_em[16]; 
  const double *ext_Bz = &ext_em[20]; 

  const double *app_curr_x = &app_current[0]; 
  const double *app_curr_y = &app_current[4]; 
  const double *app_curr_z = &app_current[8]; 

  double tot_Bx[4]; 
  double tot_By[4]; 
  double tot_Bz[4]; 
  tot_Bx[0] = Bx[0] + ext_Bx[0]; 
  tot_By[0] = By[0] + ext_By[0]; 
  tot_Bz[0] = Bz[0] + ext_Bz[0]; 
  tot_Bx[1] = Bx[1] + ext_Bx[1]; 
  tot_By[1] = By[1] + ext_By[1]; 
  tot_Bz[1] = Bz[1] + ext_Bz[1]; 
  tot_Bx[2] = Bx[2] + ext_Bx[2]; 
  tot_By[2] = By[2] + ext_By[2]; 
  tot_Bz[2] = Bz[2] + ext_Bz[2]; 
  tot_Bx[3] = Bx[3] + ext_Bx[3]; 
  tot_By[3] = By[3] + ext_By[3]; 
  tot_Bz[3] = Bz[3] + ext_Bz[3]; 

  double ext_force_x[GKYL_MAX_SPECIES][4]; 
  double ext_force_y[GKYL_MAX_SPECIES][4]; 
  double ext_force_z[GKYL_MAX_SPECIES][4]; 
  for (int i = 0; i < num_species; ++i) { 
    ext_force_x[i][0] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ex[0] + app_accel_x[i][0]); 
    ext_force_y[i][0] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ey[0] + app_accel_y[i][0]); 
    ext_force_z[i][0] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ez[0] + app_accel_z[i][0]); 

    ext_force_x[i][1] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ex[1] + app_accel_x[i][1]); 
    ext_force_y[i][1] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ey[1] + app_accel_y[i][1]); 
    ext_force_z[i][1] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ez[1] + app_accel_z[i][1]); 

    ext_force_x[i][2] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ex[2] + app_accel_x[i][2]); 
    ext_force_y[i][2] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ey[2] + app_accel_y[i][2]); 
    ext_force_z[i][2] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ez[2] + app_accel_z[i][2]); 

    ext_force_x[i][3] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ex[3] + app_accel_x[i][3]); 
    ext_force_y[i][3] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ey[3] + app_accel_y[i][3]); 
    ext_force_z[i][3] = 0.5*dt*qbym[i]*(qbym[i]*ext_Ez[3] + app_accel_z[i][3]); 

  } 

  // Set RHS for momentum equations, including solution at known time-step and external forces. 
  for (int s = 0; s < num_species; ++s) { 

    gkyl_mat_set(&rhs, 0 + s*(12), 0, qbym[s]*rhoux[s][0] + 0.5*ext_force_x[s][3]*rho[s][3]+0.5*ext_force_x[s][2]*rho[s][2]+0.5*ext_force_x[s][1]*rho[s][1]+0.5*ext_force_x[s][0]*rho[s][0]); 
    gkyl_mat_set(&rhs, 4 + s*(12), 0, qbym[s]*rhouy[s][0] + 0.5*ext_force_y[s][3]*rho[s][3]+0.5*ext_force_y[s][2]*rho[s][2]+0.5*ext_force_y[s][1]*rho[s][1]+0.5*ext_force_y[s][0]*rho[s][0]); 
    gkyl_mat_set(&rhs, 8 + s*(12), 0, qbym[s]*rhouz[s][0] + 0.5*ext_force_z[s][3]*rho[s][3]+0.5*ext_force_z[s][2]*rho[s][2]+0.5*ext_force_z[s][1]*rho[s][1]+0.5*ext_force_z[s][0]*rho[s][0]); 

    gkyl_mat_set(&rhs, 1 + s*(12), 0, qbym[s]*rhoux[s][1] + 0.5*ext_force_x[s][2]*rho[s][3]+0.5*rho[s][2]*ext_force_x[s][3]+0.5*ext_force_x[s][0]*rho[s][1]+0.5*rho[s][0]*ext_force_x[s][1]); 
    gkyl_mat_set(&rhs, 5 + s*(12), 0, qbym[s]*rhouy[s][1] + 0.5*ext_force_y[s][2]*rho[s][3]+0.5*rho[s][2]*ext_force_y[s][3]+0.5*ext_force_y[s][0]*rho[s][1]+0.5*rho[s][0]*ext_force_y[s][1]); 
    gkyl_mat_set(&rhs, 9 + s*(12), 0, qbym[s]*rhouz[s][1] + 0.5*ext_force_z[s][2]*rho[s][3]+0.5*rho[s][2]*ext_force_z[s][3]+0.5*ext_force_z[s][0]*rho[s][1]+0.5*rho[s][0]*ext_force_z[s][1]); 

    gkyl_mat_set(&rhs, 2 + s*(12), 0, qbym[s]*rhoux[s][2] + 0.5*ext_force_x[s][1]*rho[s][3]+0.5*rho[s][1]*ext_force_x[s][3]+0.5*ext_force_x[s][0]*rho[s][2]+0.5*rho[s][0]*ext_force_x[s][2]); 
    gkyl_mat_set(&rhs, 6 + s*(12), 0, qbym[s]*rhouy[s][2] + 0.5*ext_force_y[s][1]*rho[s][3]+0.5*rho[s][1]*ext_force_y[s][3]+0.5*ext_force_y[s][0]*rho[s][2]+0.5*rho[s][0]*ext_force_y[s][2]); 
    gkyl_mat_set(&rhs, 10 + s*(12), 0, qbym[s]*rhouz[s][2] + 0.5*ext_force_z[s][1]*rho[s][3]+0.5*rho[s][1]*ext_force_z[s][3]+0.5*ext_force_z[s][0]*rho[s][2]+0.5*rho[s][0]*ext_force_z[s][2]); 

    gkyl_mat_set(&rhs, 3 + s*(12), 0, qbym[s]*rhoux[s][3] + 0.5*ext_force_x[s][0]*rho[s][3]+0.5*rho[s][0]*ext_force_x[s][3]+0.5*ext_force_x[s][1]*rho[s][2]+0.5*rho[s][1]*ext_force_x[s][2]); 
    gkyl_mat_set(&rhs, 7 + s*(12), 0, qbym[s]*rhouy[s][3] + 0.5*ext_force_y[s][0]*rho[s][3]+0.5*rho[s][0]*ext_force_y[s][3]+0.5*ext_force_y[s][1]*rho[s][2]+0.5*rho[s][1]*ext_force_y[s][2]); 
    gkyl_mat_set(&rhs, 11 + s*(12), 0, qbym[s]*rhouz[s][3] + 0.5*ext_force_z[s][0]*rho[s][3]+0.5*rho[s][0]*ext_force_z[s][3]+0.5*ext_force_z[s][1]*rho[s][2]+0.5*rho[s][1]*ext_force_z[s][2]); 

  } 

  // Set RHS for Ampere's Law, including solution at known time-step and applied currents. 
  gkyl_mat_set(&rhs, 0 + num_species*(12), 0, epsilon0*Ex[0] - 0.5*dt*app_curr_x[0]); 
  gkyl_mat_set(&rhs, 4 + num_species*(12), 0, epsilon0*Ey[0] - 0.5*dt*app_curr_y[0]); 
  gkyl_mat_set(&rhs, 8 + num_species*(12), 0, epsilon0*Ez[0] - 0.5*dt*app_curr_z[0]); 

  gkyl_mat_set(&rhs, 1 + num_species*(12), 0, epsilon0*Ex[1] - 0.5*dt*app_curr_x[1]); 
  gkyl_mat_set(&rhs, 5 + num_species*(12), 0, epsilon0*Ey[1] - 0.5*dt*app_curr_y[1]); 
  gkyl_mat_set(&rhs, 9 + num_species*(12), 0, epsilon0*Ez[1] - 0.5*dt*app_curr_z[1]); 

  gkyl_mat_set(&rhs, 2 + num_species*(12), 0, epsilon0*Ex[2] - 0.5*dt*app_curr_x[2]); 
  gkyl_mat_set(&rhs, 6 + num_species*(12), 0, epsilon0*Ey[2] - 0.5*dt*app_curr_y[2]); 
  gkyl_mat_set(&rhs, 10 + num_species*(12), 0, epsilon0*Ez[2] - 0.5*dt*app_curr_z[2]); 

  gkyl_mat_set(&rhs, 3 + num_species*(12), 0, epsilon0*Ex[3] - 0.5*dt*app_curr_x[3]); 
  gkyl_mat_set(&rhs, 7 + num_species*(12), 0, epsilon0*Ey[3] - 0.5*dt*app_curr_y[3]); 
  gkyl_mat_set(&rhs, 11 + num_species*(12), 0, epsilon0*Ez[3] - 0.5*dt*app_curr_z[3]); 


  // Construct LHS. 
  // For momentum equation: J_s^{n+1} - 0.5*dt*(q_s^2/m_s^2*rho_s^n*E^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
  // For Ampere's Law: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
  for (int s = 0; s < num_species; ++s) { 
 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs, 0 + s*(12), 0 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 4 + s*(12), 4 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 8 + s*(12), 8 + s*(12), 1.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 0 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 4 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 8 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 4 + s*(12), 8 + s*(12), B_field_fac*(0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 4 + s*(12), -B_field_fac*(0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 8 + s*(12), 0 + s*(12), B_field_fac*(0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 0 + s*(12), 8 + s*(12), -B_field_fac*(0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 4 + s*(12), B_field_fac*(0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 0 + s*(12), -B_field_fac*(0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(12), 0 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(12), 4 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(12), 8 + s*(12), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 1 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 4 + s*(12), 5 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(12), 9 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 1 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 5 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 9 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 4 + s*(12), 9 + s*(12), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 5 + s*(12), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 8 + s*(12), 1 + s*(12), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 0 + s*(12), 9 + s*(12), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 5 + s*(12), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 1 + s*(12), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(12), 1 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(12), 5 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(12), 9 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 2 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 4 + s*(12), 6 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(12), 10 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 2 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 6 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 10 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 4 + s*(12), 10 + s*(12), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 6 + s*(12), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 8 + s*(12), 2 + s*(12), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 0 + s*(12), 10 + s*(12), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 6 + s*(12), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 2 + s*(12), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(12), 2 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(12), 6 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(12), 10 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 3 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 4 + s*(12), 7 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(12), 11 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 3 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 7 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 11 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 4 + s*(12), 11 + s*(12), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 8 + s*(12), 7 + s*(12), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 8 + s*(12), 3 + s*(12), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 0 + s*(12), 11 + s*(12), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 0 + s*(12), 7 + s*(12), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 4 + s*(12), 3 + s*(12), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(12), 3 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(12), 7 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(12), 11 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 0 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 5 + s*(12), 4 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(12), 8 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 0 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 4 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 8 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 5 + s*(12), 8 + s*(12), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 4 + s*(12), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 9 + s*(12), 0 + s*(12), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 1 + s*(12), 8 + s*(12), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 4 + s*(12), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 0 + s*(12), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(12), 0 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(12), 4 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(12), 8 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 1 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 5 + s*(12), 5 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 9 + s*(12), 9 + s*(12), 1.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 1 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 5 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 9 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 5 + s*(12), 9 + s*(12), B_field_fac*(0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 5 + s*(12), -B_field_fac*(0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 9 + s*(12), 1 + s*(12), B_field_fac*(0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 1 + s*(12), 9 + s*(12), -B_field_fac*(0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 5 + s*(12), B_field_fac*(0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 1 + s*(12), -B_field_fac*(0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(12), 1 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(12), 5 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(12), 9 + s*(12), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 2 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 5 + s*(12), 6 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(12), 10 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 2 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 6 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 10 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 5 + s*(12), 10 + s*(12), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 6 + s*(12), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 9 + s*(12), 2 + s*(12), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 1 + s*(12), 10 + s*(12), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 6 + s*(12), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 2 + s*(12), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(12), 2 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(12), 6 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(12), 10 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 3 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 5 + s*(12), 7 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(12), 11 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 3 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 7 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 11 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 5 + s*(12), 11 + s*(12), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 9 + s*(12), 7 + s*(12), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 9 + s*(12), 3 + s*(12), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 1 + s*(12), 11 + s*(12), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 1 + s*(12), 7 + s*(12), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 5 + s*(12), 3 + s*(12), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(12), 3 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(12), 7 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(12), 11 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 0 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 6 + s*(12), 4 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(12), 8 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 0 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 4 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 8 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 6 + s*(12), 8 + s*(12), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 4 + s*(12), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 10 + s*(12), 0 + s*(12), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 2 + s*(12), 8 + s*(12), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 4 + s*(12), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 0 + s*(12), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(12), 0 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 6 + num_species*(12), 4 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(12), 8 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 1 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 6 + s*(12), 5 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(12), 9 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 1 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 5 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 9 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 6 + s*(12), 9 + s*(12), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 5 + s*(12), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 10 + s*(12), 1 + s*(12), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 2 + s*(12), 9 + s*(12), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 5 + s*(12), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 1 + s*(12), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(12), 1 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 6 + num_species*(12), 5 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(12), 9 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 2 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 6 + s*(12), 6 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 10 + s*(12), 10 + s*(12), 1.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 2 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 6 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 10 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 6 + s*(12), 10 + s*(12), B_field_fac*(0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 6 + s*(12), -B_field_fac*(0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 10 + s*(12), 2 + s*(12), B_field_fac*(0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 2 + s*(12), 10 + s*(12), -B_field_fac*(0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 6 + s*(12), B_field_fac*(0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 2 + s*(12), -B_field_fac*(0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(12), 2 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 6 + num_species*(12), 6 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(12), 10 + s*(12), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 3 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 6 + s*(12), 7 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(12), 11 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 3 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 7 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 11 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 6 + s*(12), 11 + s*(12), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 10 + s*(12), 7 + s*(12), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 10 + s*(12), 3 + s*(12), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 2 + s*(12), 11 + s*(12), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 2 + s*(12), 7 + s*(12), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 6 + s*(12), 3 + s*(12), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(12), 3 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 6 + num_species*(12), 7 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(12), 11 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 0 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 7 + s*(12), 4 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(12), 8 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 0 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 4 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 8 + num_species*(12), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 7 + s*(12), 8 + s*(12), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 4 + s*(12), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 11 + s*(12), 0 + s*(12), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(12), 8 + s*(12), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 4 + s*(12), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 0 + s*(12), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(12), 0 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 7 + num_species*(12), 4 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(12), 8 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 1 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 7 + s*(12), 5 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(12), 9 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 1 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 5 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 9 + num_species*(12), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 7 + s*(12), 9 + s*(12), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 5 + s*(12), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 11 + s*(12), 1 + s*(12), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 3 + s*(12), 9 + s*(12), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 5 + s*(12), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 1 + s*(12), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(12), 1 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 7 + num_species*(12), 5 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(12), 9 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 2 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 7 + s*(12), 6 + s*(12), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(12), 10 + s*(12), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 2 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 6 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 10 + num_species*(12), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 7 + s*(12), 10 + s*(12), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 6 + s*(12), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 11 + s*(12), 2 + s*(12), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 3 + s*(12), 10 + s*(12), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 6 + s*(12), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 2 + s*(12), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(12), 2 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 7 + num_species*(12), 6 + s*(12), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(12), 10 + s*(12), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 3 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 7 + s*(12), 7 + s*(12), 1.0); 
    gkyl_mat_set(&lhs, 11 + s*(12), 11 + s*(12), 1.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 3 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 7 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 11 + num_species*(12), E_field_fac*(0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 7 + s*(12), 11 + s*(12), B_field_fac*(0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 11 + s*(12), 7 + s*(12), -B_field_fac*(0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 11 + s*(12), 3 + s*(12), B_field_fac*(0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 3 + s*(12), 11 + s*(12), -B_field_fac*(0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 3 + s*(12), 7 + s*(12), B_field_fac*(0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 7 + s*(12), 3 + s*(12), -B_field_fac*(0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(12), 3 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 7 + num_species*(12), 7 + s*(12), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(12), 11 + s*(12), 0.5*dt*(1.0)); 
 
  } 
  gkyl_mat_set(&lhs, 0 + num_species*(12), 0 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(12), 4 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(12), 8 + num_species*(12), 1.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(12), 1 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(12), 5 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(12), 9 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(12), 2 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(12), 6 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(12), 10 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(12), 3 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(12), 7 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(12), 11 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(12), 0 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(12), 4 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(12), 8 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(12), 1 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(12), 5 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(12), 9 + num_species*(12), 1.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(12), 2 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(12), 6 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(12), 10 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(12), 3 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(12), 7 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(12), 11 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(12), 0 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 6 + num_species*(12), 4 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(12), 8 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(12), 1 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 6 + num_species*(12), 5 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(12), 9 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(12), 2 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 6 + num_species*(12), 6 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(12), 10 + num_species*(12), 1.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(12), 3 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 6 + num_species*(12), 7 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(12), 11 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(12), 0 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 7 + num_species*(12), 4 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(12), 8 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(12), 1 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 7 + num_species*(12), 5 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(12), 9 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(12), 2 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 7 + num_species*(12), 6 + num_species*(12), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(12), 10 + num_species*(12), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(12), 3 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 7 + num_species*(12), 7 + num_species*(12), 1.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(12), 11 + num_species*(12), 1.0); 
 
} 
