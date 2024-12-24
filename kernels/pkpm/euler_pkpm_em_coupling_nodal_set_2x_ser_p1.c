#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:             integer to indicate which matrix being fetched. 
  // num_species:       number of species being evolved (number of momentum equations). 
  // qbym:              charge/mass ratio for each species. 
  // epsilon0:          permittivity of free space. 
  // pkpm_field_static: boolean for whether or not the self-consistent field is static. 
  // dt:                size of the time step. 
  // A:                 preallocated LHS matrix. 
  // rhs:               preallocated RHS vector. 
  // app_accel:         Applied accelerations (external forces).
  // ext_em:            Externally applied EM fields.
  // app_current:       Applied external currents.
  // vlasov_pkpm_moms:  [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:        [rho ux, rho uy, rho uz], PKPM Fluid momentum input state vector.
  // em:                [Ex, Ey, Ez, Bx, By, Bz], EM input state vector.

  double rho[GKYL_MAX_SPECIES][4]; 
  double rhoux[GKYL_MAX_SPECIES][4]; 
  double rhouy[GKYL_MAX_SPECIES][4]; 
  double rhouz[GKYL_MAX_SPECIES][4]; 

  double app_accel_x[GKYL_MAX_SPECIES][4]; 
  double app_accel_y[GKYL_MAX_SPECIES][4]; 
  double app_accel_z[GKYL_MAX_SPECIES][4]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_euler_pkpm = euler_pkpm[i]; 
    const double *inp_vlasov_pkpm_moms = vlasov_pkpm_moms[i]; 
    const double *inp_app_accel = app_accel[i]; 

    rho[i][0] = inp_vlasov_pkpm_moms[0]; 
    rhoux[i][0] = inp_euler_pkpm[0]; 
    rhouy[i][0] = inp_euler_pkpm[4]; 
    rhouz[i][0] = inp_euler_pkpm[8]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[4]; 
    app_accel_z[i][0] = inp_app_accel[8]; 

    rho[i][1] = inp_vlasov_pkpm_moms[1]; 
    rhoux[i][1] = inp_euler_pkpm[1]; 
    rhouy[i][1] = inp_euler_pkpm[5]; 
    rhouz[i][1] = inp_euler_pkpm[9]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[5]; 
    app_accel_z[i][1] = inp_app_accel[9]; 

    rho[i][2] = inp_vlasov_pkpm_moms[2]; 
    rhoux[i][2] = inp_euler_pkpm[2]; 
    rhouy[i][2] = inp_euler_pkpm[6]; 
    rhouz[i][2] = inp_euler_pkpm[10]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[6]; 
    app_accel_z[i][2] = inp_app_accel[10]; 

    rho[i][3] = inp_vlasov_pkpm_moms[3]; 
    rhoux[i][3] = inp_euler_pkpm[3]; 
    rhouy[i][3] = inp_euler_pkpm[7]; 
    rhouz[i][3] = inp_euler_pkpm[11]; 

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

  double rho_nodal[GKYL_MAX_SPECIES][4]; 
  double rhoux_nodal[GKYL_MAX_SPECIES][4]; 
  double rhouy_nodal[GKYL_MAX_SPECIES][4]; 
  double rhouz_nodal[GKYL_MAX_SPECIES][4]; 

  double app_accel_x_nodal[GKYL_MAX_SPECIES][4]; 
  double app_accel_y_nodal[GKYL_MAX_SPECIES][4]; 
  double app_accel_z_nodal[GKYL_MAX_SPECIES][4]; 

  double Ex_nodal[4]; 
  double Ey_nodal[4]; 
  double Ez_nodal[4]; 

  double ext_Ex_nodal[4]; 
  double ext_Ey_nodal[4]; 
  double ext_Ez_nodal[4]; 

  double app_curr_x_nodal[4]; 
  double app_curr_y_nodal[4]; 
  double app_curr_z_nodal[4]; 

  double tot_Bx_nodal[4]; 
  double tot_By_nodal[4]; 
  double tot_Bz_nodal[4]; 

  // Project modal expansions onto nodal bases. 
  for (int s = 0; s < num_species; ++s) { 
    modal_to_quad_nodal_2d_ser_p1(rho[s], rho_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(rhoux[s], rhoux_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(rhouy[s], rhouy_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(rhouz[s], rhouz_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(app_accel_x[s], app_accel_x_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(app_accel_y[s], app_accel_y_nodal[s]);
    modal_to_quad_nodal_2d_ser_p1(app_accel_z[s], app_accel_z_nodal[s]);
  } 
  modal_to_quad_nodal_2d_ser_p1(Ex, Ex_nodal);
  modal_to_quad_nodal_2d_ser_p1(Ey, Ey_nodal);
  modal_to_quad_nodal_2d_ser_p1(Ez, Ez_nodal);
  modal_to_quad_nodal_2d_ser_p1(ext_Ex, ext_Ex_nodal);
  modal_to_quad_nodal_2d_ser_p1(ext_Ey, ext_Ey_nodal);
  modal_to_quad_nodal_2d_ser_p1(ext_Ez, ext_Ez_nodal);
  modal_to_quad_nodal_2d_ser_p1(app_curr_x, app_curr_x_nodal);
  modal_to_quad_nodal_2d_ser_p1(app_curr_y, app_curr_y_nodal);
  modal_to_quad_nodal_2d_ser_p1(app_curr_z, app_curr_z_nodal);
  modal_to_quad_nodal_2d_ser_p1(tot_Bx, tot_Bx_nodal);
  modal_to_quad_nodal_2d_ser_p1(tot_By, tot_By_nodal);
  modal_to_quad_nodal_2d_ser_p1(tot_Bz, tot_Bz_nodal);
 
  struct gkyl_mat lhs0 = gkyl_nmat_get(A_n, count+0); 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs_n, count+0); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs0, 0.0); gkyl_mat_clear(&rhs0, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs0, 0 + s*3, 0, qbym[s]*rhoux_nodal[s][0] + 0.5*dt*qbym[s]*rho_nodal[s][0]*(qbym[s]*ext_Ex_nodal[0] + app_accel_x_nodal[s][0])); 
    gkyl_mat_set(&rhs0, 1 + s*3, 0, qbym[s]*rhouy_nodal[s][0] + 0.5*dt*qbym[s]*rho_nodal[s][0]*(qbym[s]*ext_Ey_nodal[0] + app_accel_y_nodal[s][0])); 
    gkyl_mat_set(&rhs0, 2 + s*3, 0, qbym[s]*rhouz_nodal[s][0] + 0.5*dt*qbym[s]*rho_nodal[s][0]*(qbym[s]*ext_Ez_nodal[0] + app_accel_z_nodal[s][0])); 
 
    // Set LHS matrix for flow velocity equation: J_s^{n+1} - 0.5*dt*((q_s/m_s)^2 rho_s^n)/epsilon0*(epsilon0*E)^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]*rho_nodal[s][0]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs0, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs0, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs0, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs0, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs0, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs0, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs0, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[0]); 
    gkyl_mat_set(&lhs0, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[0]); 
 
    gkyl_mat_set(&lhs0, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[0]); 
    gkyl_mat_set(&lhs0, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[0]); 
 
    gkyl_mat_set(&lhs0, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[0]); 
    gkyl_mat_set(&lhs0, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[0]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs0, 0 + num_species*3, 0 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs0, 1 + num_species*3, 1 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs0, 2 + num_species*3, 2 + s*3, 0.5*dt); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs0, 0 + num_species*3, 0, epsilon0*Ex_nodal[0] - 0.5*dt*app_curr_x_nodal[0]); 
  gkyl_mat_set(&rhs0, 1 + num_species*3, 0, epsilon0*Ey_nodal[0] - 0.5*dt*app_curr_y_nodal[0]); 
  gkyl_mat_set(&rhs0, 2 + num_species*3, 0, epsilon0*Ez_nodal[0] - 0.5*dt*app_curr_z_nodal[0]); 
 
  gkyl_mat_set(&lhs0, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs0, 1 + num_species*3, 1 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs0, 2 + num_species*3, 2 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs1 = gkyl_nmat_get(A_n, count+1); 
  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs_n, count+1); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs1, 0.0); gkyl_mat_clear(&rhs1, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs1, 0 + s*3, 0, qbym[s]*rhoux_nodal[s][1] + 0.5*dt*qbym[s]*rho_nodal[s][1]*(qbym[s]*ext_Ex_nodal[1] + app_accel_x_nodal[s][1])); 
    gkyl_mat_set(&rhs1, 1 + s*3, 0, qbym[s]*rhouy_nodal[s][1] + 0.5*dt*qbym[s]*rho_nodal[s][1]*(qbym[s]*ext_Ey_nodal[1] + app_accel_y_nodal[s][1])); 
    gkyl_mat_set(&rhs1, 2 + s*3, 0, qbym[s]*rhouz_nodal[s][1] + 0.5*dt*qbym[s]*rho_nodal[s][1]*(qbym[s]*ext_Ez_nodal[1] + app_accel_z_nodal[s][1])); 
 
    // Set LHS matrix for flow velocity equation: J_s^{n+1} - 0.5*dt*((q_s/m_s)^2 rho_s^n)/epsilon0*(epsilon0*E)^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]*rho_nodal[s][1]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs1, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs1, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs1, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs1, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs1, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs1, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs1, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[1]); 
    gkyl_mat_set(&lhs1, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[1]); 
 
    gkyl_mat_set(&lhs1, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[1]); 
    gkyl_mat_set(&lhs1, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[1]); 
 
    gkyl_mat_set(&lhs1, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[1]); 
    gkyl_mat_set(&lhs1, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[1]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs1, 0 + num_species*3, 0 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs1, 1 + num_species*3, 1 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs1, 2 + num_species*3, 2 + s*3, 0.5*dt); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs1, 0 + num_species*3, 0, epsilon0*Ex_nodal[1] - 0.5*dt*app_curr_x_nodal[1]); 
  gkyl_mat_set(&rhs1, 1 + num_species*3, 0, epsilon0*Ey_nodal[1] - 0.5*dt*app_curr_y_nodal[1]); 
  gkyl_mat_set(&rhs1, 2 + num_species*3, 0, epsilon0*Ez_nodal[1] - 0.5*dt*app_curr_z_nodal[1]); 
 
  gkyl_mat_set(&lhs1, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs1, 1 + num_species*3, 1 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs1, 2 + num_species*3, 2 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs2 = gkyl_nmat_get(A_n, count+2); 
  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs_n, count+2); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs2, 0.0); gkyl_mat_clear(&rhs2, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs2, 0 + s*3, 0, qbym[s]*rhoux_nodal[s][2] + 0.5*dt*qbym[s]*rho_nodal[s][2]*(qbym[s]*ext_Ex_nodal[2] + app_accel_x_nodal[s][2])); 
    gkyl_mat_set(&rhs2, 1 + s*3, 0, qbym[s]*rhouy_nodal[s][2] + 0.5*dt*qbym[s]*rho_nodal[s][2]*(qbym[s]*ext_Ey_nodal[2] + app_accel_y_nodal[s][2])); 
    gkyl_mat_set(&rhs2, 2 + s*3, 0, qbym[s]*rhouz_nodal[s][2] + 0.5*dt*qbym[s]*rho_nodal[s][2]*(qbym[s]*ext_Ez_nodal[2] + app_accel_z_nodal[s][2])); 
 
    // Set LHS matrix for flow velocity equation: J_s^{n+1} - 0.5*dt*((q_s/m_s)^2 rho_s^n)/epsilon0*(epsilon0*E)^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]*rho_nodal[s][2]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs2, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs2, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs2, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs2, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs2, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs2, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs2, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[2]); 
    gkyl_mat_set(&lhs2, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[2]); 
 
    gkyl_mat_set(&lhs2, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[2]); 
    gkyl_mat_set(&lhs2, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[2]); 
 
    gkyl_mat_set(&lhs2, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[2]); 
    gkyl_mat_set(&lhs2, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[2]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs2, 0 + num_species*3, 0 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs2, 1 + num_species*3, 1 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs2, 2 + num_species*3, 2 + s*3, 0.5*dt); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs2, 0 + num_species*3, 0, epsilon0*Ex_nodal[2] - 0.5*dt*app_curr_x_nodal[2]); 
  gkyl_mat_set(&rhs2, 1 + num_species*3, 0, epsilon0*Ey_nodal[2] - 0.5*dt*app_curr_y_nodal[2]); 
  gkyl_mat_set(&rhs2, 2 + num_species*3, 0, epsilon0*Ez_nodal[2] - 0.5*dt*app_curr_z_nodal[2]); 
 
  gkyl_mat_set(&lhs2, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs2, 1 + num_species*3, 1 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs2, 2 + num_species*3, 2 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs3 = gkyl_nmat_get(A_n, count+3); 
  struct gkyl_mat rhs3 = gkyl_nmat_get(rhs_n, count+3); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs3, 0.0); gkyl_mat_clear(&rhs3, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs3, 0 + s*3, 0, qbym[s]*rhoux_nodal[s][3] + 0.5*dt*qbym[s]*rho_nodal[s][3]*(qbym[s]*ext_Ex_nodal[3] + app_accel_x_nodal[s][3])); 
    gkyl_mat_set(&rhs3, 1 + s*3, 0, qbym[s]*rhouy_nodal[s][3] + 0.5*dt*qbym[s]*rho_nodal[s][3]*(qbym[s]*ext_Ey_nodal[3] + app_accel_y_nodal[s][3])); 
    gkyl_mat_set(&rhs3, 2 + s*3, 0, qbym[s]*rhouz_nodal[s][3] + 0.5*dt*qbym[s]*rho_nodal[s][3]*(qbym[s]*ext_Ez_nodal[3] + app_accel_z_nodal[s][3])); 
 
    // Set LHS matrix for flow velocity equation: J_s^{n+1} - 0.5*dt*((q_s/m_s)^2 rho_s^n)/epsilon0*(epsilon0*E)^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]*rho_nodal[s][3]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs3, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs3, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs3, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs3, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs3, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs3, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs3, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[3]); 
    gkyl_mat_set(&lhs3, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[3]); 
 
    gkyl_mat_set(&lhs3, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[3]); 
    gkyl_mat_set(&lhs3, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[3]); 
 
    gkyl_mat_set(&lhs3, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[3]); 
    gkyl_mat_set(&lhs3, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[3]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs3, 0 + num_species*3, 0 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs3, 1 + num_species*3, 1 + s*3, 0.5*dt); 
      gkyl_mat_set(&lhs3, 2 + num_species*3, 2 + s*3, 0.5*dt); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs3, 0 + num_species*3, 0, epsilon0*Ex_nodal[3] - 0.5*dt*app_curr_x_nodal[3]); 
  gkyl_mat_set(&rhs3, 1 + num_species*3, 0, epsilon0*Ey_nodal[3] - 0.5*dt*app_curr_y_nodal[3]); 
  gkyl_mat_set(&rhs3, 2 + num_species*3, 0, epsilon0*Ez_nodal[3] - 0.5*dt*app_curr_z_nodal[3]); 
 
  gkyl_mat_set(&lhs3, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs3, 1 + num_species*3, 1 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs3, 2 + num_species*3, 2 + num_species*3, 1.0); 
 
} 
