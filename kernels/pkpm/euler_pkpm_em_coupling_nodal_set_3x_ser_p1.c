#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em) 
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
  // pkpm_u:            [ux, uy, uz], Input flow velocity.
  // em:                [Ex, Ey, Ez, Bx, By, Bz], EM input state vector.

  double rho[GKYL_MAX_SPECIES][8]; 
  double ux[GKYL_MAX_SPECIES][8]; 
  double uy[GKYL_MAX_SPECIES][8]; 
  double uz[GKYL_MAX_SPECIES][8]; 

  double app_accel_x[GKYL_MAX_SPECIES][8]; 
  double app_accel_y[GKYL_MAX_SPECIES][8]; 
  double app_accel_z[GKYL_MAX_SPECIES][8]; 

  for (int i = 0; i < num_species; ++i) { 
    const double *inp_u = pkpm_u[i]; 
    const double *inp_app_accel = app_accel[i]; 
    const double *inp_vlasov_pkpm_moms = vlasov_pkpm_moms[i]; 

    rho[i][0] = inp_vlasov_pkpm_moms[0]; 
    ux[i][0] = inp_u[0]; 
    uy[i][0] = inp_u[8]; 
    uz[i][0] = inp_u[16]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[8]; 
    app_accel_z[i][0] = inp_app_accel[16]; 

    rho[i][1] = inp_vlasov_pkpm_moms[1]; 
    ux[i][1] = inp_u[1]; 
    uy[i][1] = inp_u[9]; 
    uz[i][1] = inp_u[17]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[9]; 
    app_accel_z[i][1] = inp_app_accel[17]; 

    rho[i][2] = inp_vlasov_pkpm_moms[2]; 
    ux[i][2] = inp_u[2]; 
    uy[i][2] = inp_u[10]; 
    uz[i][2] = inp_u[18]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[10]; 
    app_accel_z[i][2] = inp_app_accel[18]; 

    rho[i][3] = inp_vlasov_pkpm_moms[3]; 
    ux[i][3] = inp_u[3]; 
    uy[i][3] = inp_u[11]; 
    uz[i][3] = inp_u[19]; 

    app_accel_x[i][3] = inp_app_accel[3]; 
    app_accel_y[i][3] = inp_app_accel[11]; 
    app_accel_z[i][3] = inp_app_accel[19]; 

    rho[i][4] = inp_vlasov_pkpm_moms[4]; 
    ux[i][4] = inp_u[4]; 
    uy[i][4] = inp_u[12]; 
    uz[i][4] = inp_u[20]; 

    app_accel_x[i][4] = inp_app_accel[4]; 
    app_accel_y[i][4] = inp_app_accel[12]; 
    app_accel_z[i][4] = inp_app_accel[20]; 

    rho[i][5] = inp_vlasov_pkpm_moms[5]; 
    ux[i][5] = inp_u[5]; 
    uy[i][5] = inp_u[13]; 
    uz[i][5] = inp_u[21]; 

    app_accel_x[i][5] = inp_app_accel[5]; 
    app_accel_y[i][5] = inp_app_accel[13]; 
    app_accel_z[i][5] = inp_app_accel[21]; 

    rho[i][6] = inp_vlasov_pkpm_moms[6]; 
    ux[i][6] = inp_u[6]; 
    uy[i][6] = inp_u[14]; 
    uz[i][6] = inp_u[22]; 

    app_accel_x[i][6] = inp_app_accel[6]; 
    app_accel_y[i][6] = inp_app_accel[14]; 
    app_accel_z[i][6] = inp_app_accel[22]; 

    rho[i][7] = inp_vlasov_pkpm_moms[7]; 
    ux[i][7] = inp_u[7]; 
    uy[i][7] = inp_u[15]; 
    uz[i][7] = inp_u[23]; 

    app_accel_x[i][7] = inp_app_accel[7]; 
    app_accel_y[i][7] = inp_app_accel[15]; 
    app_accel_z[i][7] = inp_app_accel[23]; 

  } 

  double *Ex = &em[0]; 
  double *Ey = &em[8]; 
  double *Ez = &em[16]; 
  double *Bx = &em[24]; 
  double *By = &em[32]; 
  double *Bz = &em[40]; 

  const double *ext_Ex = &ext_em[0]; 
  const double *ext_Ey = &ext_em[8]; 
  const double *ext_Ez = &ext_em[16]; 
  const double *ext_Bx = &ext_em[24]; 
  const double *ext_By = &ext_em[32]; 
  const double *ext_Bz = &ext_em[40]; 

  const double *app_curr_x = &app_current[0]; 
  const double *app_curr_y = &app_current[8]; 
  const double *app_curr_z = &app_current[16]; 

  double tot_Bx[8]; 
  double tot_By[8]; 
  double tot_Bz[8]; 
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
  tot_Bx[4] = Bx[4] + ext_Bx[4]; 
  tot_By[4] = By[4] + ext_By[4]; 
  tot_Bz[4] = Bz[4] + ext_Bz[4]; 
  tot_Bx[5] = Bx[5] + ext_Bx[5]; 
  tot_By[5] = By[5] + ext_By[5]; 
  tot_Bz[5] = Bz[5] + ext_Bz[5]; 
  tot_Bx[6] = Bx[6] + ext_Bx[6]; 
  tot_By[6] = By[6] + ext_By[6]; 
  tot_Bz[6] = Bz[6] + ext_Bz[6]; 
  tot_Bx[7] = Bx[7] + ext_Bx[7]; 
  tot_By[7] = By[7] + ext_By[7]; 
  tot_Bz[7] = Bz[7] + ext_Bz[7]; 

  double rho_nodal[GKYL_MAX_SPECIES][8]; 
  double ux_nodal[GKYL_MAX_SPECIES][8]; 
  double uy_nodal[GKYL_MAX_SPECIES][8]; 
  double uz_nodal[GKYL_MAX_SPECIES][8]; 

  double app_accel_x_nodal[GKYL_MAX_SPECIES][8]; 
  double app_accel_y_nodal[GKYL_MAX_SPECIES][8]; 
  double app_accel_z_nodal[GKYL_MAX_SPECIES][8]; 

  double Ex_nodal[8]; 
  double Ey_nodal[8]; 
  double Ez_nodal[8]; 

  double ext_Ex_nodal[8]; 
  double ext_Ey_nodal[8]; 
  double ext_Ez_nodal[8]; 

  double app_curr_x_nodal[8]; 
  double app_curr_y_nodal[8]; 
  double app_curr_z_nodal[8]; 

  double tot_Bx_nodal[8]; 
  double tot_By_nodal[8]; 
  double tot_Bz_nodal[8]; 

  // Project modal expansions onto nodal bases. 
  for (int s = 0; s < num_species; ++s) { 
    modal_to_quad_nodal_3d_ser_p1(rho[s], rho_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(ux[s], ux_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(uy[s], uy_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(uz[s], uz_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(app_accel_x[s], app_accel_x_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(app_accel_y[s], app_accel_y_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(app_accel_z[s], app_accel_z_nodal[s]);
  } 
  modal_to_quad_nodal_3d_ser_p1(Ex, Ex_nodal);
  modal_to_quad_nodal_3d_ser_p1(Ey, Ey_nodal);
  modal_to_quad_nodal_3d_ser_p1(Ez, Ez_nodal);
  modal_to_quad_nodal_3d_ser_p1(ext_Ex, ext_Ex_nodal);
  modal_to_quad_nodal_3d_ser_p1(ext_Ey, ext_Ey_nodal);
  modal_to_quad_nodal_3d_ser_p1(ext_Ez, ext_Ez_nodal);
  modal_to_quad_nodal_3d_ser_p1(app_curr_x, app_curr_x_nodal);
  modal_to_quad_nodal_3d_ser_p1(app_curr_y, app_curr_y_nodal);
  modal_to_quad_nodal_3d_ser_p1(app_curr_z, app_curr_z_nodal);
  modal_to_quad_nodal_3d_ser_p1(tot_Bx, tot_Bx_nodal);
  modal_to_quad_nodal_3d_ser_p1(tot_By, tot_By_nodal);
  modal_to_quad_nodal_3d_ser_p1(tot_Bz, tot_Bz_nodal);
 
  struct gkyl_mat lhs0 = gkyl_nmat_get(A_n, count+0); 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs_n, count+0); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs0, 0.0); gkyl_mat_clear(&rhs0, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs0, 0 + s*3, 0, ux_nodal[s][0] + 0.5*dt*(qbym[s]*ext_Ex_nodal[0] + app_accel_x_nodal[s][0])); 
    gkyl_mat_set(&rhs0, 1 + s*3, 0, uy_nodal[s][0] + 0.5*dt*(qbym[s]*ext_Ey_nodal[0] + app_accel_y_nodal[s][0])); 
    gkyl_mat_set(&rhs0, 2 + s*3, 0, uz_nodal[s][0] + 0.5*dt*(qbym[s]*ext_Ez_nodal[0] + app_accel_z_nodal[s][0])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
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
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs0, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][0]); 
      gkyl_mat_set(&lhs0, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][0]); 
      gkyl_mat_set(&lhs0, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][0]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs0, 0 + num_species*3, 0, epsilon0*Ex_nodal[0] - 0.5*dt*app_curr_x_nodal[0]); 
  gkyl_mat_set(&rhs0, 1 + num_species*3, 0, epsilon0*Ey_nodal[0] - 0.5*dt*app_curr_y_nodal[0]); 
  gkyl_mat_set(&rhs0, 2 + num_species*3, 0, epsilon0*Ez_nodal[0] - 0.5*dt*app_curr_z_nodal[0]); 
 
  gkyl_mat_set(&lhs0, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs0, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs0, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs1 = gkyl_nmat_get(A_n, count+1); 
  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs_n, count+1); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs1, 0.0); gkyl_mat_clear(&rhs1, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs1, 0 + s*3, 0, ux_nodal[s][1] + 0.5*dt*(qbym[s]*ext_Ex_nodal[1] + app_accel_x_nodal[s][1])); 
    gkyl_mat_set(&rhs1, 1 + s*3, 0, uy_nodal[s][1] + 0.5*dt*(qbym[s]*ext_Ey_nodal[1] + app_accel_y_nodal[s][1])); 
    gkyl_mat_set(&rhs1, 2 + s*3, 0, uz_nodal[s][1] + 0.5*dt*(qbym[s]*ext_Ez_nodal[1] + app_accel_z_nodal[s][1])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
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
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs1, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][1]); 
      gkyl_mat_set(&lhs1, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][1]); 
      gkyl_mat_set(&lhs1, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][1]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs1, 0 + num_species*3, 0, epsilon0*Ex_nodal[1] - 0.5*dt*app_curr_x_nodal[1]); 
  gkyl_mat_set(&rhs1, 1 + num_species*3, 0, epsilon0*Ey_nodal[1] - 0.5*dt*app_curr_y_nodal[1]); 
  gkyl_mat_set(&rhs1, 2 + num_species*3, 0, epsilon0*Ez_nodal[1] - 0.5*dt*app_curr_z_nodal[1]); 
 
  gkyl_mat_set(&lhs1, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs1, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs1, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs2 = gkyl_nmat_get(A_n, count+2); 
  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs_n, count+2); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs2, 0.0); gkyl_mat_clear(&rhs2, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs2, 0 + s*3, 0, ux_nodal[s][2] + 0.5*dt*(qbym[s]*ext_Ex_nodal[2] + app_accel_x_nodal[s][2])); 
    gkyl_mat_set(&rhs2, 1 + s*3, 0, uy_nodal[s][2] + 0.5*dt*(qbym[s]*ext_Ey_nodal[2] + app_accel_y_nodal[s][2])); 
    gkyl_mat_set(&rhs2, 2 + s*3, 0, uz_nodal[s][2] + 0.5*dt*(qbym[s]*ext_Ez_nodal[2] + app_accel_z_nodal[s][2])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
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
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs2, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][2]); 
      gkyl_mat_set(&lhs2, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][2]); 
      gkyl_mat_set(&lhs2, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][2]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs2, 0 + num_species*3, 0, epsilon0*Ex_nodal[2] - 0.5*dt*app_curr_x_nodal[2]); 
  gkyl_mat_set(&rhs2, 1 + num_species*3, 0, epsilon0*Ey_nodal[2] - 0.5*dt*app_curr_y_nodal[2]); 
  gkyl_mat_set(&rhs2, 2 + num_species*3, 0, epsilon0*Ez_nodal[2] - 0.5*dt*app_curr_z_nodal[2]); 
 
  gkyl_mat_set(&lhs2, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs2, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs2, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs3 = gkyl_nmat_get(A_n, count+3); 
  struct gkyl_mat rhs3 = gkyl_nmat_get(rhs_n, count+3); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs3, 0.0); gkyl_mat_clear(&rhs3, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs3, 0 + s*3, 0, ux_nodal[s][3] + 0.5*dt*(qbym[s]*ext_Ex_nodal[3] + app_accel_x_nodal[s][3])); 
    gkyl_mat_set(&rhs3, 1 + s*3, 0, uy_nodal[s][3] + 0.5*dt*(qbym[s]*ext_Ey_nodal[3] + app_accel_y_nodal[s][3])); 
    gkyl_mat_set(&rhs3, 2 + s*3, 0, uz_nodal[s][3] + 0.5*dt*(qbym[s]*ext_Ez_nodal[3] + app_accel_z_nodal[s][3])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
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
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs3, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][3]); 
      gkyl_mat_set(&lhs3, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][3]); 
      gkyl_mat_set(&lhs3, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][3]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs3, 0 + num_species*3, 0, epsilon0*Ex_nodal[3] - 0.5*dt*app_curr_x_nodal[3]); 
  gkyl_mat_set(&rhs3, 1 + num_species*3, 0, epsilon0*Ey_nodal[3] - 0.5*dt*app_curr_y_nodal[3]); 
  gkyl_mat_set(&rhs3, 2 + num_species*3, 0, epsilon0*Ez_nodal[3] - 0.5*dt*app_curr_z_nodal[3]); 
 
  gkyl_mat_set(&lhs3, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs3, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs3, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs4 = gkyl_nmat_get(A_n, count+4); 
  struct gkyl_mat rhs4 = gkyl_nmat_get(rhs_n, count+4); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs4, 0.0); gkyl_mat_clear(&rhs4, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs4, 0 + s*3, 0, ux_nodal[s][4] + 0.5*dt*(qbym[s]*ext_Ex_nodal[4] + app_accel_x_nodal[s][4])); 
    gkyl_mat_set(&rhs4, 1 + s*3, 0, uy_nodal[s][4] + 0.5*dt*(qbym[s]*ext_Ey_nodal[4] + app_accel_y_nodal[s][4])); 
    gkyl_mat_set(&rhs4, 2 + s*3, 0, uz_nodal[s][4] + 0.5*dt*(qbym[s]*ext_Ez_nodal[4] + app_accel_z_nodal[s][4])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs4, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs4, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs4, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs4, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs4, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs4, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs4, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[4]); 
    gkyl_mat_set(&lhs4, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[4]); 
 
    gkyl_mat_set(&lhs4, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[4]); 
    gkyl_mat_set(&lhs4, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[4]); 
 
    gkyl_mat_set(&lhs4, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[4]); 
    gkyl_mat_set(&lhs4, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[4]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs4, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][4]); 
      gkyl_mat_set(&lhs4, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][4]); 
      gkyl_mat_set(&lhs4, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][4]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs4, 0 + num_species*3, 0, epsilon0*Ex_nodal[4] - 0.5*dt*app_curr_x_nodal[4]); 
  gkyl_mat_set(&rhs4, 1 + num_species*3, 0, epsilon0*Ey_nodal[4] - 0.5*dt*app_curr_y_nodal[4]); 
  gkyl_mat_set(&rhs4, 2 + num_species*3, 0, epsilon0*Ez_nodal[4] - 0.5*dt*app_curr_z_nodal[4]); 
 
  gkyl_mat_set(&lhs4, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs4, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs4, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs5 = gkyl_nmat_get(A_n, count+5); 
  struct gkyl_mat rhs5 = gkyl_nmat_get(rhs_n, count+5); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs5, 0.0); gkyl_mat_clear(&rhs5, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs5, 0 + s*3, 0, ux_nodal[s][5] + 0.5*dt*(qbym[s]*ext_Ex_nodal[5] + app_accel_x_nodal[s][5])); 
    gkyl_mat_set(&rhs5, 1 + s*3, 0, uy_nodal[s][5] + 0.5*dt*(qbym[s]*ext_Ey_nodal[5] + app_accel_y_nodal[s][5])); 
    gkyl_mat_set(&rhs5, 2 + s*3, 0, uz_nodal[s][5] + 0.5*dt*(qbym[s]*ext_Ez_nodal[5] + app_accel_z_nodal[s][5])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs5, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs5, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs5, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs5, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs5, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs5, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs5, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[5]); 
    gkyl_mat_set(&lhs5, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[5]); 
 
    gkyl_mat_set(&lhs5, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[5]); 
    gkyl_mat_set(&lhs5, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[5]); 
 
    gkyl_mat_set(&lhs5, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[5]); 
    gkyl_mat_set(&lhs5, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[5]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs5, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][5]); 
      gkyl_mat_set(&lhs5, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][5]); 
      gkyl_mat_set(&lhs5, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][5]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs5, 0 + num_species*3, 0, epsilon0*Ex_nodal[5] - 0.5*dt*app_curr_x_nodal[5]); 
  gkyl_mat_set(&rhs5, 1 + num_species*3, 0, epsilon0*Ey_nodal[5] - 0.5*dt*app_curr_y_nodal[5]); 
  gkyl_mat_set(&rhs5, 2 + num_species*3, 0, epsilon0*Ez_nodal[5] - 0.5*dt*app_curr_z_nodal[5]); 
 
  gkyl_mat_set(&lhs5, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs5, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs5, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs6 = gkyl_nmat_get(A_n, count+6); 
  struct gkyl_mat rhs6 = gkyl_nmat_get(rhs_n, count+6); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs6, 0.0); gkyl_mat_clear(&rhs6, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs6, 0 + s*3, 0, ux_nodal[s][6] + 0.5*dt*(qbym[s]*ext_Ex_nodal[6] + app_accel_x_nodal[s][6])); 
    gkyl_mat_set(&rhs6, 1 + s*3, 0, uy_nodal[s][6] + 0.5*dt*(qbym[s]*ext_Ey_nodal[6] + app_accel_y_nodal[s][6])); 
    gkyl_mat_set(&rhs6, 2 + s*3, 0, uz_nodal[s][6] + 0.5*dt*(qbym[s]*ext_Ez_nodal[6] + app_accel_z_nodal[s][6])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs6, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs6, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs6, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs6, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs6, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs6, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs6, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[6]); 
    gkyl_mat_set(&lhs6, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[6]); 
 
    gkyl_mat_set(&lhs6, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[6]); 
    gkyl_mat_set(&lhs6, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[6]); 
 
    gkyl_mat_set(&lhs6, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[6]); 
    gkyl_mat_set(&lhs6, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[6]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs6, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][6]); 
      gkyl_mat_set(&lhs6, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][6]); 
      gkyl_mat_set(&lhs6, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][6]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs6, 0 + num_species*3, 0, epsilon0*Ex_nodal[6] - 0.5*dt*app_curr_x_nodal[6]); 
  gkyl_mat_set(&rhs6, 1 + num_species*3, 0, epsilon0*Ey_nodal[6] - 0.5*dt*app_curr_y_nodal[6]); 
  gkyl_mat_set(&rhs6, 2 + num_species*3, 0, epsilon0*Ez_nodal[6] - 0.5*dt*app_curr_z_nodal[6]); 
 
  gkyl_mat_set(&lhs6, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs6, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs6, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
  struct gkyl_mat lhs7 = gkyl_nmat_get(A_n, count+7); 
  struct gkyl_mat rhs7 = gkyl_nmat_get(rhs_n, count+7); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs7, 0.0); gkyl_mat_clear(&rhs7, 0.0); 

  for (int s = 0; s < num_species; ++s) { 
    // Set RHS for flow velocity equations, including solution at known time-step and external forces. 
    gkyl_mat_set(&rhs7, 0 + s*3, 0, ux_nodal[s][7] + 0.5*dt*(qbym[s]*ext_Ex_nodal[7] + app_accel_x_nodal[s][7])); 
    gkyl_mat_set(&rhs7, 1 + s*3, 0, uy_nodal[s][7] + 0.5*dt*(qbym[s]*ext_Ey_nodal[7] + app_accel_y_nodal[s][7])); 
    gkyl_mat_set(&rhs7, 2 + s*3, 0, uz_nodal[s][7] + 0.5*dt*(qbym[s]*ext_Ez_nodal[7] + app_accel_z_nodal[s][7])); 
 
    // Set LHS matrix for flow velocity equation: u_s^{n+1} - 0.5*dt*(q_s/m_s*E^{n+1} + q_s/m_s*u_s^{n+1} x B^n). 
    double E_field_fac = -0.5*dt*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs7, 0 + s*3, 0 + s*3, 1.0); 
    gkyl_mat_set(&lhs7, 1 + s*3, 1 + s*3, 1.0); 
    gkyl_mat_set(&lhs7, 2 + s*3, 2 + s*3, 1.0); 
 
    gkyl_mat_set(&lhs7, 0 + s*3, 0 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs7, 1 + s*3, 1 + num_species*3, E_field_fac); 
    gkyl_mat_set(&lhs7, 2 + s*3, 2 + num_species*3, E_field_fac); 
 
    gkyl_mat_set(&lhs7, 1 + s*3, 2 + s*3, B_field_fac*tot_Bx_nodal[7]); 
    gkyl_mat_set(&lhs7, 2 + s*3, 1 + s*3, -B_field_fac*tot_Bx_nodal[7]); 
 
    gkyl_mat_set(&lhs7, 2 + s*3, 0 + s*3, B_field_fac*tot_By_nodal[7]); 
    gkyl_mat_set(&lhs7, 0 + s*3, 2 + s*3, -B_field_fac*tot_By_nodal[7]); 
 
    gkyl_mat_set(&lhs7, 0 + s*3, 1 + s*3, B_field_fac*tot_Bz_nodal[7]); 
    gkyl_mat_set(&lhs7, 1 + s*3, 0 + s*3, -B_field_fac*tot_Bz_nodal[7]); 
 
    // For Ampere's Law LHS: epsilon0*E^{n+1} + 0.5*dt*sum_s q_s/m_s*rho_s^n u_s^{n+1}. 
    if (!pkpm_field_static) { 
      gkyl_mat_set(&lhs7, 0 + num_species*3, 0 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][7]); 
      gkyl_mat_set(&lhs7, 1 + num_species*3, 1 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][7]); 
      gkyl_mat_set(&lhs7, 2 + num_species*3, 2 + s*3, 0.5*dt*qbym[s]*rho_nodal[s][7]); 
    } 
  } 
  // Set RHS for Ampere's Law, including solution at known time-step and external currents. 
  gkyl_mat_set(&rhs7, 0 + num_species*3, 0, epsilon0*Ex_nodal[7] - 0.5*dt*app_curr_x_nodal[7]); 
  gkyl_mat_set(&rhs7, 1 + num_species*3, 0, epsilon0*Ey_nodal[7] - 0.5*dt*app_curr_y_nodal[7]); 
  gkyl_mat_set(&rhs7, 2 + num_species*3, 0, epsilon0*Ez_nodal[7] - 0.5*dt*app_curr_z_nodal[7]); 
 
  gkyl_mat_set(&lhs7, 0 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs7, 1 + num_species*3, 0 + num_species*3, 1.0); 
  gkyl_mat_set(&lhs7, 2 + num_species*3, 0 + num_species*3, 1.0); 
 
} 
