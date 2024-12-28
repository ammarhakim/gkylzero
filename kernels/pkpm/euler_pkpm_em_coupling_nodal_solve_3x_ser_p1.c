#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_solve_3x_ser_p1(int num_species, 
  double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // num_species:       number of species being evolved (number of momentum equations). 
  // qbym:              charge/mass ratio for each species. 
  // epsilon0:          permittivity of free space. 
  // pkpm_field_static: boolean for whether or not the self-consistent field is static. 
  // dt:                size of the time step. 
  // app_accel:         Applied accelerations (external forces).
  // ext_em:            Externally applied EM fields.
  // app_current:       Applied external currents.
  // vlasov_pkpm_moms:  [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:        [rho ux, rho uy, rho uz], PKPM Fluid momentum input state vector.
  // em:                [Ex, Ey, Ez, Bx, By, Bz], EM input state vector.

  double rho[GKYL_MAX_SPECIES][8]; 
  double rhoux[GKYL_MAX_SPECIES][8]; 
  double rhouy[GKYL_MAX_SPECIES][8]; 
  double rhouz[GKYL_MAX_SPECIES][8]; 

  double app_accel_x[GKYL_MAX_SPECIES][8]; 
  double app_accel_y[GKYL_MAX_SPECIES][8]; 
  double app_accel_z[GKYL_MAX_SPECIES][8]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_euler_pkpm = euler_pkpm[i]; 
    const double *inp_vlasov_pkpm_moms = vlasov_pkpm_moms[i]; 
    const double *inp_app_accel = app_accel[i]; 

    rho[i][0] = inp_vlasov_pkpm_moms[0]; 
    rhoux[i][0] = inp_euler_pkpm[0]; 
    rhouy[i][0] = inp_euler_pkpm[8]; 
    rhouz[i][0] = inp_euler_pkpm[16]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[8]; 
    app_accel_z[i][0] = inp_app_accel[16]; 

    rho[i][1] = inp_vlasov_pkpm_moms[1]; 
    rhoux[i][1] = inp_euler_pkpm[1]; 
    rhouy[i][1] = inp_euler_pkpm[9]; 
    rhouz[i][1] = inp_euler_pkpm[17]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[9]; 
    app_accel_z[i][1] = inp_app_accel[17]; 

    rho[i][2] = inp_vlasov_pkpm_moms[2]; 
    rhoux[i][2] = inp_euler_pkpm[2]; 
    rhouy[i][2] = inp_euler_pkpm[10]; 
    rhouz[i][2] = inp_euler_pkpm[18]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[10]; 
    app_accel_z[i][2] = inp_app_accel[18]; 

    rho[i][3] = inp_vlasov_pkpm_moms[3]; 
    rhoux[i][3] = inp_euler_pkpm[3]; 
    rhouy[i][3] = inp_euler_pkpm[11]; 
    rhouz[i][3] = inp_euler_pkpm[19]; 

    app_accel_x[i][3] = inp_app_accel[3]; 
    app_accel_y[i][3] = inp_app_accel[11]; 
    app_accel_z[i][3] = inp_app_accel[19]; 

    rho[i][4] = inp_vlasov_pkpm_moms[4]; 
    rhoux[i][4] = inp_euler_pkpm[4]; 
    rhouy[i][4] = inp_euler_pkpm[12]; 
    rhouz[i][4] = inp_euler_pkpm[20]; 

    app_accel_x[i][4] = inp_app_accel[4]; 
    app_accel_y[i][4] = inp_app_accel[12]; 
    app_accel_z[i][4] = inp_app_accel[20]; 

    rho[i][5] = inp_vlasov_pkpm_moms[5]; 
    rhoux[i][5] = inp_euler_pkpm[5]; 
    rhouy[i][5] = inp_euler_pkpm[13]; 
    rhouz[i][5] = inp_euler_pkpm[21]; 

    app_accel_x[i][5] = inp_app_accel[5]; 
    app_accel_y[i][5] = inp_app_accel[13]; 
    app_accel_z[i][5] = inp_app_accel[21]; 

    rho[i][6] = inp_vlasov_pkpm_moms[6]; 
    rhoux[i][6] = inp_euler_pkpm[6]; 
    rhouy[i][6] = inp_euler_pkpm[14]; 
    rhouz[i][6] = inp_euler_pkpm[22]; 

    app_accel_x[i][6] = inp_app_accel[6]; 
    app_accel_y[i][6] = inp_app_accel[14]; 
    app_accel_z[i][6] = inp_app_accel[22]; 

    rho[i][7] = inp_vlasov_pkpm_moms[7]; 
    rhoux[i][7] = inp_euler_pkpm[7]; 
    rhouy[i][7] = inp_euler_pkpm[15]; 
    rhouz[i][7] = inp_euler_pkpm[23]; 

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
  double rhoux_nodal[GKYL_MAX_SPECIES][8]; 
  double rhouy_nodal[GKYL_MAX_SPECIES][8]; 
  double rhouz_nodal[GKYL_MAX_SPECIES][8]; 

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
    modal_to_quad_nodal_3d_ser_p1(rhoux[s], rhoux_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(rhouy[s], rhouy_nodal[s]);
    modal_to_quad_nodal_3d_ser_p1(rhouz[s], rhouz_nodal[s]);
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
 
  double rho_s[GKYL_MAX_SPECIES]; 
  double rhou_s[GKYL_MAX_SPECIES][3]; 
  double app_accel_s[GKYL_MAX_SPECIES][3]; 
  double E[3]; 
  double ext_E[3]; 
  double tot_B[3]; 
  double app_curr[3]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][0]; 
    rhou_s[s][0] = rhoux_nodal[s][0]; 
    rhou_s[s][1] = rhouy_nodal[s][0]; 
    rhou_s[s][2] = rhouz_nodal[s][0]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][0]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][0]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][0]; 
  } 
  E[0] = Ex_nodal[0]; 
  E[1] = Ey_nodal[0]; 
  E[2] = Ez_nodal[0]; 
  ext_E[0] = ext_Ex_nodal[0]; 
  ext_E[1] = ext_Ey_nodal[0]; 
  ext_E[2] = ext_Ez_nodal[0]; 
  tot_B[0] = tot_Bx_nodal[0]; 
  tot_B[1] = tot_By_nodal[0]; 
  tot_B[2] = tot_Bz_nodal[0]; 
  app_curr[0] = app_curr_x_nodal[0]; 
  app_curr[1] = app_curr_y_nodal[0]; 
  app_curr[2] = app_curr_z_nodal[0]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][0] = rhou_s[s][0]; 
    rhouy_nodal[s][0] = rhou_s[s][1]; 
    rhouz_nodal[s][0] = rhou_s[s][2]; 
  } 
    Ex_nodal[0] = E[0]; 
    Ey_nodal[0] = E[1]; 
    Ez_nodal[0] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][1]; 
    rhou_s[s][0] = rhoux_nodal[s][1]; 
    rhou_s[s][1] = rhouy_nodal[s][1]; 
    rhou_s[s][2] = rhouz_nodal[s][1]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][1]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][1]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][1]; 
  } 
  E[0] = Ex_nodal[1]; 
  E[1] = Ey_nodal[1]; 
  E[2] = Ez_nodal[1]; 
  ext_E[0] = ext_Ex_nodal[1]; 
  ext_E[1] = ext_Ey_nodal[1]; 
  ext_E[2] = ext_Ez_nodal[1]; 
  tot_B[0] = tot_Bx_nodal[1]; 
  tot_B[1] = tot_By_nodal[1]; 
  tot_B[2] = tot_Bz_nodal[1]; 
  app_curr[0] = app_curr_x_nodal[1]; 
  app_curr[1] = app_curr_y_nodal[1]; 
  app_curr[2] = app_curr_z_nodal[1]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][1] = rhou_s[s][0]; 
    rhouy_nodal[s][1] = rhou_s[s][1]; 
    rhouz_nodal[s][1] = rhou_s[s][2]; 
  } 
    Ex_nodal[1] = E[0]; 
    Ey_nodal[1] = E[1]; 
    Ez_nodal[1] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][2]; 
    rhou_s[s][0] = rhoux_nodal[s][2]; 
    rhou_s[s][1] = rhouy_nodal[s][2]; 
    rhou_s[s][2] = rhouz_nodal[s][2]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][2]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][2]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][2]; 
  } 
  E[0] = Ex_nodal[2]; 
  E[1] = Ey_nodal[2]; 
  E[2] = Ez_nodal[2]; 
  ext_E[0] = ext_Ex_nodal[2]; 
  ext_E[1] = ext_Ey_nodal[2]; 
  ext_E[2] = ext_Ez_nodal[2]; 
  tot_B[0] = tot_Bx_nodal[2]; 
  tot_B[1] = tot_By_nodal[2]; 
  tot_B[2] = tot_Bz_nodal[2]; 
  app_curr[0] = app_curr_x_nodal[2]; 
  app_curr[1] = app_curr_y_nodal[2]; 
  app_curr[2] = app_curr_z_nodal[2]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][2] = rhou_s[s][0]; 
    rhouy_nodal[s][2] = rhou_s[s][1]; 
    rhouz_nodal[s][2] = rhou_s[s][2]; 
  } 
    Ex_nodal[2] = E[0]; 
    Ey_nodal[2] = E[1]; 
    Ez_nodal[2] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][3]; 
    rhou_s[s][0] = rhoux_nodal[s][3]; 
    rhou_s[s][1] = rhouy_nodal[s][3]; 
    rhou_s[s][2] = rhouz_nodal[s][3]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][3]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][3]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][3]; 
  } 
  E[0] = Ex_nodal[3]; 
  E[1] = Ey_nodal[3]; 
  E[2] = Ez_nodal[3]; 
  ext_E[0] = ext_Ex_nodal[3]; 
  ext_E[1] = ext_Ey_nodal[3]; 
  ext_E[2] = ext_Ez_nodal[3]; 
  tot_B[0] = tot_Bx_nodal[3]; 
  tot_B[1] = tot_By_nodal[3]; 
  tot_B[2] = tot_Bz_nodal[3]; 
  app_curr[0] = app_curr_x_nodal[3]; 
  app_curr[1] = app_curr_y_nodal[3]; 
  app_curr[2] = app_curr_z_nodal[3]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][3] = rhou_s[s][0]; 
    rhouy_nodal[s][3] = rhou_s[s][1]; 
    rhouz_nodal[s][3] = rhou_s[s][2]; 
  } 
    Ex_nodal[3] = E[0]; 
    Ey_nodal[3] = E[1]; 
    Ez_nodal[3] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][4]; 
    rhou_s[s][0] = rhoux_nodal[s][4]; 
    rhou_s[s][1] = rhouy_nodal[s][4]; 
    rhou_s[s][2] = rhouz_nodal[s][4]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][4]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][4]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][4]; 
  } 
  E[0] = Ex_nodal[4]; 
  E[1] = Ey_nodal[4]; 
  E[2] = Ez_nodal[4]; 
  ext_E[0] = ext_Ex_nodal[4]; 
  ext_E[1] = ext_Ey_nodal[4]; 
  ext_E[2] = ext_Ez_nodal[4]; 
  tot_B[0] = tot_Bx_nodal[4]; 
  tot_B[1] = tot_By_nodal[4]; 
  tot_B[2] = tot_Bz_nodal[4]; 
  app_curr[0] = app_curr_x_nodal[4]; 
  app_curr[1] = app_curr_y_nodal[4]; 
  app_curr[2] = app_curr_z_nodal[4]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][4] = rhou_s[s][0]; 
    rhouy_nodal[s][4] = rhou_s[s][1]; 
    rhouz_nodal[s][4] = rhou_s[s][2]; 
  } 
    Ex_nodal[4] = E[0]; 
    Ey_nodal[4] = E[1]; 
    Ez_nodal[4] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][5]; 
    rhou_s[s][0] = rhoux_nodal[s][5]; 
    rhou_s[s][1] = rhouy_nodal[s][5]; 
    rhou_s[s][2] = rhouz_nodal[s][5]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][5]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][5]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][5]; 
  } 
  E[0] = Ex_nodal[5]; 
  E[1] = Ey_nodal[5]; 
  E[2] = Ez_nodal[5]; 
  ext_E[0] = ext_Ex_nodal[5]; 
  ext_E[1] = ext_Ey_nodal[5]; 
  ext_E[2] = ext_Ez_nodal[5]; 
  tot_B[0] = tot_Bx_nodal[5]; 
  tot_B[1] = tot_By_nodal[5]; 
  tot_B[2] = tot_Bz_nodal[5]; 
  app_curr[0] = app_curr_x_nodal[5]; 
  app_curr[1] = app_curr_y_nodal[5]; 
  app_curr[2] = app_curr_z_nodal[5]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][5] = rhou_s[s][0]; 
    rhouy_nodal[s][5] = rhou_s[s][1]; 
    rhouz_nodal[s][5] = rhou_s[s][2]; 
  } 
    Ex_nodal[5] = E[0]; 
    Ey_nodal[5] = E[1]; 
    Ez_nodal[5] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][6]; 
    rhou_s[s][0] = rhoux_nodal[s][6]; 
    rhou_s[s][1] = rhouy_nodal[s][6]; 
    rhou_s[s][2] = rhouz_nodal[s][6]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][6]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][6]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][6]; 
  } 
  E[0] = Ex_nodal[6]; 
  E[1] = Ey_nodal[6]; 
  E[2] = Ez_nodal[6]; 
  ext_E[0] = ext_Ex_nodal[6]; 
  ext_E[1] = ext_Ey_nodal[6]; 
  ext_E[2] = ext_Ez_nodal[6]; 
  tot_B[0] = tot_Bx_nodal[6]; 
  tot_B[1] = tot_By_nodal[6]; 
  tot_B[2] = tot_Bz_nodal[6]; 
  app_curr[0] = app_curr_x_nodal[6]; 
  app_curr[1] = app_curr_y_nodal[6]; 
  app_curr[2] = app_curr_z_nodal[6]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][6] = rhou_s[s][0]; 
    rhouy_nodal[s][6] = rhou_s[s][1]; 
    rhouz_nodal[s][6] = rhou_s[s][2]; 
  } 
    Ex_nodal[6] = E[0]; 
    Ey_nodal[6] = E[1]; 
    Ez_nodal[6] = E[2]; 

  for (int s = 0; s < num_species; ++s) { 
    rho_s[s] = rho_nodal[s][7]; 
    rhou_s[s][0] = rhoux_nodal[s][7]; 
    rhou_s[s][1] = rhouy_nodal[s][7]; 
    rhou_s[s][2] = rhouz_nodal[s][7]; 
    app_accel_s[s][0] = app_accel_x_nodal[s][7]; 
    app_accel_s[s][1] = app_accel_y_nodal[s][7]; 
    app_accel_s[s][2] = app_accel_z_nodal[s][7]; 
  } 
  E[0] = Ex_nodal[7]; 
  E[1] = Ey_nodal[7]; 
  E[2] = Ez_nodal[7]; 
  ext_E[0] = ext_Ex_nodal[7]; 
  ext_E[1] = ext_Ey_nodal[7]; 
  ext_E[2] = ext_Ez_nodal[7]; 
  tot_B[0] = tot_Bx_nodal[7]; 
  tot_B[1] = tot_By_nodal[7]; 
  tot_B[2] = tot_Bz_nodal[7]; 
  app_curr[0] = app_curr_x_nodal[7]; 
  app_curr[1] = app_curr_y_nodal[7]; 
  app_curr[2] = app_curr_z_nodal[7]; 

  implicit_nodal_pkpm_em_source_update(num_species, dt, qbym, epsilon0, 
    rho_s, rhou_s, app_accel_s, E, ext_E, tot_B, app_curr); 

  for (int s = 0; s < num_species; ++s) { 
    rhoux_nodal[s][7] = rhou_s[s][0]; 
    rhouy_nodal[s][7] = rhou_s[s][1]; 
    rhouz_nodal[s][7] = rhou_s[s][2]; 
  } 
    Ex_nodal[7] = E[0]; 
    Ey_nodal[7] = E[1]; 
    Ez_nodal[7] = E[2]; 

  // Project nodal solution onto modal bases. 
  for (int s = 0; s < num_species; ++s) { 
    double *out_rhoux = &euler_pkpm[s][0]; 
    double *out_rhouy = &euler_pkpm[s][8]; 
    double *out_rhouz = &euler_pkpm[s][16]; 
    quad_nodal_to_modal_3d_ser_p1(rhoux_nodal[s], out_rhoux);
    quad_nodal_to_modal_3d_ser_p1(rhouy_nodal[s], out_rhouy);
    quad_nodal_to_modal_3d_ser_p1(rhouz_nodal[s], out_rhouz);
  } 

  // Update electric field if PKPM field is self-consistent. 
  if (!pkpm_field_static) { 
    double *out_Ex = &em[0]; 
    double *out_Ey = &em[8]; 
    double *out_Ez = &em[16]; 
    quad_nodal_to_modal_3d_ser_p1(Ex_nodal, out_Ex);
    quad_nodal_to_modal_3d_ser_p1(Ey_nodal, out_Ey);
    quad_nodal_to_modal_3d_ser_p1(Ez_nodal, out_Ez);
  } 
} 
