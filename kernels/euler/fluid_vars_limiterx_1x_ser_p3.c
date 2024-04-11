#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_limiterx_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    double *fluid_l, double *fluid_c, double *fluid_r) 
{ 
  // limiter_fac: Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3))
  // wv_eqn:      Wave equation for computing waves for limiting characteristics
  // fluid_l/c/r: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.

  double *rho_l = &fluid_l[0]; 
  double *rhoux_l = &fluid_l[4]; 
  double *rhouy_l = &fluid_l[8]; 
  double *rhouz_l = &fluid_l[12]; 
  double *energy_l = &fluid_l[16]; 

  double *rho_c = &fluid_c[0]; 
  double *rhoux_c = &fluid_c[4]; 
  double *rhouy_c = &fluid_c[8]; 
  double *rhouz_c = &fluid_c[12]; 
  double *energy_c = &fluid_c[16]; 

  double *rho_r = &fluid_r[0]; 
  double *rhoux_r = &fluid_r[4]; 
  double *rhouy_r = &fluid_r[8]; 
  double *rhouz_r = &fluid_r[12]; 
  double *energy_r = &fluid_r[16]; 

  double fluid_avg_l[5], fluid_avg_c[5], fluid_avg_r[5] = {0.0};
  double delta_l[5], delta_c[5], delta_r[5] = {0.0};
  double waves_slope_l[15], waves_slope_c[15], waves_slope_r[15] = {0.0};
  double speeds[3];

  fluid_avg_l[0] = rho_l[0];
  fluid_avg_l[1] = rhoux_l[0];
  fluid_avg_l[2] = rhouy_l[0];
  fluid_avg_l[3] = rhouz_l[0];
  fluid_avg_l[4] = energy_l[0];

  fluid_avg_c[0] = rho_c[0];
  fluid_avg_c[1] = rhoux_c[0];
  fluid_avg_c[2] = rhouy_c[0];
  fluid_avg_c[3] = rhouz_c[0];
  fluid_avg_c[4] = energy_c[0];

  fluid_avg_r[0] = rho_r[0];
  fluid_avg_r[1] = rhoux_r[0];
  fluid_avg_r[2] = rhouy_r[0];
  fluid_avg_r[3] = rhouz_r[0];
  fluid_avg_r[4] = energy_r[0];

  delta_l[0] = limiter_fac*(fluid_avg_c[0] - fluid_avg_l[0]);
  delta_l[1] = limiter_fac*(fluid_avg_c[1] - fluid_avg_l[1]);
  delta_l[2] = limiter_fac*(fluid_avg_c[2] - fluid_avg_l[2]);
  delta_l[3] = limiter_fac*(fluid_avg_c[3] - fluid_avg_l[3]);
  delta_l[4] = limiter_fac*(fluid_avg_c[4] - fluid_avg_l[4]);

  delta_c[0] = rho_c[1];
  delta_c[1] = rhoux_c[1];
  delta_c[2] = rhouy_c[1];
  delta_c[3] = rhouz_c[1];
  delta_c[4] = energy_c[1];

  delta_r[0] = limiter_fac*(fluid_avg_r[0] - fluid_avg_c[0]);
  delta_r[1] = limiter_fac*(fluid_avg_r[1] - fluid_avg_c[1]);
  delta_r[2] = limiter_fac*(fluid_avg_r[2] - fluid_avg_c[2]);
  delta_r[3] = limiter_fac*(fluid_avg_r[3] - fluid_avg_c[3]);
  delta_r[4] = limiter_fac*(fluid_avg_r[4] - fluid_avg_c[4]);

  double my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l,
    fluid_avg_c, fluid_avg_c, waves_slope_l, speeds);
  double my_max_speed_c = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_c,
    fluid_avg_c, fluid_avg_c, waves_slope_c, speeds);
  double my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r,
    fluid_avg_c, fluid_avg_c, waves_slope_r, speeds);

  // limit wave components, accumulating it to slope
  double mm[15] = {0.0};
  double slope[5] = {0.0};
  mm[0] = minmod(waves_slope_c[0], waves_slope_l[0], waves_slope_r[0]);
  mm[1] = minmod(waves_slope_c[1], waves_slope_l[1], waves_slope_r[1]);
  mm[2] = minmod(waves_slope_c[2], waves_slope_l[2], waves_slope_r[2]);
  mm[3] = minmod(waves_slope_c[3], waves_slope_l[3], waves_slope_r[3]);
  mm[4] = minmod(waves_slope_c[4], waves_slope_l[4], waves_slope_r[4]);

  slope[0] += mm[0];
  slope[1] += mm[1];
  slope[2] += mm[2];
  slope[3] += mm[3];
  slope[4] += mm[4];

  mm[5] = minmod(waves_slope_c[5], waves_slope_l[5], waves_slope_r[5]);
  mm[6] = minmod(waves_slope_c[6], waves_slope_l[6], waves_slope_r[6]);
  mm[7] = minmod(waves_slope_c[7], waves_slope_l[7], waves_slope_r[7]);
  mm[8] = minmod(waves_slope_c[8], waves_slope_l[8], waves_slope_r[8]);
  mm[9] = minmod(waves_slope_c[9], waves_slope_l[9], waves_slope_r[9]);

  slope[0] += mm[5];
  slope[1] += mm[6];
  slope[2] += mm[7];
  slope[3] += mm[8];
  slope[4] += mm[9];

  mm[10] = minmod(waves_slope_c[10], waves_slope_l[10], waves_slope_r[10]);
  mm[11] = minmod(waves_slope_c[11], waves_slope_l[11], waves_slope_r[11]);
  mm[12] = minmod(waves_slope_c[12], waves_slope_l[12], waves_slope_r[12]);
  mm[13] = minmod(waves_slope_c[13], waves_slope_l[13], waves_slope_r[13]);
  mm[14] = minmod(waves_slope_c[14], waves_slope_l[14], waves_slope_r[14]);

  slope[0] += mm[10];
  slope[1] += mm[11];
  slope[2] += mm[12];
  slope[3] += mm[13];
  slope[4] += mm[14];

  rho_c[1] = slope[0];
  rhoux_c[1] = slope[1];
  rhouy_c[1] = slope[2];
  rhouz_c[1] = slope[3];
  energy_c[1] = slope[4];

  for (int i=0; i<3; ++i) {
    if (mm[i*5] != waves_slope_c[i*5]) {
      rho_c[2] = 0.0;
      rho_c[3] = 0.0;
    }
    if (mm[i*5+1] != waves_slope_c[i*5+1]) {
      rhoux_c[2] = 0.0;
      rhoux_c[3] = 0.0;
    }
    if (mm[i*5+2] != waves_slope_c[i*5+2]) {
      rhouy_c[2] = 0.0;
      rhouy_c[3] = 0.0;
    }
    if (mm[i*5+3] != waves_slope_c[i*5+3]) {
      rhouz_c[2] = 0.0;
      rhouz_c[3] = 0.0;
    }
    if (mm[i*5+4] != waves_slope_c[i*5+4]) {
      energy_c[2] = 0.0;
      energy_c[3] = 0.0;
    }
  }
}