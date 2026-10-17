#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_wv_ten_moment_priv.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv,
  const double *moms, const double *bvar, const double *u, const double *euler,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_x_1x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfx_1x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_2x_ser_p1(const double *w, const double *dxv,
  const double *moms, const double *bvar, const double *u, const double *euler,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_x_2x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_y_2x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_3x_ser_p1(const double *w, const double *dxv,
  const double *moms, const double *bvar, const double *u, const double *euler,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_x_3x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfx_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_y_3x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfy_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_z_3x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH double pkpm_fluid_surfz_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

EXTERN_C_END 
