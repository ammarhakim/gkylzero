#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_wv_ten_moment_priv.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_solve_1x_ser_p1(int num_species, 
  double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_1x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_1x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_1x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_1x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_solve_2x_ser_p1(int num_species, 
  double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_2x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_2x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_2x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_2x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_2x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_solve_3x_ser_p1(int num_species, 
  double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_3x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_3x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_3x_ser_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_3x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_3x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_tensor_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_tensor_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_tensor_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_set_2x_tensor_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_nodal_copy_2x_tensor_p2(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_source_2x_tensor_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_tensor_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_tensor_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_tensor_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_tensor_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_tensor_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH 
static inline void 
modal_to_quad_nodal_1d_ser_p1(const double* fmodal, double* fquad) { 
  fquad[0] = 0.7071067811865475*fmodal[0]-0.7071067811865475*fmodal[1]; 
  fquad[1] = 0.7071067811865475*fmodal[1]+0.7071067811865475*fmodal[0]; 
} 
GKYL_CU_DH 
static inline void 
modal_to_quad_nodal_2d_ser_p1(const double* fmodal, double* fquad) { 
  fquad[0] = 0.5*fmodal[3]-0.5*fmodal[2]-0.5*fmodal[1]+0.5*fmodal[0]; 
  fquad[1] = (-0.5*fmodal[3])+0.5*fmodal[2]-0.5*fmodal[1]+0.5*fmodal[0]; 
  fquad[2] = (-0.5*fmodal[3])-0.5*fmodal[2]+0.5*fmodal[1]+0.5*fmodal[0]; 
  fquad[3] = 0.5*fmodal[3]+0.5*fmodal[2]+0.5*fmodal[1]+0.5*fmodal[0]; 
} 
GKYL_CU_DH 
static inline void 
modal_to_quad_nodal_3d_ser_p1(const double* fmodal, double* fquad) { 
  fquad[0] = (-0.3535533905932737*fmodal[7])+0.3535533905932737*fmodal[6]+0.3535533905932737*fmodal[5]+0.3535533905932737*fmodal[4]-0.3535533905932737*fmodal[3]-0.3535533905932737*fmodal[2]-0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[1] = 0.3535533905932737*fmodal[7]-0.3535533905932737*fmodal[6]-0.3535533905932737*fmodal[5]+0.3535533905932737*fmodal[4]+0.3535533905932737*fmodal[3]-0.3535533905932737*fmodal[2]-0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[2] = 0.3535533905932737*fmodal[7]-0.3535533905932737*fmodal[6]+0.3535533905932737*fmodal[5]-0.3535533905932737*fmodal[4]-0.3535533905932737*fmodal[3]+0.3535533905932737*fmodal[2]-0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[3] = (-0.3535533905932737*fmodal[7])+0.3535533905932737*fmodal[6]-0.3535533905932737*fmodal[5]-0.3535533905932737*fmodal[4]+0.3535533905932737*fmodal[3]+0.3535533905932737*fmodal[2]-0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[4] = 0.3535533905932737*fmodal[7]+0.3535533905932737*fmodal[6]-0.3535533905932737*fmodal[5]-0.3535533905932737*fmodal[4]-0.3535533905932737*fmodal[3]-0.3535533905932737*fmodal[2]+0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[5] = (-0.3535533905932737*fmodal[7])-0.3535533905932737*fmodal[6]+0.3535533905932737*fmodal[5]-0.3535533905932737*fmodal[4]+0.3535533905932737*fmodal[3]-0.3535533905932737*fmodal[2]+0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[6] = (-0.3535533905932737*fmodal[7])-0.3535533905932737*fmodal[6]-0.3535533905932737*fmodal[5]+0.3535533905932737*fmodal[4]-0.3535533905932737*fmodal[3]+0.3535533905932737*fmodal[2]+0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
  fquad[7] = 0.3535533905932737*fmodal[7]+0.3535533905932737*fmodal[6]+0.3535533905932737*fmodal[5]+0.3535533905932737*fmodal[4]+0.3535533905932737*fmodal[3]+0.3535533905932737*fmodal[2]+0.3535533905932737*fmodal[1]+0.3535533905932737*fmodal[0]; 
} 
GKYL_CU_DH 
static inline void 
quad_nodal_to_modal_1d_ser_p1(const double* fquad, double* fmodal) { 
  fmodal[0] = 0.7071067811865475*fquad[1]+0.7071067811865475*fquad[0]; 
  fmodal[1] = 0.7071067811865475*fquad[1]-0.7071067811865475*fquad[0]; 
} 
GKYL_CU_DH 
static inline void 
quad_nodal_to_modal_2d_ser_p1(const double* fquad, double* fmodal) { 
  fmodal[0] = 0.5*fquad[3]+0.5*fquad[2]+0.5*fquad[1]+0.5*fquad[0]; 
  fmodal[1] = 0.5*fquad[3]+0.5*fquad[2]-0.5*fquad[1]-0.5*fquad[0]; 
  fmodal[2] = 0.5*fquad[3]-0.5*fquad[2]+0.5*fquad[1]-0.5*fquad[0]; 
  fmodal[3] = 0.5*fquad[3]-0.5*fquad[2]-0.5*fquad[1]+0.5*fquad[0]; 
} 
GKYL_CU_DH 
static inline void 
quad_nodal_to_modal_3d_ser_p1(const double* fquad, double* fmodal) { 
  fmodal[0] = 0.3535533905932737*fquad[7]+0.3535533905932737*fquad[6]+0.3535533905932737*fquad[5]+0.3535533905932737*fquad[4]+0.3535533905932737*fquad[3]+0.3535533905932737*fquad[2]+0.3535533905932737*fquad[1]+0.3535533905932737*fquad[0]; 
  fmodal[1] = 0.3535533905932737*fquad[7]+0.3535533905932737*fquad[6]+0.3535533905932737*fquad[5]+0.3535533905932737*fquad[4]-0.3535533905932737*fquad[3]-0.3535533905932737*fquad[2]-0.3535533905932737*fquad[1]-0.3535533905932737*fquad[0]; 
  fmodal[2] = 0.3535533905932737*fquad[7]+0.3535533905932737*fquad[6]-0.3535533905932737*fquad[5]-0.3535533905932737*fquad[4]+0.3535533905932737*fquad[3]+0.3535533905932737*fquad[2]-0.3535533905932737*fquad[1]-0.3535533905932737*fquad[0]; 
  fmodal[3] = 0.3535533905932737*fquad[7]-0.3535533905932737*fquad[6]+0.3535533905932737*fquad[5]-0.3535533905932737*fquad[4]+0.3535533905932737*fquad[3]-0.3535533905932737*fquad[2]+0.3535533905932737*fquad[1]-0.3535533905932737*fquad[0]; 
  fmodal[4] = 0.3535533905932737*fquad[7]+0.3535533905932737*fquad[6]-0.3535533905932737*fquad[5]-0.3535533905932737*fquad[4]-0.3535533905932737*fquad[3]-0.3535533905932737*fquad[2]+0.3535533905932737*fquad[1]+0.3535533905932737*fquad[0]; 
  fmodal[5] = 0.3535533905932737*fquad[7]-0.3535533905932737*fquad[6]+0.3535533905932737*fquad[5]-0.3535533905932737*fquad[4]-0.3535533905932737*fquad[3]+0.3535533905932737*fquad[2]-0.3535533905932737*fquad[1]+0.3535533905932737*fquad[0]; 
  fmodal[6] = 0.3535533905932737*fquad[7]-0.3535533905932737*fquad[6]-0.3535533905932737*fquad[5]+0.3535533905932737*fquad[4]+0.3535533905932737*fquad[3]-0.3535533905932737*fquad[2]-0.3535533905932737*fquad[1]+0.3535533905932737*fquad[0]; 
  fmodal[7] = 0.3535533905932737*fquad[7]-0.3535533905932737*fquad[6]-0.3535533905932737*fquad[5]+0.3535533905932737*fquad[4]-0.3535533905932737*fquad[3]+0.3535533905932737*fquad[2]+0.3535533905932737*fquad[1]-0.3535533905932737*fquad[0]; 
} 

GKYL_CU_DH 
static void
implicit_nodal_pkpm_em_source_update(int nfluids, double dt, double q_over_m[GKYL_MAX_SPECIES], double epsilon0, 
  double rho_s[GKYL_MAX_SPECIES], double fluid_s[GKYL_MAX_SPECIES][3],
  double app_accel_s[GKYL_MAX_SPECIES][3], 
  double E[3], double ext_E[3], double tot_B[3], double app_current[3])
{
  double Bx = tot_B[0];
  double By = tot_B[1];
  double Bz = tot_B[2];
  double B_mag = sqrt((Bx * Bx) + (By * By) + (Bz * Bz));
  
  double bx = 0.0, by = 0.0, bz = 0.0;
  if (B_mag > 0.0) {
    bx = Bx / B_mag;
    by = By / B_mag;
    bz = Bz / B_mag;
  }

  double wc_dt[GKYL_MAX_SPECIES];
  double wp_dt_sq[GKYL_MAX_SPECIES];
  double J_old[GKYL_MAX_SPECIES][3];
  double J[GKYL_MAX_SPECIES][3];

  double w0_sq = 0.0;
  double gam_sq = 0.0;
  double delta = 0.0;
  double Kx = 0.0, Ky = 0.0, Kz = 0.0;

  for (int i = 0; i < nfluids; i++) {
    const double *app_accel = app_accel_s[i];

    double rho = rho_s[i];
    double mom_x = fluid_s[i][0], mom_y = fluid_s[i][1], mom_z = fluid_s[i][2];

    J_old[i][0] = mom_x * q_over_m[i];
    J_old[i][1] = mom_y * q_over_m[i];
    J_old[i][2] = mom_z * q_over_m[i];

    J[i][0] = J_old[i][0] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_E[0]) + app_accel[0]));
    J[i][1] = J_old[i][1] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_E[1]) + app_accel[1]));
    J[i][2] = J_old[i][2] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_E[2]) + app_accel[2]));

    wc_dt[i] = q_over_m[i] * B_mag * dt;
    wp_dt_sq[i] = (rho * (q_over_m[i] * q_over_m[i]) * (dt * dt)) / epsilon0;

    double denom = 1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0);
    w0_sq += wp_dt_sq[i] / denom;
    gam_sq += (wp_dt_sq[i] * (wc_dt[i] * wc_dt[i])) / denom;
    delta += (wp_dt_sq[i] * wc_dt[i]) / denom;

    Kx -= (dt / denom) * (J[i][0] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((by * J[i][2]) - (bz * J[i][1]))));
    Ky -= (dt / denom) * (J[i][1] + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bz * J[i][0]) - (bx * J[i][2]))));
    Kz -= (dt / denom) * (J[i][2] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bx * J[i][1]) - (by * J[i][0]))));
  }

  double Delta_sq = (delta * delta) / (1.0 + (w0_sq / 4.0));

  double Fx_old = E[0] * epsilon0;
  double Fy_old = E[1] * epsilon0;
  double Fz_old = E[2] * epsilon0;

  double Fx = Fx_old - (0.5 * dt * app_current[0]);
  double Fy = Fy_old - (0.5 * dt * app_current[1]);
  double Fz = Fz_old - (0.5 * dt * app_current[2]);

  double Fx_K = Fx + (0.5 * Kx);
  double Fy_K = Fy + (0.5 * Ky);
  double Fz_K = Fz + (0.5 * Kz);

  double Fx_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fx_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bx * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((by * Fz_K) - (bz * Fy_K))));
  double Fy_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fy_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * by * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bz * Fx_K) - (bx * Fz_K))));
  double Fz_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fz_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bz * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bx * Fy_K) - (by * Fx_K))));
  
  E[0] = ((2.0 * Fx_bar) - Fx_old) / epsilon0;
  E[1] = ((2.0 * Fy_bar) - Fy_old) / epsilon0;
  E[2] = ((2.0 * Fz_bar) - Fz_old) / epsilon0;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double Jx_star = J[i][0] + (Fx_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jy_star = J[i][1] + (Fy_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jz_star = J[i][2] + (Fz_bar * ((wp_dt_sq[i] / dt) / 2.0));

    double Jx_new = ((2.0 * (Jx_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((by * Jz_star) - (bz * Jy_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][0];
    double Jy_new = ((2.0 * (Jy_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bz * Jx_star) - (bx * Jz_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][1];
    double Jz_new = ((2.0 * (Jz_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bx * Jy_star) - (by * Jx_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][2];
    
    f[0] = Jx_new / q_over_m[i];
    f[1] = Jy_new / q_over_m[i];
    f[2] = Jz_new / q_over_m[i];
  }
}

EXTERN_C_END 
