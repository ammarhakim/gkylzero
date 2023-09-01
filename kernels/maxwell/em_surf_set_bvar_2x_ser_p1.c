#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_surf_set_bvar_2x_ser_p1(const double* bvar, double* GKYL_RESTRICT bvar_surf) 
{ 
  // bvar:      Input volume expansion of b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
  // bvar_surf: Output surface expansion of magnetic field unit vector and tensor in each direction. 
 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[2]; 
  double *bxbx_xl = &bvar_surf[4]; 
  double *bxbx_xr = &bvar_surf[6]; 
  double *bxby_xl = &bvar_surf[8]; 
  double *bxby_xr = &bvar_surf[10]; 
  double *bxbz_xl = &bvar_surf[12]; 
  double *bxbz_xr = &bvar_surf[14]; 
 
  bx_xl[0] = 0.7071067811865475*bx[0]-1.224744871391589*bx[1]; 
  bx_xl[1] = 0.7071067811865475*bx[2]-1.224744871391589*bx[3]; 
  bxbx_xl[0] = 0.7071067811865475*bxbx[0]-1.224744871391589*bxbx[1]; 
  bxbx_xl[1] = 0.7071067811865475*bxbx[2]-1.224744871391589*bxbx[3]; 
  bxby_xl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[1]; 
  bxby_xl[1] = 0.7071067811865475*bxby[2]-1.224744871391589*bxby[3]; 
  bxbz_xl[0] = 0.7071067811865475*bxbz[0]-1.224744871391589*bxbz[1]; 
  bxbz_xl[1] = 0.7071067811865475*bxbz[2]-1.224744871391589*bxbz[3]; 
 
  bx_xr[0] = 1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[1] = 1.224744871391589*bx[3]+0.7071067811865475*bx[2]; 
  bxbx_xr[0] = 1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xr[1] = 1.224744871391589*bxbx[3]+0.7071067811865475*bxbx[2]; 
  bxby_xr[0] = 1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xr[1] = 1.224744871391589*bxby[3]+0.7071067811865475*bxby[2]; 
  bxbz_xr[0] = 1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xr[1] = 1.224744871391589*bxbz[3]+0.7071067811865475*bxbz[2]; 
 
  double *by_yl = &bvar_surf[16]; 
  double *by_yr = &bvar_surf[18]; 
  double *bxby_yl = &bvar_surf[20]; 
  double *bxby_yr = &bvar_surf[22]; 
  double *byby_yl = &bvar_surf[24]; 
  double *byby_yr = &bvar_surf[26]; 
  double *bybz_yl = &bvar_surf[28]; 
  double *bybz_yr = &bvar_surf[30]; 
 
  by_yl[0] = 0.7071067811865475*by[0]-1.224744871391589*by[2]; 
  by_yl[1] = 0.7071067811865475*by[1]-1.224744871391589*by[3]; 
  bxby_yl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[2]; 
  bxby_yl[1] = 0.7071067811865475*bxby[1]-1.224744871391589*bxby[3]; 
  byby_yl[0] = 0.7071067811865475*byby[0]-1.224744871391589*byby[2]; 
  byby_yl[1] = 0.7071067811865475*byby[1]-1.224744871391589*byby[3]; 
  bybz_yl[0] = 0.7071067811865475*bybz[0]-1.224744871391589*bybz[2]; 
  bybz_yl[1] = 0.7071067811865475*bybz[1]-1.224744871391589*bybz[3]; 
 
  by_yr[0] = 1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yr[1] = 1.224744871391589*by[3]+0.7071067811865475*by[1]; 
  bxby_yr[0] = 1.224744871391589*bxby[2]+0.7071067811865475*bxby[0]; 
  bxby_yr[1] = 1.224744871391589*bxby[3]+0.7071067811865475*bxby[1]; 
  byby_yr[0] = 1.224744871391589*byby[2]+0.7071067811865475*byby[0]; 
  byby_yr[1] = 1.224744871391589*byby[3]+0.7071067811865475*byby[1]; 
  bybz_yr[0] = 1.224744871391589*bybz[2]+0.7071067811865475*bybz[0]; 
  bybz_yr[1] = 1.224744871391589*bybz[3]+0.7071067811865475*bybz[1]; 
 
} 
 
