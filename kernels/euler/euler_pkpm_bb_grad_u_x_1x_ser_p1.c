#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_bb_grad_u_x_1x_ser_p1(const double *dxv, const double *bvar, const double *u_il, const double *u_ic, const double *u_ir, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]:      Cell spacing.
  // bvar:           magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_il/u_ic/u_ir: flow velocity in left/center/right cells.
  // out:            Increment to volume expansion of bb : grad_u in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[2]; 
  const double *uz_l = &u_il[4]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[2]; 
  const double *uz_c = &u_ic[4]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[2]; 
  const double *uz_r = &u_ir[4]; 

  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 

  double grad_u_x[2] = {0.0}; 
  double grad_u_y[2] = {0.0}; 
  double grad_u_z[2] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[1])-0.2886751345948129*ux_l[1]+0.5773502691896258*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = (-0.5*ux_r[1])+0.5*ux_l[1]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[1])-0.2886751345948129*uy_l[1]+0.5773502691896258*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = (-0.5*uy_r[1])+0.5*uy_l[1]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[1])-0.2886751345948129*uz_l[1]+0.5773502691896258*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = (-0.5*uz_r[1])+0.5*uz_l[1]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 

  out[0] += 0.7071067811865475*(bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  out[1] += 0.7071067811865475*(bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1])*dx1; 
} 
