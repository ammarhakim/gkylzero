#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_bb_grad_u_x_1x_ser_p2(const double *dxv, const double *bvar, const double *u_il, const double *u_ic, const double *u_ir, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]:      Cell spacing.
  // bvar:           magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_il/u_ic/u_ir: flow velocity in left/center/right cells.
  // out:            Increment to volume expansion of bb : grad_u in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[3]; 
  const double *uz_l = &u_il[6]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[3]; 
  const double *uz_c = &u_ic[6]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[3]; 
  const double *uz_r = &u_ir[6]; 

  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 

  double grad_u_x[3] = {0.0}; 
  double grad_u_y[3] = {0.0}; 
  double grad_u_z[3] = {0.0}; 
  grad_u_x[0] = 0.2445699350390395*ux_r[2]-0.2445699350390395*ux_l[2]-0.3518228202874282*ux_r[1]-0.3518228202874282*ux_l[1]+0.7036456405748563*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = 0.4236075534914363*ux_r[2]+0.4236075534914363*ux_l[2]+0.8472151069828725*ux_c[2]-0.609375*ux_r[1]+0.609375*ux_l[1]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[2] = 0.546875*ux_r[2]-0.546875*ux_l[2]-0.7866997421983816*ux_r[1]-0.7866997421983816*ux_l[1]-2.299583861810654*ux_c[1]+0.5590169943749475*ux_r[0]-0.5590169943749475*ux_l[0]; 

  grad_u_y[0] = 0.2445699350390395*uy_r[2]-0.2445699350390395*uy_l[2]-0.3518228202874282*uy_r[1]-0.3518228202874282*uy_l[1]+0.7036456405748563*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = 0.4236075534914363*uy_r[2]+0.4236075534914363*uy_l[2]+0.8472151069828725*uy_c[2]-0.609375*uy_r[1]+0.609375*uy_l[1]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[2] = 0.546875*uy_r[2]-0.546875*uy_l[2]-0.7866997421983816*uy_r[1]-0.7866997421983816*uy_l[1]-2.299583861810654*uy_c[1]+0.5590169943749475*uy_r[0]-0.5590169943749475*uy_l[0]; 

  grad_u_z[0] = 0.2445699350390395*uz_r[2]-0.2445699350390395*uz_l[2]-0.3518228202874282*uz_r[1]-0.3518228202874282*uz_l[1]+0.7036456405748563*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = 0.4236075534914363*uz_r[2]+0.4236075534914363*uz_l[2]+0.8472151069828725*uz_c[2]-0.609375*uz_r[1]+0.609375*uz_l[1]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[2] = 0.546875*uz_r[2]-0.546875*uz_l[2]-0.7866997421983816*uz_r[1]-0.7866997421983816*uz_l[1]-2.299583861810654*uz_c[1]+0.5590169943749475*uz_r[0]-0.5590169943749475*uz_l[0]; 

  out[0] += 0.7071067811865475*(bxbz[2]*grad_u_z[2]+bxby[2]*grad_u_y[2]+bxbx[2]*grad_u_x[2]+bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  out[1] += (0.6324555320336759*(bxbz[1]*grad_u_z[2]+bxby[1]*grad_u_y[2]+bxbx[1]*grad_u_x[2]+grad_u_z[1]*bxbz[2]+grad_u_y[1]*bxby[2]+grad_u_x[1]*bxbx[2])+0.7071067811865475*(bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1]))*dx1; 
  out[2] += ((0.4517539514526256*bxbz[2]+0.7071067811865475*bxbz[0])*grad_u_z[2]+(0.4517539514526256*bxby[2]+0.7071067811865475*bxby[0])*grad_u_y[2]+0.4517539514526256*bxbx[2]*grad_u_x[2]+0.7071067811865475*(bxbx[0]*grad_u_x[2]+grad_u_z[0]*bxbz[2]+grad_u_y[0]*bxby[2]+grad_u_x[0]*bxbx[2])+0.6324555320336759*(bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]))*dx1; 
} 
