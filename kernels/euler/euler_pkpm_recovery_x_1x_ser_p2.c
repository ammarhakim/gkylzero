#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p2(const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, const double *u_il, const double *u_ic, const double *u_ir, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* GKYL_RESTRICT div_b, double* GKYL_RESTRICT bb_grad_u, double* GKYL_RESTRICT div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[3]; 
  const double *uz_l = &u_il[6]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[3]; 
  const double *uz_c = &u_ic[6]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[3]; 
  const double *uz_r = &u_ir[6]; 

  const double *bxbx = &bvarc[9]; 
  const double *bxby = &bvarc[12]; 
  const double *bxbz = &bvarc[15]; 
  const double *byby = &bvarc[18]; 
  const double *bybz = &bvarc[21]; 
  const double *bzbz = &bvarc[24]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[3]; 
  const double *Pxz_l = &p_ijl[6]; 
  const double *Pyy_l = &p_ijl[9]; 
  const double *Pyz_l = &p_ijl[12]; 
  const double *Pzz_l = &p_ijl[15]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[3]; 
  const double *Pxz_c = &p_ijc[6]; 
  const double *Pyy_c = &p_ijc[9]; 
  const double *Pyz_c = &p_ijc[12]; 
  const double *Pzz_c = &p_ijc[15]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[3]; 
  const double *Pxz_r = &p_ijr[6]; 
  const double *Pyy_r = &p_ijr[9]; 
  const double *Pyz_r = &p_ijr[12]; 
  const double *Pzz_r = &p_ijr[15]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[3]; 
  double *div_p_z = &div_p[6]; 
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

  div_b[0] += (0.2445699350390395*b_r[2]-0.2445699350390395*b_l[2]-0.3518228202874282*(b_r[1]+b_l[1])+0.7036456405748563*b_c[1]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += (0.4236075534914363*(b_r[2]+b_l[2])+0.8472151069828725*b_c[2]-0.609375*b_r[1]+0.609375*b_l[1]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[2] += (0.546875*b_r[2]-0.546875*b_l[2]-0.7866997421983816*(b_r[1]+b_l[1])-2.299583861810654*b_c[1]+0.5590169943749475*b_r[0]-0.5590169943749475*b_l[0])*dx1; 

  bb_grad_u[0] += 0.7071067811865475*(bxbz[2]*grad_u_z[2]+bxby[2]*grad_u_y[2]+bxbx[2]*grad_u_x[2]+bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += (0.6324555320336759*(bxbz[1]*grad_u_z[2]+bxby[1]*grad_u_y[2]+bxbx[1]*grad_u_x[2]+grad_u_z[1]*bxbz[2]+grad_u_y[1]*bxby[2]+grad_u_x[1]*bxbx[2])+0.7071067811865475*(bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1]))*dx1; 
  bb_grad_u[2] += ((0.4517539514526256*bxbz[2]+0.7071067811865475*bxbz[0])*grad_u_z[2]+(0.4517539514526256*bxby[2]+0.7071067811865475*bxby[0])*grad_u_y[2]+0.4517539514526256*bxbx[2]*grad_u_x[2]+0.7071067811865475*(bxbx[0]*grad_u_x[2]+grad_u_z[0]*bxbz[2]+grad_u_y[0]*bxby[2]+grad_u_x[0]*bxbx[2])+0.6324555320336759*(bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]))*dx1; 

  div_p_x[0] += (0.2445699350390395*Pxx_r[2]-0.2445699350390395*Pxx_l[2]-0.3518228202874282*(Pxx_r[1]+Pxx_l[1])+0.7036456405748563*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += (0.4236075534914363*(Pxx_r[2]+Pxx_l[2])+0.8472151069828725*Pxx_c[2]-0.609375*Pxx_r[1]+0.609375*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 
  div_p_x[2] += (0.546875*Pxx_r[2]-0.546875*Pxx_l[2]-0.7866997421983816*(Pxx_r[1]+Pxx_l[1])-2.299583861810654*Pxx_c[1]+0.5590169943749475*Pxx_r[0]-0.5590169943749475*Pxx_l[0])*dx1; 

  div_p_y[0] += (0.2445699350390395*Pxy_r[2]-0.2445699350390395*Pxy_l[2]-0.3518228202874282*(Pxy_r[1]+Pxy_l[1])+0.7036456405748563*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += (0.4236075534914363*(Pxy_r[2]+Pxy_l[2])+0.8472151069828725*Pxy_c[2]-0.609375*Pxy_r[1]+0.609375*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_y[2] += (0.546875*Pxy_r[2]-0.546875*Pxy_l[2]-0.7866997421983816*(Pxy_r[1]+Pxy_l[1])-2.299583861810654*Pxy_c[1]+0.5590169943749475*Pxy_r[0]-0.5590169943749475*Pxy_l[0])*dx1; 

  div_p_z[0] += (0.2445699350390395*Pxz_r[2]-0.2445699350390395*Pxz_l[2]-0.3518228202874282*(Pxz_r[1]+Pxz_l[1])+0.7036456405748563*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += (0.4236075534914363*(Pxz_r[2]+Pxz_l[2])+0.8472151069828725*Pxz_c[2]-0.609375*Pxz_r[1]+0.609375*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_z[2] += (0.546875*Pxz_r[2]-0.546875*Pxz_l[2]-0.7866997421983816*(Pxz_r[1]+Pxz_l[1])-2.299583861810654*Pxz_c[1]+0.5590169943749475*Pxz_r[0]-0.5590169943749475*Pxz_l[0])*dx1; 

} 
