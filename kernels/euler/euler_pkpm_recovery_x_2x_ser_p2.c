#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p2(const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, const double *u_il, const double *u_ic, const double *u_ir, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* div_b, double* bb_grad_u, double* div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[8]; 
  const double *uz_l = &u_il[16]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[8]; 
  const double *uz_c = &u_ic[16]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[8]; 
  const double *uz_r = &u_ir[16]; 

  const double *bxbx = &bvarc[24]; 
  const double *bxby = &bvarc[32]; 
  const double *bxbz = &bvarc[40]; 
  const double *byby = &bvarc[48]; 
  const double *bybz = &bvarc[56]; 
  const double *bzbz = &bvarc[64]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[8]; 
  const double *Pxz_l = &p_ijl[16]; 
  const double *Pyy_l = &p_ijl[24]; 
  const double *Pyz_l = &p_ijl[32]; 
  const double *Pzz_l = &p_ijl[40]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[8]; 
  const double *Pxz_c = &p_ijc[16]; 
  const double *Pyy_c = &p_ijc[24]; 
  const double *Pyz_c = &p_ijc[32]; 
  const double *Pzz_c = &p_ijc[40]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[8]; 
  const double *Pxz_r = &p_ijr[16]; 
  const double *Pyy_r = &p_ijr[24]; 
  const double *Pyz_r = &p_ijr[32]; 
  const double *Pzz_r = &p_ijr[40]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[8]; 
  double *div_p_z = &div_p[16]; 
  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = 0.2445699350390395*ux_r[4]-0.2445699350390395*ux_l[4]-0.3518228202874282*ux_r[1]-0.3518228202874282*ux_l[1]+0.7036456405748563*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = 0.4236075534914363*ux_r[4]+0.4236075534914363*ux_l[4]+0.8472151069828725*ux_c[4]-0.609375*ux_r[1]+0.609375*ux_l[1]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[2] = 0.2445699350390395*ux_r[6]-0.2445699350390395*ux_l[6]-0.3518228202874282*ux_r[3]-0.3518228202874282*ux_l[3]+0.7036456405748563*ux_c[3]+0.25*ux_r[2]-0.25*ux_l[2]; 
  grad_u_x[3] = 0.4236075534914363*ux_r[6]+0.4236075534914363*ux_l[6]+0.8472151069828725*ux_c[6]-0.609375*ux_r[3]+0.609375*ux_l[3]+0.4330127018922193*ux_r[2]+0.4330127018922193*ux_l[2]-0.8660254037844386*ux_c[2]; 
  grad_u_x[4] = 0.546875*ux_r[4]-0.546875*ux_l[4]-0.7866997421983816*ux_r[1]-0.7866997421983816*ux_l[1]-2.299583861810654*ux_c[1]+0.5590169943749475*ux_r[0]-0.5590169943749475*ux_l[0]; 
  grad_u_x[5] = (-0.3518228202874282*ux_r[7])-0.3518228202874282*ux_l[7]+0.7036456405748563*ux_c[7]+0.25*ux_r[5]-0.25*ux_l[5]; 
  grad_u_x[6] = 0.546875*ux_r[6]-0.546875*ux_l[6]-0.7866997421983816*ux_r[3]-0.7866997421983816*ux_l[3]-2.299583861810654*ux_c[3]+0.5590169943749476*ux_r[2]-0.5590169943749476*ux_l[2]; 
  grad_u_x[7] = (-0.609375*ux_r[7])+0.609375*ux_l[7]+0.4330127018922194*ux_r[5]+0.4330127018922194*ux_l[5]-0.8660254037844387*ux_c[5]; 

  grad_u_y[0] = 0.2445699350390395*uy_r[4]-0.2445699350390395*uy_l[4]-0.3518228202874282*uy_r[1]-0.3518228202874282*uy_l[1]+0.7036456405748563*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = 0.4236075534914363*uy_r[4]+0.4236075534914363*uy_l[4]+0.8472151069828725*uy_c[4]-0.609375*uy_r[1]+0.609375*uy_l[1]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[2] = 0.2445699350390395*uy_r[6]-0.2445699350390395*uy_l[6]-0.3518228202874282*uy_r[3]-0.3518228202874282*uy_l[3]+0.7036456405748563*uy_c[3]+0.25*uy_r[2]-0.25*uy_l[2]; 
  grad_u_y[3] = 0.4236075534914363*uy_r[6]+0.4236075534914363*uy_l[6]+0.8472151069828725*uy_c[6]-0.609375*uy_r[3]+0.609375*uy_l[3]+0.4330127018922193*uy_r[2]+0.4330127018922193*uy_l[2]-0.8660254037844386*uy_c[2]; 
  grad_u_y[4] = 0.546875*uy_r[4]-0.546875*uy_l[4]-0.7866997421983816*uy_r[1]-0.7866997421983816*uy_l[1]-2.299583861810654*uy_c[1]+0.5590169943749475*uy_r[0]-0.5590169943749475*uy_l[0]; 
  grad_u_y[5] = (-0.3518228202874282*uy_r[7])-0.3518228202874282*uy_l[7]+0.7036456405748563*uy_c[7]+0.25*uy_r[5]-0.25*uy_l[5]; 
  grad_u_y[6] = 0.546875*uy_r[6]-0.546875*uy_l[6]-0.7866997421983816*uy_r[3]-0.7866997421983816*uy_l[3]-2.299583861810654*uy_c[3]+0.5590169943749476*uy_r[2]-0.5590169943749476*uy_l[2]; 
  grad_u_y[7] = (-0.609375*uy_r[7])+0.609375*uy_l[7]+0.4330127018922194*uy_r[5]+0.4330127018922194*uy_l[5]-0.8660254037844387*uy_c[5]; 

  grad_u_z[0] = 0.2445699350390395*uz_r[4]-0.2445699350390395*uz_l[4]-0.3518228202874282*uz_r[1]-0.3518228202874282*uz_l[1]+0.7036456405748563*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = 0.4236075534914363*uz_r[4]+0.4236075534914363*uz_l[4]+0.8472151069828725*uz_c[4]-0.609375*uz_r[1]+0.609375*uz_l[1]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[2] = 0.2445699350390395*uz_r[6]-0.2445699350390395*uz_l[6]-0.3518228202874282*uz_r[3]-0.3518228202874282*uz_l[3]+0.7036456405748563*uz_c[3]+0.25*uz_r[2]-0.25*uz_l[2]; 
  grad_u_z[3] = 0.4236075534914363*uz_r[6]+0.4236075534914363*uz_l[6]+0.8472151069828725*uz_c[6]-0.609375*uz_r[3]+0.609375*uz_l[3]+0.4330127018922193*uz_r[2]+0.4330127018922193*uz_l[2]-0.8660254037844386*uz_c[2]; 
  grad_u_z[4] = 0.546875*uz_r[4]-0.546875*uz_l[4]-0.7866997421983816*uz_r[1]-0.7866997421983816*uz_l[1]-2.299583861810654*uz_c[1]+0.5590169943749475*uz_r[0]-0.5590169943749475*uz_l[0]; 
  grad_u_z[5] = (-0.3518228202874282*uz_r[7])-0.3518228202874282*uz_l[7]+0.7036456405748563*uz_c[7]+0.25*uz_r[5]-0.25*uz_l[5]; 
  grad_u_z[6] = 0.546875*uz_r[6]-0.546875*uz_l[6]-0.7866997421983816*uz_r[3]-0.7866997421983816*uz_l[3]-2.299583861810654*uz_c[3]+0.5590169943749476*uz_r[2]-0.5590169943749476*uz_l[2]; 
  grad_u_z[7] = (-0.609375*uz_r[7])+0.609375*uz_l[7]+0.4330127018922194*uz_r[5]+0.4330127018922194*uz_l[5]-0.8660254037844387*uz_c[5]; 

  div_b[0] += (0.2445699350390395*b_r[4]-0.2445699350390395*b_l[4]-0.3518228202874282*(b_r[1]+b_l[1])+0.7036456405748563*b_c[1]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += (0.4236075534914363*(b_r[4]+b_l[4])+0.8472151069828725*b_c[4]-0.609375*b_r[1]+0.609375*b_l[1]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[2] += (0.2445699350390395*b_r[6]-0.2445699350390395*b_l[6]-0.3518228202874282*(b_r[3]+b_l[3])+0.7036456405748563*b_c[3]+0.25*b_r[2]-0.25*b_l[2])*dx1; 
  div_b[3] += (0.4236075534914363*(b_r[6]+b_l[6])+0.8472151069828725*b_c[6]-0.609375*b_r[3]+0.609375*b_l[3]+0.4330127018922193*(b_r[2]+b_l[2])-0.8660254037844386*b_c[2])*dx1; 
  div_b[4] += (0.546875*b_r[4]-0.546875*b_l[4]-0.7866997421983816*(b_r[1]+b_l[1])-2.299583861810654*b_c[1]+0.5590169943749475*b_r[0]-0.5590169943749475*b_l[0])*dx1; 
  div_b[5] += ((-0.3518228202874282*(b_r[7]+b_l[7]))+0.7036456405748563*b_c[7]+0.25*b_r[5]-0.25*b_l[5])*dx1; 
  div_b[6] += (0.546875*b_r[6]-0.546875*b_l[6]-0.7866997421983816*(b_r[3]+b_l[3])-2.299583861810654*b_c[3]+0.5590169943749476*b_r[2]-0.5590169943749476*b_l[2])*dx1; 
  div_b[7] += ((-0.609375*b_r[7])+0.609375*b_l[7]+0.4330127018922194*(b_r[5]+b_l[5])-0.8660254037844387*b_c[5])*dx1; 

  bb_grad_u[0] += 0.5*(bxbz[7]*grad_u_z[7]+bxby[7]*grad_u_y[7]+bxbx[7]*grad_u_x[7]+bxbz[6]*grad_u_z[6]+bxby[6]*grad_u_y[6]+bxbx[6]*grad_u_x[6]+bxbz[5]*grad_u_z[5]+bxby[5]*grad_u_y[5]+bxbx[5]*grad_u_x[5]+bxbz[4]*grad_u_z[4]+bxby[4]*grad_u_y[4]+bxbx[4]*grad_u_x[4]+bxbz[3]*grad_u_z[3]+bxby[3]*grad_u_y[3]+bxbx[3]*grad_u_x[3]+bxbz[2]*grad_u_z[2]+bxby[2]*grad_u_y[2]+bxbx[2]*grad_u_x[2]+bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += (0.5000000000000001*(bxbz[5]*grad_u_z[7]+bxby[5]*grad_u_y[7]+bxbx[5]*grad_u_x[7]+grad_u_z[5]*bxbz[7]+grad_u_y[5]*bxby[7]+grad_u_x[5]*bxbx[7])+0.447213595499958*(bxbz[3]*grad_u_z[6]+bxby[3]*grad_u_y[6]+bxbx[3]*grad_u_x[6]+grad_u_z[3]*bxbz[6]+grad_u_y[3]*bxby[6]+grad_u_x[3]*bxbx[6])+0.4472135954999579*(bxbz[1]*grad_u_z[4]+bxby[1]*grad_u_y[4]+bxbx[1]*grad_u_x[4]+grad_u_z[1]*bxbz[4]+grad_u_y[1]*bxby[4]+grad_u_x[1]*bxbx[4])+0.5*(bxbz[2]*grad_u_z[3]+bxby[2]*grad_u_y[3]+bxbx[2]*grad_u_x[3]+grad_u_z[2]*bxbz[3]+grad_u_y[2]*bxby[3]+grad_u_x[2]*bxbx[3]+bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1]))*dx1; 
  bb_grad_u[2] += (0.447213595499958*(bxbz[3]*grad_u_z[7]+bxby[3]*grad_u_y[7]+bxbx[3]*grad_u_x[7]+grad_u_z[3]*bxbz[7]+grad_u_y[3]*bxby[7]+grad_u_x[3]*bxbx[7])+0.5000000000000001*(bxbz[4]*grad_u_z[6]+bxby[4]*grad_u_y[6]+bxbx[4]*grad_u_x[6]+grad_u_z[4]*bxbz[6]+grad_u_y[4]*bxby[6]+grad_u_x[4]*bxbx[6])+0.4472135954999579*(bxbz[2]*grad_u_z[5]+bxby[2]*grad_u_y[5]+bxbx[2]*grad_u_x[5]+grad_u_z[2]*bxbz[5]+grad_u_y[2]*bxby[5]+grad_u_x[2]*bxbx[5])+0.5*(bxbz[1]*grad_u_z[3]+bxby[1]*grad_u_y[3]+bxbx[1]*grad_u_x[3]+grad_u_z[1]*bxbz[3]+grad_u_y[1]*bxby[3]+grad_u_x[1]*bxbx[3]+bxbz[0]*grad_u_z[2]+bxby[0]*grad_u_y[2]+bxbx[0]*grad_u_x[2]+grad_u_z[0]*bxbz[2]+grad_u_y[0]*bxby[2]+grad_u_x[0]*bxbx[2]))*dx1; 
  bb_grad_u[3] += ((0.4*bxbz[6]+0.447213595499958*bxbz[2])*grad_u_z[7]+(0.4*bxby[6]+0.447213595499958*bxby[2])*grad_u_y[7]+(0.4*bxbx[6]+0.447213595499958*bxbx[2])*grad_u_x[7]+(0.4*grad_u_z[6]+0.447213595499958*grad_u_z[2])*bxbz[7]+(0.4*grad_u_y[6]+0.447213595499958*grad_u_y[2])*bxby[7]+0.4*grad_u_x[6]*bxbx[7]+0.447213595499958*(grad_u_x[2]*bxbx[7]+bxbz[1]*grad_u_z[6]+bxby[1]*grad_u_y[6]+bxbx[1]*grad_u_x[6]+grad_u_z[1]*bxbz[6]+grad_u_y[1]*bxby[6]+grad_u_x[1]*bxbx[6])+0.4472135954999579*(bxbz[3]*grad_u_z[5]+bxby[3]*grad_u_y[5]+bxbx[3]*grad_u_x[5]+grad_u_z[3]*bxbz[5]+grad_u_y[3]*bxby[5]+grad_u_x[3]*bxbx[5]+bxbz[3]*grad_u_z[4]+bxby[3]*grad_u_y[4]+bxbx[3]*grad_u_x[4]+grad_u_z[3]*bxbz[4]+grad_u_y[3]*bxby[4]+grad_u_x[3]*bxbx[4])+0.5*(bxbz[0]*grad_u_z[3]+bxby[0]*grad_u_y[3]+bxbx[0]*grad_u_x[3]+grad_u_z[0]*bxbz[3]+grad_u_y[0]*bxby[3]+grad_u_x[0]*bxbx[3]+bxbz[1]*grad_u_z[2]+bxby[1]*grad_u_y[2]+bxbx[1]*grad_u_x[2]+grad_u_z[1]*bxbz[2]+grad_u_y[1]*bxby[2]+grad_u_x[1]*bxbx[2]))*dx1; 
  bb_grad_u[4] += (0.4472135954999579*(bxbz[7]*grad_u_z[7]+bxby[7]*grad_u_y[7]+bxbx[7]*grad_u_x[7])+(0.31943828249997*bxbz[6]+0.5000000000000001*bxbz[2])*grad_u_z[6]+(0.31943828249997*bxby[6]+0.5000000000000001*bxby[2])*grad_u_y[6]+0.31943828249997*bxbx[6]*grad_u_x[6]+0.5000000000000001*(bxbx[2]*grad_u_x[6]+grad_u_z[2]*bxbz[6]+grad_u_y[2]*bxby[6]+grad_u_x[2]*bxbx[6])+(0.31943828249997*bxbz[4]+0.5*bxbz[0])*grad_u_z[4]+(0.31943828249997*bxby[4]+0.5*bxby[0])*grad_u_y[4]+0.31943828249997*bxbx[4]*grad_u_x[4]+0.5*(bxbx[0]*grad_u_x[4]+grad_u_z[0]*bxbz[4]+grad_u_y[0]*bxby[4]+grad_u_x[0]*bxbx[4])+0.4472135954999579*(bxbz[3]*grad_u_z[3]+bxby[3]*grad_u_y[3]+bxbx[3]*grad_u_x[3]+bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]))*dx1; 
  bb_grad_u[5] += ((0.31943828249997*bxbz[7]+0.5000000000000001*bxbz[1])*grad_u_z[7]+(0.31943828249997*bxby[7]+0.5000000000000001*bxby[1])*grad_u_y[7]+0.31943828249997*bxbx[7]*grad_u_x[7]+0.5000000000000001*(bxbx[1]*grad_u_x[7]+grad_u_z[1]*bxbz[7]+grad_u_y[1]*bxby[7]+grad_u_x[1]*bxbx[7])+0.4472135954999579*(bxbz[6]*grad_u_z[6]+bxby[6]*grad_u_y[6]+bxbx[6]*grad_u_x[6])+(0.31943828249997*bxbz[5]+0.5*bxbz[0])*grad_u_z[5]+(0.31943828249997*bxby[5]+0.5*bxby[0])*grad_u_y[5]+0.31943828249997*bxbx[5]*grad_u_x[5]+0.5*(bxbx[0]*grad_u_x[5]+grad_u_z[0]*bxbz[5]+grad_u_y[0]*bxby[5]+grad_u_x[0]*bxbx[5])+0.4472135954999579*(bxbz[3]*grad_u_z[3]+bxby[3]*grad_u_y[3]+bxbx[3]*grad_u_x[3]+bxbz[2]*grad_u_z[2]+bxby[2]*grad_u_y[2]+bxbx[2]*grad_u_x[2]))*dx1; 
  bb_grad_u[6] += (0.4*(bxbz[3]*grad_u_z[7]+bxby[3]*grad_u_y[7]+bxbx[3]*grad_u_x[7]+grad_u_z[3]*bxbz[7]+grad_u_y[3]*bxby[7]+grad_u_x[3]*bxbx[7])+(0.4472135954999579*bxbz[5]+0.31943828249997*bxbz[4]+0.5*bxbz[0])*grad_u_z[6]+(0.4472135954999579*bxby[5]+0.31943828249997*bxby[4]+0.5*bxby[0])*grad_u_y[6]+(0.4472135954999579*bxbx[5]+0.31943828249997*bxbx[4]+0.5*bxbx[0])*grad_u_x[6]+(0.4472135954999579*grad_u_z[5]+0.31943828249997*grad_u_z[4]+0.5*grad_u_z[0])*bxbz[6]+(0.4472135954999579*grad_u_y[5]+0.31943828249997*grad_u_y[4]+0.5*grad_u_y[0])*bxby[6]+(0.4472135954999579*grad_u_x[5]+0.31943828249997*grad_u_x[4]+0.5*grad_u_x[0])*bxbx[6]+0.5000000000000001*(bxbz[2]*grad_u_z[4]+bxby[2]*grad_u_y[4]+bxbx[2]*grad_u_x[4]+grad_u_z[2]*bxbz[4]+grad_u_y[2]*bxby[4]+grad_u_x[2]*bxbx[4])+0.447213595499958*(bxbz[1]*grad_u_z[3]+bxby[1]*grad_u_y[3]+bxbx[1]*grad_u_x[3]+grad_u_z[1]*bxbz[3]+grad_u_y[1]*bxby[3]+grad_u_x[1]*bxbx[3]))*dx1; 
  bb_grad_u[7] += ((0.31943828249997*bxbz[5]+0.4472135954999579*bxbz[4]+0.5*bxbz[0])*grad_u_z[7]+(0.31943828249997*bxby[5]+0.4472135954999579*bxby[4]+0.5*bxby[0])*grad_u_y[7]+(0.31943828249997*bxbx[5]+0.4472135954999579*bxbx[4]+0.5*bxbx[0])*grad_u_x[7]+(0.31943828249997*grad_u_z[5]+0.4472135954999579*grad_u_z[4]+0.5*grad_u_z[0])*bxbz[7]+(0.31943828249997*grad_u_y[5]+0.4472135954999579*grad_u_y[4]+0.5*grad_u_y[0])*bxby[7]+(0.31943828249997*grad_u_x[5]+0.4472135954999579*grad_u_x[4]+0.5*grad_u_x[0])*bxbx[7]+0.4*(bxbz[3]*grad_u_z[6]+bxby[3]*grad_u_y[6]+bxbx[3]*grad_u_x[6]+grad_u_z[3]*bxbz[6]+grad_u_y[3]*bxby[6]+grad_u_x[3]*bxbx[6])+0.5000000000000001*(bxbz[1]*grad_u_z[5]+bxby[1]*grad_u_y[5]+bxbx[1]*grad_u_x[5]+grad_u_z[1]*bxbz[5]+grad_u_y[1]*bxby[5]+grad_u_x[1]*bxbx[5])+0.447213595499958*(bxbz[2]*grad_u_z[3]+bxby[2]*grad_u_y[3]+bxbx[2]*grad_u_x[3]+grad_u_z[2]*bxbz[3]+grad_u_y[2]*bxby[3]+grad_u_x[2]*bxbx[3]))*dx1; 

  div_p_x[0] += (0.2445699350390395*Pxx_r[4]-0.2445699350390395*Pxx_l[4]-0.3518228202874282*(Pxx_r[1]+Pxx_l[1])+0.7036456405748563*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += (0.4236075534914363*(Pxx_r[4]+Pxx_l[4])+0.8472151069828725*Pxx_c[4]-0.609375*Pxx_r[1]+0.609375*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 
  div_p_x[2] += (0.2445699350390395*Pxx_r[6]-0.2445699350390395*Pxx_l[6]-0.3518228202874282*(Pxx_r[3]+Pxx_l[3])+0.7036456405748563*Pxx_c[3]+0.25*Pxx_r[2]-0.25*Pxx_l[2])*dx1; 
  div_p_x[3] += (0.4236075534914363*(Pxx_r[6]+Pxx_l[6])+0.8472151069828725*Pxx_c[6]-0.609375*Pxx_r[3]+0.609375*Pxx_l[3]+0.4330127018922193*(Pxx_r[2]+Pxx_l[2])-0.8660254037844386*Pxx_c[2])*dx1; 
  div_p_x[4] += (0.546875*Pxx_r[4]-0.546875*Pxx_l[4]-0.7866997421983816*(Pxx_r[1]+Pxx_l[1])-2.299583861810654*Pxx_c[1]+0.5590169943749475*Pxx_r[0]-0.5590169943749475*Pxx_l[0])*dx1; 
  div_p_x[5] += ((-0.3518228202874282*(Pxx_r[7]+Pxx_l[7]))+0.7036456405748563*Pxx_c[7]+0.25*Pxx_r[5]-0.25*Pxx_l[5])*dx1; 
  div_p_x[6] += (0.546875*Pxx_r[6]-0.546875*Pxx_l[6]-0.7866997421983816*(Pxx_r[3]+Pxx_l[3])-2.299583861810654*Pxx_c[3]+0.5590169943749476*Pxx_r[2]-0.5590169943749476*Pxx_l[2])*dx1; 
  div_p_x[7] += ((-0.609375*Pxx_r[7])+0.609375*Pxx_l[7]+0.4330127018922194*(Pxx_r[5]+Pxx_l[5])-0.8660254037844387*Pxx_c[5])*dx1; 

  div_p_y[0] += (0.2445699350390395*Pxy_r[4]-0.2445699350390395*Pxy_l[4]-0.3518228202874282*(Pxy_r[1]+Pxy_l[1])+0.7036456405748563*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += (0.4236075534914363*(Pxy_r[4]+Pxy_l[4])+0.8472151069828725*Pxy_c[4]-0.609375*Pxy_r[1]+0.609375*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_y[2] += (0.2445699350390395*Pxy_r[6]-0.2445699350390395*Pxy_l[6]-0.3518228202874282*(Pxy_r[3]+Pxy_l[3])+0.7036456405748563*Pxy_c[3]+0.25*Pxy_r[2]-0.25*Pxy_l[2])*dx1; 
  div_p_y[3] += (0.4236075534914363*(Pxy_r[6]+Pxy_l[6])+0.8472151069828725*Pxy_c[6]-0.609375*Pxy_r[3]+0.609375*Pxy_l[3]+0.4330127018922193*(Pxy_r[2]+Pxy_l[2])-0.8660254037844386*Pxy_c[2])*dx1; 
  div_p_y[4] += (0.546875*Pxy_r[4]-0.546875*Pxy_l[4]-0.7866997421983816*(Pxy_r[1]+Pxy_l[1])-2.299583861810654*Pxy_c[1]+0.5590169943749475*Pxy_r[0]-0.5590169943749475*Pxy_l[0])*dx1; 
  div_p_y[5] += ((-0.3518228202874282*(Pxy_r[7]+Pxy_l[7]))+0.7036456405748563*Pxy_c[7]+0.25*Pxy_r[5]-0.25*Pxy_l[5])*dx1; 
  div_p_y[6] += (0.546875*Pxy_r[6]-0.546875*Pxy_l[6]-0.7866997421983816*(Pxy_r[3]+Pxy_l[3])-2.299583861810654*Pxy_c[3]+0.5590169943749476*Pxy_r[2]-0.5590169943749476*Pxy_l[2])*dx1; 
  div_p_y[7] += ((-0.609375*Pxy_r[7])+0.609375*Pxy_l[7]+0.4330127018922194*(Pxy_r[5]+Pxy_l[5])-0.8660254037844387*Pxy_c[5])*dx1; 

  div_p_z[0] += (0.2445699350390395*Pxz_r[4]-0.2445699350390395*Pxz_l[4]-0.3518228202874282*(Pxz_r[1]+Pxz_l[1])+0.7036456405748563*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += (0.4236075534914363*(Pxz_r[4]+Pxz_l[4])+0.8472151069828725*Pxz_c[4]-0.609375*Pxz_r[1]+0.609375*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_z[2] += (0.2445699350390395*Pxz_r[6]-0.2445699350390395*Pxz_l[6]-0.3518228202874282*(Pxz_r[3]+Pxz_l[3])+0.7036456405748563*Pxz_c[3]+0.25*Pxz_r[2]-0.25*Pxz_l[2])*dx1; 
  div_p_z[3] += (0.4236075534914363*(Pxz_r[6]+Pxz_l[6])+0.8472151069828725*Pxz_c[6]-0.609375*Pxz_r[3]+0.609375*Pxz_l[3]+0.4330127018922193*(Pxz_r[2]+Pxz_l[2])-0.8660254037844386*Pxz_c[2])*dx1; 
  div_p_z[4] += (0.546875*Pxz_r[4]-0.546875*Pxz_l[4]-0.7866997421983816*(Pxz_r[1]+Pxz_l[1])-2.299583861810654*Pxz_c[1]+0.5590169943749475*Pxz_r[0]-0.5590169943749475*Pxz_l[0])*dx1; 
  div_p_z[5] += ((-0.3518228202874282*(Pxz_r[7]+Pxz_l[7]))+0.7036456405748563*Pxz_c[7]+0.25*Pxz_r[5]-0.25*Pxz_l[5])*dx1; 
  div_p_z[6] += (0.546875*Pxz_r[6]-0.546875*Pxz_l[6]-0.7866997421983816*(Pxz_r[3]+Pxz_l[3])-2.299583861810654*Pxz_c[3]+0.5590169943749476*Pxz_r[2]-0.5590169943749476*Pxz_l[2])*dx1; 
  div_p_z[7] += ((-0.609375*Pxz_r[7])+0.609375*Pxz_l[7]+0.4330127018922194*(Pxz_r[5]+Pxz_l[5])-0.8660254037844387*Pxz_c[5])*dx1; 

} 
