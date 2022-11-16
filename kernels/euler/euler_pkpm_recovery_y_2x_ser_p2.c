#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p2(const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, const double *u_il, const double *u_ic, const double *u_ir, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* GKYL_RESTRICT div_b, double* GKYL_RESTRICT bb_grad_u, double* GKYL_RESTRICT div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[1]; 

  const double *b_l = &bvarl[8]; 
  const double *b_c = &bvarc[8]; 
  const double *b_r = &bvarr[8]; 

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
  grad_u_x[0] = 0.2445699350390395*ux_r[5]-0.2445699350390395*ux_l[5]-0.3518228202874282*ux_r[2]-0.3518228202874282*ux_l[2]+0.7036456405748563*ux_c[2]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = 0.2445699350390395*ux_r[7]-0.2445699350390395*ux_l[7]-0.3518228202874282*ux_r[3]-0.3518228202874282*ux_l[3]+0.7036456405748563*ux_c[3]+0.25*ux_r[1]-0.25*ux_l[1]; 
  grad_u_x[2] = 0.4236075534914363*ux_r[5]+0.4236075534914363*ux_l[5]+0.8472151069828725*ux_c[5]-0.609375*ux_r[2]+0.609375*ux_l[2]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[3] = 0.4236075534914363*ux_r[7]+0.4236075534914363*ux_l[7]+0.8472151069828725*ux_c[7]-0.609375*ux_r[3]+0.609375*ux_l[3]+0.4330127018922193*ux_r[1]+0.4330127018922193*ux_l[1]-0.8660254037844386*ux_c[1]; 
  grad_u_x[4] = (-0.3518228202874282*ux_r[6])-0.3518228202874282*ux_l[6]+0.7036456405748563*ux_c[6]+0.25*ux_r[4]-0.25*ux_l[4]; 
  grad_u_x[5] = 0.546875*ux_r[5]-0.546875*ux_l[5]-0.7866997421983816*ux_r[2]-0.7866997421983816*ux_l[2]-2.299583861810654*ux_c[2]+0.5590169943749475*ux_r[0]-0.5590169943749475*ux_l[0]; 
  grad_u_x[6] = (-0.609375*ux_r[6])+0.609375*ux_l[6]+0.4330127018922194*ux_r[4]+0.4330127018922194*ux_l[4]-0.8660254037844387*ux_c[4]; 
  grad_u_x[7] = 0.546875*ux_r[7]-0.546875*ux_l[7]-0.7866997421983816*ux_r[3]-0.7866997421983816*ux_l[3]-2.299583861810654*ux_c[3]+0.5590169943749476*ux_r[1]-0.5590169943749476*ux_l[1]; 

  grad_u_y[0] = 0.2445699350390395*uy_r[5]-0.2445699350390395*uy_l[5]-0.3518228202874282*uy_r[2]-0.3518228202874282*uy_l[2]+0.7036456405748563*uy_c[2]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = 0.2445699350390395*uy_r[7]-0.2445699350390395*uy_l[7]-0.3518228202874282*uy_r[3]-0.3518228202874282*uy_l[3]+0.7036456405748563*uy_c[3]+0.25*uy_r[1]-0.25*uy_l[1]; 
  grad_u_y[2] = 0.4236075534914363*uy_r[5]+0.4236075534914363*uy_l[5]+0.8472151069828725*uy_c[5]-0.609375*uy_r[2]+0.609375*uy_l[2]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[3] = 0.4236075534914363*uy_r[7]+0.4236075534914363*uy_l[7]+0.8472151069828725*uy_c[7]-0.609375*uy_r[3]+0.609375*uy_l[3]+0.4330127018922193*uy_r[1]+0.4330127018922193*uy_l[1]-0.8660254037844386*uy_c[1]; 
  grad_u_y[4] = (-0.3518228202874282*uy_r[6])-0.3518228202874282*uy_l[6]+0.7036456405748563*uy_c[6]+0.25*uy_r[4]-0.25*uy_l[4]; 
  grad_u_y[5] = 0.546875*uy_r[5]-0.546875*uy_l[5]-0.7866997421983816*uy_r[2]-0.7866997421983816*uy_l[2]-2.299583861810654*uy_c[2]+0.5590169943749475*uy_r[0]-0.5590169943749475*uy_l[0]; 
  grad_u_y[6] = (-0.609375*uy_r[6])+0.609375*uy_l[6]+0.4330127018922194*uy_r[4]+0.4330127018922194*uy_l[4]-0.8660254037844387*uy_c[4]; 
  grad_u_y[7] = 0.546875*uy_r[7]-0.546875*uy_l[7]-0.7866997421983816*uy_r[3]-0.7866997421983816*uy_l[3]-2.299583861810654*uy_c[3]+0.5590169943749476*uy_r[1]-0.5590169943749476*uy_l[1]; 

  grad_u_z[0] = 0.2445699350390395*uz_r[5]-0.2445699350390395*uz_l[5]-0.3518228202874282*uz_r[2]-0.3518228202874282*uz_l[2]+0.7036456405748563*uz_c[2]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = 0.2445699350390395*uz_r[7]-0.2445699350390395*uz_l[7]-0.3518228202874282*uz_r[3]-0.3518228202874282*uz_l[3]+0.7036456405748563*uz_c[3]+0.25*uz_r[1]-0.25*uz_l[1]; 
  grad_u_z[2] = 0.4236075534914363*uz_r[5]+0.4236075534914363*uz_l[5]+0.8472151069828725*uz_c[5]-0.609375*uz_r[2]+0.609375*uz_l[2]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[3] = 0.4236075534914363*uz_r[7]+0.4236075534914363*uz_l[7]+0.8472151069828725*uz_c[7]-0.609375*uz_r[3]+0.609375*uz_l[3]+0.4330127018922193*uz_r[1]+0.4330127018922193*uz_l[1]-0.8660254037844386*uz_c[1]; 
  grad_u_z[4] = (-0.3518228202874282*uz_r[6])-0.3518228202874282*uz_l[6]+0.7036456405748563*uz_c[6]+0.25*uz_r[4]-0.25*uz_l[4]; 
  grad_u_z[5] = 0.546875*uz_r[5]-0.546875*uz_l[5]-0.7866997421983816*uz_r[2]-0.7866997421983816*uz_l[2]-2.299583861810654*uz_c[2]+0.5590169943749475*uz_r[0]-0.5590169943749475*uz_l[0]; 
  grad_u_z[6] = (-0.609375*uz_r[6])+0.609375*uz_l[6]+0.4330127018922194*uz_r[4]+0.4330127018922194*uz_l[4]-0.8660254037844387*uz_c[4]; 
  grad_u_z[7] = 0.546875*uz_r[7]-0.546875*uz_l[7]-0.7866997421983816*uz_r[3]-0.7866997421983816*uz_l[3]-2.299583861810654*uz_c[3]+0.5590169943749476*uz_r[1]-0.5590169943749476*uz_l[1]; 

  div_b[0] += (0.2445699350390395*b_r[5]-0.2445699350390395*b_l[5]-0.3518228202874282*(b_r[2]+b_l[2])+0.7036456405748563*b_c[2]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += (0.2445699350390395*b_r[7]-0.2445699350390395*b_l[7]-0.3518228202874282*(b_r[3]+b_l[3])+0.7036456405748563*b_c[3]+0.25*b_r[1]-0.25*b_l[1])*dx1; 
  div_b[2] += (0.4236075534914363*(b_r[5]+b_l[5])+0.8472151069828725*b_c[5]-0.609375*b_r[2]+0.609375*b_l[2]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[3] += (0.4236075534914363*(b_r[7]+b_l[7])+0.8472151069828725*b_c[7]-0.609375*b_r[3]+0.609375*b_l[3]+0.4330127018922193*(b_r[1]+b_l[1])-0.8660254037844386*b_c[1])*dx1; 
  div_b[4] += ((-0.3518228202874282*(b_r[6]+b_l[6]))+0.7036456405748563*b_c[6]+0.25*b_r[4]-0.25*b_l[4])*dx1; 
  div_b[5] += (0.546875*b_r[5]-0.546875*b_l[5]-0.7866997421983816*(b_r[2]+b_l[2])-2.299583861810654*b_c[2]+0.5590169943749475*b_r[0]-0.5590169943749475*b_l[0])*dx1; 
  div_b[6] += ((-0.609375*b_r[6])+0.609375*b_l[6]+0.4330127018922194*(b_r[4]+b_l[4])-0.8660254037844387*b_c[4])*dx1; 
  div_b[7] += (0.546875*b_r[7]-0.546875*b_l[7]-0.7866997421983816*(b_r[3]+b_l[3])-2.299583861810654*b_c[3]+0.5590169943749476*b_r[1]-0.5590169943749476*b_l[1])*dx1; 

  bb_grad_u[0] += 0.5*(bybz[7]*grad_u_z[7]+byby[7]*grad_u_y[7]+bxby[7]*grad_u_x[7]+bybz[6]*grad_u_z[6]+byby[6]*grad_u_y[6]+bxby[6]*grad_u_x[6]+bybz[5]*grad_u_z[5]+byby[5]*grad_u_y[5]+bxby[5]*grad_u_x[5]+bybz[4]*grad_u_z[4]+byby[4]*grad_u_y[4]+bxby[4]*grad_u_x[4]+bybz[3]*grad_u_z[3]+byby[3]*grad_u_y[3]+bxby[3]*grad_u_x[3]+bybz[2]*grad_u_z[2]+byby[2]*grad_u_y[2]+bxby[2]*grad_u_x[2]+bybz[1]*grad_u_z[1]+byby[1]*grad_u_y[1]+bxby[1]*grad_u_x[1]+bybz[0]*grad_u_z[0]+byby[0]*grad_u_y[0]+bxby[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += (0.5000000000000001*(bybz[5]*grad_u_z[7]+byby[5]*grad_u_y[7]+bxby[5]*grad_u_x[7]+grad_u_z[5]*bybz[7]+grad_u_y[5]*byby[7]+grad_u_x[5]*bxby[7])+0.447213595499958*(bybz[3]*grad_u_z[6]+byby[3]*grad_u_y[6]+bxby[3]*grad_u_x[6]+grad_u_z[3]*bybz[6]+grad_u_y[3]*byby[6]+grad_u_x[3]*bxby[6])+0.4472135954999579*(bybz[1]*grad_u_z[4]+byby[1]*grad_u_y[4]+bxby[1]*grad_u_x[4]+grad_u_z[1]*bybz[4]+grad_u_y[1]*byby[4]+grad_u_x[1]*bxby[4])+0.5*(bybz[2]*grad_u_z[3]+byby[2]*grad_u_y[3]+bxby[2]*grad_u_x[3]+grad_u_z[2]*bybz[3]+grad_u_y[2]*byby[3]+grad_u_x[2]*bxby[3]+bybz[0]*grad_u_z[1]+byby[0]*grad_u_y[1]+bxby[0]*grad_u_x[1]+grad_u_z[0]*bybz[1]+grad_u_y[0]*byby[1]+grad_u_x[0]*bxby[1]))*dx1; 
  bb_grad_u[2] += (0.447213595499958*(bybz[3]*grad_u_z[7]+byby[3]*grad_u_y[7]+bxby[3]*grad_u_x[7]+grad_u_z[3]*bybz[7]+grad_u_y[3]*byby[7]+grad_u_x[3]*bxby[7])+0.5000000000000001*(bybz[4]*grad_u_z[6]+byby[4]*grad_u_y[6]+bxby[4]*grad_u_x[6]+grad_u_z[4]*bybz[6]+grad_u_y[4]*byby[6]+grad_u_x[4]*bxby[6])+0.4472135954999579*(bybz[2]*grad_u_z[5]+byby[2]*grad_u_y[5]+bxby[2]*grad_u_x[5]+grad_u_z[2]*bybz[5]+grad_u_y[2]*byby[5]+grad_u_x[2]*bxby[5])+0.5*(bybz[1]*grad_u_z[3]+byby[1]*grad_u_y[3]+bxby[1]*grad_u_x[3]+grad_u_z[1]*bybz[3]+grad_u_y[1]*byby[3]+grad_u_x[1]*bxby[3]+bybz[0]*grad_u_z[2]+byby[0]*grad_u_y[2]+bxby[0]*grad_u_x[2]+grad_u_z[0]*bybz[2]+grad_u_y[0]*byby[2]+grad_u_x[0]*bxby[2]))*dx1; 
  bb_grad_u[3] += ((0.4*bybz[6]+0.447213595499958*bybz[2])*grad_u_z[7]+(0.4*byby[6]+0.447213595499958*byby[2])*grad_u_y[7]+(0.4*bxby[6]+0.447213595499958*bxby[2])*grad_u_x[7]+(0.4*grad_u_z[6]+0.447213595499958*grad_u_z[2])*bybz[7]+(0.4*grad_u_y[6]+0.447213595499958*grad_u_y[2])*byby[7]+0.4*grad_u_x[6]*bxby[7]+0.447213595499958*(grad_u_x[2]*bxby[7]+bybz[1]*grad_u_z[6]+byby[1]*grad_u_y[6]+bxby[1]*grad_u_x[6]+grad_u_z[1]*bybz[6]+grad_u_y[1]*byby[6]+grad_u_x[1]*bxby[6])+0.4472135954999579*(bybz[3]*grad_u_z[5]+byby[3]*grad_u_y[5]+bxby[3]*grad_u_x[5]+grad_u_z[3]*bybz[5]+grad_u_y[3]*byby[5]+grad_u_x[3]*bxby[5]+bybz[3]*grad_u_z[4]+byby[3]*grad_u_y[4]+bxby[3]*grad_u_x[4]+grad_u_z[3]*bybz[4]+grad_u_y[3]*byby[4]+grad_u_x[3]*bxby[4])+0.5*(bybz[0]*grad_u_z[3]+byby[0]*grad_u_y[3]+bxby[0]*grad_u_x[3]+grad_u_z[0]*bybz[3]+grad_u_y[0]*byby[3]+grad_u_x[0]*bxby[3]+bybz[1]*grad_u_z[2]+byby[1]*grad_u_y[2]+bxby[1]*grad_u_x[2]+grad_u_z[1]*bybz[2]+grad_u_y[1]*byby[2]+grad_u_x[1]*bxby[2]))*dx1; 
  bb_grad_u[4] += (0.4472135954999579*(bybz[7]*grad_u_z[7]+byby[7]*grad_u_y[7]+bxby[7]*grad_u_x[7])+(0.31943828249997*bybz[6]+0.5000000000000001*bybz[2])*grad_u_z[6]+(0.31943828249997*byby[6]+0.5000000000000001*byby[2])*grad_u_y[6]+0.31943828249997*bxby[6]*grad_u_x[6]+0.5000000000000001*(bxby[2]*grad_u_x[6]+grad_u_z[2]*bybz[6]+grad_u_y[2]*byby[6]+grad_u_x[2]*bxby[6])+(0.31943828249997*bybz[4]+0.5*bybz[0])*grad_u_z[4]+(0.31943828249997*byby[4]+0.5*byby[0])*grad_u_y[4]+0.31943828249997*bxby[4]*grad_u_x[4]+0.5*(bxby[0]*grad_u_x[4]+grad_u_z[0]*bybz[4]+grad_u_y[0]*byby[4]+grad_u_x[0]*bxby[4])+0.4472135954999579*(bybz[3]*grad_u_z[3]+byby[3]*grad_u_y[3]+bxby[3]*grad_u_x[3]+bybz[1]*grad_u_z[1]+byby[1]*grad_u_y[1]+bxby[1]*grad_u_x[1]))*dx1; 
  bb_grad_u[5] += ((0.31943828249997*bybz[7]+0.5000000000000001*bybz[1])*grad_u_z[7]+(0.31943828249997*byby[7]+0.5000000000000001*byby[1])*grad_u_y[7]+0.31943828249997*bxby[7]*grad_u_x[7]+0.5000000000000001*(bxby[1]*grad_u_x[7]+grad_u_z[1]*bybz[7]+grad_u_y[1]*byby[7]+grad_u_x[1]*bxby[7])+0.4472135954999579*(bybz[6]*grad_u_z[6]+byby[6]*grad_u_y[6]+bxby[6]*grad_u_x[6])+(0.31943828249997*bybz[5]+0.5*bybz[0])*grad_u_z[5]+(0.31943828249997*byby[5]+0.5*byby[0])*grad_u_y[5]+0.31943828249997*bxby[5]*grad_u_x[5]+0.5*(bxby[0]*grad_u_x[5]+grad_u_z[0]*bybz[5]+grad_u_y[0]*byby[5]+grad_u_x[0]*bxby[5])+0.4472135954999579*(bybz[3]*grad_u_z[3]+byby[3]*grad_u_y[3]+bxby[3]*grad_u_x[3]+bybz[2]*grad_u_z[2]+byby[2]*grad_u_y[2]+bxby[2]*grad_u_x[2]))*dx1; 
  bb_grad_u[6] += (0.4*(bybz[3]*grad_u_z[7]+byby[3]*grad_u_y[7]+bxby[3]*grad_u_x[7]+grad_u_z[3]*bybz[7]+grad_u_y[3]*byby[7]+grad_u_x[3]*bxby[7])+(0.4472135954999579*bybz[5]+0.31943828249997*bybz[4]+0.5*bybz[0])*grad_u_z[6]+(0.4472135954999579*byby[5]+0.31943828249997*byby[4]+0.5*byby[0])*grad_u_y[6]+(0.4472135954999579*bxby[5]+0.31943828249997*bxby[4]+0.5*bxby[0])*grad_u_x[6]+(0.4472135954999579*grad_u_z[5]+0.31943828249997*grad_u_z[4]+0.5*grad_u_z[0])*bybz[6]+(0.4472135954999579*grad_u_y[5]+0.31943828249997*grad_u_y[4]+0.5*grad_u_y[0])*byby[6]+(0.4472135954999579*grad_u_x[5]+0.31943828249997*grad_u_x[4]+0.5*grad_u_x[0])*bxby[6]+0.5000000000000001*(bybz[2]*grad_u_z[4]+byby[2]*grad_u_y[4]+bxby[2]*grad_u_x[4]+grad_u_z[2]*bybz[4]+grad_u_y[2]*byby[4]+grad_u_x[2]*bxby[4])+0.447213595499958*(bybz[1]*grad_u_z[3]+byby[1]*grad_u_y[3]+bxby[1]*grad_u_x[3]+grad_u_z[1]*bybz[3]+grad_u_y[1]*byby[3]+grad_u_x[1]*bxby[3]))*dx1; 
  bb_grad_u[7] += ((0.31943828249997*bybz[5]+0.4472135954999579*bybz[4]+0.5*bybz[0])*grad_u_z[7]+(0.31943828249997*byby[5]+0.4472135954999579*byby[4]+0.5*byby[0])*grad_u_y[7]+(0.31943828249997*bxby[5]+0.4472135954999579*bxby[4]+0.5*bxby[0])*grad_u_x[7]+(0.31943828249997*grad_u_z[5]+0.4472135954999579*grad_u_z[4]+0.5*grad_u_z[0])*bybz[7]+(0.31943828249997*grad_u_y[5]+0.4472135954999579*grad_u_y[4]+0.5*grad_u_y[0])*byby[7]+(0.31943828249997*grad_u_x[5]+0.4472135954999579*grad_u_x[4]+0.5*grad_u_x[0])*bxby[7]+0.4*(bybz[3]*grad_u_z[6]+byby[3]*grad_u_y[6]+bxby[3]*grad_u_x[6]+grad_u_z[3]*bybz[6]+grad_u_y[3]*byby[6]+grad_u_x[3]*bxby[6])+0.5000000000000001*(bybz[1]*grad_u_z[5]+byby[1]*grad_u_y[5]+bxby[1]*grad_u_x[5]+grad_u_z[1]*bybz[5]+grad_u_y[1]*byby[5]+grad_u_x[1]*bxby[5])+0.447213595499958*(bybz[2]*grad_u_z[3]+byby[2]*grad_u_y[3]+bxby[2]*grad_u_x[3]+grad_u_z[2]*bybz[3]+grad_u_y[2]*byby[3]+grad_u_x[2]*bxby[3]))*dx1; 

  div_p_x[0] += (0.2445699350390395*Pxy_r[5]-0.2445699350390395*Pxy_l[5]-0.3518228202874282*(Pxy_r[2]+Pxy_l[2])+0.7036456405748563*Pxy_c[2]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_x[1] += (0.2445699350390395*Pxy_r[7]-0.2445699350390395*Pxy_l[7]-0.3518228202874282*(Pxy_r[3]+Pxy_l[3])+0.7036456405748563*Pxy_c[3]+0.25*Pxy_r[1]-0.25*Pxy_l[1])*dx1; 
  div_p_x[2] += (0.4236075534914363*(Pxy_r[5]+Pxy_l[5])+0.8472151069828725*Pxy_c[5]-0.609375*Pxy_r[2]+0.609375*Pxy_l[2]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_x[3] += (0.4236075534914363*(Pxy_r[7]+Pxy_l[7])+0.8472151069828725*Pxy_c[7]-0.609375*Pxy_r[3]+0.609375*Pxy_l[3]+0.4330127018922193*(Pxy_r[1]+Pxy_l[1])-0.8660254037844386*Pxy_c[1])*dx1; 
  div_p_x[4] += ((-0.3518228202874282*(Pxy_r[6]+Pxy_l[6]))+0.7036456405748563*Pxy_c[6]+0.25*Pxy_r[4]-0.25*Pxy_l[4])*dx1; 
  div_p_x[5] += (0.546875*Pxy_r[5]-0.546875*Pxy_l[5]-0.7866997421983816*(Pxy_r[2]+Pxy_l[2])-2.299583861810654*Pxy_c[2]+0.5590169943749475*Pxy_r[0]-0.5590169943749475*Pxy_l[0])*dx1; 
  div_p_x[6] += ((-0.609375*Pxy_r[6])+0.609375*Pxy_l[6]+0.4330127018922194*(Pxy_r[4]+Pxy_l[4])-0.8660254037844387*Pxy_c[4])*dx1; 
  div_p_x[7] += (0.546875*Pxy_r[7]-0.546875*Pxy_l[7]-0.7866997421983816*(Pxy_r[3]+Pxy_l[3])-2.299583861810654*Pxy_c[3]+0.5590169943749476*Pxy_r[1]-0.5590169943749476*Pxy_l[1])*dx1; 

  div_p_y[0] += (0.2445699350390395*Pyy_r[5]-0.2445699350390395*Pyy_l[5]-0.3518228202874282*(Pyy_r[2]+Pyy_l[2])+0.7036456405748563*Pyy_c[2]+0.25*Pyy_r[0]-0.25*Pyy_l[0])*dx1; 
  div_p_y[1] += (0.2445699350390395*Pyy_r[7]-0.2445699350390395*Pyy_l[7]-0.3518228202874282*(Pyy_r[3]+Pyy_l[3])+0.7036456405748563*Pyy_c[3]+0.25*Pyy_r[1]-0.25*Pyy_l[1])*dx1; 
  div_p_y[2] += (0.4236075534914363*(Pyy_r[5]+Pyy_l[5])+0.8472151069828725*Pyy_c[5]-0.609375*Pyy_r[2]+0.609375*Pyy_l[2]+0.4330127018922193*(Pyy_r[0]+Pyy_l[0])-0.8660254037844386*Pyy_c[0])*dx1; 
  div_p_y[3] += (0.4236075534914363*(Pyy_r[7]+Pyy_l[7])+0.8472151069828725*Pyy_c[7]-0.609375*Pyy_r[3]+0.609375*Pyy_l[3]+0.4330127018922193*(Pyy_r[1]+Pyy_l[1])-0.8660254037844386*Pyy_c[1])*dx1; 
  div_p_y[4] += ((-0.3518228202874282*(Pyy_r[6]+Pyy_l[6]))+0.7036456405748563*Pyy_c[6]+0.25*Pyy_r[4]-0.25*Pyy_l[4])*dx1; 
  div_p_y[5] += (0.546875*Pyy_r[5]-0.546875*Pyy_l[5]-0.7866997421983816*(Pyy_r[2]+Pyy_l[2])-2.299583861810654*Pyy_c[2]+0.5590169943749475*Pyy_r[0]-0.5590169943749475*Pyy_l[0])*dx1; 
  div_p_y[6] += ((-0.609375*Pyy_r[6])+0.609375*Pyy_l[6]+0.4330127018922194*(Pyy_r[4]+Pyy_l[4])-0.8660254037844387*Pyy_c[4])*dx1; 
  div_p_y[7] += (0.546875*Pyy_r[7]-0.546875*Pyy_l[7]-0.7866997421983816*(Pyy_r[3]+Pyy_l[3])-2.299583861810654*Pyy_c[3]+0.5590169943749476*Pyy_r[1]-0.5590169943749476*Pyy_l[1])*dx1; 

  div_p_z[0] += (0.2445699350390395*Pyz_r[5]-0.2445699350390395*Pyz_l[5]-0.3518228202874282*(Pyz_r[2]+Pyz_l[2])+0.7036456405748563*Pyz_c[2]+0.25*Pyz_r[0]-0.25*Pyz_l[0])*dx1; 
  div_p_z[1] += (0.2445699350390395*Pyz_r[7]-0.2445699350390395*Pyz_l[7]-0.3518228202874282*(Pyz_r[3]+Pyz_l[3])+0.7036456405748563*Pyz_c[3]+0.25*Pyz_r[1]-0.25*Pyz_l[1])*dx1; 
  div_p_z[2] += (0.4236075534914363*(Pyz_r[5]+Pyz_l[5])+0.8472151069828725*Pyz_c[5]-0.609375*Pyz_r[2]+0.609375*Pyz_l[2]+0.4330127018922193*(Pyz_r[0]+Pyz_l[0])-0.8660254037844386*Pyz_c[0])*dx1; 
  div_p_z[3] += (0.4236075534914363*(Pyz_r[7]+Pyz_l[7])+0.8472151069828725*Pyz_c[7]-0.609375*Pyz_r[3]+0.609375*Pyz_l[3]+0.4330127018922193*(Pyz_r[1]+Pyz_l[1])-0.8660254037844386*Pyz_c[1])*dx1; 
  div_p_z[4] += ((-0.3518228202874282*(Pyz_r[6]+Pyz_l[6]))+0.7036456405748563*Pyz_c[6]+0.25*Pyz_r[4]-0.25*Pyz_l[4])*dx1; 
  div_p_z[5] += (0.546875*Pyz_r[5]-0.546875*Pyz_l[5]-0.7866997421983816*(Pyz_r[2]+Pyz_l[2])-2.299583861810654*Pyz_c[2]+0.5590169943749475*Pyz_r[0]-0.5590169943749475*Pyz_l[0])*dx1; 
  div_p_z[6] += ((-0.609375*Pyz_r[6])+0.609375*Pyz_l[6]+0.4330127018922194*(Pyz_r[4]+Pyz_l[4])-0.8660254037844387*Pyz_c[4])*dx1; 
  div_p_z[7] += (0.546875*Pyz_r[7]-0.546875*Pyz_l[7]-0.7866997421983816*(Pyz_r[3]+Pyz_l[3])-2.299583861810654*Pyz_c[3]+0.5590169943749476*Pyz_r[1]-0.5590169943749476*Pyz_l[1])*dx1; 

} 
