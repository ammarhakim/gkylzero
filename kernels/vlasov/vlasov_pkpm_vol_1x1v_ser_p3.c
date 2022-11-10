#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *bvar, const double *rho_inv_b, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // div_p:     divergence of the pressure tensor.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[4]; 
  const double *div_p_z = &div_p[8]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[4]; 
  const double *rho_inv_bz = &rho_inv_b[8]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[12] = {0.0}; 
  double alpha_vdim[12] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 
  alpha_cdim[4] = 1.414213562373095*dx0*(bx[2]*wvpar+ux[2]); 
  alpha_cdim[6] = 0.408248290463863*bx[2]*dvpar*dx0; 
  alpha_cdim[8] = 1.414213562373095*dx0*(bx[3]*wvpar+ux[3]); 
  alpha_cdim[10] = 0.408248290463863*bx[3]*dvpar*dx0; 
  cflFreq_mid += 7.0*fabs(0.25*alpha_cdim[0]-0.2795084971874737*alpha_cdim[4]); 

  alpha_vdim[0] = ((-5.916079783099617*bxbz[2]*uz[3])-2.645751311064591*bxbz[0]*uz[3]-5.916079783099617*bxby[2]*uy[3]-2.645751311064591*bxby[0]*uy[3]-5.916079783099617*bxbx[2]*ux[3]-2.645751311064591*bxbx[0]*ux[3]-3.872983346207417*bxbz[1]*uz[2]-3.872983346207417*bxby[1]*uy[2]-3.872983346207417*bxbx[1]*ux[2]-1.732050807568877*bxbz[0]*uz[1]-1.732050807568877*bxby[0]*uy[1]-1.732050807568877*bxbx[0]*ux[1])*dv1par*dx0*wvpar+(div_p_z[3]*rho_inv_bz[3]+div_p_y[3]*rho_inv_by[3]+div_p_x[3]*rho_inv_bx[3]+div_p_z[2]*rho_inv_bz[2]+div_p_y[2]*rho_inv_by[2]+div_p_x[2]*rho_inv_bx[2]+div_p_z[1]*rho_inv_bz[1]+div_p_y[1]*rho_inv_by[1]+div_p_x[1]*rho_inv_bx[1]+div_p_z[0]*rho_inv_bz[0]+div_p_y[0]*rho_inv_by[0]+div_p_x[0]*rho_inv_bx[0])*dv1par; 
  alpha_vdim[1] = ((-5.196152422706631*bxbz[3]*uz[3])-7.937253933193772*bxbz[1]*uz[3]-5.196152422706631*bxby[3]*uy[3]-7.937253933193772*bxby[1]*uy[3]-5.196152422706631*bxbx[3]*ux[3]-7.937253933193772*bxbx[1]*ux[3]-3.464101615137754*bxbz[2]*uz[2]-3.872983346207417*bxbz[0]*uz[2]-3.464101615137754*bxby[2]*uy[2]-3.872983346207417*bxby[0]*uy[2]-3.464101615137754*bxbx[2]*ux[2]-3.872983346207417*bxbx[0]*ux[2]-1.732050807568877*bxbz[1]*uz[1]-1.732050807568877*bxby[1]*uy[1]-1.732050807568877*bxbx[1]*ux[1])*dv1par*dx0*wvpar+(0.8783100656536796*div_p_z[2]*rho_inv_bz[3]+0.8783100656536796*div_p_y[2]*rho_inv_by[3]+0.8783100656536796*div_p_x[2]*rho_inv_bx[3]+0.8783100656536796*rho_inv_bz[2]*div_p_z[3]+0.8783100656536796*rho_inv_by[2]*div_p_y[3]+0.8783100656536796*rho_inv_bx[2]*div_p_x[3]+0.8944271909999159*div_p_z[1]*rho_inv_bz[2]+0.8944271909999159*div_p_y[1]*rho_inv_by[2]+0.8944271909999159*div_p_x[1]*rho_inv_bx[2]+0.8944271909999159*rho_inv_bz[1]*div_p_z[2]+0.8944271909999159*rho_inv_by[1]*div_p_y[2]+0.8944271909999159*rho_inv_bx[1]*div_p_x[2]+div_p_z[0]*rho_inv_bz[1]+div_p_y[0]*rho_inv_by[1]+div_p_x[0]*rho_inv_bx[1]+rho_inv_bz[0]*div_p_z[1]+rho_inv_by[0]*div_p_y[1]+rho_inv_bx[0]*div_p_x[1])*dv1par; 
  alpha_vdim[2] = ((-1.707825127659933*bxbz[2]*uz[3])-0.7637626158259735*bxbz[0]*uz[3]-1.707825127659933*bxby[2]*uy[3]-0.7637626158259735*bxby[0]*uy[3]-1.707825127659933*bxbx[2]*ux[3]-0.7637626158259735*bxbx[0]*ux[3]-1.118033988749895*bxbz[1]*uz[2]-1.118033988749895*bxby[1]*uy[2]-1.118033988749895*bxbx[1]*ux[2]-0.5*bxbz[0]*uz[1]-0.5*bxby[0]*uy[1]-0.5*bxbx[0]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[3] = ((-1.5*bxbz[3]*uz[3])-2.29128784747792*bxbz[1]*uz[3]-1.5*bxby[3]*uy[3]-2.29128784747792*bxby[1]*uy[3]-1.5*bxbx[3]*ux[3]-2.29128784747792*bxbx[1]*ux[3]-1.0*bxbz[2]*uz[2]-1.118033988749895*bxbz[0]*uz[2]-1.0*bxby[2]*uy[2]-1.118033988749895*bxby[0]*uy[2]-1.0*bxbx[2]*ux[2]-1.118033988749895*bxbx[0]*ux[2]-0.5*bxbz[1]*uz[1]-0.5*bxby[1]*uy[1]-0.5*bxbx[1]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[4] = ((-6.425396041156862*bxbz[2]*uz[3])-5.916079783099617*bxbz[0]*uz[3]-6.425396041156862*bxby[2]*uy[3]-5.916079783099617*bxby[0]*uy[3]-6.425396041156862*bxbx[2]*ux[3]-5.916079783099617*bxbx[0]*ux[3]-3.401680257083045*uz[2]*bxbz[3]-3.401680257083045*uy[2]*bxby[3]-3.401680257083045*ux[2]*bxbx[3]-3.464101615137754*bxbz[1]*uz[2]-3.464101615137754*bxby[1]*uy[2]-3.464101615137754*bxbx[1]*ux[2]-1.732050807568877*uz[1]*bxbz[2]-1.732050807568877*uy[1]*bxby[2]-1.732050807568877*ux[1]*bxbx[2])*dv1par*dx0*wvpar+(0.5962847939999438*div_p_z[3]*rho_inv_bz[3]+0.8783100656536796*div_p_z[1]*rho_inv_bz[3]+0.5962847939999438*div_p_y[3]*rho_inv_by[3]+0.8783100656536796*div_p_y[1]*rho_inv_by[3]+0.5962847939999438*div_p_x[3]*rho_inv_bx[3]+0.8783100656536796*div_p_x[1]*rho_inv_bx[3]+0.8783100656536796*rho_inv_bz[1]*div_p_z[3]+0.8783100656536796*rho_inv_by[1]*div_p_y[3]+0.8783100656536796*rho_inv_bx[1]*div_p_x[3]+0.6388765649999399*div_p_z[2]*rho_inv_bz[2]+div_p_z[0]*rho_inv_bz[2]+0.6388765649999399*div_p_y[2]*rho_inv_by[2]+div_p_y[0]*rho_inv_by[2]+0.6388765649999399*div_p_x[2]*rho_inv_bx[2]+div_p_x[0]*rho_inv_bx[2]+rho_inv_bz[0]*div_p_z[2]+rho_inv_by[0]*div_p_y[2]+rho_inv_bx[0]*div_p_x[2]+0.8944271909999159*div_p_z[1]*rho_inv_bz[1]+0.8944271909999159*div_p_y[1]*rho_inv_by[1]+0.8944271909999159*div_p_x[1]*rho_inv_bx[1])*dv1par; 
  alpha_vdim[6] = ((-1.854852067005935*bxbz[2]*uz[3])-1.707825127659933*bxbz[0]*uz[3]-1.854852067005935*bxby[2]*uy[3]-1.707825127659933*bxby[0]*uy[3]-1.854852067005935*bxbx[2]*ux[3]-1.707825127659933*bxbx[0]*ux[3]-0.9819805060619657*uz[2]*bxbz[3]-0.9819805060619657*uy[2]*bxby[3]-0.9819805060619657*ux[2]*bxbx[3]-1.0*bxbz[1]*uz[2]-1.0*bxby[1]*uy[2]-1.0*bxbx[1]*ux[2]-0.5000000000000001*uz[1]*bxbz[2]-0.5000000000000001*uy[1]*bxby[2]-0.5000000000000001*ux[1]*bxbx[2])*dv1par*dvpar*dx0; 
  alpha_vdim[8] = ((-6.173419725817379*bxbz[3]*uz[3])-5.196152422706631*bxbz[1]*uz[3]-6.173419725817379*bxby[3]*uy[3]-5.196152422706631*bxby[1]*uy[3]-6.173419725817379*bxbx[3]*ux[3]-5.196152422706631*bxbx[1]*ux[3]-1.732050807568877*uz[1]*bxbz[3]-1.732050807568877*uy[1]*bxby[3]-1.732050807568877*ux[1]*bxbx[3]-3.401680257083045*bxbz[2]*uz[2]-3.401680257083045*bxby[2]*uy[2]-3.401680257083045*bxbx[2]*ux[2])*dv1par*dx0*wvpar+(0.5962847939999438*div_p_z[2]*rho_inv_bz[3]+div_p_z[0]*rho_inv_bz[3]+0.5962847939999438*div_p_y[2]*rho_inv_by[3]+div_p_y[0]*rho_inv_by[3]+0.5962847939999438*div_p_x[2]*rho_inv_bx[3]+div_p_x[0]*rho_inv_bx[3]+0.5962847939999438*rho_inv_bz[2]*div_p_z[3]+rho_inv_bz[0]*div_p_z[3]+0.5962847939999438*rho_inv_by[2]*div_p_y[3]+rho_inv_by[0]*div_p_y[3]+0.5962847939999438*rho_inv_bx[2]*div_p_x[3]+rho_inv_bx[0]*div_p_x[3]+0.8783100656536796*div_p_z[1]*rho_inv_bz[2]+0.8783100656536796*div_p_y[1]*rho_inv_by[2]+0.8783100656536796*div_p_x[1]*rho_inv_bx[2]+0.8783100656536796*rho_inv_bz[1]*div_p_z[2]+0.8783100656536796*rho_inv_by[1]*div_p_y[2]+0.8783100656536796*rho_inv_bx[1]*div_p_x[2])*dv1par; 
  alpha_vdim[10] = ((-1.782112770260604*bxbz[3]*uz[3])-1.5*bxbz[1]*uz[3]-1.782112770260604*bxby[3]*uy[3]-1.5*bxby[1]*uy[3]-1.782112770260604*bxbx[3]*ux[3]-1.5*bxbx[1]*ux[3]-0.5*uz[1]*bxbz[3]-0.5*uy[1]*bxby[3]-0.5*ux[1]*bxbx[3]-0.9819805060619656*bxbz[2]*uz[2]-0.9819805060619656*bxby[2]*uy[2]-0.9819805060619656*bxbx[2]*ux[2])*dv1par*dvpar*dx0; 
  cflFreq_mid += 7.0*fabs(0.25*alpha_vdim[0]-0.2795084971874737*alpha_vdim[4]); 

  out[1] += 0.8660254037844386*(alpha_cdim[10]*f[10]+alpha_cdim[8]*f[8]+alpha_cdim[6]*f[6]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[1]*f[1]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[10]*f[10]+alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.8660254037844386*alpha_cdim[8]*f[10]+0.7606388292556648*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.8660254037844386*f[8]*alpha_cdim[10]+0.7606388292556648*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+0.7745966692414833*alpha_cdim[3]*f[7]+0.8660254037844386*alpha_cdim[4]*f[6]+0.7745966692414833*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+0.8660254037844386*f[4]*alpha_cdim[6]+0.7745966692414833*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*((alpha_vdim[2]+alpha_cdim[1])*f[3]+f[2]*alpha_vdim[3]+f[1]*alpha_cdim[3]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.700840128541522*(alpha_cdim[6]*f[10]+f[6]*alpha_cdim[10]+alpha_cdim[4]*f[8]+f[4]*alpha_cdim[8])+1.732050807568877*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_cdim[1]*f[4]+f[1]*alpha_cdim[4])+1.936491673103709*(alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]+alpha_cdim[0]*f[1]+f[0]*alpha_cdim[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[8]*f[10]+f[8]*alpha_vdim[10])+1.732050807568877*alpha_vdim[3]*f[7]+1.936491673103709*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6])+1.732050807568877*alpha_vdim[2]*f[5]+1.936491673103709*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[6] += (0.5163977794943223*alpha_vdim[10]+1.700840128541522*alpha_cdim[4])*f[10]+0.7606388292556648*(alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10])+1.700840128541522*f[4]*alpha_cdim[10]+(0.5163977794943223*alpha_vdim[8]+1.700840128541522*alpha_cdim[6])*f[8]+0.7606388292556648*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8])+1.700840128541522*f[6]*alpha_cdim[8]+(1.549193338482967*alpha_cdim[6]+1.732050807568877*alpha_cdim[2])*f[7]+(0.5532833351724881*alpha_vdim[6]+0.8660254037844386*alpha_vdim[2]+1.732050807568877*alpha_cdim[1])*f[6]+0.8660254037844386*f[2]*alpha_vdim[6]+1.732050807568877*(f[1]*alpha_cdim[6]+alpha_cdim[3]*f[5])+(0.5532833351724881*alpha_vdim[4]+1.732050807568877*alpha_cdim[3])*f[4]+0.8660254037844386*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+f[3]*(1.732050807568877*alpha_cdim[4]+0.7745966692414833*alpha_vdim[3])+1.936491673103709*(alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_cdim[1]*f[2])+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 0.7606388292556648*alpha_cdim[3]*f[11]+0.7745966692414833*alpha_cdim[10]*f[10]+1.700840128541522*(alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10])+0.7606388292556648*alpha_cdim[2]*f[9]+1.700840128541522*(alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+(1.549193338482967*alpha_vdim[6]+1.732050807568877*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[7]+0.7745966692414833*alpha_cdim[6]*f[6]+1.732050807568877*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+(1.732050807568877*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[5]+1.732050807568877*alpha_vdim[3]*f[4]+f[3]*(1.732050807568877*alpha_vdim[4]+0.7745966692414833*alpha_cdim[3])+1.936491673103709*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3])+0.7745966692414833*alpha_cdim[2]*f[2]+1.936491673103709*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[8] += 3.086709862908689*alpha_cdim[10]*f[10]+2.598076211353316*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10])+3.086709862908689*alpha_cdim[8]*f[8]+2.598076211353316*(alpha_cdim[1]*f[8]+f[1]*alpha_cdim[8])+3.212698020578431*alpha_cdim[6]*f[6]+2.958039891549809*(alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6])+3.212698020578431*alpha_cdim[4]*f[4]+2.958039891549809*(alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4])+3.968626966596886*alpha_cdim[3]*f[3]+1.322875655532295*alpha_cdim[2]*f[2]+3.968626966596886*alpha_cdim[1]*f[1]+1.322875655532295*alpha_cdim[0]*f[0]; 
  out[9] += 2.598076211353316*alpha_vdim[3]*f[11]+3.968626966596886*alpha_vdim[10]*f[10]+2.598076211353316*alpha_vdim[2]*f[9]+1.322875655532295*alpha_vdim[8]*f[8]+2.958039891549809*alpha_vdim[1]*f[7]+3.968626966596886*alpha_vdim[6]*f[6]+2.958039891549809*alpha_vdim[0]*f[5]+1.322875655532295*alpha_vdim[4]*f[4]+3.968626966596886*(alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2])+1.322875655532295*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[10] += (3.086709862908689*alpha_cdim[8]+0.5163977794943223*alpha_vdim[6]+0.8660254037844386*alpha_vdim[2]+2.598076211353316*alpha_cdim[1])*f[10]+(0.5163977794943223*f[6]+0.8660254037844386*f[2])*alpha_vdim[10]+(3.086709862908689*f[8]+2.32379000772445*f[7]+2.598076211353316*f[1])*alpha_cdim[10]+(0.5163977794943223*alpha_vdim[4]+2.598076211353316*alpha_cdim[3]+0.8660254037844386*alpha_vdim[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alpha_vdim[8]+2.598076211353316*f[3]*alpha_cdim[8]+3.54964786985977*alpha_cdim[3]*f[7]+(3.212698020578431*alpha_cdim[4]+0.7606388292556648*alpha_vdim[3]+2.958039891549809*alpha_cdim[0])*f[6]+0.7606388292556648*f[3]*alpha_vdim[6]+(2.645751311064591*f[5]+3.212698020578431*f[4]+2.958039891549809*f[0])*alpha_cdim[6]+alpha_cdim[2]*(1.183215956619923*f[5]+2.958039891549809*f[4])+0.7606388292556648*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+2.958039891549809*f[2]*alpha_cdim[4]+3.968626966596886*(alpha_cdim[1]*f[3]+f[1]*alpha_cdim[3])+1.322875655532295*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[11] += (2.32379000772445*alpha_vdim[6]+2.598076211353316*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[11]+3.485685011586674*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+(2.598076211353316*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[9]+1.161895003862225*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+(2.645751311064591*alpha_vdim[4]+0.7606388292556648*alpha_cdim[3]+2.958039891549809*alpha_vdim[0])*f[7]+3.54964786985977*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+(0.7606388292556648*alpha_cdim[2]+2.958039891549809*alpha_vdim[1])*f[5]+1.183215956619923*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+3.968626966596886*(alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3])+1.322875655532295*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 

  return cflFreq_mid; 
} 
