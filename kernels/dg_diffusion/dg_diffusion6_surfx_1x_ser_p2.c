#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion6_surfx_1x_ser_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[0]; 
  const double J = pow(dx1, 6.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*(35.21807064562169*q0r[2]+35.21807064562169*q0l[2]-70.43614129124337*q0c[2]-34.09975027401226*q0r[1]+34.09975027401226*q0l[1]+19.6875*q0r[0]+19.6875*q0l[0]-39.375*q0c[0]); 
  out0[1] += J*D*(51.46831774920947*q0r[2]-51.46831774920947*q0l[2]-56.6015625*q0r[1]-56.6015625*q0l[1]-123.046875*q0c[1]+34.09975027401226*q0r[0]-34.09975027401226*q0l[0]); 
  out0[2] += J*D*((-3.1640625*q0r[2])-3.1640625*q0l[2]-141.328125*q0c[2]-12.2543613688594*q0r[1]+12.2543613688594*q0l[1]+12.57788237343632*q0r[0]+12.57788237343632*q0l[0]-25.15576474687264*q0c[0]); 

  return 0.;

} 
