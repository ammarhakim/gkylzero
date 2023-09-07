#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[0]*qr[2]+107.3312629199899*coeff[0]*ql[2]-214.6625258399798*coeff[0]*qc[2]-129.9038105676658*coeff[0]*qr[1]+129.9038105676658*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*qr[2]-290.4737509655563*coeff[0]*ql[2]-405.0*coeff[0]*qr[1]-405.0*coeff[0]*ql[1]-990.0*coeff[0]*qc[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.03125*(21.0*coeff[0]*qr[2]+21.0*coeff[0]*ql[2]-1302.0*coeff[0]*qc[2]-151.0463505020892*coeff[0]*qr[1]+151.0463505020892*coeff[0]*ql[1]+134.1640786499874*coeff[0]*qr[0]+134.1640786499874*coeff[0]*ql[0]-268.3281572999748*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.02209708691207959*((795.0*coeff[2]+453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*qr[2]+(795.0*coeff[2]-453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*ql[2]+((-330.0*coeff[2])-429.3250516799596*coeff[0])*qc[2]+((-755.2317525104463*qr[1])+755.2317525104463*ql[1]+335.4101966249685*qr[0]+335.4101966249685*ql[0]-670.8203932499371*qc[0])*coeff[2]+((-495.0*coeff[1])-259.8076211353315*coeff[0])*qr[1]+(259.8076211353315*coeff[0]-495.0*coeff[1])*ql[1]-810.0*coeff[1]*qc[1]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[1]+150.0*coeff[0]*qr[0]+150.0*coeff[0]*ql[0]-300.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.02209708691207959*((1195.115057222525*coeff[2]+643.9875775199395*coeff[1]+290.4737509655563*coeff[0])*qr[2]+((-1195.115057222525*coeff[2])+643.9875775199395*coeff[1]-290.4737509655563*coeff[0])*ql[2]-1287.975155039879*coeff[1]*qc[2]+((-1207.476707849887*qr[1])-1207.476707849887*ql[1]-1609.968943799849*qc[1]+580.9475019311126*qr[0]-580.9475019311126*ql[0])*coeff[2]+((-779.4228634059946*coeff[1])-405.0*coeff[0])*qr[1]+(779.4228634059946*coeff[1]-405.0*coeff[0])*ql[1]-990.0*coeff[0]*qc[1]+(450.0*qr[0]+450.0*ql[0]-900.0*qc[0])*coeff[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0110485434560398*((2206.999093792293*coeff[2]+618.3421383020891*coeff[1]+42.0*coeff[0])*qr[2]+(2206.999093792293*coeff[2]-618.3421383020891*coeff[1]+42.0*coeff[0])*ql[2]+((-1596.55253593485*coeff[2])-2604.0*coeff[0])*qc[2]+((-2468.17240078565*qr[1])+2468.17240078565*ql[1]+1320.0*qr[0]+1320.0*ql[0]-2640.0*qc[0])*coeff[2]+((-986.1059780774073*coeff[1])-302.0927010041785*coeff[0])*qr[1]+(302.0927010041785*coeff[0]-986.1059780774073*coeff[1])*ql[1]-2535.701086484762*coeff[1]*qc[1]+(650.6612021628459*qr[0]-650.6612021628459*ql[0])*coeff[1]+268.3281572999748*coeff[0]*qr[0]+268.3281572999748*coeff[0]*ql[0]-536.6563145999496*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

