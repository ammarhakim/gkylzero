#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[2]+563.489130329947*coeff[0]*ql[2]-1126.978260659894*coeff[0]*qc[2]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[2]-6587.944671898812*coeff[0]*ql[2]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0078125*(405.0*coeff[0]*qr[2]+405.0*coeff[0]*ql[2]+18090.0*coeff[0]*qc[2]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_tensor_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0110485434560398*((6885.0*coeff[2]+5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*qr[2]+(6885.0*coeff[2]-5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*ql[2]+(5130.0*coeff[2]-4507.913042639576*coeff[0])*qc[2]+((-2614.263758690006*qr[1])+2614.263758690006*ql[1]+804.9844718999244*qr[0]+804.9844718999244*ql[0]-1609.968943799849*qc[0])*coeff[2]+((-4095.0*coeff[1])-2182.384017536785*coeff[0])*qr[1]+(2182.384017536785*coeff[0]-4095.0*coeff[1])*ql[1]-6930.0*coeff[1]*qc[1]+(2182.384017536785*qr[0]-2182.384017536785*ql[0])*coeff[1]+1260.0*coeff[0]*qr[0]+1260.0*coeff[0]*ql[0]-2520.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.005524271728019897*((31098.97224989918*coeff[2]+18212.77367673579*coeff[1]+6587.944671898812*coeff[0])*qr[2]+((-31098.97224989918*coeff[2])+18212.77367673579*coeff[1]-6587.944671898812*coeff[0])*ql[2]-27973.21039852238*coeff[1]*qc[2]+((-20426.48097446058*qr[1])-20426.48097446058*ql[1]-26765.73369067249*qc[1]+9759.918032442689*qr[0]-9759.918032442689*ql[0])*coeff[2]+((-16757.59156322888*coeff[1])-7245.0*coeff[0])*qr[1]+(16757.59156322888*coeff[1]-7245.0*coeff[0])*ql[1]-15750.0*coeff[0]*qc[1]+(9360.0*qr[0]+9360.0*ql[0]-18720.0*qc[0])*coeff[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.005524271728019897*((45984.73795728319*coeff[2]+14731.09211837329*coeff[1]-405.0*coeff[0])*qr[2]+(45984.73795728319*coeff[2]-14731.09211837329*coeff[1]-405.0*coeff[0])*ql[2]+((-49707.79113982034*coeff[2])-18090.0*coeff[0])*qc[2]+((-40140.27746540872*qr[1])+40140.27746540872*ql[1]+21600.0*qr[0]+21600.0*ql[0]-43200.0*qc[0])*coeff[2]+((-16200.31249698598*coeff[1])-1568.558255214003*coeff[0])*qr[1]+(1568.558255214003*coeff[0]-16200.31249698598*coeff[1])*ql[1]-35218.07064562169*coeff[1]*qc[1]+(9759.918032442689*qr[0]-9759.918032442689*ql[0])*coeff[1]+1609.968943799849*coeff[0]*qr[0]+1609.968943799849*coeff[0]*ql[0]-3219.937887599698*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

