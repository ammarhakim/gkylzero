#include <gkyl_diffusion_kernels.h>

GKYL_CU_DH void
diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
  const double* Dl, const double* Dc, const double* Dr,
  const double* ql, const double* qc, const double* qr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dl: Diffusion coefficient in the left cell
  // Dc: Diffusion coefficient in the center cell
  // Dr: Diffusion coefficient in the right cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (-0.06629126073623882*Dr[1]*qr[1]-0.110485434560398*Dc[1]*qr[1]-0.140335349846953*Dr[0]*qr[1]-0.2423974224629187*Dc[0]*qr[1]-0.06629126073623882*Dl[1]*ql[1]-0.110485434560398*Dc[1]*ql[1]+0.140335349846953*Dl[0]*ql[1]+0.2423974224629187*Dc[0]*ql[1]+0.110485434560398*Dr[1]*qc[1]+0.110485434560398*Dl[1]*qc[1]+0.1325825214724776*Dc[1]*qc[1]-0.2423974224629187*Dr[0]*qc[1]+0.2423974224629187*Dl[0]*qc[1]-0.1148198316929614*qr[0]*Dr[1]+0.1148198316929614*qc[0]*Dr[1]+0.1148198316929614*ql[0]*Dl[1]-0.1148198316929614*qc[0]*Dl[1]+0.1148198316929614*qr[0]*Dc[1]-0.1148198316929614*ql[0]*Dc[1]+0.1988737822087165*Dr[0]*qr[0]+0.1988737822087165*Dc[0]*qr[0]+0.1988737822087165*Dl[0]*ql[0]+0.1988737822087165*Dc[0]*ql[0]-0.1988737822087165*Dr[0]*qc[0]-0.1988737822087165*Dl[0]*qc[0]-0.397747564417433*Dc[0]*qc[0])*J;
  out[1] += (-0.1148198316929615*Dr[1]*qr[1]-0.1913663861549358*Dc[1]*qr[1]-0.2430679560328759*Dr[0]*qr[1]-0.4198446513295127*Dc[0]*qr[1]+0.1148198316929615*Dl[1]*ql[1]+0.1913663861549358*Dc[1]*ql[1]-0.2430679560328759*Dl[0]*ql[1]-0.4198446513295127*Dc[0]*ql[1]+0.1913663861549358*Dr[1]*qc[1]-0.1913663861549358*Dl[1]*qc[1]-0.4198446513295127*Dr[0]*qc[1]-0.4198446513295127*Dl[0]*qc[1]-0.4861359120657519*Dc[0]*qc[1]-0.1988737822087164*qr[0]*Dr[1]+0.1988737822087164*qc[0]*Dr[1]-0.1988737822087164*ql[0]*Dl[1]+0.1988737822087164*qc[0]*Dl[1]+0.1988737822087164*qr[0]*Dc[1]+0.1988737822087164*ql[0]*Dc[1]-0.3977475644174329*qc[0]*Dc[1]+0.3444594950788844*Dr[0]*qr[0]+0.3444594950788844*Dc[0]*qr[0]-0.3444594950788844*Dl[0]*ql[0]-0.3444594950788844*Dc[0]*ql[0]-0.3444594950788844*Dr[0]*qc[0]+0.3444594950788844*Dl[0]*qc[0])*J;
}
