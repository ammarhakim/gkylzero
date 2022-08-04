#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_2x2v_ser_p3(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[12]; 
  const double *rhou1l = &statevecl[24]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[12]; 
  const double *rhou1r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[12]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[12]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[12]; 
  double *outlrhouy = &out[24]; 
  double *outrrho = &out[36]; 
  double *outrrhoux = &out[48]; 
  double *outrrhouy = &out[60]; 
  double incr[12]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.6614378277661477*rhou0r[8])+0.6614378277661477*rhou0l[8]+0.5590169943749475*rhou0r[4]+0.5590169943749475*rhou0l[4]-0.4330127018922193*rhou0r[1]+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = 1.14564392373896*rhou0r[8]-1.14564392373896*rhou0l[8]-0.9682458365518543*rhou0r[4]-0.9682458365518543*rhou0l[4]+0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = (-0.6614378277661477*rhou0r[10])+0.6614378277661477*rhou0l[10]+0.5590169943749476*rhou0r[6]+0.5590169943749476*rhou0l[6]-0.4330127018922193*rhou0r[3]+0.4330127018922193*rhou0l[3]+0.25*rhou0r[2]+0.25*rhou0l[2]; 
  incr[3] = 1.14564392373896*rhou0r[10]-1.14564392373896*rhou0l[10]-0.9682458365518543*rhou0r[6]-0.9682458365518543*rhou0l[6]+0.75*rhou0r[3]-0.75*rhou0l[3]-0.4330127018922194*rhou0r[2]-0.4330127018922194*rhou0l[2]; 
  incr[4] = (-1.479019945774904*rhou0r[8])+1.479019945774904*rhou0l[8]+1.25*rhou0r[4]+1.25*rhou0l[4]-0.9682458365518543*rhou0r[1]+0.9682458365518543*rhou0l[1]+0.5590169943749475*rhou0r[0]+0.5590169943749475*rhou0l[0]; 
  incr[5] = (-0.4330127018922193*rhou0r[7])+0.4330127018922193*rhou0l[7]+0.25*rhou0r[5]+0.25*rhou0l[5]; 
  incr[6] = (-1.479019945774904*rhou0r[10])+1.479019945774904*rhou0l[10]+1.25*rhou0r[6]+1.25*rhou0l[6]-0.9682458365518543*rhou0r[3]+0.9682458365518543*rhou0l[3]+0.5590169943749476*rhou0r[2]+0.5590169943749476*rhou0l[2]; 
  incr[7] = 0.75*rhou0r[7]-0.75*rhou0l[7]-0.4330127018922194*rhou0r[5]-0.4330127018922194*rhou0l[5]; 
  incr[8] = 1.75*rhou0r[8]-1.75*rhou0l[8]-1.479019945774904*rhou0r[4]-1.479019945774904*rhou0l[4]+1.14564392373896*rhou0r[1]-1.14564392373896*rhou0l[1]-0.6614378277661477*rhou0r[0]-0.6614378277661477*rhou0l[0]; 
  incr[9] = (-0.4330127018922193*rhou0r[11])+0.4330127018922193*rhou0l[11]+0.25*rhou0r[9]+0.25*rhou0l[9]; 
  incr[10] = 1.75*rhou0r[10]-1.75*rhou0l[10]-1.479019945774904*rhou0r[6]-1.479019945774904*rhou0l[6]+1.14564392373896*rhou0r[3]-1.14564392373896*rhou0l[3]-0.6614378277661477*rhou0r[2]-0.6614378277661477*rhou0l[2]; 
  incr[11] = 0.75*rhou0r[11]-0.75*rhou0l[11]-0.4330127018922193*rhou0r[9]-0.4330127018922193*rhou0l[9]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 
  outrrho[8] += incr[8]*dxr1; 
  outrrho[9] += incr[9]*dxr1; 
  outrrho[10] += incr[10]*dxr1; 
  outrrho[11] += incr[11]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += -1.0*incr[6]*dxl1; 
  outlrho[7] += incr[7]*dxl1; 
  outlrho[8] += incr[8]*dxl1; 
  outlrho[9] += -1.0*incr[9]*dxl1; 
  outlrho[10] += incr[10]*dxl1; 
  outlrho[11] += incr[11]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = 3.968626966596886*rhor[8]*vthsq+3.968626966596886*rhol[8]*vthsq-1.677050983124842*rhor[4]*vthsq+1.677050983124842*rhol[4]*vthsq+0.4330127018922193*rhor[1]*vthsq+0.4330127018922193*rhol[1]*vthsq+0.75*rhor[11]*uvar0r[11]-0.4330127018922193*rhor[9]*uvar0r[11]+0.75*rhol[11]*uvar0l[11]+0.4330127018922193*rhol[9]*uvar0l[11]-0.4330127018922193*uvar0r[9]*rhor[11]+0.4330127018922193*uvar0l[9]*rhol[11]+1.75*rhor[10]*uvar0r[10]-1.479019945774903*rhor[6]*uvar0r[10]+1.14564392373896*rhor[3]*uvar0r[10]-0.6614378277661479*rhor[2]*uvar0r[10]+1.75*rhol[10]*uvar0l[10]+1.479019945774903*rhol[6]*uvar0l[10]+1.14564392373896*rhol[3]*uvar0l[10]+0.6614378277661479*rhol[2]*uvar0l[10]-1.479019945774903*uvar0r[6]*rhor[10]+1.14564392373896*uvar0r[3]*rhor[10]-0.6614378277661479*uvar0r[2]*rhor[10]+1.479019945774903*uvar0l[6]*rhol[10]+1.14564392373896*uvar0l[3]*rhol[10]+0.6614378277661479*uvar0l[2]*rhol[10]+0.25*rhor[9]*uvar0r[9]+0.25*rhol[9]*uvar0l[9]+1.75*rhor[8]*uvar0r[8]-1.479019945774905*rhor[4]*uvar0r[8]+1.14564392373896*rhor[1]*uvar0r[8]-0.6614378277661477*rhor[0]*uvar0r[8]+1.75*rhol[8]*uvar0l[8]+1.479019945774905*rhol[4]*uvar0l[8]+1.14564392373896*rhol[1]*uvar0l[8]+0.6614378277661477*rhol[0]*uvar0l[8]-1.479019945774905*uvar0r[4]*rhor[8]+1.14564392373896*uvar0r[1]*rhor[8]-0.6614378277661477*uvar0r[0]*rhor[8]+1.479019945774905*uvar0l[4]*rhol[8]+1.14564392373896*uvar0l[1]*rhol[8]+0.6614378277661477*uvar0l[0]*rhol[8]+0.75*rhor[7]*uvar0r[7]-0.4330127018922194*rhor[5]*uvar0r[7]+0.75*rhol[7]*uvar0l[7]+0.4330127018922194*rhol[5]*uvar0l[7]-0.4330127018922194*uvar0r[5]*rhor[7]+0.4330127018922194*uvar0l[5]*rhol[7]+1.25*rhor[6]*uvar0r[6]-0.9682458365518543*rhor[3]*uvar0r[6]+0.5590169943749472*rhor[2]*uvar0r[6]+1.25*rhol[6]*uvar0l[6]+0.9682458365518543*rhol[3]*uvar0l[6]+0.5590169943749472*rhol[2]*uvar0l[6]-0.9682458365518543*uvar0r[3]*rhor[6]+0.5590169943749472*uvar0r[2]*rhor[6]+0.9682458365518543*uvar0l[3]*rhol[6]+0.5590169943749472*uvar0l[2]*rhol[6]+0.25*rhor[5]*uvar0r[5]+0.25*rhol[5]*uvar0l[5]+1.25*rhor[4]*uvar0r[4]-0.9682458365518541*rhor[1]*uvar0r[4]+0.5590169943749475*rhor[0]*uvar0r[4]+1.25*rhol[4]*uvar0l[4]+0.9682458365518541*rhol[1]*uvar0l[4]+0.5590169943749475*rhol[0]*uvar0l[4]-0.9682458365518541*uvar0r[1]*rhor[4]+0.5590169943749475*uvar0r[0]*rhor[4]+0.9682458365518541*uvar0l[1]*rhol[4]+0.5590169943749475*uvar0l[0]*rhol[4]+0.75*rhor[3]*uvar0r[3]-0.4330127018922193*rhor[2]*uvar0r[3]+0.75*rhol[3]*uvar0l[3]+0.4330127018922193*rhol[2]*uvar0l[3]-0.4330127018922193*uvar0r[2]*rhor[3]+0.4330127018922193*uvar0l[2]*rhol[3]+0.25*rhor[2]*uvar0r[2]+0.25*rhol[2]*uvar0l[2]+0.75*rhor[1]*uvar0r[1]-0.4330127018922193*rhor[0]*uvar0r[1]+0.75*rhol[1]*uvar0l[1]+0.4330127018922193*rhol[0]*uvar0l[1]-0.4330127018922193*uvar0r[0]*rhor[1]+0.4330127018922193*uvar0l[0]*rhol[1]+0.25*rhor[0]*uvar0r[0]+0.25*rhol[0]*uvar0l[0]; 
  incr[1] = 11.90588089979066*rhor[8]*vthsq-11.90588089979066*rhol[8]*vthsq-5.031152949374527*rhor[4]*vthsq-5.031152949374527*rhol[4]*vthsq+1.299038105676658*rhor[1]*vthsq-1.299038105676658*rhol[1]*vthsq-1.299038105676658*rhor[11]*uvar0r[11]+0.7499999999999999*rhor[9]*uvar0r[11]-1.299038105676658*rhol[11]*uvar0l[11]-0.7499999999999999*rhol[9]*uvar0l[11]+0.7499999999999999*uvar0r[9]*rhor[11]-0.7499999999999999*uvar0l[9]*rhol[11]-3.031088913245535*rhor[10]*uvar0r[10]+2.561737691489897*rhor[6]*uvar0r[10]-1.984313483298443*rhor[3]*uvar0r[10]+1.14564392373896*rhor[2]*uvar0r[10]-3.031088913245535*rhol[10]*uvar0l[10]-2.561737691489897*rhol[6]*uvar0l[10]-1.984313483298443*rhol[3]*uvar0l[10]-1.14564392373896*rhol[2]*uvar0l[10]+2.561737691489897*uvar0r[6]*rhor[10]-1.984313483298443*uvar0r[3]*rhor[10]+1.14564392373896*uvar0r[2]*rhor[10]-2.561737691489897*uvar0l[6]*rhol[10]-1.984313483298443*uvar0l[3]*rhol[10]-1.14564392373896*uvar0l[2]*rhol[10]-0.4330127018922193*rhor[9]*uvar0r[9]-0.4330127018922193*rhol[9]*uvar0l[9]-3.031088913245535*rhor[8]*uvar0r[8]+2.561737691489901*rhor[4]*uvar0r[8]-1.984313483298443*rhor[1]*uvar0r[8]+1.14564392373896*rhor[0]*uvar0r[8]-3.031088913245535*rhol[8]*uvar0l[8]-2.561737691489901*rhol[4]*uvar0l[8]-1.984313483298443*rhol[1]*uvar0l[8]-1.14564392373896*rhol[0]*uvar0l[8]+2.561737691489901*uvar0r[4]*rhor[8]-1.984313483298443*uvar0r[1]*rhor[8]+1.14564392373896*uvar0r[0]*rhor[8]-2.561737691489901*uvar0l[4]*rhol[8]-1.984313483298443*uvar0l[1]*rhol[8]-1.14564392373896*uvar0l[0]*rhol[8]-1.299038105676658*rhor[7]*uvar0r[7]+0.75*rhor[5]*uvar0r[7]-1.299038105676658*rhol[7]*uvar0l[7]-0.75*rhol[5]*uvar0l[7]+0.75*uvar0r[5]*rhor[7]-0.75*uvar0l[5]*rhol[7]-2.165063509461096*rhor[6]*uvar0r[6]+1.677050983124842*rhor[3]*uvar0r[6]-0.9682458365518543*rhor[2]*uvar0r[6]-2.165063509461096*rhol[6]*uvar0l[6]-1.677050983124842*rhol[3]*uvar0l[6]-0.9682458365518543*rhol[2]*uvar0l[6]+1.677050983124842*uvar0r[3]*rhor[6]-0.9682458365518543*uvar0r[2]*rhor[6]-1.677050983124842*uvar0l[3]*rhol[6]-0.9682458365518543*uvar0l[2]*rhol[6]-0.4330127018922193*rhor[5]*uvar0r[5]-0.4330127018922193*rhol[5]*uvar0l[5]-2.165063509461096*rhor[4]*uvar0r[4]+1.677050983124842*rhor[1]*uvar0r[4]-0.9682458365518543*rhor[0]*uvar0r[4]-2.165063509461096*rhol[4]*uvar0l[4]-1.677050983124842*rhol[1]*uvar0l[4]-0.9682458365518543*rhol[0]*uvar0l[4]+1.677050983124842*uvar0r[1]*rhor[4]-0.9682458365518543*uvar0r[0]*rhor[4]-1.677050983124842*uvar0l[1]*rhol[4]-0.9682458365518543*uvar0l[0]*rhol[4]-1.299038105676658*rhor[3]*uvar0r[3]+0.75*rhor[2]*uvar0r[3]-1.299038105676658*rhol[3]*uvar0l[3]-0.75*rhol[2]*uvar0l[3]+0.75*uvar0r[2]*rhor[3]-0.75*uvar0l[2]*rhol[3]-0.4330127018922193*rhor[2]*uvar0r[2]-0.4330127018922193*rhol[2]*uvar0l[2]-1.299038105676658*rhor[1]*uvar0r[1]+0.75*rhor[0]*uvar0r[1]-1.299038105676658*rhol[1]*uvar0l[1]-0.75*rhol[0]*uvar0l[1]+0.75*uvar0r[0]*rhor[1]-0.75*uvar0l[0]*rhol[1]-0.4330127018922193*rhor[0]*uvar0r[0]-0.4330127018922193*rhol[0]*uvar0l[0]; 
  incr[4] = 19.84313483298443*rhor[8]*vthsq+19.84313483298443*rhol[8]*vthsq-8.385254915624213*rhor[4]*vthsq+8.385254915624213*rhol[4]*vthsq+2.165063509461096*rhor[1]*vthsq+2.165063509461096*rhol[1]*vthsq+1.677050983124842*rhor[11]*uvar0r[11]-0.9682458365518541*rhor[9]*uvar0r[11]+1.677050983124842*rhol[11]*uvar0l[11]+0.9682458365518541*rhol[9]*uvar0l[11]-0.9682458365518541*uvar0r[9]*rhor[11]+0.9682458365518541*uvar0l[9]*rhol[11]+3.913118960624632*rhor[10]*uvar0r[10]-3.307189138830736*rhor[6]*uvar0r[10]+2.561737691489899*rhor[3]*uvar0r[10]-1.479019945774905*rhor[2]*uvar0r[10]+3.913118960624632*rhol[10]*uvar0l[10]+3.307189138830736*rhol[6]*uvar0l[10]+2.561737691489899*rhol[3]*uvar0l[10]+1.479019945774905*rhol[2]*uvar0l[10]-3.307189138830736*uvar0r[6]*rhor[10]+2.561737691489899*uvar0r[3]*rhor[10]-1.479019945774905*uvar0r[2]*rhor[10]+3.307189138830736*uvar0l[6]*rhol[10]+2.561737691489899*uvar0l[3]*rhol[10]+1.479019945774905*uvar0l[2]*rhol[10]+0.5590169943749475*rhor[9]*uvar0r[9]+0.5590169943749475*rhol[9]*uvar0l[9]+3.913118960624632*rhor[8]*uvar0r[8]-3.307189138830738*rhor[4]*uvar0r[8]+2.561737691489901*rhor[1]*uvar0r[8]-1.479019945774904*rhor[0]*uvar0r[8]+3.913118960624632*rhol[8]*uvar0l[8]+3.307189138830738*rhol[4]*uvar0l[8]+2.561737691489901*rhol[1]*uvar0l[8]+1.479019945774904*rhol[0]*uvar0l[8]-3.307189138830738*uvar0r[4]*rhor[8]+2.561737691489901*uvar0r[1]*rhor[8]-1.479019945774904*uvar0r[0]*rhor[8]+3.307189138830738*uvar0l[4]*rhol[8]+2.561737691489901*uvar0l[1]*rhol[8]+1.479019945774904*uvar0l[0]*rhol[8]+1.677050983124842*rhor[7]*uvar0r[7]-0.9682458365518543*rhor[5]*uvar0r[7]+1.677050983124842*rhol[7]*uvar0l[7]+0.9682458365518543*rhol[5]*uvar0l[7]-0.9682458365518543*uvar0r[5]*rhor[7]+0.9682458365518543*uvar0l[5]*rhol[7]+2.795084971874738*rhor[6]*uvar0r[6]-2.165063509461097*rhor[3]*uvar0r[6]+1.25*rhor[2]*uvar0r[6]+2.795084971874738*rhol[6]*uvar0l[6]+2.165063509461097*rhol[3]*uvar0l[6]+1.25*rhol[2]*uvar0l[6]-2.165063509461097*uvar0r[3]*rhor[6]+1.25*uvar0r[2]*rhor[6]+2.165063509461097*uvar0l[3]*rhol[6]+1.25*uvar0l[2]*rhol[6]+0.5590169943749475*rhor[5]*uvar0r[5]+0.5590169943749475*rhol[5]*uvar0l[5]+2.795084971874738*rhor[4]*uvar0r[4]-2.165063509461096*rhor[1]*uvar0r[4]+1.25*rhor[0]*uvar0r[4]+2.795084971874738*rhol[4]*uvar0l[4]+2.165063509461096*rhol[1]*uvar0l[4]+1.25*rhol[0]*uvar0l[4]-2.165063509461096*uvar0r[1]*rhor[4]+1.25*uvar0r[0]*rhor[4]+2.165063509461096*uvar0l[1]*rhol[4]+1.25*uvar0l[0]*rhol[4]+1.677050983124842*rhor[3]*uvar0r[3]-0.9682458365518543*rhor[2]*uvar0r[3]+1.677050983124842*rhol[3]*uvar0l[3]+0.9682458365518543*rhol[2]*uvar0l[3]-0.9682458365518543*uvar0r[2]*rhor[3]+0.9682458365518543*uvar0l[2]*rhol[3]+0.5590169943749475*rhor[2]*uvar0r[2]+0.5590169943749475*rhol[2]*uvar0l[2]+1.677050983124842*rhor[1]*uvar0r[1]-0.9682458365518543*rhor[0]*uvar0r[1]+1.677050983124842*rhol[1]*uvar0l[1]+0.9682458365518543*rhol[0]*uvar0l[1]-0.9682458365518543*uvar0r[0]*rhor[1]+0.9682458365518543*uvar0l[0]*rhol[1]+0.5590169943749475*rhor[0]*uvar0r[0]+0.5590169943749475*rhol[0]*uvar0l[0]; 
  incr[8] = 27.78038876617821*rhor[8]*vthsq-27.78038876617821*rhol[8]*vthsq-11.7393568818739*rhor[4]*vthsq-11.7393568818739*rhol[4]*vthsq+3.031088913245535*rhor[1]*vthsq-3.031088913245535*rhol[1]*vthsq-1.984313483298443*rhor[11]*uvar0r[11]+1.14564392373896*rhor[9]*uvar0r[11]-1.984313483298443*rhol[11]*uvar0l[11]-1.14564392373896*rhol[9]*uvar0l[11]+1.14564392373896*uvar0r[9]*rhor[11]-1.14564392373896*uvar0l[9]*rhol[11]-4.630064794363034*rhor[10]*uvar0r[10]+3.913118960624629*rhor[6]*uvar0r[10]-3.031088913245535*rhor[3]*uvar0r[10]+1.75*rhor[2]*uvar0r[10]-4.630064794363034*rhol[10]*uvar0l[10]-3.913118960624629*rhol[6]*uvar0l[10]-3.031088913245535*rhol[3]*uvar0l[10]-1.75*rhol[2]*uvar0l[10]+3.913118960624629*uvar0r[6]*rhor[10]-3.031088913245535*uvar0r[3]*rhor[10]+1.75*uvar0r[2]*rhor[10]-3.913118960624629*uvar0l[6]*rhol[10]-3.031088913245535*uvar0l[3]*rhol[10]-1.75*uvar0l[2]*rhol[10]-0.6614378277661477*rhor[9]*uvar0r[9]-0.6614378277661477*rhol[9]*uvar0l[9]-4.630064794363034*rhor[8]*uvar0r[8]+3.913118960624635*rhor[4]*uvar0r[8]-3.031088913245536*rhor[1]*uvar0r[8]+1.75*rhor[0]*uvar0r[8]-4.630064794363034*rhol[8]*uvar0l[8]-3.913118960624635*rhol[4]*uvar0l[8]-3.031088913245536*rhol[1]*uvar0l[8]-1.75*rhol[0]*uvar0l[8]+3.913118960624635*uvar0r[4]*rhor[8]-3.031088913245536*uvar0r[1]*rhor[8]+1.75*uvar0r[0]*rhor[8]-3.913118960624635*uvar0l[4]*rhol[8]-3.031088913245536*uvar0l[1]*rhol[8]-1.75*uvar0l[0]*rhol[8]-1.984313483298443*rhor[7]*uvar0r[7]+1.14564392373896*rhor[5]*uvar0r[7]-1.984313483298443*rhol[7]*uvar0l[7]-1.14564392373896*rhol[5]*uvar0l[7]+1.14564392373896*uvar0r[5]*rhor[7]-1.14564392373896*uvar0l[5]*rhol[7]-3.307189138830738*rhor[6]*uvar0r[6]+2.5617376914899*rhor[3]*uvar0r[6]-1.479019945774904*rhor[2]*uvar0r[6]-3.307189138830738*rhol[6]*uvar0l[6]-2.5617376914899*rhol[3]*uvar0l[6]-1.479019945774904*rhol[2]*uvar0l[6]+2.5617376914899*uvar0r[3]*rhor[6]-1.479019945774904*uvar0r[2]*rhor[6]-2.5617376914899*uvar0l[3]*rhol[6]-1.479019945774904*uvar0l[2]*rhol[6]-0.6614378277661477*rhor[5]*uvar0r[5]-0.6614378277661477*rhol[5]*uvar0l[5]-3.307189138830738*rhor[4]*uvar0r[4]+2.561737691489899*rhor[1]*uvar0r[4]-1.479019945774904*rhor[0]*uvar0r[4]-3.307189138830738*rhol[4]*uvar0l[4]-2.561737691489899*rhol[1]*uvar0l[4]-1.479019945774904*rhol[0]*uvar0l[4]+2.561737691489899*uvar0r[1]*rhor[4]-1.479019945774904*uvar0r[0]*rhor[4]-2.561737691489899*uvar0l[1]*rhol[4]-1.479019945774904*uvar0l[0]*rhol[4]-1.984313483298443*rhor[3]*uvar0r[3]+1.14564392373896*rhor[2]*uvar0r[3]-1.984313483298443*rhol[3]*uvar0l[3]-1.14564392373896*rhol[2]*uvar0l[3]+1.14564392373896*uvar0r[2]*rhor[3]-1.14564392373896*uvar0l[2]*rhol[3]-0.6614378277661477*rhor[2]*uvar0r[2]-0.6614378277661477*rhol[2]*uvar0l[2]-1.984313483298443*rhor[1]*uvar0r[1]+1.14564392373896*rhor[0]*uvar0r[1]-1.984313483298443*rhol[1]*uvar0l[1]-1.14564392373896*rhol[0]*uvar0l[1]+1.14564392373896*uvar0r[0]*rhor[1]-1.14564392373896*uvar0l[0]*rhol[1]-0.6614378277661477*rhor[0]*uvar0r[0]-0.6614378277661477*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 
  outrrhoux[8] += incr[8]*dxr1; 
  outrrhoux[9] += incr[9]*dxr1; 
  outrrhoux[10] += incr[10]*dxr1; 
  outrrhoux[11] += incr[11]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 
  outlrhoux[4] += -1.0*incr[4]*dxl1; 
  outlrhoux[8] += incr[8]*dxl1; 

 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 
  outrrhouy[8] += incr[8]*dxr1; 
  outrrhouy[9] += incr[9]*dxr1; 
  outrrhouy[10] += incr[10]*dxr1; 
  outrrhouy[11] += incr[11]*dxr1; 


 
} 
GKYL_CU_DH void isoeuler_surfy_2x2v_ser_p3(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[1]; 
  const double dxr1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[12]; 
  const double *rhou1l = &statevecl[24]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[12]; 
  const double *rhou1r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[12]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[12]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[12]; 
  double *outlrhouy = &out[24]; 
  double *outrrho = &out[36]; 
  double *outrrhoux = &out[48]; 
  double *outrrhouy = &out[60]; 
  double incr[12]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.6614378277661477*rhou1r[9])+0.6614378277661477*rhou1l[9]+0.5590169943749475*rhou1r[5]+0.5590169943749475*rhou1l[5]-0.4330127018922193*rhou1r[2]+0.4330127018922193*rhou1l[2]+0.25*rhou1r[0]+0.25*rhou1l[0]; 
  incr[1] = (-0.6614378277661477*rhou1r[11])+0.6614378277661477*rhou1l[11]+0.5590169943749476*rhou1r[7]+0.5590169943749476*rhou1l[7]-0.4330127018922193*rhou1r[3]+0.4330127018922193*rhou1l[3]+0.25*rhou1r[1]+0.25*rhou1l[1]; 
  incr[2] = 1.14564392373896*rhou1r[9]-1.14564392373896*rhou1l[9]-0.9682458365518543*rhou1r[5]-0.9682458365518543*rhou1l[5]+0.75*rhou1r[2]-0.75*rhou1l[2]-0.4330127018922193*rhou1r[0]-0.4330127018922193*rhou1l[0]; 
  incr[3] = 1.14564392373896*rhou1r[11]-1.14564392373896*rhou1l[11]-0.9682458365518543*rhou1r[7]-0.9682458365518543*rhou1l[7]+0.75*rhou1r[3]-0.75*rhou1l[3]-0.4330127018922194*rhou1r[1]-0.4330127018922194*rhou1l[1]; 
  incr[4] = (-0.4330127018922193*rhou1r[6])+0.4330127018922193*rhou1l[6]+0.25*rhou1r[4]+0.25*rhou1l[4]; 
  incr[5] = (-1.479019945774904*rhou1r[9])+1.479019945774904*rhou1l[9]+1.25*rhou1r[5]+1.25*rhou1l[5]-0.9682458365518543*rhou1r[2]+0.9682458365518543*rhou1l[2]+0.5590169943749475*rhou1r[0]+0.5590169943749475*rhou1l[0]; 
  incr[6] = 0.75*rhou1r[6]-0.75*rhou1l[6]-0.4330127018922194*rhou1r[4]-0.4330127018922194*rhou1l[4]; 
  incr[7] = (-1.479019945774904*rhou1r[11])+1.479019945774904*rhou1l[11]+1.25*rhou1r[7]+1.25*rhou1l[7]-0.9682458365518543*rhou1r[3]+0.9682458365518543*rhou1l[3]+0.5590169943749476*rhou1r[1]+0.5590169943749476*rhou1l[1]; 
  incr[8] = (-0.4330127018922193*rhou1r[10])+0.4330127018922193*rhou1l[10]+0.25*rhou1r[8]+0.25*rhou1l[8]; 
  incr[9] = 1.75*rhou1r[9]-1.75*rhou1l[9]-1.479019945774904*rhou1r[5]-1.479019945774904*rhou1l[5]+1.14564392373896*rhou1r[2]-1.14564392373896*rhou1l[2]-0.6614378277661477*rhou1r[0]-0.6614378277661477*rhou1l[0]; 
  incr[10] = 0.75*rhou1r[10]-0.75*rhou1l[10]-0.4330127018922193*rhou1r[8]-0.4330127018922193*rhou1l[8]; 
  incr[11] = 1.75*rhou1r[11]-1.75*rhou1l[11]-1.479019945774904*rhou1r[7]-1.479019945774904*rhou1l[7]+1.14564392373896*rhou1r[3]-1.14564392373896*rhou1l[3]-0.6614378277661477*rhou1r[1]-0.6614378277661477*rhou1l[1]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 
  outrrho[8] += incr[8]*dxr1; 
  outrrho[9] += incr[9]*dxr1; 
  outrrho[10] += incr[10]*dxr1; 
  outrrho[11] += incr[11]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += -1.0*incr[7]*dxl1; 
  outlrho[8] += -1.0*incr[8]*dxl1; 
  outlrho[9] += incr[9]*dxl1; 
  outlrho[10] += incr[10]*dxl1; 
  outlrho[11] += incr[11]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 
  outrrhoux[8] += incr[8]*dxr1; 
  outrrhoux[9] += incr[9]*dxr1; 
  outrrhoux[10] += incr[10]*dxr1; 
  outrrhoux[11] += incr[11]*dxr1; 


 
  //FluxRhoUy; 
  incr[0] = 3.968626966596886*rhor[9]*vthsq+3.968626966596886*rhol[9]*vthsq-1.677050983124842*rhor[5]*vthsq+1.677050983124842*rhol[5]*vthsq+0.4330127018922193*rhor[2]*vthsq+0.4330127018922193*rhol[2]*vthsq+1.75*rhor[11]*uvar1r[11]-1.479019945774903*rhor[7]*uvar1r[11]+1.14564392373896*rhor[3]*uvar1r[11]-0.6614378277661479*rhor[1]*uvar1r[11]+1.75*rhol[11]*uvar1l[11]+1.479019945774903*rhol[7]*uvar1l[11]+1.14564392373896*rhol[3]*uvar1l[11]+0.6614378277661479*rhol[1]*uvar1l[11]-1.479019945774903*uvar1r[7]*rhor[11]+1.14564392373896*uvar1r[3]*rhor[11]-0.6614378277661479*uvar1r[1]*rhor[11]+1.479019945774903*uvar1l[7]*rhol[11]+1.14564392373896*uvar1l[3]*rhol[11]+0.6614378277661479*uvar1l[1]*rhol[11]+0.75*rhor[10]*uvar1r[10]-0.4330127018922193*rhor[8]*uvar1r[10]+0.75*rhol[10]*uvar1l[10]+0.4330127018922193*rhol[8]*uvar1l[10]-0.4330127018922193*uvar1r[8]*rhor[10]+0.4330127018922193*uvar1l[8]*rhol[10]+1.75*rhor[9]*uvar1r[9]-1.479019945774905*rhor[5]*uvar1r[9]+1.14564392373896*rhor[2]*uvar1r[9]-0.6614378277661477*rhor[0]*uvar1r[9]+1.75*rhol[9]*uvar1l[9]+1.479019945774905*rhol[5]*uvar1l[9]+1.14564392373896*rhol[2]*uvar1l[9]+0.6614378277661477*rhol[0]*uvar1l[9]-1.479019945774905*uvar1r[5]*rhor[9]+1.14564392373896*uvar1r[2]*rhor[9]-0.6614378277661477*uvar1r[0]*rhor[9]+1.479019945774905*uvar1l[5]*rhol[9]+1.14564392373896*uvar1l[2]*rhol[9]+0.6614378277661477*uvar1l[0]*rhol[9]+0.25*rhor[8]*uvar1r[8]+0.25*rhol[8]*uvar1l[8]+1.25*rhor[7]*uvar1r[7]-0.9682458365518543*rhor[3]*uvar1r[7]+0.5590169943749472*rhor[1]*uvar1r[7]+1.25*rhol[7]*uvar1l[7]+0.9682458365518543*rhol[3]*uvar1l[7]+0.5590169943749472*rhol[1]*uvar1l[7]-0.9682458365518543*uvar1r[3]*rhor[7]+0.5590169943749472*uvar1r[1]*rhor[7]+0.9682458365518543*uvar1l[3]*rhol[7]+0.5590169943749472*uvar1l[1]*rhol[7]+0.75*rhor[6]*uvar1r[6]-0.4330127018922194*rhor[4]*uvar1r[6]+0.75*rhol[6]*uvar1l[6]+0.4330127018922194*rhol[4]*uvar1l[6]-0.4330127018922194*uvar1r[4]*rhor[6]+0.4330127018922194*uvar1l[4]*rhol[6]+1.25*rhor[5]*uvar1r[5]-0.9682458365518541*rhor[2]*uvar1r[5]+0.5590169943749475*rhor[0]*uvar1r[5]+1.25*rhol[5]*uvar1l[5]+0.9682458365518541*rhol[2]*uvar1l[5]+0.5590169943749475*rhol[0]*uvar1l[5]-0.9682458365518541*uvar1r[2]*rhor[5]+0.5590169943749475*uvar1r[0]*rhor[5]+0.9682458365518541*uvar1l[2]*rhol[5]+0.5590169943749475*uvar1l[0]*rhol[5]+0.25*rhor[4]*uvar1r[4]+0.25*rhol[4]*uvar1l[4]+0.75*rhor[3]*uvar1r[3]-0.4330127018922193*rhor[1]*uvar1r[3]+0.75*rhol[3]*uvar1l[3]+0.4330127018922193*rhol[1]*uvar1l[3]-0.4330127018922193*uvar1r[1]*rhor[3]+0.4330127018922193*uvar1l[1]*rhol[3]+0.75*rhor[2]*uvar1r[2]-0.4330127018922193*rhor[0]*uvar1r[2]+0.75*rhol[2]*uvar1l[2]+0.4330127018922193*rhol[0]*uvar1l[2]-0.4330127018922193*uvar1r[0]*rhor[2]+0.4330127018922193*uvar1l[0]*rhol[2]+0.25*rhor[1]*uvar1r[1]+0.25*rhol[1]*uvar1l[1]+0.25*rhor[0]*uvar1r[0]+0.25*rhol[0]*uvar1l[0]; 
  incr[2] = 11.90588089979066*rhor[9]*vthsq-11.90588089979066*rhol[9]*vthsq-5.031152949374527*rhor[5]*vthsq-5.031152949374527*rhol[5]*vthsq+1.299038105676658*rhor[2]*vthsq-1.299038105676658*rhol[2]*vthsq-3.031088913245535*rhor[11]*uvar1r[11]+2.561737691489897*rhor[7]*uvar1r[11]-1.984313483298443*rhor[3]*uvar1r[11]+1.14564392373896*rhor[1]*uvar1r[11]-3.031088913245535*rhol[11]*uvar1l[11]-2.561737691489897*rhol[7]*uvar1l[11]-1.984313483298443*rhol[3]*uvar1l[11]-1.14564392373896*rhol[1]*uvar1l[11]+2.561737691489897*uvar1r[7]*rhor[11]-1.984313483298443*uvar1r[3]*rhor[11]+1.14564392373896*uvar1r[1]*rhor[11]-2.561737691489897*uvar1l[7]*rhol[11]-1.984313483298443*uvar1l[3]*rhol[11]-1.14564392373896*uvar1l[1]*rhol[11]-1.299038105676658*rhor[10]*uvar1r[10]+0.7499999999999999*rhor[8]*uvar1r[10]-1.299038105676658*rhol[10]*uvar1l[10]-0.7499999999999999*rhol[8]*uvar1l[10]+0.7499999999999999*uvar1r[8]*rhor[10]-0.7499999999999999*uvar1l[8]*rhol[10]-3.031088913245535*rhor[9]*uvar1r[9]+2.561737691489901*rhor[5]*uvar1r[9]-1.984313483298443*rhor[2]*uvar1r[9]+1.14564392373896*rhor[0]*uvar1r[9]-3.031088913245535*rhol[9]*uvar1l[9]-2.561737691489901*rhol[5]*uvar1l[9]-1.984313483298443*rhol[2]*uvar1l[9]-1.14564392373896*rhol[0]*uvar1l[9]+2.561737691489901*uvar1r[5]*rhor[9]-1.984313483298443*uvar1r[2]*rhor[9]+1.14564392373896*uvar1r[0]*rhor[9]-2.561737691489901*uvar1l[5]*rhol[9]-1.984313483298443*uvar1l[2]*rhol[9]-1.14564392373896*uvar1l[0]*rhol[9]-0.4330127018922193*rhor[8]*uvar1r[8]-0.4330127018922193*rhol[8]*uvar1l[8]-2.165063509461096*rhor[7]*uvar1r[7]+1.677050983124842*rhor[3]*uvar1r[7]-0.9682458365518543*rhor[1]*uvar1r[7]-2.165063509461096*rhol[7]*uvar1l[7]-1.677050983124842*rhol[3]*uvar1l[7]-0.9682458365518543*rhol[1]*uvar1l[7]+1.677050983124842*uvar1r[3]*rhor[7]-0.9682458365518543*uvar1r[1]*rhor[7]-1.677050983124842*uvar1l[3]*rhol[7]-0.9682458365518543*uvar1l[1]*rhol[7]-1.299038105676658*rhor[6]*uvar1r[6]+0.75*rhor[4]*uvar1r[6]-1.299038105676658*rhol[6]*uvar1l[6]-0.75*rhol[4]*uvar1l[6]+0.75*uvar1r[4]*rhor[6]-0.75*uvar1l[4]*rhol[6]-2.165063509461096*rhor[5]*uvar1r[5]+1.677050983124842*rhor[2]*uvar1r[5]-0.9682458365518543*rhor[0]*uvar1r[5]-2.165063509461096*rhol[5]*uvar1l[5]-1.677050983124842*rhol[2]*uvar1l[5]-0.9682458365518543*rhol[0]*uvar1l[5]+1.677050983124842*uvar1r[2]*rhor[5]-0.9682458365518543*uvar1r[0]*rhor[5]-1.677050983124842*uvar1l[2]*rhol[5]-0.9682458365518543*uvar1l[0]*rhol[5]-0.4330127018922193*rhor[4]*uvar1r[4]-0.4330127018922193*rhol[4]*uvar1l[4]-1.299038105676658*rhor[3]*uvar1r[3]+0.75*rhor[1]*uvar1r[3]-1.299038105676658*rhol[3]*uvar1l[3]-0.75*rhol[1]*uvar1l[3]+0.75*uvar1r[1]*rhor[3]-0.75*uvar1l[1]*rhol[3]-1.299038105676658*rhor[2]*uvar1r[2]+0.75*rhor[0]*uvar1r[2]-1.299038105676658*rhol[2]*uvar1l[2]-0.75*rhol[0]*uvar1l[2]+0.75*uvar1r[0]*rhor[2]-0.75*uvar1l[0]*rhol[2]-0.4330127018922193*rhor[1]*uvar1r[1]-0.4330127018922193*rhol[1]*uvar1l[1]-0.4330127018922193*rhor[0]*uvar1r[0]-0.4330127018922193*rhol[0]*uvar1l[0]; 
  incr[5] = 19.84313483298443*rhor[9]*vthsq+19.84313483298443*rhol[9]*vthsq-8.385254915624213*rhor[5]*vthsq+8.385254915624213*rhol[5]*vthsq+2.165063509461096*rhor[2]*vthsq+2.165063509461096*rhol[2]*vthsq+3.913118960624632*rhor[11]*uvar1r[11]-3.307189138830736*rhor[7]*uvar1r[11]+2.561737691489899*rhor[3]*uvar1r[11]-1.479019945774905*rhor[1]*uvar1r[11]+3.913118960624632*rhol[11]*uvar1l[11]+3.307189138830736*rhol[7]*uvar1l[11]+2.561737691489899*rhol[3]*uvar1l[11]+1.479019945774905*rhol[1]*uvar1l[11]-3.307189138830736*uvar1r[7]*rhor[11]+2.561737691489899*uvar1r[3]*rhor[11]-1.479019945774905*uvar1r[1]*rhor[11]+3.307189138830736*uvar1l[7]*rhol[11]+2.561737691489899*uvar1l[3]*rhol[11]+1.479019945774905*uvar1l[1]*rhol[11]+1.677050983124842*rhor[10]*uvar1r[10]-0.9682458365518541*rhor[8]*uvar1r[10]+1.677050983124842*rhol[10]*uvar1l[10]+0.9682458365518541*rhol[8]*uvar1l[10]-0.9682458365518541*uvar1r[8]*rhor[10]+0.9682458365518541*uvar1l[8]*rhol[10]+3.913118960624632*rhor[9]*uvar1r[9]-3.307189138830738*rhor[5]*uvar1r[9]+2.561737691489901*rhor[2]*uvar1r[9]-1.479019945774904*rhor[0]*uvar1r[9]+3.913118960624632*rhol[9]*uvar1l[9]+3.307189138830738*rhol[5]*uvar1l[9]+2.561737691489901*rhol[2]*uvar1l[9]+1.479019945774904*rhol[0]*uvar1l[9]-3.307189138830738*uvar1r[5]*rhor[9]+2.561737691489901*uvar1r[2]*rhor[9]-1.479019945774904*uvar1r[0]*rhor[9]+3.307189138830738*uvar1l[5]*rhol[9]+2.561737691489901*uvar1l[2]*rhol[9]+1.479019945774904*uvar1l[0]*rhol[9]+0.5590169943749475*rhor[8]*uvar1r[8]+0.5590169943749475*rhol[8]*uvar1l[8]+2.795084971874738*rhor[7]*uvar1r[7]-2.165063509461097*rhor[3]*uvar1r[7]+1.25*rhor[1]*uvar1r[7]+2.795084971874738*rhol[7]*uvar1l[7]+2.165063509461097*rhol[3]*uvar1l[7]+1.25*rhol[1]*uvar1l[7]-2.165063509461097*uvar1r[3]*rhor[7]+1.25*uvar1r[1]*rhor[7]+2.165063509461097*uvar1l[3]*rhol[7]+1.25*uvar1l[1]*rhol[7]+1.677050983124842*rhor[6]*uvar1r[6]-0.9682458365518543*rhor[4]*uvar1r[6]+1.677050983124842*rhol[6]*uvar1l[6]+0.9682458365518543*rhol[4]*uvar1l[6]-0.9682458365518543*uvar1r[4]*rhor[6]+0.9682458365518543*uvar1l[4]*rhol[6]+2.795084971874738*rhor[5]*uvar1r[5]-2.165063509461096*rhor[2]*uvar1r[5]+1.25*rhor[0]*uvar1r[5]+2.795084971874738*rhol[5]*uvar1l[5]+2.165063509461096*rhol[2]*uvar1l[5]+1.25*rhol[0]*uvar1l[5]-2.165063509461096*uvar1r[2]*rhor[5]+1.25*uvar1r[0]*rhor[5]+2.165063509461096*uvar1l[2]*rhol[5]+1.25*uvar1l[0]*rhol[5]+0.5590169943749475*rhor[4]*uvar1r[4]+0.5590169943749475*rhol[4]*uvar1l[4]+1.677050983124842*rhor[3]*uvar1r[3]-0.9682458365518543*rhor[1]*uvar1r[3]+1.677050983124842*rhol[3]*uvar1l[3]+0.9682458365518543*rhol[1]*uvar1l[3]-0.9682458365518543*uvar1r[1]*rhor[3]+0.9682458365518543*uvar1l[1]*rhol[3]+1.677050983124842*rhor[2]*uvar1r[2]-0.9682458365518543*rhor[0]*uvar1r[2]+1.677050983124842*rhol[2]*uvar1l[2]+0.9682458365518543*rhol[0]*uvar1l[2]-0.9682458365518543*uvar1r[0]*rhor[2]+0.9682458365518543*uvar1l[0]*rhol[2]+0.5590169943749475*rhor[1]*uvar1r[1]+0.5590169943749475*rhol[1]*uvar1l[1]+0.5590169943749475*rhor[0]*uvar1r[0]+0.5590169943749475*rhol[0]*uvar1l[0]; 
  incr[9] = 27.78038876617821*rhor[9]*vthsq-27.78038876617821*rhol[9]*vthsq-11.7393568818739*rhor[5]*vthsq-11.7393568818739*rhol[5]*vthsq+3.031088913245535*rhor[2]*vthsq-3.031088913245535*rhol[2]*vthsq-4.630064794363034*rhor[11]*uvar1r[11]+3.913118960624629*rhor[7]*uvar1r[11]-3.031088913245535*rhor[3]*uvar1r[11]+1.75*rhor[1]*uvar1r[11]-4.630064794363034*rhol[11]*uvar1l[11]-3.913118960624629*rhol[7]*uvar1l[11]-3.031088913245535*rhol[3]*uvar1l[11]-1.75*rhol[1]*uvar1l[11]+3.913118960624629*uvar1r[7]*rhor[11]-3.031088913245535*uvar1r[3]*rhor[11]+1.75*uvar1r[1]*rhor[11]-3.913118960624629*uvar1l[7]*rhol[11]-3.031088913245535*uvar1l[3]*rhol[11]-1.75*uvar1l[1]*rhol[11]-1.984313483298443*rhor[10]*uvar1r[10]+1.14564392373896*rhor[8]*uvar1r[10]-1.984313483298443*rhol[10]*uvar1l[10]-1.14564392373896*rhol[8]*uvar1l[10]+1.14564392373896*uvar1r[8]*rhor[10]-1.14564392373896*uvar1l[8]*rhol[10]-4.630064794363034*rhor[9]*uvar1r[9]+3.913118960624635*rhor[5]*uvar1r[9]-3.031088913245536*rhor[2]*uvar1r[9]+1.75*rhor[0]*uvar1r[9]-4.630064794363034*rhol[9]*uvar1l[9]-3.913118960624635*rhol[5]*uvar1l[9]-3.031088913245536*rhol[2]*uvar1l[9]-1.75*rhol[0]*uvar1l[9]+3.913118960624635*uvar1r[5]*rhor[9]-3.031088913245536*uvar1r[2]*rhor[9]+1.75*uvar1r[0]*rhor[9]-3.913118960624635*uvar1l[5]*rhol[9]-3.031088913245536*uvar1l[2]*rhol[9]-1.75*uvar1l[0]*rhol[9]-0.6614378277661477*rhor[8]*uvar1r[8]-0.6614378277661477*rhol[8]*uvar1l[8]-3.307189138830738*rhor[7]*uvar1r[7]+2.5617376914899*rhor[3]*uvar1r[7]-1.479019945774904*rhor[1]*uvar1r[7]-3.307189138830738*rhol[7]*uvar1l[7]-2.5617376914899*rhol[3]*uvar1l[7]-1.479019945774904*rhol[1]*uvar1l[7]+2.5617376914899*uvar1r[3]*rhor[7]-1.479019945774904*uvar1r[1]*rhor[7]-2.5617376914899*uvar1l[3]*rhol[7]-1.479019945774904*uvar1l[1]*rhol[7]-1.984313483298443*rhor[6]*uvar1r[6]+1.14564392373896*rhor[4]*uvar1r[6]-1.984313483298443*rhol[6]*uvar1l[6]-1.14564392373896*rhol[4]*uvar1l[6]+1.14564392373896*uvar1r[4]*rhor[6]-1.14564392373896*uvar1l[4]*rhol[6]-3.307189138830738*rhor[5]*uvar1r[5]+2.561737691489899*rhor[2]*uvar1r[5]-1.479019945774904*rhor[0]*uvar1r[5]-3.307189138830738*rhol[5]*uvar1l[5]-2.561737691489899*rhol[2]*uvar1l[5]-1.479019945774904*rhol[0]*uvar1l[5]+2.561737691489899*uvar1r[2]*rhor[5]-1.479019945774904*uvar1r[0]*rhor[5]-2.561737691489899*uvar1l[2]*rhol[5]-1.479019945774904*uvar1l[0]*rhol[5]-0.6614378277661477*rhor[4]*uvar1r[4]-0.6614378277661477*rhol[4]*uvar1l[4]-1.984313483298443*rhor[3]*uvar1r[3]+1.14564392373896*rhor[1]*uvar1r[3]-1.984313483298443*rhol[3]*uvar1l[3]-1.14564392373896*rhol[1]*uvar1l[3]+1.14564392373896*uvar1r[1]*rhor[3]-1.14564392373896*uvar1l[1]*rhol[3]-1.984313483298443*rhor[2]*uvar1r[2]+1.14564392373896*rhor[0]*uvar1r[2]-1.984313483298443*rhol[2]*uvar1l[2]-1.14564392373896*rhol[0]*uvar1l[2]+1.14564392373896*uvar1r[0]*rhor[2]-1.14564392373896*uvar1l[0]*rhol[2]-0.6614378277661477*rhor[1]*uvar1r[1]-0.6614378277661477*rhol[1]*uvar1l[1]-0.6614378277661477*rhor[0]*uvar1r[0]-0.6614378277661477*rhol[0]*uvar1l[0]; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 
  outrrhouy[8] += incr[8]*dxr1; 
  outrrhouy[9] += incr[9]*dxr1; 
  outrrhouy[10] += incr[10]*dxr1; 
  outrrhouy[11] += incr[11]*dxr1; 

  outlrhouy[0] += -1.0*incr[0]*dxl1; 
  outlrhouy[2] += incr[2]*dxl1; 
  outlrhouy[5] += -1.0*incr[5]*dxl1; 
  outlrhouy[9] += incr[9]*dxl1; 

 
} 
