#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_1x1v_ser_p1(const double *w, const double *dxv, double mass, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented pressure force in center cell (involves integral over velocity space).
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  const double volFact = dxv[1]/2.0; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  double alpha_l[6] = {0.0}; 
  double alpha_c[6] = {0.0}; 
  double alpha_r[6] = {0.0}; 
  alpha_l[0] = 0.7071067811865475*bl[1]*fl[1]*wvpar_sq+0.7071067811865475*bl[0]*fl[0]*wvpar_sq+0.408248290463863*bl[1]*fl[3]*dvpar*wvpar+0.408248290463863*bl[0]*fl[2]*dvpar*wvpar+0.05270462766947299*bl[1]*fl[5]*dvpar_sq+0.05270462766947297*bl[0]*fl[4]*dvpar_sq+0.05892556509887893*bl[1]*fl[1]*dvpar_sq+0.05892556509887893*bl[0]*fl[0]*dvpar_sq; 
  alpha_l[1] = 0.7071067811865475*bl[0]*fl[1]*wvpar_sq+0.7071067811865475*fl[0]*bl[1]*wvpar_sq+0.408248290463863*bl[0]*fl[3]*dvpar*wvpar+0.408248290463863*bl[1]*fl[2]*dvpar*wvpar+0.05270462766947299*bl[0]*fl[5]*dvpar_sq+0.05270462766947297*bl[1]*fl[4]*dvpar_sq+0.05892556509887893*bl[0]*fl[1]*dvpar_sq+0.05892556509887893*fl[0]*bl[1]*dvpar_sq; 
  alpha_l[2] = 0.7071067811865475*bl[1]*fl[3]*wvpar_sq+0.7071067811865475*bl[0]*fl[2]*wvpar_sq+0.3651483716701107*bl[1]*fl[5]*dvpar*wvpar+0.3651483716701108*bl[0]*fl[4]*dvpar*wvpar+0.408248290463863*bl[1]*fl[1]*dvpar*wvpar+0.408248290463863*bl[0]*fl[0]*dvpar*wvpar+0.1060660171779821*bl[1]*fl[3]*dvpar_sq+0.1060660171779821*bl[0]*fl[2]*dvpar_sq; 
  alpha_l[3] = 0.7071067811865475*bl[0]*fl[3]*wvpar_sq+0.7071067811865475*bl[1]*fl[2]*wvpar_sq+0.3651483716701107*bl[0]*fl[5]*dvpar*wvpar+0.3651483716701108*bl[1]*fl[4]*dvpar*wvpar+0.408248290463863*bl[0]*fl[1]*dvpar*wvpar+0.408248290463863*fl[0]*bl[1]*dvpar*wvpar+0.1060660171779821*bl[0]*fl[3]*dvpar_sq+0.1060660171779821*bl[1]*fl[2]*dvpar_sq; 
  alpha_l[4] = 0.7071067811865475*bl[1]*fl[5]*wvpar_sq+0.7071067811865475*bl[0]*fl[4]*wvpar_sq+0.3651483716701108*bl[1]*fl[3]*dvpar*wvpar+0.3651483716701108*bl[0]*fl[2]*dvpar*wvpar+0.09259731658395264*bl[1]*fl[5]*dvpar_sq+0.09259731658395262*bl[0]*fl[4]*dvpar_sq+0.05270462766947297*bl[1]*fl[1]*dvpar_sq+0.05270462766947297*bl[0]*fl[0]*dvpar_sq; 
  alpha_l[5] = 0.7071067811865475*bl[0]*fl[5]*wvpar_sq+0.7071067811865475*bl[1]*fl[4]*wvpar_sq+0.3651483716701107*bl[0]*fl[3]*dvpar*wvpar+0.3651483716701107*bl[1]*fl[2]*dvpar*wvpar+0.09259731658395262*bl[0]*fl[5]*dvpar_sq+0.09259731658395264*bl[1]*fl[4]*dvpar_sq+0.05270462766947299*bl[0]*fl[1]*dvpar_sq+0.05270462766947299*fl[0]*bl[1]*dvpar_sq; 

  alpha_c[0] = 0.7071067811865475*bc[1]*fc[1]*wvpar_sq+0.7071067811865475*bc[0]*fc[0]*wvpar_sq+0.408248290463863*bc[1]*fc[3]*dvpar*wvpar+0.408248290463863*bc[0]*fc[2]*dvpar*wvpar+0.05270462766947299*bc[1]*fc[5]*dvpar_sq+0.05270462766947297*bc[0]*fc[4]*dvpar_sq+0.05892556509887893*bc[1]*fc[1]*dvpar_sq+0.05892556509887893*bc[0]*fc[0]*dvpar_sq; 
  alpha_c[1] = 0.7071067811865475*bc[0]*fc[1]*wvpar_sq+0.7071067811865475*fc[0]*bc[1]*wvpar_sq+0.408248290463863*bc[0]*fc[3]*dvpar*wvpar+0.408248290463863*bc[1]*fc[2]*dvpar*wvpar+0.05270462766947299*bc[0]*fc[5]*dvpar_sq+0.05270462766947297*bc[1]*fc[4]*dvpar_sq+0.05892556509887893*bc[0]*fc[1]*dvpar_sq+0.05892556509887893*fc[0]*bc[1]*dvpar_sq; 
  alpha_c[2] = 0.7071067811865475*bc[1]*fc[3]*wvpar_sq+0.7071067811865475*bc[0]*fc[2]*wvpar_sq+0.3651483716701107*bc[1]*fc[5]*dvpar*wvpar+0.3651483716701108*bc[0]*fc[4]*dvpar*wvpar+0.408248290463863*bc[1]*fc[1]*dvpar*wvpar+0.408248290463863*bc[0]*fc[0]*dvpar*wvpar+0.1060660171779821*bc[1]*fc[3]*dvpar_sq+0.1060660171779821*bc[0]*fc[2]*dvpar_sq; 
  alpha_c[3] = 0.7071067811865475*bc[0]*fc[3]*wvpar_sq+0.7071067811865475*bc[1]*fc[2]*wvpar_sq+0.3651483716701107*bc[0]*fc[5]*dvpar*wvpar+0.3651483716701108*bc[1]*fc[4]*dvpar*wvpar+0.408248290463863*bc[0]*fc[1]*dvpar*wvpar+0.408248290463863*fc[0]*bc[1]*dvpar*wvpar+0.1060660171779821*bc[0]*fc[3]*dvpar_sq+0.1060660171779821*bc[1]*fc[2]*dvpar_sq; 
  alpha_c[4] = 0.7071067811865475*bc[1]*fc[5]*wvpar_sq+0.7071067811865475*bc[0]*fc[4]*wvpar_sq+0.3651483716701108*bc[1]*fc[3]*dvpar*wvpar+0.3651483716701108*bc[0]*fc[2]*dvpar*wvpar+0.09259731658395264*bc[1]*fc[5]*dvpar_sq+0.09259731658395262*bc[0]*fc[4]*dvpar_sq+0.05270462766947297*bc[1]*fc[1]*dvpar_sq+0.05270462766947297*bc[0]*fc[0]*dvpar_sq; 
  alpha_c[5] = 0.7071067811865475*bc[0]*fc[5]*wvpar_sq+0.7071067811865475*bc[1]*fc[4]*wvpar_sq+0.3651483716701107*bc[0]*fc[3]*dvpar*wvpar+0.3651483716701107*bc[1]*fc[2]*dvpar*wvpar+0.09259731658395262*bc[0]*fc[5]*dvpar_sq+0.09259731658395264*bc[1]*fc[4]*dvpar_sq+0.05270462766947299*bc[0]*fc[1]*dvpar_sq+0.05270462766947299*fc[0]*bc[1]*dvpar_sq; 

  alpha_r[0] = 0.7071067811865475*br[1]*fr[1]*wvpar_sq+0.7071067811865475*br[0]*fr[0]*wvpar_sq+0.408248290463863*br[1]*fr[3]*dvpar*wvpar+0.408248290463863*br[0]*fr[2]*dvpar*wvpar+0.05270462766947299*br[1]*fr[5]*dvpar_sq+0.05270462766947297*br[0]*fr[4]*dvpar_sq+0.05892556509887893*br[1]*fr[1]*dvpar_sq+0.05892556509887893*br[0]*fr[0]*dvpar_sq; 
  alpha_r[1] = 0.7071067811865475*br[0]*fr[1]*wvpar_sq+0.7071067811865475*fr[0]*br[1]*wvpar_sq+0.408248290463863*br[0]*fr[3]*dvpar*wvpar+0.408248290463863*br[1]*fr[2]*dvpar*wvpar+0.05270462766947299*br[0]*fr[5]*dvpar_sq+0.05270462766947297*br[1]*fr[4]*dvpar_sq+0.05892556509887893*br[0]*fr[1]*dvpar_sq+0.05892556509887893*fr[0]*br[1]*dvpar_sq; 
  alpha_r[2] = 0.7071067811865475*br[1]*fr[3]*wvpar_sq+0.7071067811865475*br[0]*fr[2]*wvpar_sq+0.3651483716701107*br[1]*fr[5]*dvpar*wvpar+0.3651483716701108*br[0]*fr[4]*dvpar*wvpar+0.408248290463863*br[1]*fr[1]*dvpar*wvpar+0.408248290463863*br[0]*fr[0]*dvpar*wvpar+0.1060660171779821*br[1]*fr[3]*dvpar_sq+0.1060660171779821*br[0]*fr[2]*dvpar_sq; 
  alpha_r[3] = 0.7071067811865475*br[0]*fr[3]*wvpar_sq+0.7071067811865475*br[1]*fr[2]*wvpar_sq+0.3651483716701107*br[0]*fr[5]*dvpar*wvpar+0.3651483716701108*br[1]*fr[4]*dvpar*wvpar+0.408248290463863*br[0]*fr[1]*dvpar*wvpar+0.408248290463863*fr[0]*br[1]*dvpar*wvpar+0.1060660171779821*br[0]*fr[3]*dvpar_sq+0.1060660171779821*br[1]*fr[2]*dvpar_sq; 
  alpha_r[4] = 0.7071067811865475*br[1]*fr[5]*wvpar_sq+0.7071067811865475*br[0]*fr[4]*wvpar_sq+0.3651483716701108*br[1]*fr[3]*dvpar*wvpar+0.3651483716701108*br[0]*fr[2]*dvpar*wvpar+0.09259731658395264*br[1]*fr[5]*dvpar_sq+0.09259731658395262*br[0]*fr[4]*dvpar_sq+0.05270462766947297*br[1]*fr[1]*dvpar_sq+0.05270462766947297*br[0]*fr[0]*dvpar_sq; 
  alpha_r[5] = 0.7071067811865475*br[0]*fr[5]*wvpar_sq+0.7071067811865475*br[1]*fr[4]*wvpar_sq+0.3651483716701107*br[0]*fr[3]*dvpar*wvpar+0.3651483716701107*br[1]*fr[2]*dvpar*wvpar+0.09259731658395262*br[0]*fr[5]*dvpar_sq+0.09259731658395264*br[1]*fr[4]*dvpar_sq+0.05270462766947299*br[0]*fr[1]*dvpar_sq+0.05270462766947299*fr[0]*br[1]*dvpar_sq; 

  out[0] += ((-0.408248290463863*(alpha_r[1]+alpha_l[1]))+0.8164965809277261*alpha_c[1]+0.3535533905932737*alpha_r[0]-0.3535533905932737*alpha_l[0])*dx1*mass*volFact; 
  out[1] += ((-0.7071067811865475*alpha_r[1])+0.7071067811865475*alpha_l[1]+0.6123724356957944*(alpha_r[0]+alpha_l[0])-1.224744871391589*alpha_c[0])*dx1*mass*volFact; 

} 
