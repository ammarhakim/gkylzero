#include <gkyl_bc_sheath_gyrokinetic_kernels.h> 


GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZQuad[2][6]; 
  

  vcutSq = -0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[3]-1.0*f[2]))/(f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.7071067811865475*(f[1]-1.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[3]+f[2]))/(f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.7071067811865475*(f[1]+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } 
  } 
  fRefl[0] = 0.7071067811865475*(fReflZQuad[1][0]+fReflZQuad[0][0]); 
  fRefl[1] = 0.7071067811865475*(fReflZQuad[1][0]-1.0*fReflZQuad[0][0]); 
  fRefl[2] = 0.7071067811865475*(fReflZQuad[1][1]+fReflZQuad[0][1]); 
  fRefl[3] = 0.7071067811865475*(fReflZQuad[1][1]-1.0*fReflZQuad[0][1]); 
  fRefl[4] = 0.7071067811865475*(fReflZQuad[1][2]+fReflZQuad[0][2]); 
  fRefl[5] = 0.7071067811865475*(fReflZQuad[1][2]-1.0*fReflZQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZQuad[2][6]; 
  

  vcutSq = 0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[3]-1.0*f[2]))/(f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.7071067811865475*(f[1]-1.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[3]+f[2]))/(f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.7071067811865475*(f[1]+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } 
  } 
  fRefl[0] = 0.7071067811865475*(fReflZQuad[1][0]+fReflZQuad[0][0]); 
  fRefl[1] = 0.7071067811865475*(fReflZQuad[1][0]-1.0*fReflZQuad[0][0]); 
  fRefl[2] = 0.7071067811865475*(fReflZQuad[1][1]+fReflZQuad[0][1]); 
  fRefl[3] = 0.7071067811865475*(fReflZQuad[1][1]-1.0*fReflZQuad[0][1]); 
  fRefl[4] = 0.7071067811865475*(fReflZQuad[1][2]+fReflZQuad[0][2]); 
  fRefl[5] = 0.7071067811865475*(fReflZQuad[1][2]-1.0*fReflZQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZMuQuad[4][6]; 
  

  vcutSq = -0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.0; 
  fRefl[10] = 0.0; 
  fRefl[11] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  fRefl[8] = f[8]; 
  fRefl[9] = f[9]; 
  fRefl[10] = f[10]; 
  fRefl[11] = f[11]; 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[7]-1.0*(f[6]+f[4])+f[2]))/(f[5]-1.0*(f[3]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))/(f[5]-1.0*f[3]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]-1.0*(f[4]+f[2])))/(f[5]+f[3]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]+f[4]+f[2]))/(f[5]+f[3]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]+f[3]+f[1]+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]+fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]-1.0*(fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[2] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]+fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 0.5*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[2][0]+fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[4] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]-1.0*(fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[5] = 0.5*(fReflZMuQuad[3][0]-1.0*(fReflZMuQuad[2][0]+fReflZMuQuad[1][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.5*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[2][1]+fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[7] = 0.5*(fReflZMuQuad[3][1]-1.0*(fReflZMuQuad[2][1]+fReflZMuQuad[1][1])+fReflZMuQuad[0][1]); 
  fRefl[8] = 0.5*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]+fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.5000000000000001*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]-1.0*(fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[10] = 0.5000000000000001*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[2][2]+fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[11] = 0.5*(fReflZMuQuad[3][2]-1.0*(fReflZMuQuad[2][2]+fReflZMuQuad[1][2])+fReflZMuQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZMuQuad[4][6]; 
  

  vcutSq = 0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.0; 
  fRefl[10] = 0.0; 
  fRefl[11] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  fRefl[8] = f[8]; 
  fRefl[9] = f[9]; 
  fRefl[10] = f[10]; 
  fRefl[11] = f[11]; 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[7]-1.0*(f[6]+f[4])+f[2]))/(f[5]-1.0*(f[3]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))/(f[5]-1.0*f[3]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]-1.0*(f[4]+f[2])))/(f[5]+f[3]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]+f[4]+f[2]))/(f[5]+f[3]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]+f[3]+f[1]+f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<2e-8? (1.+xc)/2. : (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<2e-8? (1.-xc)/2. : ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]+fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]-1.0*(fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[2] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]+fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 0.5*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[2][0]+fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[4] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]-1.0*(fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[5] = 0.5*(fReflZMuQuad[3][0]-1.0*(fReflZMuQuad[2][0]+fReflZMuQuad[1][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.5*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[2][1]+fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[7] = 0.5*(fReflZMuQuad[3][1]-1.0*(fReflZMuQuad[2][1]+fReflZMuQuad[1][1])+fReflZMuQuad[0][1]); 
  fRefl[8] = 0.5*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]+fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.5000000000000001*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]-1.0*(fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[10] = 0.5000000000000001*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[2][2]+fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[11] = 0.5*(fReflZMuQuad[3][2]-1.0*(fReflZMuQuad[2][2]+fReflZMuQuad[1][2])+fReflZMuQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1(double wv, double dv, double vlowerSq, double vupperSq, double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq;
  double fReflSurfXY[4][6] = {0.}; 
  // node (x,y)_0 
  vcutSq = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*phi[6]-2.449489742783178*phiWall[5]+2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*phi[2]+1.414213562373095*phiWall[1]-1.414213562373095*(phi[1]+phiWall[0])+1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[0][0] = 0.0; 
  fReflSurfXY[0][1] = 0.0; 
  fReflSurfXY[0][2] = 0.0; 
  fReflSurfXY[0][3] = 0.0; 
  fReflSurfXY[0][4] = 0.0; 
  fReflSurfXY[0][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[0][0] = -0.3535533905932737*(1.732050807568877*f[16]-1.0*(1.732050807568877*(f[8]+f[7])+f[6])+1.732050807568877*f[3]+f[2]+f[1]-1.0*f[0]); 
  fReflSurfXY[0][1] = -0.3535533905932737*(1.732050807568877*f[26]-1.0*(1.732050807568877*(f[19]+f[18])+f[17])+1.732050807568877*f[11]+f[10]+f[9]-1.0*f[4]); 
  fReflSurfXY[0][2] = -0.3535533905932737*(1.732050807568877*f[27]-1.0*(1.732050807568877*(f[22]+f[21])+f[20])+1.732050807568877*f[14]+f[13]+f[12]-1.0*f[5]); 
  fReflSurfXY[0][3] = -0.3535533905932737*(1.732050807568877*f[31]-1.0*(1.732050807568877*(f[30]+f[29])+f[28])+1.732050807568877*f[25]+f[24]+f[23]-1.0*f[15]); 
  fReflSurfXY[0][4] = -0.02357022603955158*(25.98076211353316*f[43]-5.0*(5.196152422706631*(f[39]+f[38])+3.0*f[37])+8.660254037844387*(3.0*f[35]+1.732050807568877*(f[34]+f[33]))-15.0*f[32]); 
  fReflSurfXY[0][5] = -0.02357022603955158*(25.98076211353316*f[47]-5.0*(5.196152422706631*(f[46]+f[45])+3.0*f[44])+8.660254037844387*(3.0*f[42]+1.732050807568877*(f[41]+f[40]))-15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]-8.0*(f[10]+f[9])+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33])-7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])-6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]+6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33])-3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])+12.0*(f[22]+f[21])+6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])+8.0*f[17]+12.0*f[16]-8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]-8.0*(f[10]+f[9])-12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])-174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]+174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9])+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])+116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]-116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9])-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])-150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]-150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])-120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])+120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]-120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])+58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]-58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9])-58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]-8.0*(f[10]+f[9])+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33])-7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])-6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]+6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33])-3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*(f[30]+f[29])+8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*(f[24]+f[23])+12.0*(f[22]+f[21])+6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])-8.0*f[17]+12.0*f[16]+8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]+8.0*(f[10]+f[9])-12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])-174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]+174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9])+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])+116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]-116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9])-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])-150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]-150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])+120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])-120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]+120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])+58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]-58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9])-58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])+13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])-8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*(f[10]+f[9])-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33])+7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])-6.928203230275509*f[20]+12.0*(f[16]+f[14])+6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33])+3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])+12.0*(f[22]+f[21])+6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])-8.0*f[17]-12.0*f[16]-8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]+8.0*(f[10]+f[9])+12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])-174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])+174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9])-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])+116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])-116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9])+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])-150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]+150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])-120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])-120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]+120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])+58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])-58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9])+58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])+13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])-8.0*f[28]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])-8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*(f[10]+f[9])-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33])+7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])-6.928203230275509*f[20]+12.0*(f[16]+f[14])+6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33])+3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*(f[30]+f[29])+8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])+12.0*(f[22]+f[21])+6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])+8.0*f[17]-12.0*f[16]+8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]-8.0*(f[10]+f[9])+12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])-174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])+174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9])-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])+116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])-116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9])+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])-150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]+150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])+120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])+120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]-120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])+58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])-58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9])+58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[0][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[0][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[0][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[0][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[0][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[0][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_1 
  vcutSq = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*(phi[6]+phiWall[5])-2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*(phi[2]+phiWall[1])+1.414213562373095*(phi[1]+phiWall[0])-1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[1][0] = 0.0; 
  fReflSurfXY[1][1] = 0.0; 
  fReflSurfXY[1][2] = 0.0; 
  fReflSurfXY[1][3] = 0.0; 
  fReflSurfXY[1][4] = 0.0; 
  fReflSurfXY[1][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[1][0] = 0.3535533905932737*(1.732050807568877*(f[16]-1.0*f[8]+f[7])-1.0*(f[6]+1.732050807568877*f[3])+f[2]-1.0*f[1]+f[0]); 
  fReflSurfXY[1][1] = 0.3535533905932737*(1.732050807568877*(f[26]-1.0*f[19]+f[18])-1.0*(f[17]+1.732050807568877*f[11])+f[10]-1.0*f[9]+f[4]); 
  fReflSurfXY[1][2] = 0.3535533905932737*(1.732050807568877*(f[27]-1.0*f[22]+f[21])-1.0*(f[20]+1.732050807568877*f[14])+f[13]-1.0*f[12]+f[5]); 
  fReflSurfXY[1][3] = 0.3535533905932737*(1.732050807568877*(f[31]-1.0*f[30]+f[29])-1.0*(f[28]+1.732050807568877*f[25])+f[24]-1.0*f[23]+f[15]); 
  fReflSurfXY[1][4] = 0.02357022603955158*(25.98076211353316*f[43]+5.0*(5.196152422706631*(f[38]-1.0*f[39])-3.0*f[37])+8.660254037844387*(1.732050807568877*(f[34]-1.0*f[33])-3.0*f[35])+15.0*f[32]); 
  fReflSurfXY[1][5] = 0.02357022603955158*(25.98076211353316*f[47]+5.0*(5.196152422706631*(f[45]-1.0*f[46])-3.0*f[44])+8.660254037844387*(1.732050807568877*(f[41]-1.0*f[40])-3.0*f[42])+15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]-11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])+13.41640786499874*f[13]-13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]+11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]-11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]+7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*f[33]+7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]-12.0*(f[16]+f[14])+6.928203230275509*f[13]-6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]+6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]+3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*f[33]+3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]+12.0*f[22]-12.0*f[21]+6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]+12.0*f[16]+8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-12.0*f[8]+12.0*f[7]-6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]-269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]-174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])+174.2842505793337*f[24]-174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*f[9]-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]+300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]+300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]+116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])-116.1895003862225*f[24]+116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*f[9]+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]-150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]-150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]-120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]+120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]+58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])-58.09475019311126*f[24]+58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*f[9]+58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]-11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])+13.41640786499874*f[13]-13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]+11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]+13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]-11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]+7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*f[33]+7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]-12.0*(f[16]+f[14])+6.928203230275509*f[13]-6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]+6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]+3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*f[33]+3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*f[30]-13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]+12.0*f[22]-12.0*f[21]+6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]+12.0*f[16]-8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]-12.0*f[8]+12.0*f[7]-6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]-269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]-174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])+174.2842505793337*f[24]-174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*f[9]-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]+300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]+300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]+116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])-116.1895003862225*f[24]+116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*f[9]+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]-150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]-150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]+120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]-120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]+58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])-58.09475019311126*f[24]+58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*f[9]+58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]-11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]+11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]-11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]+7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*f[33]-7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]+3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*f[33]-3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]+12.0*f[22]-12.0*f[21]+6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]-12.0*f[16]+8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+12.0*f[8]-12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]-269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]-174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]+174.2842505793337*f[24]-174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*f[9]+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]+300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]-300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]+116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]-116.1895003862225*f[24]+116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*f[9]-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]-150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]+150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]-120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]-120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]+58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]-58.09475019311126*f[24]+58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*f[9]-58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]-11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]+11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]+13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]-11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]-8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]+7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*f[33]-7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]+3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*f[33]-3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*f[30]-13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]+12.0*f[22]-12.0*f[21]+6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]-12.0*f[16]-8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]+12.0*f[8]-12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]-269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]-174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]+174.2842505793337*f[24]-174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*f[9]+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]+300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]-300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]+116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]-116.1895003862225*f[24]+116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*f[9]-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]-150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]+150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]+120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]+120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]+58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]-58.09475019311126*f[24]+58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*f[9]-58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[1][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[1][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[1][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[1][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[1][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[1][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_2 
  vcutSq = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*(phi[6]+phiWall[5])+2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*(phi[2]+phiWall[1])-1.414213562373095*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[2][0] = 0.0; 
  fReflSurfXY[2][1] = 0.0; 
  fReflSurfXY[2][2] = 0.0; 
  fReflSurfXY[2][3] = 0.0; 
  fReflSurfXY[2][4] = 0.0; 
  fReflSurfXY[2][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[2][0] = 0.3535533905932737*(1.732050807568877*(f[16]+f[8])-1.0*(1.732050807568877*f[7]+f[6]+1.732050807568877*f[3]+f[2])+f[1]+f[0]); 
  fReflSurfXY[2][1] = 0.3535533905932737*(1.732050807568877*(f[26]+f[19])-1.0*(1.732050807568877*f[18]+f[17]+1.732050807568877*f[11]+f[10])+f[9]+f[4]); 
  fReflSurfXY[2][2] = 0.3535533905932737*(1.732050807568877*(f[27]+f[22])-1.0*(1.732050807568877*f[21]+f[20]+1.732050807568877*f[14]+f[13])+f[12]+f[5]); 
  fReflSurfXY[2][3] = 0.3535533905932737*(1.732050807568877*(f[31]+f[30])-1.0*(1.732050807568877*f[29]+f[28]+1.732050807568877*f[25]+f[24])+f[23]+f[15]); 
  fReflSurfXY[2][4] = 0.02357022603955158*(25.98076211353316*f[43]+5.0*(5.196152422706631*f[39]-1.0*(5.196152422706631*f[38]+3.0*f[37]))+8.660254037844387*(1.732050807568877*(f[33]-1.0*f[34])-3.0*f[35])+15.0*f[32]); 
  fReflSurfXY[2][5] = 0.02357022603955158*(25.98076211353316*f[47]+5.0*(5.196152422706631*f[46]-1.0*(5.196152422706631*f[45]+3.0*f[44]))+8.660254037844387*(1.732050807568877*(f[40]-1.0*f[41])-3.0*f[42])+15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]-11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]-13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])-13.41640786499874*f[13]+13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]+11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]-13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]-11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*f[10]-8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]+7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]-6.928203230275509*f[20]-12.0*(f[16]+f[14])-6.928203230275509*f[13]+6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]+6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]+3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*(f[33]+f[32])+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]+12.0*f[16]+8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+12.0*f[8]-12.0*f[7]-6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]-269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]-174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])-174.2842505793337*f[24]+174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]+300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]+300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]+116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])+116.1895003862225*f[24]-116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]-150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]-150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]-120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]+120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]-120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]+58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])+58.09475019311126*f[24]-58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*(f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]-11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]-13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])-13.41640786499874*f[13]+13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]+11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]-11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*f[10]-8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]+7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]-6.928203230275509*f[20]-12.0*(f[16]+f[14])-6.928203230275509*f[13]+6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]+6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]+3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*(f[33]+f[32])-13.85640646055102*(f[31]+f[30])+13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]+12.0*f[16]-8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]+12.0*f[8]-12.0*f[7]-6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]-269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]-174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])-174.2842505793337*f[24]+174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]+300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]+300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]+116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])+116.1895003862225*f[24]-116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]-150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]-150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]+120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]-120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]+120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]+58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])+58.09475019311126*f[24]-58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*(f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]-11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]-13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]+11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]-13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]-11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]-8.0*f[10]+8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]+7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]-6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]+3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*(f[33]+f[32])+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]-12.0*f[16]+8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-12.0*f[8]+12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]-269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]-174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]-174.2842505793337*f[24]+174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]+300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]-300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]+116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]+116.1895003862225*f[24]-116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]-150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]+150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]-120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]-120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]+120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]+58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]+58.09475019311126*f[24]-58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*(f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]-11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]-13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]+11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]-11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]-8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]-8.0*f[10]+8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]+7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]-6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]+3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*(f[33]+f[32])-13.85640646055102*(f[31]+f[30])+13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]-12.0*f[16]-8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]-12.0*f[8]+12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]-269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]-174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]-174.2842505793337*f[24]+174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]+300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]-300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]+116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]+116.1895003862225*f[24]-116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]-150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]+150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]+120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]+120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]-120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]+58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]+58.09475019311126*f[24]-58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*(f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[2][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[2][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[2][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[2][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[2][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[2][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_3 
  vcutSq = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*phi[6]+2.449489742783178*phiWall[5]-2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*phi[2]-1.414213562373095*phiWall[1]+1.414213562373095*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[3][0] = 0.0; 
  fReflSurfXY[3][1] = 0.0; 
  fReflSurfXY[3][2] = 0.0; 
  fReflSurfXY[3][3] = 0.0; 
  fReflSurfXY[3][4] = 0.0; 
  fReflSurfXY[3][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[3][0] = -0.3535533905932737*(1.732050807568877*(f[16]+f[8]+f[7])-1.0*f[6]+1.732050807568877*f[3]-1.0*(f[2]+f[1]+f[0])); 
  fReflSurfXY[3][1] = -0.3535533905932737*(1.732050807568877*(f[26]+f[19]+f[18])-1.0*f[17]+1.732050807568877*f[11]-1.0*(f[10]+f[9]+f[4])); 
  fReflSurfXY[3][2] = -0.3535533905932737*(1.732050807568877*(f[27]+f[22]+f[21])-1.0*f[20]+1.732050807568877*f[14]-1.0*(f[13]+f[12]+f[5])); 
  fReflSurfXY[3][3] = -0.3535533905932737*(1.732050807568877*(f[31]+f[30]+f[29])-1.0*f[28]+1.732050807568877*f[25]-1.0*(f[24]+f[23]+f[15])); 
  fReflSurfXY[3][4] = -0.02357022603955158*(25.98076211353316*f[43]+25.98076211353316*(f[39]+f[38])+3.0*(8.660254037844387*f[35]-5.0*f[37])-1.0*(15.0*(f[34]+f[33])+15.0*f[32])); 
  fReflSurfXY[3][5] = -0.02357022603955158*(25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+3.0*(8.660254037844387*f[42]-5.0*f[44])-1.0*(15.0*(f[41]+f[40])+15.0*f[36])); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])-13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]+8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])-6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]-6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33]+f[32])+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])-12.0*(f[22]+f[21])+6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])+8.0*f[17]+12.0*f[16]-8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]+8.0*(f[10]+f[9])+12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])-174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]-174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])+116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]+116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])-150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]+150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])-120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])+120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]+120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])+58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]+58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])-13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]+8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])-6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]-6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33]+f[32])-13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*(f[24]+f[23])-12.0*(f[22]+f[21])+6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])-8.0*f[17]+12.0*f[16]+8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]-8.0*(f[10]+f[9])+12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])-174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]-174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])+116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]+116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])-150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]+150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])+120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])-120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]-120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])+58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]+58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])-13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])-13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])-8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])-6.928203230275509*f[20]+12.0*(f[16]+f[14])-6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33]+f[32])+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])-12.0*(f[22]+f[21])+6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])-8.0*f[17]-12.0*f[16]-8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]-8.0*(f[10]+f[9])-12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])-174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])-174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])+116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])+116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])-150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]-150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])-120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])-120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]-120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])+58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])+58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])-13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])-13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])-11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])-8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])+7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])-6.928203230275509*f[20]+12.0*(f[16]+f[14])-6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])+3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33]+f[32])-13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])-12.0*(f[22]+f[21])+6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])+8.0*f[17]-12.0*f[16]+8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]+8.0*(f[10]+f[9])-12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])-269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])-174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])-174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])+300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])+116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])+116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])-150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]-150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])+120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])+120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]+120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])+58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])+58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[3][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[3][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[3][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[3][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[3][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[3][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  fRefl[0] = 0.7071067811865475*(fReflSurfXY[3][0]+fReflSurfXY[2][0]+fReflSurfXY[1][0]+fReflSurfXY[0][0]); 
  fRefl[1] = 0.7071067811865475*(fReflSurfXY[3][0]+fReflSurfXY[2][0]-1.0*(fReflSurfXY[1][0]+fReflSurfXY[0][0])); 
  fRefl[2] = 0.7071067811865475*(fReflSurfXY[3][0]-1.0*fReflSurfXY[2][0]+fReflSurfXY[1][0]-1.0*fReflSurfXY[0][0]); 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.7071067811865475*(fReflSurfXY[3][1]+fReflSurfXY[2][1]+fReflSurfXY[1][1]+fReflSurfXY[0][1]); 
  fRefl[5] = 0.7071067811865475*(fReflSurfXY[3][2]+fReflSurfXY[2][2]+fReflSurfXY[1][2]+fReflSurfXY[0][2]); 
  fRefl[6] = 0.7071067811865475*(fReflSurfXY[3][0]-1.0*(fReflSurfXY[2][0]+fReflSurfXY[1][0])+fReflSurfXY[0][0]); 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.7071067811865475*(fReflSurfXY[3][1]+fReflSurfXY[2][1]-1.0*(fReflSurfXY[1][1]+fReflSurfXY[0][1])); 
  fRefl[10] = 0.7071067811865475*(fReflSurfXY[3][1]-1.0*fReflSurfXY[2][1]+fReflSurfXY[1][1]-1.0*fReflSurfXY[0][1]); 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.7071067811865475*(fReflSurfXY[3][2]+fReflSurfXY[2][2]-1.0*(fReflSurfXY[1][2]+fReflSurfXY[0][2])); 
  fRefl[13] = 0.7071067811865475*(fReflSurfXY[3][2]-1.0*fReflSurfXY[2][2]+fReflSurfXY[1][2]-1.0*fReflSurfXY[0][2]); 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.7071067811865475*(fReflSurfXY[3][3]+fReflSurfXY[2][3]+fReflSurfXY[1][3]+fReflSurfXY[0][3]); 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.7071067811865475*(fReflSurfXY[3][1]-1.0*(fReflSurfXY[2][1]+fReflSurfXY[1][1])+fReflSurfXY[0][1]); 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  fRefl[20] = 0.7071067811865475*(fReflSurfXY[3][2]-1.0*(fReflSurfXY[2][2]+fReflSurfXY[1][2])+fReflSurfXY[0][2]); 
  fRefl[21] = 0.0; 
  fRefl[22] = 0.0; 
  fRefl[23] = 0.7071067811865475*(fReflSurfXY[3][3]+fReflSurfXY[2][3]-1.0*(fReflSurfXY[1][3]+fReflSurfXY[0][3])); 
  fRefl[24] = 0.7071067811865475*(fReflSurfXY[3][3]-1.0*fReflSurfXY[2][3]+fReflSurfXY[1][3]-1.0*fReflSurfXY[0][3]); 
  fRefl[25] = 0.0; 
  fRefl[26] = 0.0; 
  fRefl[27] = 0.0; 
  fRefl[28] = 0.7071067811865475*(fReflSurfXY[3][3]-1.0*(fReflSurfXY[2][3]+fReflSurfXY[1][3])+fReflSurfXY[0][3]); 
  fRefl[29] = 0.0; 
  fRefl[30] = 0.0; 
  fRefl[31] = 0.0; 
  fRefl[32] = 0.7071067811865475*(fReflSurfXY[3][4]+fReflSurfXY[2][4]+fReflSurfXY[1][4]+fReflSurfXY[0][4]); 
  fRefl[33] = 0.7071067811865475*(fReflSurfXY[3][4]+fReflSurfXY[2][4]-1.0*(fReflSurfXY[1][4]+fReflSurfXY[0][4])); 
  fRefl[34] = 0.7071067811865475*(fReflSurfXY[3][4]-1.0*fReflSurfXY[2][4]+fReflSurfXY[1][4]-1.0*fReflSurfXY[0][4]); 
  fRefl[35] = 0.0; 
  fRefl[36] = 0.7071067811865475*(fReflSurfXY[3][5]+fReflSurfXY[2][5]+fReflSurfXY[1][5]+fReflSurfXY[0][5]); 
  fRefl[37] = 0.7071067811865475*(fReflSurfXY[3][4]-1.0*(fReflSurfXY[2][4]+fReflSurfXY[1][4])+fReflSurfXY[0][4]); 
  fRefl[38] = 0.0; 
  fRefl[39] = 0.0; 
  fRefl[40] = 0.7071067811865475*(fReflSurfXY[3][5]+fReflSurfXY[2][5]-1.0*(fReflSurfXY[1][5]+fReflSurfXY[0][5])); 
  fRefl[41] = 0.7071067811865475*(fReflSurfXY[3][5]-1.0*fReflSurfXY[2][5]+fReflSurfXY[1][5]-1.0*fReflSurfXY[0][5]); 
  fRefl[42] = 0.0; 
  fRefl[43] = 0.0; 
  fRefl[44] = 0.7071067811865475*(fReflSurfXY[3][5]-1.0*(fReflSurfXY[2][5]+fReflSurfXY[1][5])+fReflSurfXY[0][5]); 
  fRefl[45] = 0.0; 
  fRefl[46] = 0.0; 
  fRefl[47] = 0.0; 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1(double wv, double dv, double vlowerSq, double vupperSq, double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq;
  double fReflSurfXY[4][6] = {0.}; 
  // node (x,y)_0 
  vcutSq = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*phi[6]-2.449489742783178*phiWall[5]+2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*phi[2]-1.414213562373095*phiWall[1]+1.414213562373095*(phi[1]+phiWall[0])-1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[0][0] = 0.0; 
  fReflSurfXY[0][1] = 0.0; 
  fReflSurfXY[0][2] = 0.0; 
  fReflSurfXY[0][3] = 0.0; 
  fReflSurfXY[0][4] = 0.0; 
  fReflSurfXY[0][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[0][0] = 0.3535533905932737*(1.732050807568877*(f[16]-1.0*(f[8]+f[7]))+f[6]+1.732050807568877*f[3]-1.0*(f[2]+f[1])+f[0]); 
  fReflSurfXY[0][1] = 0.3535533905932737*(1.732050807568877*(f[26]-1.0*(f[19]+f[18]))+f[17]+1.732050807568877*f[11]-1.0*(f[10]+f[9])+f[4]); 
  fReflSurfXY[0][2] = 0.3535533905932737*(1.732050807568877*(f[27]-1.0*(f[22]+f[21]))+f[20]+1.732050807568877*f[14]-1.0*(f[13]+f[12])+f[5]); 
  fReflSurfXY[0][3] = 0.3535533905932737*(1.732050807568877*(f[31]-1.0*(f[30]+f[29]))+f[28]+1.732050807568877*f[25]-1.0*(f[24]+f[23])+f[15]); 
  fReflSurfXY[0][4] = 0.02357022603955158*(25.98076211353316*f[43]+5.0*(3.0*f[37]-5.196152422706631*(f[39]+f[38]))+8.660254037844387*(3.0*f[35]-1.732050807568877*(f[34]+f[33]))+15.0*f[32]); 
  fReflSurfXY[0][5] = 0.02357022603955158*(25.98076211353316*f[47]+5.0*(3.0*f[44]-5.196152422706631*(f[46]+f[45]))+8.660254037844387*(3.0*f[42]-1.732050807568877*(f[41]+f[40]))+15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]+11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]+8.0*(f[10]+f[9])-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33])+7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])+6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]-6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33])+3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])+12.0*(f[22]+f[21])-6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])-8.0*f[17]+12.0*f[16]+8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]+8.0*(f[10]+f[9])-12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])+174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]-174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9])-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])-116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]+116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9])+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])+150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]+150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])+120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])-120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]+120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])-58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]+58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9])+58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]+11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]+8.0*(f[10]+f[9])-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33])+7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])+6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]-6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33])+3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*(f[30]+f[29])-8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*(f[24]+f[23])+12.0*(f[22]+f[21])-6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])+8.0*f[17]+12.0*f[16]-8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]-8.0*(f[10]+f[9])-12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])+174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]-174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9])-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])-116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]+116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9])+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])+150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]+150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])-120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])+120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]-120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])-58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]+58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9])+58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])-13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])+8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*(f[10]+f[9])+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33])-7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])+6.928203230275509*f[20]+12.0*(f[16]+f[14])-6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33])-3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])+12.0*(f[22]+f[21])-6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])+8.0*f[17]-12.0*f[16]+8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]-8.0*(f[10]+f[9])+12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])+174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])-174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9])+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])-116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])+116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9])-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])+150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]-150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])+120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])+120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]-120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])-58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])+58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9])-58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*(f[30]+f[29])+11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])-13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1])+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*(f[30]+f[29])-11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1])-13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]-11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*(f[30]+f[29])+8.0*f[28]+13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])+8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*(f[10]+f[9])+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]+7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33])-7.745966692414834*f[32]+12.0*f[27]-12.0*(f[22]+f[21])+6.928203230275509*f[20]+12.0*(f[16]+f[14])-6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*(f[2]+f[1])+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]+3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33])-3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*(f[30]+f[29])-8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])+12.0*(f[22]+f[21])-6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])-8.0*f[17]-12.0*f[16]-8.0*f[15]-12.0*f[14]+6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]+8.0*(f[10]+f[9])+12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]+6.928203230275509*(f[2]+f[1])-6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]-270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*(f[30]+f[29])+174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])-174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9])+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]+300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*f[27]-232.379000772445*(f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*(f[2]+f[1])+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*(f[30]+f[29])-116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])+116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9])-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*(f[46]+f[45])+150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]-150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]-150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*f[27]+232.379000772445*(f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*(f[2]+f[1])-134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])-120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])-120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]+120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*(f[30]+f[29])-58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])+58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9])-58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[0][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[0][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[0][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[0][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[0][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[0][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_1 
  vcutSq = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*(phi[6]+phiWall[5])-2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*(phi[2]+phiWall[1])-1.414213562373095*(phi[1]+phiWall[0])+1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[1][0] = 0.0; 
  fReflSurfXY[1][1] = 0.0; 
  fReflSurfXY[1][2] = 0.0; 
  fReflSurfXY[1][3] = 0.0; 
  fReflSurfXY[1][4] = 0.0; 
  fReflSurfXY[1][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[1][0] = -0.3535533905932737*(1.732050807568877*(f[16]-1.0*f[8]+f[7])+f[6]-1.0*(1.732050807568877*f[3]+f[2])+f[1]-1.0*f[0]); 
  fReflSurfXY[1][1] = -0.3535533905932737*(1.732050807568877*(f[26]-1.0*f[19]+f[18])+f[17]-1.0*(1.732050807568877*f[11]+f[10])+f[9]-1.0*f[4]); 
  fReflSurfXY[1][2] = -0.3535533905932737*(1.732050807568877*(f[27]-1.0*f[22]+f[21])+f[20]-1.0*(1.732050807568877*f[14]+f[13])+f[12]-1.0*f[5]); 
  fReflSurfXY[1][3] = -0.3535533905932737*(1.732050807568877*(f[31]-1.0*f[30]+f[29])+f[28]-1.0*(1.732050807568877*f[25]+f[24])+f[23]-1.0*f[15]); 
  fReflSurfXY[1][4] = -0.02357022603955158*(25.98076211353316*f[43]+5.0*(5.196152422706631*(f[38]-1.0*f[39])+3.0*f[37])+8.660254037844387*(1.732050807568877*(f[33]-1.0*f[34])-3.0*f[35])-15.0*f[32]); 
  fReflSurfXY[1][5] = -0.02357022603955158*(25.98076211353316*f[47]+5.0*(5.196152422706631*(f[45]-1.0*f[46])+3.0*f[44])+8.660254037844387*(1.732050807568877*(f[40]-1.0*f[41])-3.0*f[42])-15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]+11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])-13.41640786499874*f[13]+13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]-11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]+11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]-7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*f[33]-7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]-12.0*(f[16]+f[14])-6.928203230275509*f[13]+6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]-6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]-3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*f[33]-3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]+12.0*f[22]-12.0*f[21]-6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]+12.0*f[16]-8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]-12.0*f[8]+12.0*f[7]+6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]+269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]+174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])-174.2842505793337*f[24]+174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*f[9]+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]-300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]-300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]-116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])+116.1895003862225*f[24]-116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*f[9]-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]+150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]+150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]+120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]-120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]-58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])+58.09475019311126*f[24]-58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*f[9]-58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]+11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])-13.41640786499874*f[13]+13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]-11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]-13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0]); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]+11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]-13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]-7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*f[33]-7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]-12.0*(f[16]+f[14])-6.928203230275509*f[13]+6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]-6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]-3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*f[33]-3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*f[30]-13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]+12.0*f[22]-12.0*f[21]-6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]+12.0*f[16]+8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-12.0*f[8]+12.0*f[7]+6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0]); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]+269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]+174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])-174.2842505793337*f[24]+174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*f[9]+174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]-300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]-300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]-116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])+116.1895003862225*f[24]-116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*f[9]-116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]+150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]+150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]-120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]+120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]-58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])+58.09475019311126*f[24]-58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*f[9]-58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]+11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]-11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]+11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]+11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]-7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*f[33]+7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]-3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*f[33]+3.872983346207417*f[32]+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]+12.0*f[22]-12.0*f[21]-6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]-12.0*f[16]-8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]+12.0*f[8]-12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]+269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]+174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]-174.2842505793337*f[24]+174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*f[9]-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]-300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]+300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]-116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]+116.1895003862225*f[24]-116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*f[9]+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]+150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]-150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])+207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]+120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]+120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]-120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]-58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]+58.09475019311126*f[24]-58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*f[9]+58.09475019311126*f[4]); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*f[31]-20.12461179749811*f[30]+20.12461179749811*f[29]+11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-11.61895003862225*f[4])-25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]+15.0*f[32]+23.2379000772445*f[27]-23.2379000772445*f[22]+23.2379000772445*f[21]+13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*f[1]-13.41640786499874*f[0])-20.12461179749811*f[31]+20.12461179749811*f[30]-20.12461179749811*f[29]-11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]-13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]-11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]+11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*f[1]+13.41640786499874*f[0]); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]-20.1246117974981*f[46]+20.1246117974981*f[45]+11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]-11.61895003862225*f[41]+11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*f[31]-13.85640646055102*f[30]+13.85640646055102*f[29]+8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-8.0*f[4])-13.41640786499874*f[47]+13.41640786499874*f[46]-13.41640786499874*f[45]-7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]+7.745966692414834*f[41]-7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*f[33]+7.745966692414834*f[32]+12.0*f[27]-12.0*f[22]+12.0*f[21]+6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*f[1]-6.928203230275509*f[0])-6.708203932499369*f[47]+6.708203932499369*f[46]-6.708203932499369*f[45]-3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]+3.872983346207417*f[41]-3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*f[33]+3.872983346207417*f[32]-13.85640646055102*f[31]+13.85640646055102*f[30]-13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]+12.0*f[22]-12.0*f[21]-6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]-8.0*f[17]-12.0*f[16]+8.0*f[15]+12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+12.0*f[8]-12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*f[1]+6.928203230275509*f[0]); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]-467.6537180435969*f[46]+467.6537180435969*f[45]+269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]-270.0*f[41]+270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*f[31]-301.8691769624717*f[30]+301.8691769624717*f[29]+174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]-174.2842505793337*f[24]+174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*f[9]-174.2842505793337*f[4])-519.6152422706631*f[47]+519.6152422706633*f[46]-519.6152422706633*f[45]-300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]+300.0*f[41]-300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]+300.0*f[32]+232.379000772445*f[27]-232.379000772445*f[22]+232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*f[1]-134.1640786499874*f[0])-201.2461179749811*f[31]+201.2461179749811*f[30]-201.2461179749811*f[29]-116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]+116.1895003862225*f[24]-116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*f[9]+116.1895003862225*f[4])+259.8076211353315*f[47]-259.8076211353317*f[46]+259.8076211353317*f[45]+150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]-150.0*f[41]+150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]-150.0*f[32]-232.379000772445*f[27]+232.379000772445*f[22]-232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*f[1]+134.1640786499874*f[0])-207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]-120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]-120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]+120.0*f[32]-100.6230589874906*f[31]+100.6230589874906*f[30]-100.6230589874906*f[29]-58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]+58.09475019311126*f[24]-58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*f[9]+58.09475019311126*f[4]); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[1][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[1][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[1][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[1][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[1][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[1][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_2 
  vcutSq = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*(phi[6]+phiWall[5])+2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*(phi[2]+phiWall[1])+1.414213562373095*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[2][0] = 0.0; 
  fReflSurfXY[2][1] = 0.0; 
  fReflSurfXY[2][2] = 0.0; 
  fReflSurfXY[2][3] = 0.0; 
  fReflSurfXY[2][4] = 0.0; 
  fReflSurfXY[2][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[2][0] = -0.3535533905932737*(1.732050807568877*(f[16]+f[8]-1.0*f[7])+f[6]-1.732050807568877*f[3]+f[2]-1.0*(f[1]+f[0])); 
  fReflSurfXY[2][1] = -0.3535533905932737*(1.732050807568877*(f[26]+f[19]-1.0*f[18])+f[17]-1.732050807568877*f[11]+f[10]-1.0*(f[9]+f[4])); 
  fReflSurfXY[2][2] = -0.3535533905932737*(1.732050807568877*(f[27]+f[22]-1.0*f[21])+f[20]-1.732050807568877*f[14]+f[13]-1.0*(f[12]+f[5])); 
  fReflSurfXY[2][3] = -0.3535533905932737*(1.732050807568877*(f[31]+f[30]-1.0*f[29])+f[28]-1.732050807568877*f[25]+f[24]-1.0*(f[23]+f[15])); 
  fReflSurfXY[2][4] = -0.02357022603955158*(25.98076211353316*f[43]+5.0*(5.196152422706631*(f[39]-1.0*f[38])+3.0*f[37])+8.660254037844387*(1.732050807568877*f[34]-3.0*f[35])-1.0*(15.0*f[33]+15.0*f[32])); 
  fReflSurfXY[2][5] = -0.02357022603955158*(25.98076211353316*f[47]+5.0*(5.196152422706631*(f[46]-1.0*f[45])+3.0*f[44])+8.660254037844387*(1.732050807568877*f[41]-3.0*f[42])-1.0*(15.0*f[40]+15.0*f[36])); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]+11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]+13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])+13.41640786499874*f[13]-13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]-11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]+13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]+11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*f[10]+8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]-7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]+6.928203230275509*f[20]-12.0*(f[16]+f[14])+6.928203230275509*f[13]-6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]-6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]-3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*(f[33]+f[32])+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]+12.0*f[16]-8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]+12.0*f[8]-12.0*f[7]+6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]+269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]+174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])+174.2842505793337*f[24]-174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]-300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]-300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]-116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])-116.1895003862225*f[24]+116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]+150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]+150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]+120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]-120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]+120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]-58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])-58.09475019311126*f[24]+58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*(f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]+11.61895003862225*f[28]-20.12461179749811*(f[26]+f[25])+11.61895003862225*f[24]-11.61895003862225*f[23]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]+13.41640786499874*f[20]-23.2379000772445*(f[16]+f[14])+13.41640786499874*f[13]-13.41640786499874*f[12]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]-11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*(f[26]+f[25])-11.61895003862225*f[24]+11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]+11.61895003862225*f[44]-20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]-20.12461179749811*f[39]+20.12461179749811*f[38]-11.61895003862225*f[37]-11.61895003862225*f[36]+20.1246117974981*f[35]-11.61895003862225*f[34]+11.61895003862225*f[33]+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]-13.85640646055102*(f[26]+f[25])+8.0*f[24]-8.0*f[23]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*(f[17]+f[15])+13.85640646055102*f[11]-8.0*f[10]+8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]-7.745966692414834*f[44]+13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]+13.41640786499874*f[39]-13.41640786499874*f[38]+7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]+7.745966692414834*f[34]-7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]+6.928203230275509*f[20]-12.0*(f[16]+f[14])+6.928203230275509*f[13]-6.928203230275509*f[12]-12.0*f[8]+12.0*f[7]-6.928203230275509*(f[6]+f[5])+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]-3.872983346207417*f[44]+6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]+6.708203932499369*f[39]-6.708203932499369*f[38]+3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]+3.872983346207417*f[34]-3.872983346207417*(f[33]+f[32])-13.85640646055102*(f[31]+f[30])+13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])-8.0*f[24]+8.0*f[23]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]+12.0*f[16]+8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]+12.0*f[8]-12.0*f[7]+6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]+269.9999999999999*f[44]-467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]-467.6537180435967*f[39]+467.6537180435967*f[38]-270.0*f[37]-269.9999999999999*f[36]+467.6537180435969*f[35]-269.9999999999999*f[34]+269.9999999999999*f[33]+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]+174.2842505793337*f[28]-301.8691769624717*(f[26]+f[25])+174.2842505793337*f[24]-174.2842505793337*f[23]-301.8691769624717*f[19]+301.8691769624717*f[18]-174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]-174.2842505793337*f[10]+174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]-300.0000000000001*f[44]+519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]+519.6152422706631*f[39]-519.6152422706631*f[38]+300.0*f[37]+300.0000000000001*f[36]-519.6152422706633*f[35]+300.0000000000001*f[34]-300.0000000000001*f[33]-300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]+134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])+134.1640786499874*f[13]-134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]-116.1895003862225*f[28]+201.2461179749811*(f[26]+f[25])-116.1895003862225*f[24]+116.1895003862225*f[23]+201.2461179749811*f[19]-201.2461179749811*f[18]+116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]+116.1895003862225*f[10]-116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]+150.0*f[44]-259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]-259.8076211353315*f[39]+259.8076211353315*f[38]-150.0*f[37]-150.0*f[36]+259.8076211353317*f[35]-150.0*f[34]+150.0*f[33]+150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]-134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])-134.1640786499874*f[13]+134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]-120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]+120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]-120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]-58.09475019311126*f[28]+100.6230589874906*(f[26]+f[25])-58.09475019311126*f[24]+58.09475019311126*f[23]+100.6230589874906*f[19]-100.6230589874906*f[18]+58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]+58.09475019311126*f[10]-58.09475019311126*(f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]+11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]+13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]-11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]+23.2379000772445*f[22]-23.2379000772445*f[21]+13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]+23.2379000772445*f[16]+11.61895003862225*f[15]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]+11.61895003862225*f[4]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]+11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]+8.0*f[10]-8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]-7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]+6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]-3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*(f[33]+f[32])+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]-12.0*f[16]-8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]-13.85640646055102*f[11]+8.0*f[10]-8.0*f[9]-12.0*f[8]+12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]+269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]+174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]+174.2842505793337*f[24]-174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]-300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]+300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]-116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]-116.1895003862225*f[24]+116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]+150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]-150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*f[46]-207.8460969082653*f[45]+120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]+120.0*f[41]-120.0*f[40]+207.8460969082653*f[39]-207.8460969082653*f[38]+120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]+120.0*f[34]-120.0*f[33]-120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]-58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]-58.09475019311126*f[24]+58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*(f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*f[46]-25.98076211353316*f[45]+15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]+15.0*f[41]-15.0*f[40]+25.98076211353316*f[39]-25.98076211353316*f[38]+15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]+15.0*f[34]-15.0*f[33]-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30])-20.12461179749811*f[29]+11.61895003862225*f[28]+20.12461179749811*f[26]-20.12461179749811*f[25]+11.61895003862225*f[24]-11.61895003862225*f[23]+20.12461179749811*f[19]-20.12461179749811*f[18]+11.61895003862225*f[17]-11.61895003862225*f[15]-20.12461179749811*f[11]+11.61895003862225*f[10]-11.61895003862225*(f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*f[46]+25.98076211353316*f[45]-15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]-15.0*f[41]+15.0*f[40]-25.98076211353316*f[39]+25.98076211353316*f[38]-15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]-15.0*f[34]+15.0*f[33]+15.0*f[32]+23.2379000772445*(f[27]+f[22])-23.2379000772445*f[21]+13.41640786499874*f[20]+23.2379000772445*f[16]-23.2379000772445*f[14]+13.41640786499874*f[13]-13.41640786499874*f[12]+23.2379000772445*f[8]-23.2379000772445*f[7]+13.41640786499874*f[6]-13.41640786499874*f[5]-23.2379000772445*f[3]+13.41640786499874*f[2]-13.41640786499874*(f[1]+f[0]))-20.12461179749811*(f[31]+f[30])+20.12461179749811*f[29]-11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*f[26]+20.12461179749811*f[25]-11.61895003862225*f[24]+11.61895003862225*f[23]-23.2379000772445*f[22]+23.2379000772445*f[21]-13.41640786499874*f[20]-20.12461179749811*f[19]+20.12461179749811*f[18]-11.61895003862225*f[17]-23.2379000772445*f[16]+11.61895003862225*f[15]+23.2379000772445*f[14]-13.41640786499874*f[13]+13.41640786499874*f[12]+20.12461179749811*f[11]-11.61895003862225*f[10]+11.61895003862225*f[9]-23.2379000772445*f[8]+23.2379000772445*f[7]-13.41640786499874*f[6]+13.41640786499874*f[5]+11.61895003862225*f[4]+23.2379000772445*f[3]-13.41640786499874*f[2]+13.41640786499874*(f[1]+f[0])); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*f[46]-20.1246117974981*f[45]+11.61895003862225*f[44]+20.1246117974981*f[43]-20.12461179749811*f[42]+11.61895003862225*f[41]-11.61895003862225*f[40]+20.12461179749811*f[39]-20.12461179749811*f[38]+11.61895003862225*f[37]-11.61895003862225*f[36]-20.1246117974981*f[35]+11.61895003862225*f[34]-11.61895003862225*f[33]-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30])-13.85640646055102*f[29]+8.0*f[28]+13.85640646055102*f[26]-13.85640646055102*f[25]+8.0*f[24]-8.0*f[23]+13.85640646055102*f[19]-13.85640646055102*f[18]+8.0*f[17]-8.0*f[15]-13.85640646055102*f[11]+8.0*f[10]-8.0*(f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*f[46]+13.41640786499874*f[45]-7.745966692414834*f[44]-13.41640786499874*f[43]+13.41640786499874*f[42]-7.745966692414834*f[41]+7.745966692414834*f[40]-13.41640786499874*f[39]+13.41640786499874*f[38]-7.745966692414834*f[37]+7.745966692414834*f[36]+13.41640786499874*f[35]-7.745966692414834*f[34]+7.745966692414834*(f[33]+f[32])+12.0*(f[27]+f[22])-12.0*f[21]+6.928203230275509*f[20]+12.0*f[16]-12.0*f[14]+6.928203230275509*f[13]-6.928203230275509*f[12]+12.0*f[8]-12.0*f[7]+6.928203230275509*f[6]-6.928203230275509*f[5]-12.0*f[3]+6.928203230275509*f[2]-6.928203230275509*(f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*f[46]+6.708203932499369*f[45]-3.872983346207417*f[44]-6.708203932499369*f[43]+6.708203932499369*f[42]-3.872983346207417*f[41]+3.872983346207417*f[40]-6.708203932499369*f[39]+6.708203932499369*f[38]-3.872983346207417*f[37]+3.872983346207417*f[36]+6.708203932499369*f[35]-3.872983346207417*f[34]+3.872983346207417*(f[33]+f[32])-13.85640646055102*(f[31]+f[30])+13.85640646055102*f[29]-8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]-8.0*f[24]+8.0*f[23]-12.0*f[22]+12.0*f[21]-6.928203230275509*f[20]-13.85640646055102*f[19]+13.85640646055102*f[18]-8.0*f[17]-12.0*f[16]+8.0*f[15]+12.0*f[14]-6.928203230275509*f[13]+6.928203230275509*f[12]+13.85640646055102*f[11]-8.0*f[10]+8.0*f[9]-12.0*f[8]+12.0*f[7]-6.928203230275509*f[6]+6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]-6.928203230275509*f[2]+6.928203230275509*(f[1]+f[0])); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*f[46]-467.6537180435969*f[45]+269.9999999999999*f[44]+467.6537180435969*f[43]-467.6537180435967*f[42]+270.0*f[41]-270.0*f[40]+467.6537180435967*f[39]-467.6537180435967*f[38]+270.0*f[37]-269.9999999999999*f[36]-467.6537180435969*f[35]+269.9999999999999*f[34]-269.9999999999999*f[33]-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30])-301.8691769624717*f[29]+174.2842505793337*f[28]+301.8691769624717*f[26]-301.8691769624717*f[25]+174.2842505793337*f[24]-174.2842505793337*f[23]+301.8691769624717*f[19]-301.8691769624717*f[18]+174.2842505793337*f[17]-174.2842505793337*f[15]-301.8691769624717*f[11]+174.2842505793337*f[10]-174.2842505793337*(f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*f[46]+519.6152422706633*f[45]-300.0000000000001*f[44]-519.6152422706633*f[43]+519.6152422706631*f[42]-300.0*f[41]+300.0*f[40]-519.6152422706631*f[39]+519.6152422706631*f[38]-300.0*f[37]+300.0000000000001*f[36]+519.6152422706633*f[35]-300.0000000000001*f[34]+300.0000000000001*f[33]+300.0*f[32]+232.379000772445*(f[27]+f[22])-232.379000772445*f[21]+134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]+134.1640786499874*f[13]-134.1640786499874*f[12]+232.379000772445*f[8]-232.379000772445*f[7]+134.1640786499874*f[6]-134.1640786499874*f[5]-232.379000772445*f[3]+134.1640786499874*f[2]-134.1640786499874*(f[1]+f[0]))-201.2461179749811*(f[31]+f[30])+201.2461179749811*f[29]-116.1895003862225*f[28]-201.2461179749811*f[26]+201.2461179749811*f[25]-116.1895003862225*f[24]+116.1895003862225*f[23]-201.2461179749811*f[19]+201.2461179749811*f[18]-116.1895003862225*f[17]+116.1895003862225*f[15]+201.2461179749811*f[11]-116.1895003862225*f[10]+116.1895003862225*(f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*f[46]-259.8076211353317*f[45]+150.0*f[44]+259.8076211353317*f[43]-259.8076211353315*f[42]+150.0*f[41]-150.0*f[40]+259.8076211353315*f[39]-259.8076211353315*f[38]+150.0*f[37]-150.0*f[36]-259.8076211353317*f[35]+150.0*f[34]-150.0*f[33]-150.0*f[32]-232.379000772445*(f[27]+f[22])+232.379000772445*f[21]-134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]-134.1640786499874*f[13]+134.1640786499874*f[12]-232.379000772445*f[8]+232.379000772445*f[7]-134.1640786499874*f[6]+134.1640786499874*f[5]+232.379000772445*f[3]-134.1640786499874*f[2]+134.1640786499874*(f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*f[46]+207.8460969082653*f[45]-120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]-120.0*f[41]+120.0*f[40]-207.8460969082653*f[39]+207.8460969082653*f[38]-120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]-120.0*f[34]+120.0*f[33]+120.0*f[32]-100.6230589874906*(f[31]+f[30])+100.6230589874906*f[29]-58.09475019311126*f[28]-100.6230589874906*f[26]+100.6230589874906*f[25]-58.09475019311126*f[24]+58.09475019311126*f[23]-100.6230589874906*f[19]+100.6230589874906*f[18]-58.09475019311126*f[17]+58.09475019311126*f[15]+100.6230589874906*f[11]-58.09475019311126*f[10]+58.09475019311126*(f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[2][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[2][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[2][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[2][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[2][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[2][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  // node (x,y)_3 
  vcutSq = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*phi[6]+2.449489742783178*phiWall[5]-2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*phi[2]+1.414213562373095*phiWall[1]-1.414213562373095*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm;

  if (vcutSq <= vlowerSq) { // absorb (no reflection)

  fReflSurfXY[3][0] = 0.0; 
  fReflSurfXY[3][1] = 0.0; 
  fReflSurfXY[3][2] = 0.0; 
  fReflSurfXY[3][3] = 0.0; 
  fReflSurfXY[3][4] = 0.0; 
  fReflSurfXY[3][5] = 0.0; 

  } else if (vcutSq > vupperSq) { // full reflection

  fReflSurfXY[3][0] = 0.3535533905932737*(1.732050807568877*(f[16]+f[8]+f[7])+f[6]+1.732050807568877*f[3]+f[2]+f[1]+f[0]); 
  fReflSurfXY[3][1] = 0.3535533905932737*(1.732050807568877*(f[26]+f[19]+f[18])+f[17]+1.732050807568877*f[11]+f[10]+f[9]+f[4]); 
  fReflSurfXY[3][2] = 0.3535533905932737*(1.732050807568877*(f[27]+f[22]+f[21])+f[20]+1.732050807568877*f[14]+f[13]+f[12]+f[5]); 
  fReflSurfXY[3][3] = 0.3535533905932737*(1.732050807568877*(f[31]+f[30]+f[29])+f[28]+1.732050807568877*f[25]+f[24]+f[23]+f[15]); 
  fReflSurfXY[3][4] = 0.02357022603955158*(25.98076211353316*f[43]+5.0*(5.196152422706631*(f[39]+f[38])+3.0*f[37])+8.660254037844387*(3.0*f[35]+1.732050807568877*(f[34]+f[33]))+15.0*f[32]); 
  fReflSurfXY[3][5] = 0.02357022603955158*(25.98076211353316*f[47]+5.0*(5.196152422706631*(f[46]+f[45])+3.0*f[44])+8.660254037844387*(3.0*f[42]+1.732050807568877*(f[41]+f[40]))+15.0*f[36]); 

  } else { // partial reflection

    double xBar, xSqBar;
    double fReflSurfXYMu[2][6] = {0.}; 

    // node (mu)_0 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])+13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]+23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]+11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]-8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])+6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]+6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33]+f[32])+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-12.0*f[27]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])-12.0*(f[22]+f[21])-6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])-8.0*f[17]+12.0*f[16]+8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]-8.0*(f[10]+f[9])+12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]-8.0*f[4]+12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])+174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]+174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])-116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]-116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])+150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]-150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])+120.0*f[44]-207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])-120.0*f[37]+120.0*f[36]-207.8460969082653*f[35]-120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])-58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]-58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[0][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]-25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]+15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]-20.12461179749811*f[26]+20.12461179749811*f[25]+11.61895003862225*(f[24]+f[23])-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+11.61895003862225*f[15]-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]+25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]-15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])+13.41640786499874*f[20]-23.2379000772445*f[16]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-23.2379000772445*(f[8]+f[7])-13.41640786499874*f[6]+13.41640786499874*f[5]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]-23.2379000772445*f[27]+20.12461179749811*f[26]-20.12461179749811*f[25]-11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]+20.12461179749811*(f[19]+f[18])+11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])+13.41640786499874*f[6]-13.41640786499874*f[5]+11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]-20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])-20.12461179749811*(f[39]+f[38])-11.61895003862225*f[37]+11.61895003862225*f[36]-20.1246117974981*f[35]-11.61895003862225*(f[34]+f[33])-11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-13.85640646055102*f[26]+13.85640646055102*f[25]+8.0*(f[24]+f[23])-13.85640646055102*(f[19]+f[18])-8.0*f[17]+8.0*f[15]-13.85640646055102*f[11]-8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]+13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])+13.41640786499874*(f[39]+f[38])+7.745966692414834*f[37]-7.745966692414834*f[36]+13.41640786499874*f[35]+7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])+6.928203230275509*f[20]-12.0*f[16]+12.0*f[14]+6.928203230275509*(f[13]+f[12])-12.0*(f[8]+f[7])-6.928203230275509*f[6]+6.928203230275509*f[5]-12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]+6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])+6.708203932499369*(f[39]+f[38])+3.872983346207417*f[37]-3.872983346207417*f[36]+6.708203932499369*f[35]+3.872983346207417*(f[34]+f[33]+f[32])-13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-12.0*f[27]+13.85640646055102*f[26]-13.85640646055102*f[25]-8.0*(f[24]+f[23])-12.0*(f[22]+f[21])-6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])+8.0*f[17]+12.0*f[16]-8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]+8.0*(f[10]+f[9])+12.0*(f[8]+f[7])+6.928203230275509*f[6]-6.928203230275509*f[5]+8.0*f[4]+12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[0][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]-467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])-467.6537180435967*(f[39]+f[38])-270.0*f[37]+269.9999999999999*f[36]-467.6537180435969*f[35]-269.9999999999999*(f[34]+f[33])-270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])+174.2842505793337*f[28]-301.8691769624717*f[26]+301.8691769624717*f[25]+174.2842505793337*(f[24]+f[23])-301.8691769624717*(f[19]+f[18])-174.2842505793337*f[17]+174.2842505793337*f[15]-301.8691769624717*f[11]-174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]+519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])+519.6152422706631*(f[39]+f[38])+300.0*f[37]-300.0000000000001*f[36]+519.6152422706633*f[35]+300.0000000000001*(f[34]+f[33])+300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]-232.379000772445*f[16]+232.379000772445*f[14]+134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*f[6]+134.1640786499874*f[5]-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])-116.1895003862225*f[28]+201.2461179749811*f[26]-201.2461179749811*f[25]-116.1895003862225*(f[24]+f[23])+201.2461179749811*(f[19]+f[18])+116.1895003862225*f[17]-116.1895003862225*f[15]+201.2461179749811*f[11]+116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])+150.0*f[44]-259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])-259.8076211353315*(f[39]+f[38])-150.0*f[37]+150.0*f[36]-259.8076211353317*f[35]-150.0*(f[34]+f[33])-150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]+232.379000772445*f[16]-232.379000772445*f[14]-134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*f[6]-134.1640786499874*f[5]+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])-120.0*f[44]+207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])+120.0*f[37]-120.0*f[36]+207.8460969082653*f[35]+120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])-58.09475019311126*f[28]+100.6230589874906*f[26]-100.6230589874906*f[25]-58.09475019311126*(f[24]+f[23])+100.6230589874906*(f[19]+f[18])+58.09475019311126*f[17]-58.09475019311126*f[15]+100.6230589874906*f[11]+58.09475019311126*(f[10]+f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[0][1])/fReflSurfXYMu[0][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[0][2]+5.0*fReflSurfXYMu[0][0]))/fReflSurfXYMu[0][0];
  if (fReflSurfXYMu[0][0]<0.) {
  fReflSurfXYMu[0][0] = 0.0; 
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[0][1] = 0.0; 
  fReflSurfXYMu[0][2] = 0.0; 
  }

    // node (mu)_1 
    if (wv > 0.) {
      // vcut in logical space.
      double vcut_l = 2.*(sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = 0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])+13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])+13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]+23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])+23.2379000772445*(f[22]+f[21])+13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]+23.2379000772445*f[16]-11.61895003862225*f[15]+23.2379000772445*f[14]+13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])+23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][1] = 0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])+8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])+6.928203230275509*f[20]+12.0*(f[16]+f[14])+6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33]+f[32])+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]-12.0*f[27]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])-12.0*(f[22]+f[21])-6.928203230275509*f[20]+13.85640646055102*(f[19]+f[18])+8.0*f[17]-12.0*f[16]+8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])+13.85640646055102*f[11]+8.0*(f[10]+f[9])-12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])+8.0*f[4]-12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][2] = 0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])+174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])+174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])-116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])-116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])+150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]+150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))+207.8460969082653*f[47]+207.8460969082653*(f[46]+f[45])+120.0*f[44]+207.8460969082653*f[43]+207.8460969082653*f[42]+120.0*(f[41]+f[40])+207.8460969082653*(f[39]+f[38])+120.0*f[37]+120.0*f[36]+207.8460969082653*f[35]+120.0*(f[34]+f[33])+120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])-58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])-58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9]+f[4])); 

    } else {
      // vcut in logical space.
      double vcut_l = 2.*(-sqrt(vcutSq)-wv)/dv;

  fReflSurfXYMu[1][0] = -0.009316949906249122*(vcut_l*(vcut_l*((25.98076211353316*f[47]+25.98076211353316*(f[46]+f[45])+15.0*f[44]+25.98076211353316*f[43]+25.98076211353316*f[42]+15.0*(f[41]+f[40])+25.98076211353316*(f[39]+f[38])+15.0*f[37]+15.0*f[36]+25.98076211353316*f[35]+15.0*(f[34]+f[33])+15.0*f[32])*vcut_l+20.12461179749811*(f[31]+f[30]+f[29])+11.61895003862225*f[28]+20.12461179749811*(f[26]+f[25])+11.61895003862225*(f[24]+f[23])+20.12461179749811*(f[19]+f[18])+11.61895003862225*(f[17]+f[15])+20.12461179749811*f[11]+11.61895003862225*(f[10]+f[9]+f[4]))-25.98076211353316*f[47]-25.98076211353316*(f[46]+f[45])-15.0*f[44]-25.98076211353316*f[43]-25.98076211353316*f[42]-15.0*(f[41]+f[40])-25.98076211353316*(f[39]+f[38])-15.0*f[37]-15.0*f[36]-25.98076211353316*f[35]-15.0*(f[34]+f[33])-15.0*f[32]+23.2379000772445*(f[27]+f[22]+f[21])+13.41640786499874*f[20]+23.2379000772445*(f[16]+f[14])+13.41640786499874*(f[13]+f[12])+23.2379000772445*(f[8]+f[7])+13.41640786499874*(f[6]+f[5])+23.2379000772445*f[3]+13.41640786499874*(f[2]+f[1]+f[0]))-20.12461179749811*(f[31]+f[30]+f[29])-11.61895003862225*f[28]-23.2379000772445*f[27]-20.12461179749811*(f[26]+f[25])-11.61895003862225*(f[24]+f[23])-23.2379000772445*(f[22]+f[21])-13.41640786499874*f[20]-20.12461179749811*(f[19]+f[18])-11.61895003862225*f[17]-23.2379000772445*f[16]-11.61895003862225*f[15]-23.2379000772445*f[14]-13.41640786499874*(f[13]+f[12])-20.12461179749811*f[11]-11.61895003862225*(f[10]+f[9])-23.2379000772445*(f[8]+f[7])-13.41640786499874*(f[6]+f[5])-11.61895003862225*f[4]-23.2379000772445*f[3]-13.41640786499874*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][1] = -0.015625*(vcut_l*vcut_l*(vcut_l*((20.12461179749811*f[47]+20.1246117974981*(f[46]+f[45])+11.61895003862225*f[44]+20.1246117974981*f[43]+20.12461179749811*f[42]+11.61895003862225*(f[41]+f[40])+20.12461179749811*(f[39]+f[38])+11.61895003862225*f[37]+11.61895003862225*f[36]+20.1246117974981*f[35]+11.61895003862225*(f[34]+f[33])+11.61895003862225*f[32])*vcut_l+13.85640646055102*(f[31]+f[30]+f[29])+8.0*f[28]+13.85640646055102*(f[26]+f[25])+8.0*(f[24]+f[23])+13.85640646055102*(f[19]+f[18])+8.0*(f[17]+f[15])+13.85640646055102*f[11]+8.0*(f[10]+f[9]+f[4]))-13.41640786499874*f[47]-13.41640786499874*(f[46]+f[45])-7.745966692414834*f[44]-13.41640786499874*f[43]-13.41640786499874*f[42]-7.745966692414834*(f[41]+f[40])-13.41640786499874*(f[39]+f[38])-7.745966692414834*(f[37]+f[36])-13.41640786499874*f[35]-7.745966692414834*(f[34]+f[33]+f[32])+12.0*(f[27]+f[22]+f[21])+6.928203230275509*f[20]+12.0*(f[16]+f[14])+6.928203230275509*(f[13]+f[12])+12.0*(f[8]+f[7])+6.928203230275509*(f[6]+f[5])+12.0*f[3]+6.928203230275509*(f[2]+f[1]+f[0]))-6.708203932499369*f[47]-6.708203932499369*(f[46]+f[45])-3.872983346207417*f[44]-6.708203932499369*f[43]-6.708203932499369*f[42]-3.872983346207417*(f[41]+f[40])-6.708203932499369*(f[39]+f[38])-3.872983346207417*(f[37]+f[36])-6.708203932499369*f[35]-3.872983346207417*(f[34]+f[33]+f[32])-13.85640646055102*(f[31]+f[30]+f[29])-8.0*f[28]-12.0*f[27]-13.85640646055102*(f[26]+f[25])-8.0*(f[24]+f[23])-12.0*(f[22]+f[21])-6.928203230275509*f[20]-13.85640646055102*(f[19]+f[18])-8.0*f[17]-12.0*f[16]-8.0*f[15]-12.0*f[14]-6.928203230275509*(f[13]+f[12])-13.85640646055102*f[11]-8.0*(f[10]+f[9])-12.0*(f[8]+f[7])-6.928203230275509*(f[6]+f[5])-8.0*f[4]-12.0*f[3]-6.928203230275509*(f[2]+f[1]+f[0])); 
  fReflSurfXYMu[1][2] = -0.001041666666666667*(vcut_l*(vcut_l*(vcut_l*(vcut_l*((467.6537180435967*f[47]+467.6537180435969*(f[46]+f[45])+269.9999999999999*f[44]+467.6537180435969*f[43]+467.6537180435967*f[42]+270.0*(f[41]+f[40])+467.6537180435967*(f[39]+f[38])+270.0*f[37]+269.9999999999999*f[36]+467.6537180435969*f[35]+269.9999999999999*(f[34]+f[33])+270.0*f[32])*vcut_l+301.8691769624717*(f[31]+f[30]+f[29])+174.2842505793337*f[28]+301.8691769624717*(f[26]+f[25])+174.2842505793337*(f[24]+f[23])+301.8691769624717*(f[19]+f[18])+174.2842505793337*(f[17]+f[15])+301.8691769624717*f[11]+174.2842505793337*(f[10]+f[9]+f[4]))-519.6152422706631*f[47]-519.6152422706633*(f[46]+f[45])-300.0000000000001*f[44]-519.6152422706633*f[43]-519.6152422706631*f[42]-300.0*(f[41]+f[40])-519.6152422706631*(f[39]+f[38])-300.0*f[37]-300.0000000000001*f[36]-519.6152422706633*f[35]-300.0000000000001*(f[34]+f[33])-300.0*f[32]+232.379000772445*(f[27]+f[22]+f[21])+134.1640786499874*f[20]+232.379000772445*(f[16]+f[14])+134.1640786499874*(f[13]+f[12])+232.379000772445*(f[8]+f[7])+134.1640786499874*(f[6]+f[5])+232.379000772445*f[3]+134.1640786499874*(f[2]+f[1]+f[0]))-201.2461179749811*(f[31]+f[30]+f[29])-116.1895003862225*f[28]-201.2461179749811*(f[26]+f[25])-116.1895003862225*(f[24]+f[23])-201.2461179749811*(f[19]+f[18])-116.1895003862225*(f[17]+f[15])-201.2461179749811*f[11]-116.1895003862225*(f[10]+f[9]+f[4]))+259.8076211353315*f[47]+259.8076211353317*(f[46]+f[45])+150.0*f[44]+259.8076211353317*f[43]+259.8076211353315*f[42]+150.0*(f[41]+f[40])+259.8076211353315*(f[39]+f[38])+150.0*f[37]+150.0*f[36]+259.8076211353317*f[35]+150.0*(f[34]+f[33])+150.0*f[32]-232.379000772445*(f[27]+f[22]+f[21])-134.1640786499874*f[20]-232.379000772445*(f[16]+f[14])-134.1640786499874*(f[13]+f[12])-232.379000772445*(f[8]+f[7])-134.1640786499874*(f[6]+f[5])-232.379000772445*f[3]-134.1640786499874*(f[2]+f[1]+f[0]))-207.8460969082653*f[47]-207.8460969082653*(f[46]+f[45])-120.0*f[44]-207.8460969082653*f[43]-207.8460969082653*f[42]-120.0*(f[41]+f[40])-207.8460969082653*(f[39]+f[38])-120.0*f[37]-120.0*f[36]-207.8460969082653*f[35]-120.0*(f[34]+f[33])-120.0*f[32]-100.6230589874906*(f[31]+f[30]+f[29])-58.09475019311126*f[28]-100.6230589874906*(f[26]+f[25])-58.09475019311126*(f[24]+f[23])-100.6230589874906*(f[19]+f[18])-58.09475019311126*(f[17]+f[15])-100.6230589874906*f[11]-58.09475019311126*(f[10]+f[9]+f[4])); 

    }

  // If the cut distribution f(vpar), where vpar in [-1,vcut_l] or [vcut_l/1],
  // has a cell average < 0, set to 0. If it's >0 but not realizable, set to p=0.
  xBar = (0.5773502691896258*fReflSurfXYMu[1][1])/fReflSurfXYMu[1][0];
  xSqBar = (0.06666666666666667*(4.47213595499958*fReflSurfXYMu[1][2]+5.0*fReflSurfXYMu[1][0]))/fReflSurfXYMu[1][0];
  if (fReflSurfXYMu[1][0]<0.) {
  fReflSurfXYMu[1][0] = 0.0; 
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  } else if (fabs(xBar)>=1. || fabs(xSqBar)>=1.) {
  fReflSurfXYMu[1][1] = 0.0; 
  fReflSurfXYMu[1][2] = 0.0; 
  }

  fReflSurfXY[3][0] = 0.7071067811865475*(fReflSurfXYMu[1][0]+fReflSurfXYMu[0][0]); 
  fReflSurfXY[3][1] = 0.7071067811865475*(fReflSurfXYMu[1][1]+fReflSurfXYMu[0][1]); 
  fReflSurfXY[3][2] = 0.7071067811865475*(fReflSurfXYMu[1][0]-1.0*fReflSurfXYMu[0][0]); 
  fReflSurfXY[3][3] = 0.7071067811865475*(fReflSurfXYMu[1][1]-1.0*fReflSurfXYMu[0][1]); 
  fReflSurfXY[3][4] = 0.7071067811865475*(fReflSurfXYMu[1][2]+fReflSurfXYMu[0][2]); 
  fReflSurfXY[3][5] = 0.7071067811865475*(fReflSurfXYMu[1][2]-1.0*fReflSurfXYMu[0][2]); 

  }

  fRefl[0] = 0.7071067811865475*(fReflSurfXY[3][0]+fReflSurfXY[2][0]+fReflSurfXY[1][0]+fReflSurfXY[0][0]); 
  fRefl[1] = 0.7071067811865475*(fReflSurfXY[3][0]+fReflSurfXY[2][0]-1.0*(fReflSurfXY[1][0]+fReflSurfXY[0][0])); 
  fRefl[2] = 0.7071067811865475*(fReflSurfXY[3][0]-1.0*fReflSurfXY[2][0]+fReflSurfXY[1][0]-1.0*fReflSurfXY[0][0]); 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.7071067811865475*(fReflSurfXY[3][1]+fReflSurfXY[2][1]+fReflSurfXY[1][1]+fReflSurfXY[0][1]); 
  fRefl[5] = 0.7071067811865475*(fReflSurfXY[3][2]+fReflSurfXY[2][2]+fReflSurfXY[1][2]+fReflSurfXY[0][2]); 
  fRefl[6] = 0.7071067811865475*(fReflSurfXY[3][0]-1.0*(fReflSurfXY[2][0]+fReflSurfXY[1][0])+fReflSurfXY[0][0]); 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.7071067811865475*(fReflSurfXY[3][1]+fReflSurfXY[2][1]-1.0*(fReflSurfXY[1][1]+fReflSurfXY[0][1])); 
  fRefl[10] = 0.7071067811865475*(fReflSurfXY[3][1]-1.0*fReflSurfXY[2][1]+fReflSurfXY[1][1]-1.0*fReflSurfXY[0][1]); 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.7071067811865475*(fReflSurfXY[3][2]+fReflSurfXY[2][2]-1.0*(fReflSurfXY[1][2]+fReflSurfXY[0][2])); 
  fRefl[13] = 0.7071067811865475*(fReflSurfXY[3][2]-1.0*fReflSurfXY[2][2]+fReflSurfXY[1][2]-1.0*fReflSurfXY[0][2]); 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.7071067811865475*(fReflSurfXY[3][3]+fReflSurfXY[2][3]+fReflSurfXY[1][3]+fReflSurfXY[0][3]); 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.7071067811865475*(fReflSurfXY[3][1]-1.0*(fReflSurfXY[2][1]+fReflSurfXY[1][1])+fReflSurfXY[0][1]); 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  fRefl[20] = 0.7071067811865475*(fReflSurfXY[3][2]-1.0*(fReflSurfXY[2][2]+fReflSurfXY[1][2])+fReflSurfXY[0][2]); 
  fRefl[21] = 0.0; 
  fRefl[22] = 0.0; 
  fRefl[23] = 0.7071067811865475*(fReflSurfXY[3][3]+fReflSurfXY[2][3]-1.0*(fReflSurfXY[1][3]+fReflSurfXY[0][3])); 
  fRefl[24] = 0.7071067811865475*(fReflSurfXY[3][3]-1.0*fReflSurfXY[2][3]+fReflSurfXY[1][3]-1.0*fReflSurfXY[0][3]); 
  fRefl[25] = 0.0; 
  fRefl[26] = 0.0; 
  fRefl[27] = 0.0; 
  fRefl[28] = 0.7071067811865475*(fReflSurfXY[3][3]-1.0*(fReflSurfXY[2][3]+fReflSurfXY[1][3])+fReflSurfXY[0][3]); 
  fRefl[29] = 0.0; 
  fRefl[30] = 0.0; 
  fRefl[31] = 0.0; 
  fRefl[32] = 0.7071067811865475*(fReflSurfXY[3][4]+fReflSurfXY[2][4]+fReflSurfXY[1][4]+fReflSurfXY[0][4]); 
  fRefl[33] = 0.7071067811865475*(fReflSurfXY[3][4]+fReflSurfXY[2][4]-1.0*(fReflSurfXY[1][4]+fReflSurfXY[0][4])); 
  fRefl[34] = 0.7071067811865475*(fReflSurfXY[3][4]-1.0*fReflSurfXY[2][4]+fReflSurfXY[1][4]-1.0*fReflSurfXY[0][4]); 
  fRefl[35] = 0.0; 
  fRefl[36] = 0.7071067811865475*(fReflSurfXY[3][5]+fReflSurfXY[2][5]+fReflSurfXY[1][5]+fReflSurfXY[0][5]); 
  fRefl[37] = 0.7071067811865475*(fReflSurfXY[3][4]-1.0*(fReflSurfXY[2][4]+fReflSurfXY[1][4])+fReflSurfXY[0][4]); 
  fRefl[38] = 0.0; 
  fRefl[39] = 0.0; 
  fRefl[40] = 0.7071067811865475*(fReflSurfXY[3][5]+fReflSurfXY[2][5]-1.0*(fReflSurfXY[1][5]+fReflSurfXY[0][5])); 
  fRefl[41] = 0.7071067811865475*(fReflSurfXY[3][5]-1.0*fReflSurfXY[2][5]+fReflSurfXY[1][5]-1.0*fReflSurfXY[0][5]); 
  fRefl[42] = 0.0; 
  fRefl[43] = 0.0; 
  fRefl[44] = 0.7071067811865475*(fReflSurfXY[3][5]-1.0*(fReflSurfXY[2][5]+fReflSurfXY[1][5])+fReflSurfXY[0][5]); 
  fRefl[45] = 0.0; 
  fRefl[46] = 0.0; 
  fRefl[47] = 0.0; 
}
