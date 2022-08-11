#include <gkyl_bc_sheath_gyrokinetic_kernels.h> 


GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZQuad[3][8]; 
  

  vcutSq = 0.5*(3.16227766016838*phiWall[2]-3.16227766016838*phi[2]-2.449489742783178*phiWall[1]+2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))/(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*f[6]-6.0*f[2]))/(2.23606797749979*f[4]-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))/(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[2][0] = 0.0; 
  fReflZQuad[2][1] = 0.0; 
  fReflZQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } 
  } 
  fRefl[0] = 0.07856742013183861*(5.0*fReflZQuad[2][0]+8.0*fReflZQuad[1][0]+5.0*fReflZQuad[0][0]); 
  fRefl[1] = 0.5270462766947298*(fReflZQuad[2][0]-1.0*fReflZQuad[0][0]); 
  fRefl[2] = 0.07856742013183861*(5.0*fReflZQuad[2][1]+8.0*fReflZQuad[1][1]+5.0*fReflZQuad[0][1]); 
  fRefl[3] = 0.5270462766947298*(fReflZQuad[2][1]-1.0*fReflZQuad[0][1]); 
  fRefl[4] = 0.3513641844631533*(fReflZQuad[2][0]-2.0*fReflZQuad[1][0]+fReflZQuad[0][0]); 
  fRefl[5] = 0.07856742013183861*(5.0*fReflZQuad[2][2]+8.0*fReflZQuad[1][2]+5.0*fReflZQuad[0][2]); 
  fRefl[6] = 0.3513641844631534*(fReflZQuad[2][1]-2.0*fReflZQuad[1][1]+fReflZQuad[0][1]); 
  fRefl[7] = 0.5270462766947299*(fReflZQuad[2][2]-1.0*fReflZQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZQuad[3][8]; 
  

  vcutSq = 0.5*(3.16227766016838*phiWall[2]-3.16227766016838*phi[2]+2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))/(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*f[6]-6.0*f[2]))/(2.23606797749979*f[4]-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))/(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZQuad[2][0] = 0.0; 
  fReflZQuad[2][1] = 0.0; 
  fReflZQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } 
  } 
  fRefl[0] = 0.07856742013183861*(5.0*fReflZQuad[2][0]+8.0*fReflZQuad[1][0]+5.0*fReflZQuad[0][0]); 
  fRefl[1] = 0.5270462766947298*(fReflZQuad[2][0]-1.0*fReflZQuad[0][0]); 
  fRefl[2] = 0.07856742013183861*(5.0*fReflZQuad[2][1]+8.0*fReflZQuad[1][1]+5.0*fReflZQuad[0][1]); 
  fRefl[3] = 0.5270462766947298*(fReflZQuad[2][1]-1.0*fReflZQuad[0][1]); 
  fRefl[4] = 0.3513641844631533*(fReflZQuad[2][0]-2.0*fReflZQuad[1][0]+fReflZQuad[0][0]); 
  fRefl[5] = 0.07856742013183861*(5.0*fReflZQuad[2][2]+8.0*fReflZQuad[1][2]+5.0*fReflZQuad[0][2]); 
  fRefl[6] = 0.3513641844631534*(fReflZQuad[2][1]-2.0*fReflZQuad[1][1]+fReflZQuad[0][1]); 
  fRefl[7] = 0.5270462766947299*(fReflZQuad[2][2]-1.0*fReflZQuad[0][2]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZMuQuad[9][8]; 
  

  vcutSq = 0.5*(2.23606797749979*(1.414213562373095*phiWall[2]-1.414213562373095*phi[2])+1.732050807568877*(1.414213562373095*phi[1]-1.414213562373095*phiWall[1])+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
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
  fRefl[12] = 0.0; 
  fRefl[13] = 0.0; 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.0; 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.0; 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
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
  fRefl[12] = f[12]; 
  fRefl[13] = f[13]; 
  fRefl[14] = f[14]; 
  fRefl[15] = f[15]; 
  fRefl[16] = f[16]; 
  fRefl[17] = f[17]; 
  fRefl[18] = f[18]; 
  fRefl[19] = f[19]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))/(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))/(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))/(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]-1.0*(6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]-1.0*(6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))/(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))/(2.23606797749979*(f[9]+f[7])-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.08333333333333333*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.5*f[8])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.08333333333333333*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.5*f[8])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))/(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))/(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[8][0] = 0.0; 
  fReflZMuQuad[8][1] = 0.0; 
  fReflZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0])+8.0*(5.0*fReflZMuQuad[5][0]+8.0*fReflZMuQuad[4][0])+5.0*(8.0*fReflZMuQuad[3][0]+5.0*fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0])); 
  fRefl[1] = 0.0414086662499961*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0]-1.0*(5.0*fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0])); 
  fRefl[2] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1])+8.0*(5.0*fReflZMuQuad[5][1]+8.0*fReflZMuQuad[4][1])+5.0*(8.0*fReflZMuQuad[3][1]+5.0*fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1])); 
  fRefl[3] = 0.0414086662499961*(5.0*(fReflZMuQuad[8][0]-1.0*fReflZMuQuad[6][0])+8.0*(fReflZMuQuad[5][0]-1.0*fReflZMuQuad[3][0])+5.0*(fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[4] = 0.0414086662499961*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1]-1.0*(5.0*fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1])); 
  fRefl[5] = 0.2777777777777778*(fReflZMuQuad[8][0]-1.0*(fReflZMuQuad[6][0]+fReflZMuQuad[2][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.0414086662499961*(5.0*(fReflZMuQuad[8][1]-1.0*fReflZMuQuad[6][1])+8.0*(fReflZMuQuad[5][1]-1.0*fReflZMuQuad[3][1])+5.0*(fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1])); 
  fRefl[7] = 0.0276057774999974*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0]-2.0*(5.0*fReflZMuQuad[5][0]+8.0*fReflZMuQuad[4][0])+5.0*(fReflZMuQuad[2][0]-2.0*fReflZMuQuad[3][0])+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0]); 
  fRefl[8] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2])+8.0*(5.0*fReflZMuQuad[5][2]+8.0*fReflZMuQuad[4][2])+5.0*(8.0*fReflZMuQuad[3][2]+5.0*fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2])); 
  fRefl[9] = 0.0276057774999974*(5.0*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0])+8.0*(fReflZMuQuad[5][0]-2.0*fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+5.0*(fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[10] = 0.2777777777777778*(fReflZMuQuad[8][1]-1.0*(fReflZMuQuad[6][1]+fReflZMuQuad[2][1])+fReflZMuQuad[0][1]); 
  fRefl[11] = 0.02760577749999742*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1]-2.0*(5.0*fReflZMuQuad[5][1]+8.0*fReflZMuQuad[4][1])+5.0*(fReflZMuQuad[2][1]-2.0*fReflZMuQuad[3][1])+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1]); 
  fRefl[12] = 0.04140866624999612*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2]-1.0*(5.0*fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2])); 
  fRefl[13] = 0.1851851851851853*(fReflZMuQuad[8][0]-1.0*fReflZMuQuad[6][0]+2.0*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[5][0])+fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[14] = 0.04140866624999612*(5.0*(fReflZMuQuad[8][2]-1.0*fReflZMuQuad[6][2])+8.0*(fReflZMuQuad[5][2]-1.0*fReflZMuQuad[3][2])+5.0*(fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2])); 
  fRefl[15] = 0.1851851851851853*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0]-1.0*fReflZMuQuad[2][0]+2.0*fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[16] = 0.02760577749999742*(5.0*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1])+8.0*(fReflZMuQuad[5][1]-2.0*fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+5.0*(fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[17] = 0.1851851851851852*(fReflZMuQuad[8][1]-1.0*fReflZMuQuad[6][1]+2.0*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[5][1])+fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[18] = 0.2777777777777778*(fReflZMuQuad[8][2]-1.0*(fReflZMuQuad[6][2]+fReflZMuQuad[2][2])+fReflZMuQuad[0][2]); 
  fRefl[19] = 0.1851851851851852*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1]-1.0*fReflZMuQuad[2][1]+2.0*fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; double xc, b, xbarVal, fac; 
  double fReflZMuQuad[9][8]; 
  

  vcutSq = 0.5*(2.23606797749979*(1.414213562373095*phiWall[2]-1.414213562373095*phi[2])+1.732050807568877*(1.414213562373095*phiWall[1]-1.414213562373095*phi[1])+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
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
  fRefl[12] = 0.0; 
  fRefl[13] = 0.0; 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.0; 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.0; 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
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
  fRefl[12] = f[12]; 
  fRefl[13] = f[13]; 
  fRefl[14] = f[14]; 
  fRefl[15] = f[15]; 
  fRefl[16] = f[16]; 
  fRefl[17] = f[17]; 
  fRefl[18] = f[18]; 
  fRefl[19] = f[19]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))/(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(30.0*(f[15]+f[13])-1.0*(11.18033988749895*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))+5.0*(9.0*f[5]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))/(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-1.0*(5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7])))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[19]-6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))/(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]-1.0*(6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])-5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1])))+5.0*(9.0*f[5]-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(9.0*f[10]-1.0*(6.708203932499369*(f[6]-1.0*f[4])+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]-1.0*(6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))/(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))/(2.23606797749979*(f[9]+f[7])-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.08333333333333333*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.5*f[8])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.25*(2.23606797749979*(f[9]+f[7])-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.08333333333333333*(6.708203932499369*(f[16]+f[11])-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.5*f[8])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))/(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(45.0*f[17]-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))/(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(45.0*f[19]+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflZMuQuad[8][0] = 0.0; 
  fReflZMuQuad[8][1] = 0.0; 
  fReflZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0])+8.0*(5.0*fReflZMuQuad[5][0]+8.0*fReflZMuQuad[4][0])+5.0*(8.0*fReflZMuQuad[3][0]+5.0*fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0])); 
  fRefl[1] = 0.0414086662499961*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0]-1.0*(5.0*fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0])); 
  fRefl[2] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1])+8.0*(5.0*fReflZMuQuad[5][1]+8.0*fReflZMuQuad[4][1])+5.0*(8.0*fReflZMuQuad[3][1]+5.0*fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1])); 
  fRefl[3] = 0.0414086662499961*(5.0*(fReflZMuQuad[8][0]-1.0*fReflZMuQuad[6][0])+8.0*(fReflZMuQuad[5][0]-1.0*fReflZMuQuad[3][0])+5.0*(fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[4] = 0.0414086662499961*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1]-1.0*(5.0*fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1])); 
  fRefl[5] = 0.2777777777777778*(fReflZMuQuad[8][0]-1.0*(fReflZMuQuad[6][0]+fReflZMuQuad[2][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.0414086662499961*(5.0*(fReflZMuQuad[8][1]-1.0*fReflZMuQuad[6][1])+8.0*(fReflZMuQuad[5][1]-1.0*fReflZMuQuad[3][1])+5.0*(fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1])); 
  fRefl[7] = 0.0276057774999974*(5.0*fReflZMuQuad[8][0]+8.0*fReflZMuQuad[7][0]+5.0*fReflZMuQuad[6][0]-2.0*(5.0*fReflZMuQuad[5][0]+8.0*fReflZMuQuad[4][0])+5.0*(fReflZMuQuad[2][0]-2.0*fReflZMuQuad[3][0])+8.0*fReflZMuQuad[1][0]+5.0*fReflZMuQuad[0][0]); 
  fRefl[8] = 0.006172839506172839*(5.0*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2])+8.0*(5.0*fReflZMuQuad[5][2]+8.0*fReflZMuQuad[4][2])+5.0*(8.0*fReflZMuQuad[3][2]+5.0*fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2])); 
  fRefl[9] = 0.0276057774999974*(5.0*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0])+8.0*(fReflZMuQuad[5][0]-2.0*fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+5.0*(fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[10] = 0.2777777777777778*(fReflZMuQuad[8][1]-1.0*(fReflZMuQuad[6][1]+fReflZMuQuad[2][1])+fReflZMuQuad[0][1]); 
  fRefl[11] = 0.02760577749999742*(5.0*fReflZMuQuad[8][1]+8.0*fReflZMuQuad[7][1]+5.0*fReflZMuQuad[6][1]-2.0*(5.0*fReflZMuQuad[5][1]+8.0*fReflZMuQuad[4][1])+5.0*(fReflZMuQuad[2][1]-2.0*fReflZMuQuad[3][1])+8.0*fReflZMuQuad[1][1]+5.0*fReflZMuQuad[0][1]); 
  fRefl[12] = 0.04140866624999612*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2]-1.0*(5.0*fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2])); 
  fRefl[13] = 0.1851851851851853*(fReflZMuQuad[8][0]-1.0*fReflZMuQuad[6][0]+2.0*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[5][0])+fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[14] = 0.04140866624999612*(5.0*(fReflZMuQuad[8][2]-1.0*fReflZMuQuad[6][2])+8.0*(fReflZMuQuad[5][2]-1.0*fReflZMuQuad[3][2])+5.0*(fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2])); 
  fRefl[15] = 0.1851851851851853*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0]-1.0*fReflZMuQuad[2][0]+2.0*fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[16] = 0.02760577749999742*(5.0*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1])+8.0*(fReflZMuQuad[5][1]-2.0*fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+5.0*(fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[17] = 0.1851851851851852*(fReflZMuQuad[8][1]-1.0*fReflZMuQuad[6][1]+2.0*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[5][1])+fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[18] = 0.2777777777777778*(fReflZMuQuad[8][2]-1.0*(fReflZMuQuad[6][2]+fReflZMuQuad[2][2])+fReflZMuQuad[0][2]); 
  fRefl[19] = 0.1851851851851852*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1]-1.0*fReflZMuQuad[2][1]+2.0*fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  } 

 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; double xc, b, xbarVal, fac; 
  double fReflXYQuad[9][20]; 
  double fReflXYZMuQuad[9][8]; 
  

// node (x,y)_1 
  vcutSq_i = -0.01*(2.23606797749979*((-63.63961030678928*phiWall[19])+63.63961030678928*phi[19]-35.35533905932738*phiWall[9]+35.35533905932738*phi[9]-14.14213562373095*phiWall[8]+14.14213562373095*phi[8]-14.14213562373095*phiWall[7]+14.14213562373095*phi[7]+1.732050807568877*((-21.21320343559643*phiWall[6])+21.21320343559643*phi[6]-21.21320343559643*phiWall[5]+21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*phi[2]+21.21320343559643*phiWall[1]-21.21320343559643*phi[1])+1.732050807568877*((-42.42640687119286*phiWall[18])+42.42640687119286*phi[18]-42.42640687119286*phiWall[17]+42.42640687119286*phi[17]+63.63961030678928*phiWall[10]-63.63961030678928*phi[10]+35.35533905932738*phiWall[3]-35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*phi[16]+7.071067811865476*phiWall[15]-7.071067811865476*phi[15]+2.828427124746191*phiWall[12]-2.828427124746191*phi[12]+2.828427124746191*phiWall[11]-2.828427124746191*phi[11])+14.14213562373095*phiWall[14]-14.14213562373095*phi[14]+14.14213562373095*phiWall[13]-14.14213562373095*phi[13])-63.63961030678928*phiWall[4]+63.63961030678928*phi[4]-35.35533905932738*phiWall[0]+35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[0][0] = 0.0; 
  fReflXYQuad[0][1] = 0.0; 
  fReflXYQuad[0][2] = 0.0; 
  fReflXYQuad[0][3] = 0.0; 
  fReflXYQuad[0][4] = 0.0; 
  fReflXYQuad[0][5] = 0.0; 
  fReflXYQuad[0][6] = 0.0; 
  fReflXYQuad[0][7] = 0.0; 
  fReflXYQuad[0][8] = 0.0; 
  fReflXYQuad[0][9] = 0.0; 
  fReflXYQuad[0][10] = 0.0; 
  fReflXYQuad[0][11] = 0.0; 
  fReflXYQuad[0][12] = 0.0; 
  fReflXYQuad[0][13] = 0.0; 
  fReflXYQuad[0][14] = 0.0; 
  fReflXYQuad[0][15] = 0.0; 
  fReflXYQuad[0][16] = 0.0; 
  fReflXYQuad[0][17] = 0.0; 
  fReflXYQuad[0][18] = 0.0; 
  fReflXYQuad[0][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = -0.02*(30.0*(f[32]+f[31])-1.0*(11.18033988749895*(2.0*(f[17]+f[16])-3.0*(f[2]+f[1]))+5.0*(9.0*f[6]+5.0*f[0]))); 
  fReflXYQuad[0][1] = -0.03333333333333333*(2.0*(9.0*(f[57]+f[56])-6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[8]+f[7])-3.0*f[21])-5.0*f[3])); 
  fReflXYQuad[0][2] = -0.03333333333333333*(2.0*(9.0*(f[60]+f[59])-6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[10]+f[9])-3.0*f[22])-5.0*f[4])); 
  fReflXYQuad[0][3] = -0.03333333333333333*(2.0*(9.0*(f[69]+f[68])-6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[13]+f[12])-3.0*f[25])-5.0*f[5])); 
  fReflXYQuad[0][4] = -0.02*(30.0*(f[88]+f[87])-1.0*(11.18033988749895*(2.0*(f[62]+f[61])-3.0*(f[24]+f[23]))+5.0*(9.0*f[51]+5.0*f[11]))); 
  fReflXYQuad[0][5] = -0.02*(30.0*(f[92]+f[91])-1.0*(11.18033988749895*(2.0*(f[71]+f[70])-3.0*(f[27]+f[26]))+5.0*(9.0*f[52]+5.0*f[14]))); 
  fReflXYQuad[0][6] = -0.02*(30.0*(f[95]+f[94])-1.0*(11.18033988749895*(2.0*(f[75]+f[74])-3.0*(f[29]+f[28]))+5.0*(9.0*f[53]+5.0*f[15]))); 
  fReflXYQuad[0][7] = 0.1*(9.0*f[58]-6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[0][8] = 0.1*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[0][9] = 0.1*(9.0*f[80]-6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[0][10] = -0.03333333333333333*(2.0*(9.0*(f[108]+f[107])-6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[55]+f[54])-3.0*f[86])-5.0*f[30])); 
  fReflXYQuad[0][11] = 0.1*(9.0*f[89]-6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[0][12] = 0.1*(9.0*f[90]-6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[0][13] = 0.1*(9.0*f[93]-6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[0][14] = 0.1*(9.0*f[100]-6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[0][15] = 0.1*(9.0*f[103]-6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[0][16] = 0.1*(9.0*f[104]-6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[0][17] = 0.1*(9.0*f[109]-6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[0][18] = 0.1*(9.0*f[110]-6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[0][19] = 0.1*(9.0*f[111]-6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))/(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))-5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))-5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))-5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))-5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))/(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_2 
  vcutSq_i = -0.025*(1.732050807568877*(21.21320343559643*phiWall[18]-21.21320343559643*phi[18]+14.14213562373095*phiWall[3]-14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[15]-2.828427124746191*phi[15]-1.414213562373095*phiWall[12]+1.414213562373095*phi[12])-7.071067811865476*phiWall[14]+7.071067811865476*phi[14]+5.656854249492382*phiWall[13]-5.656854249492382*phi[13])+2.23606797749979*((-14.14213562373095*phiWall[9])+14.14213562373095*phi[9]+7.071067811865476*phiWall[8]-7.071067811865476*phi[8]-5.656854249492382*phiWall[7]+5.656854249492382*phi[7]+1.732050807568877*(8.485281374238571*phi[5]-8.485281374238571*phiWall[5])+8.485281374238571*phiWall[1]-8.485281374238571*phi[1])-14.14213562373095*phiWall[0]+14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[1][0] = 0.0; 
  fReflXYQuad[1][1] = 0.0; 
  fReflXYQuad[1][2] = 0.0; 
  fReflXYQuad[1][3] = 0.0; 
  fReflXYQuad[1][4] = 0.0; 
  fReflXYQuad[1][5] = 0.0; 
  fReflXYQuad[1][6] = 0.0; 
  fReflXYQuad[1][7] = 0.0; 
  fReflXYQuad[1][8] = 0.0; 
  fReflXYQuad[1][9] = 0.0; 
  fReflXYQuad[1][10] = 0.0; 
  fReflXYQuad[1][11] = 0.0; 
  fReflXYQuad[1][12] = 0.0; 
  fReflXYQuad[1][13] = 0.0; 
  fReflXYQuad[1][14] = 0.0; 
  fReflXYQuad[1][15] = 0.0; 
  fReflXYQuad[1][16] = 0.0; 
  fReflXYQuad[1][17] = 0.0; 
  fReflXYQuad[1][18] = 0.0; 
  fReflXYQuad[1][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = 0.05*(2.23606797749979*(6.708203932499369*f[32]-1.0*(5.0*f[17]+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]); 
  fReflXYQuad[1][1] = 0.01666666666666667*(45.0*f[57]-6.708203932499369*(5.0*f[34]-4.0*f[33])+6.0*(5.0*f[3]-6.708203932499369*f[7])); 
  fReflXYQuad[1][2] = 0.01666666666666667*(45.0*f[60]-6.708203932499369*(5.0*f[38]-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])); 
  fReflXYQuad[1][3] = 0.01666666666666667*(45.0*f[69]-6.708203932499369*(5.0*f[44]-4.0*f[43])+6.0*(5.0*f[5]-6.708203932499369*f[12])); 
  fReflXYQuad[1][4] = 0.05*(2.23606797749979*(6.708203932499369*f[88]-1.0*(5.0*f[62]+2.0*(3.0*f[23]-2.0*f[61])))+10.0*f[11]); 
  fReflXYQuad[1][5] = 0.05*(2.23606797749979*(6.708203932499369*f[92]-1.0*(5.0*f[71]+2.0*(3.0*f[26]-2.0*f[70])))+10.0*f[14]); 
  fReflXYQuad[1][6] = 0.05*(2.23606797749979*(6.708203932499369*f[95]-1.0*(5.0*f[75]+2.0*(3.0*f[28]-2.0*f[74])))+10.0*f[15]); 
  fReflXYQuad[1][7] = -0.1*(6.708203932499369*f[35]-5.0*f[18]); 
  fReflXYQuad[1][8] = -0.1*(6.708203932499369*f[40]-5.0*f[19]); 
  fReflXYQuad[1][9] = -0.1*(6.708203932499369*f[47]-5.0*f[20]); 
  fReflXYQuad[1][10] = 0.01666666666666667*(45.0*f[108]-6.708203932499369*(5.0*f[97]-4.0*f[96])+6.0*(5.0*f[30]-6.708203932499369*f[54])); 
  fReflXYQuad[1][11] = -0.1*(6.708203932499369*f[63]-5.0*f[39]); 
  fReflXYQuad[1][12] = -0.1*(6.708203932499369*f[66]-5.0*f[42]); 
  fReflXYQuad[1][13] = -0.1*(6.708203932499369*f[72]-5.0*f[45]); 
  fReflXYQuad[1][14] = -0.1*(6.708203932499369*f[77]-5.0*f[46]); 
  fReflXYQuad[1][15] = -0.1*(6.708203932499369*f[81]-5.0*f[49]); 
  fReflXYQuad[1][16] = -0.1*(6.708203932499369*f[83]-5.0*f[50]); 
  fReflXYQuad[1][17] = -0.1*(6.708203932499369*f[98]-5.0*f[76]); 
  fReflXYQuad[1][18] = -0.1*(6.708203932499369*f[101]-5.0*f[79]); 
  fReflXYQuad[1][19] = -0.1*(6.708203932499369*f[105]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))-5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))-5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))/(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))/(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[40]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[40]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[77]-1.0*(6.708203932499369*(f[46]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[77]-1.0*(6.708203932499369*(f[46]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[66]-1.0*(6.708203932499369*(f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[66]-1.0*(6.708203932499369*(f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_3 
  vcutSq_i = 0.01*(2.23606797749979*((-63.63961030678928*phiWall[19])+63.63961030678928*phi[19]+35.35533905932738*phiWall[9]-35.35533905932738*phi[9]+14.14213562373095*phiWall[8]-14.14213562373095*phi[8]+14.14213562373095*phiWall[7]-14.14213562373095*phi[7]+1.732050807568877*((-21.21320343559643*phiWall[6])+21.21320343559643*(phi[6]+phiWall[5])-21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*(phi[2]+phiWall[1])+21.21320343559643*phi[1])+1.732050807568877*(42.42640687119286*phiWall[18]-42.42640687119286*(phi[18]+phiWall[17])+42.42640687119286*phi[17]+63.63961030678928*phiWall[10]-63.63961030678928*phi[10]-35.35533905932738*phiWall[3]+35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*(phi[16]+phiWall[15])+7.071067811865476*phi[15]-2.828427124746191*phiWall[12]+2.828427124746191*(phi[12]+phiWall[11])-2.828427124746191*phi[11])-14.14213562373095*phiWall[14]+14.14213562373095*phi[14]-14.14213562373095*phiWall[13]+14.14213562373095*phi[13])-63.63961030678928*phiWall[4]+63.63961030678928*phi[4]+35.35533905932738*phiWall[0]-35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[2][0] = 0.0; 
  fReflXYQuad[2][1] = 0.0; 
  fReflXYQuad[2][2] = 0.0; 
  fReflXYQuad[2][3] = 0.0; 
  fReflXYQuad[2][4] = 0.0; 
  fReflXYQuad[2][5] = 0.0; 
  fReflXYQuad[2][6] = 0.0; 
  fReflXYQuad[2][7] = 0.0; 
  fReflXYQuad[2][8] = 0.0; 
  fReflXYQuad[2][9] = 0.0; 
  fReflXYQuad[2][10] = 0.0; 
  fReflXYQuad[2][11] = 0.0; 
  fReflXYQuad[2][12] = 0.0; 
  fReflXYQuad[2][13] = 0.0; 
  fReflXYQuad[2][14] = 0.0; 
  fReflXYQuad[2][15] = 0.0; 
  fReflXYQuad[2][16] = 0.0; 
  fReflXYQuad[2][17] = 0.0; 
  fReflXYQuad[2][18] = 0.0; 
  fReflXYQuad[2][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.02*(2.23606797749979*(13.41640786499874*(f[32]-1.0*f[31])-5.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))+5.0*(9.0*f[6]-5.0*f[0])); 
  fReflXYQuad[2][1] = -0.03333333333333333*(2.0*(9.0*f[57]-1.0*(9.0*f[56]+6.708203932499369*(f[34]+f[33])))+3.0*(9.0*f[21]-1.0*(6.708203932499369*(f[8]-1.0*f[7])+5.0*f[3]))); 
  fReflXYQuad[2][2] = -0.03333333333333333*(2.0*(9.0*f[60]-1.0*(9.0*f[59]+6.708203932499369*(f[38]+f[37])))+3.0*(9.0*f[22]-1.0*(6.708203932499369*(f[10]-1.0*f[9])+5.0*f[4]))); 
  fReflXYQuad[2][3] = -0.03333333333333333*(2.0*(9.0*f[69]-1.0*(9.0*f[68]+6.708203932499369*(f[44]+f[43])))+3.0*(9.0*f[25]-1.0*(6.708203932499369*(f[13]-1.0*f[12])+5.0*f[5]))); 
  fReflXYQuad[2][4] = -0.02*(2.23606797749979*(13.41640786499874*(f[88]-1.0*f[87])-5.0*(2.0*(f[62]+f[61])+3.0*(f[24]-1.0*f[23])))+5.0*(9.0*f[51]-5.0*f[11])); 
  fReflXYQuad[2][5] = -0.02*(2.23606797749979*(13.41640786499874*(f[92]-1.0*f[91])-5.0*(2.0*(f[71]+f[70])+3.0*(f[27]-1.0*f[26])))+5.0*(9.0*f[52]-5.0*f[14])); 
  fReflXYQuad[2][6] = -0.02*(2.23606797749979*(13.41640786499874*(f[95]-1.0*f[94])-5.0*(2.0*(f[75]+f[74])+3.0*(f[29]-1.0*f[28])))+5.0*(9.0*f[53]-5.0*f[15])); 
  fReflXYQuad[2][7] = -0.1*(9.0*f[58]-1.0*(6.708203932499369*(f[36]-1.0*f[35])+5.0*f[18])); 
  fReflXYQuad[2][8] = -0.1*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])); 
  fReflXYQuad[2][9] = -0.1*(9.0*f[80]-1.0*(6.708203932499369*(f[48]-1.0*f[47])+5.0*f[20])); 
  fReflXYQuad[2][10] = -0.03333333333333333*(2.0*(9.0*f[108]-1.0*(9.0*f[107]+6.708203932499369*(f[97]+f[96])))+3.0*(9.0*f[86]-1.0*(6.708203932499369*(f[55]-1.0*f[54])+5.0*f[30]))); 
  fReflXYQuad[2][11] = -0.1*(9.0*f[89]-1.0*(6.708203932499369*(f[64]-1.0*f[63])+5.0*f[39])); 
  fReflXYQuad[2][12] = -0.1*(9.0*f[90]-1.0*(6.708203932499369*(f[67]-1.0*f[66])+5.0*f[42])); 
  fReflXYQuad[2][13] = -0.1*(9.0*f[93]-1.0*(6.708203932499369*(f[73]-1.0*f[72])+5.0*f[45])); 
  fReflXYQuad[2][14] = -0.1*(9.0*f[100]-1.0*(6.708203932499369*(f[78]-1.0*f[77])+5.0*f[46])); 
  fReflXYQuad[2][15] = -0.1*(9.0*f[103]-1.0*(6.708203932499369*(f[82]-1.0*f[81])+5.0*f[49])); 
  fReflXYQuad[2][16] = -0.1*(9.0*f[104]-1.0*(6.708203932499369*(f[84]-1.0*f[83])+5.0*f[50])); 
  fReflXYQuad[2][17] = -0.1*(9.0*f[109]-1.0*(6.708203932499369*(f[99]-1.0*f[98])+5.0*f[76])); 
  fReflXYQuad[2][18] = -0.1*(9.0*f[110]-1.0*(6.708203932499369*(f[102]-1.0*f[101])+5.0*f[79])); 
  fReflXYQuad[2][19] = -0.1*(9.0*f[111]-1.0*(6.708203932499369*(f[106]-1.0*f[105])+5.0*f[85])); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_4 
  vcutSq_i = -0.025*(1.732050807568877*(21.21320343559643*phiWall[17]-21.21320343559643*phi[17]+14.14213562373095*phiWall[3]-14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[16]-2.828427124746191*phi[16]-1.414213562373095*phiWall[11]+1.414213562373095*phi[11])+5.656854249492382*phiWall[14]-5.656854249492382*phi[14]-7.071067811865476*phiWall[13]+7.071067811865476*phi[13])+2.23606797749979*((-14.14213562373095*phiWall[9])+14.14213562373095*phi[9]-5.656854249492382*phiWall[8]+5.656854249492382*phi[8]+7.071067811865476*phiWall[7]-7.071067811865476*phi[7]+1.732050807568877*(8.485281374238571*phi[6]-8.485281374238571*phiWall[6])+8.485281374238571*phiWall[2]-8.485281374238571*phi[2])-14.14213562373095*phiWall[0]+14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[3][0] = 0.0; 
  fReflXYQuad[3][1] = 0.0; 
  fReflXYQuad[3][2] = 0.0; 
  fReflXYQuad[3][3] = 0.0; 
  fReflXYQuad[3][4] = 0.0; 
  fReflXYQuad[3][5] = 0.0; 
  fReflXYQuad[3][6] = 0.0; 
  fReflXYQuad[3][7] = 0.0; 
  fReflXYQuad[3][8] = 0.0; 
  fReflXYQuad[3][9] = 0.0; 
  fReflXYQuad[3][10] = 0.0; 
  fReflXYQuad[3][11] = 0.0; 
  fReflXYQuad[3][12] = 0.0; 
  fReflXYQuad[3][13] = 0.0; 
  fReflXYQuad[3][14] = 0.0; 
  fReflXYQuad[3][15] = 0.0; 
  fReflXYQuad[3][16] = 0.0; 
  fReflXYQuad[3][17] = 0.0; 
  fReflXYQuad[3][18] = 0.0; 
  fReflXYQuad[3][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.05*(2.23606797749979*(6.708203932499369*f[31]+4.0*f[17]-1.0*(5.0*f[16]+6.0*f[2]))+10.0*f[0]); 
  fReflXYQuad[3][1] = 0.01666666666666667*(45.0*f[56]+6.708203932499369*(4.0*f[34]-5.0*f[33])+6.0*(5.0*f[3]-6.708203932499369*f[8])); 
  fReflXYQuad[3][2] = 0.01666666666666667*(45.0*f[59]+6.708203932499369*(4.0*f[38]-5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])); 
  fReflXYQuad[3][3] = 0.01666666666666667*(45.0*f[68]+6.708203932499369*(4.0*f[44]-5.0*f[43])+6.0*(5.0*f[5]-6.708203932499369*f[13])); 
  fReflXYQuad[3][4] = 0.05*(2.23606797749979*(6.708203932499369*f[87]+4.0*f[62]-1.0*(5.0*f[61]+6.0*f[24]))+10.0*f[11]); 
  fReflXYQuad[3][5] = 0.05*(2.23606797749979*(6.708203932499369*f[91]+4.0*f[71]-1.0*(5.0*f[70]+6.0*f[27]))+10.0*f[14]); 
  fReflXYQuad[3][6] = 0.05*(2.23606797749979*(6.708203932499369*f[94]+4.0*f[75]-1.0*(5.0*f[74]+6.0*f[29]))+10.0*f[15]); 
  fReflXYQuad[3][7] = -0.1*(6.708203932499369*f[36]-5.0*f[18]); 
  fReflXYQuad[3][8] = -0.1*(6.708203932499369*f[41]-5.0*f[19]); 
  fReflXYQuad[3][9] = -0.1*(6.708203932499369*f[48]-5.0*f[20]); 
  fReflXYQuad[3][10] = 0.01666666666666667*(45.0*f[107]+6.708203932499369*(4.0*f[97]-5.0*f[96])+6.0*(5.0*f[30]-6.708203932499369*f[55])); 
  fReflXYQuad[3][11] = -0.1*(6.708203932499369*f[64]-5.0*f[39]); 
  fReflXYQuad[3][12] = -0.1*(6.708203932499369*f[67]-5.0*f[42]); 
  fReflXYQuad[3][13] = -0.1*(6.708203932499369*f[73]-5.0*f[45]); 
  fReflXYQuad[3][14] = -0.1*(6.708203932499369*f[78]-5.0*f[46]); 
  fReflXYQuad[3][15] = -0.1*(6.708203932499369*f[82]-5.0*f[49]); 
  fReflXYQuad[3][16] = -0.1*(6.708203932499369*f[84]-5.0*f[50]); 
  fReflXYQuad[3][17] = -0.1*(6.708203932499369*f[99]-5.0*f[76]); 
  fReflXYQuad[3][18] = -0.1*(6.708203932499369*f[102]-5.0*f[79]); 
  fReflXYQuad[3][19] = -0.1*(6.708203932499369*f[106]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))-5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))-5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))/(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))/(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[41]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[41]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[78]-1.0*(6.708203932499369*(f[46]-1.0*f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[78]-1.0*(6.708203932499369*(f[46]-1.0*f[41])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[67]-1.0*(6.708203932499369*(f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[67]-1.0*(6.708203932499369*(f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_5 
  vcutSq_i = 0.125*(3.872983346207417*(1.414213562373095*phiWall[14]-1.414213562373095*phi[14]+1.414213562373095*phiWall[13]-1.414213562373095*phi[13])+2.23606797749979*(2.828427124746191*phiWall[9]-2.828427124746191*phi[9]-1.414213562373095*phiWall[8]+1.414213562373095*phi[8]-1.414213562373095*phiWall[7]+1.414213562373095*phi[7])+1.732050807568877*(2.828427124746191*phi[3]-2.828427124746191*phiWall[3])+2.828427124746191*phiWall[0]-2.828427124746191*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[4][0] = 0.0; 
  fReflXYQuad[4][1] = 0.0; 
  fReflXYQuad[4][2] = 0.0; 
  fReflXYQuad[4][3] = 0.0; 
  fReflXYQuad[4][4] = 0.0; 
  fReflXYQuad[4][5] = 0.0; 
  fReflXYQuad[4][6] = 0.0; 
  fReflXYQuad[4][7] = 0.0; 
  fReflXYQuad[4][8] = 0.0; 
  fReflXYQuad[4][9] = 0.0; 
  fReflXYQuad[4][10] = 0.0; 
  fReflXYQuad[4][11] = 0.0; 
  fReflXYQuad[4][12] = 0.0; 
  fReflXYQuad[4][13] = 0.0; 
  fReflXYQuad[4][14] = 0.0; 
  fReflXYQuad[4][15] = 0.0; 
  fReflXYQuad[4][16] = 0.0; 
  fReflXYQuad[4][17] = 0.0; 
  fReflXYQuad[4][18] = 0.0; 
  fReflXYQuad[4][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[4][0] = -0.25*(2.23606797749979*(f[17]+f[16])-2.0*f[0]); 
  fReflXYQuad[4][1] = -0.08333333333333333*(6.708203932499369*(f[34]+f[33])-6.0*f[3]); 
  fReflXYQuad[4][2] = -0.08333333333333333*(6.708203932499369*(f[38]+f[37])-6.0*f[4]); 
  fReflXYQuad[4][3] = -0.08333333333333333*(6.708203932499369*(f[44]+f[43])-6.0*f[5]); 
  fReflXYQuad[4][4] = -0.25*(2.23606797749979*(f[62]+f[61])-2.0*f[11]); 
  fReflXYQuad[4][5] = -0.25*(2.23606797749979*(f[71]+f[70])-2.0*f[14]); 
  fReflXYQuad[4][6] = -0.25*(2.23606797749979*(f[75]+f[74])-2.0*f[15]); 
  fReflXYQuad[4][7] = 0.5*f[18]; 
  fReflXYQuad[4][8] = 0.5*f[19]; 
  fReflXYQuad[4][9] = 0.5*f[20]; 
  fReflXYQuad[4][10] = -0.08333333333333333*(6.708203932499369*(f[97]+f[96])-6.0*f[30]); 
  fReflXYQuad[4][11] = 0.5*f[39]; 
  fReflXYQuad[4][12] = 0.5*f[42]; 
  fReflXYQuad[4][13] = 0.5*f[45]; 
  fReflXYQuad[4][14] = 0.5*f[46]; 
  fReflXYQuad[4][15] = 0.5*f[49]; 
  fReflXYQuad[4][16] = 0.5*f[50]; 
  fReflXYQuad[4][17] = 0.5*f[76]; 
  fReflXYQuad[4][18] = 0.5*f[79]; 
  fReflXYQuad[4][19] = 0.5*f[85]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.05*(9.0*f[79]-6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.05*(9.0*f[79]-6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))/(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(6.708203932499369*f[42]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(6.708203932499369*f[42]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.05*(9.0*f[79]-1.0*(6.708203932499369*(f[46]-1.0*f[42])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.05*(9.0*f[79]-1.0*(6.708203932499369*(f[46]-1.0*f[42])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))/(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(6.708203932499369*f[46]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(6.708203932499369*f[46]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))/(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.04166666666666666*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.25*f[19])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.04166666666666666*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.25*f[19])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))/(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(6.708203932499369*f[46]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(6.708203932499369*f[46]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[79]+6.708203932499369*(f[46]-1.0*f[42])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[79]+6.708203932499369*(f[46]-1.0*f[42])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(6.708203932499369*f[42]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(6.708203932499369*f[42]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.05*(9.0*f[79]+6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.05*(9.0*f[79]+6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } 
  } 
  fReflXYQuad[4][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[4][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_6 
  vcutSq_i = 0.025*(1.732050807568877*(21.21320343559643*phiWall[17]-21.21320343559643*phi[17]-14.14213562373095*phiWall[3]+14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[16]-2.828427124746191*phi[16]-1.414213562373095*phiWall[11]+1.414213562373095*phi[11])-5.656854249492382*phiWall[14]+5.656854249492382*phi[14]+7.071067811865476*phiWall[13]-7.071067811865476*phi[13])+2.23606797749979*(14.14213562373095*phiWall[9]-14.14213562373095*phi[9]+5.656854249492382*phiWall[8]-5.656854249492382*phi[8]-7.071067811865476*phiWall[7]+7.071067811865476*phi[7]+1.732050807568877*(8.485281374238571*phi[6]-8.485281374238571*phiWall[6])+8.485281374238571*phiWall[2]-8.485281374238571*phi[2])+14.14213562373095*phiWall[0]-14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[5][0] = 0.0; 
  fReflXYQuad[5][1] = 0.0; 
  fReflXYQuad[5][2] = 0.0; 
  fReflXYQuad[5][3] = 0.0; 
  fReflXYQuad[5][4] = 0.0; 
  fReflXYQuad[5][5] = 0.0; 
  fReflXYQuad[5][6] = 0.0; 
  fReflXYQuad[5][7] = 0.0; 
  fReflXYQuad[5][8] = 0.0; 
  fReflXYQuad[5][9] = 0.0; 
  fReflXYQuad[5][10] = 0.0; 
  fReflXYQuad[5][11] = 0.0; 
  fReflXYQuad[5][12] = 0.0; 
  fReflXYQuad[5][13] = 0.0; 
  fReflXYQuad[5][14] = 0.0; 
  fReflXYQuad[5][15] = 0.0; 
  fReflXYQuad[5][16] = 0.0; 
  fReflXYQuad[5][17] = 0.0; 
  fReflXYQuad[5][18] = 0.0; 
  fReflXYQuad[5][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[5][0] = -0.05*(15.0*f[31]-1.0*(2.23606797749979*(4.0*f[17]-5.0*f[16]+6.0*f[2])+10.0*f[0])); 
  fReflXYQuad[5][1] = -0.01666666666666667*(45.0*f[56]-1.0*(6.708203932499369*(4.0*f[34]-5.0*f[33])+6.0*(6.708203932499369*f[8]+5.0*f[3]))); 
  fReflXYQuad[5][2] = -0.01666666666666667*(45.0*f[59]-1.0*(6.708203932499369*(4.0*f[38]-5.0*f[37])+6.0*(6.708203932499369*f[10]+5.0*f[4]))); 
  fReflXYQuad[5][3] = -0.01666666666666667*(45.0*f[68]-1.0*(6.708203932499369*(4.0*f[44]-5.0*f[43])+6.0*(6.708203932499369*f[13]+5.0*f[5]))); 
  fReflXYQuad[5][4] = -0.05*(15.0*f[87]-1.0*(2.23606797749979*(4.0*f[62]-5.0*f[61]+6.0*f[24])+10.0*f[11])); 
  fReflXYQuad[5][5] = -0.05*(15.0*f[91]-1.0*(2.23606797749979*(4.0*f[71]-5.0*f[70]+6.0*f[27])+10.0*f[14])); 
  fReflXYQuad[5][6] = -0.05*(15.0*f[94]-1.0*(2.23606797749979*(4.0*f[75]-5.0*f[74]+6.0*f[29])+10.0*f[15])); 
  fReflXYQuad[5][7] = 0.1*(6.708203932499369*f[36]+5.0*f[18]); 
  fReflXYQuad[5][8] = 0.1*(6.708203932499369*f[41]+5.0*f[19]); 
  fReflXYQuad[5][9] = 0.1*(6.708203932499369*f[48]+5.0*f[20]); 
  fReflXYQuad[5][10] = -0.01666666666666667*(45.0*f[107]-1.0*(6.708203932499369*(4.0*f[97]-5.0*f[96])+6.0*(6.708203932499369*f[55]+5.0*f[30]))); 
  fReflXYQuad[5][11] = 0.1*(6.708203932499369*f[64]+5.0*f[39]); 
  fReflXYQuad[5][12] = 0.1*(6.708203932499369*f[67]+5.0*f[42]); 
  fReflXYQuad[5][13] = 0.1*(6.708203932499369*f[73]+5.0*f[45]); 
  fReflXYQuad[5][14] = 0.1*(6.708203932499369*f[78]+5.0*f[46]); 
  fReflXYQuad[5][15] = 0.1*(6.708203932499369*f[82]+5.0*f[49]); 
  fReflXYQuad[5][16] = 0.1*(6.708203932499369*f[84]+5.0*f[50]); 
  fReflXYQuad[5][17] = 0.1*(6.708203932499369*f[99]+5.0*f[76]); 
  fReflXYQuad[5][18] = 0.1*(6.708203932499369*f[102]+5.0*f[79]); 
  fReflXYQuad[5][19] = 0.1*(6.708203932499369*f[106]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[42]-1.0*f[41])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[42]-1.0*f[41])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*(f[79]-1.0*f[78]+f[67])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*(f[79]-1.0*f[78]+f[67])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[46]-1.0*f[41])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[46]-1.0*f[41])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[41]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[41]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[79]+f[78])-1.0*(9.0*f[67]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[79]+f[78])-1.0*(9.0*f[67]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))/(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))+5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))+5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[5][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[5][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_7 
  vcutSq_i = -0.01*(2.23606797749979*(63.63961030678928*phiWall[19]-63.63961030678928*phi[19]-35.35533905932738*phiWall[9]+35.35533905932738*phi[9]-14.14213562373095*phiWall[8]+14.14213562373095*phi[8]-14.14213562373095*phiWall[7]+14.14213562373095*phi[7]+1.732050807568877*((-21.21320343559643*phiWall[6])+21.21320343559643*(phi[6]+phiWall[5])-21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*(phi[2]+phiWall[1])+21.21320343559643*phi[1])+1.732050807568877*(42.42640687119286*phiWall[18]-42.42640687119286*(phi[18]+phiWall[17])+42.42640687119286*phi[17]-63.63961030678928*phiWall[10]+63.63961030678928*phi[10]+35.35533905932738*phiWall[3]-35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*(phi[16]+phiWall[15])+7.071067811865476*phi[15]-2.828427124746191*phiWall[12]+2.828427124746191*(phi[12]+phiWall[11])-2.828427124746191*phi[11])+14.14213562373095*phiWall[14]-14.14213562373095*phi[14]+14.14213562373095*phiWall[13]-14.14213562373095*phi[13])+63.63961030678928*phiWall[4]-63.63961030678928*phi[4]-35.35533905932738*phiWall[0]+35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[6][0] = 0.0; 
  fReflXYQuad[6][1] = 0.0; 
  fReflXYQuad[6][2] = 0.0; 
  fReflXYQuad[6][3] = 0.0; 
  fReflXYQuad[6][4] = 0.0; 
  fReflXYQuad[6][5] = 0.0; 
  fReflXYQuad[6][6] = 0.0; 
  fReflXYQuad[6][7] = 0.0; 
  fReflXYQuad[6][8] = 0.0; 
  fReflXYQuad[6][9] = 0.0; 
  fReflXYQuad[6][10] = 0.0; 
  fReflXYQuad[6][11] = 0.0; 
  fReflXYQuad[6][12] = 0.0; 
  fReflXYQuad[6][13] = 0.0; 
  fReflXYQuad[6][14] = 0.0; 
  fReflXYQuad[6][15] = 0.0; 
  fReflXYQuad[6][16] = 0.0; 
  fReflXYQuad[6][17] = 0.0; 
  fReflXYQuad[6][18] = 0.0; 
  fReflXYQuad[6][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[6][0] = 0.02*(2.23606797749979*(13.41640786499874*(f[32]-1.0*f[31])+5.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+5.0*(5.0*f[0]-9.0*f[6])); 
  fReflXYQuad[6][1] = 0.03333333333333333*(2.0*(9.0*(f[57]-1.0*f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[7]-1.0*f[8])-3.0*f[21])+5.0*f[3])); 
  fReflXYQuad[6][2] = 0.03333333333333333*(2.0*(9.0*(f[60]-1.0*f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[10])-3.0*f[22])+5.0*f[4])); 
  fReflXYQuad[6][3] = 0.03333333333333333*(2.0*(9.0*(f[69]-1.0*f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[12]-1.0*f[13])-3.0*f[25])+5.0*f[5])); 
  fReflXYQuad[6][4] = 0.02*(2.23606797749979*(13.41640786499874*(f[88]-1.0*f[87])+5.0*(2.0*(f[62]+f[61])+3.0*(f[23]-1.0*f[24])))+5.0*(5.0*f[11]-9.0*f[51])); 
  fReflXYQuad[6][5] = 0.02*(2.23606797749979*(13.41640786499874*(f[92]-1.0*f[91])+5.0*(2.0*(f[71]+f[70])+3.0*(f[26]-1.0*f[27])))+5.0*(5.0*f[14]-9.0*f[52])); 
  fReflXYQuad[6][6] = 0.02*(2.23606797749979*(13.41640786499874*(f[95]-1.0*f[94])+5.0*(2.0*(f[75]+f[74])+3.0*(f[28]-1.0*f[29])))+5.0*(5.0*f[15]-9.0*f[53])); 
  fReflXYQuad[6][7] = -0.1*(9.0*f[58]+6.708203932499369*(f[36]-1.0*f[35])-5.0*f[18]); 
  fReflXYQuad[6][8] = -0.1*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]); 
  fReflXYQuad[6][9] = -0.1*(9.0*f[80]+6.708203932499369*(f[48]-1.0*f[47])-5.0*f[20]); 
  fReflXYQuad[6][10] = 0.03333333333333333*(2.0*(9.0*(f[108]-1.0*f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[54]-1.0*f[55])-3.0*f[86])+5.0*f[30])); 
  fReflXYQuad[6][11] = -0.1*(9.0*f[89]+6.708203932499369*(f[64]-1.0*f[63])-5.0*f[39]); 
  fReflXYQuad[6][12] = -0.1*(9.0*f[90]+6.708203932499369*(f[67]-1.0*f[66])-5.0*f[42]); 
  fReflXYQuad[6][13] = -0.1*(9.0*f[93]+6.708203932499369*(f[73]-1.0*f[72])-5.0*f[45]); 
  fReflXYQuad[6][14] = -0.1*(9.0*f[100]+6.708203932499369*(f[78]-1.0*f[77])-5.0*f[46]); 
  fReflXYQuad[6][15] = -0.1*(9.0*f[103]+6.708203932499369*(f[82]-1.0*f[81])-5.0*f[49]); 
  fReflXYQuad[6][16] = -0.1*(9.0*f[104]+6.708203932499369*(f[84]-1.0*f[83])-5.0*f[50]); 
  fReflXYQuad[6][17] = -0.1*(9.0*f[109]+6.708203932499369*(f[99]-1.0*f[98])-5.0*f[76]); 
  fReflXYQuad[6][18] = -0.1*(9.0*f[110]+6.708203932499369*(f[102]-1.0*f[101])-5.0*f[79]); 
  fReflXYQuad[6][19] = -0.1*(9.0*f[111]+6.708203932499369*(f[106]-1.0*f[105])-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*(f[78]-1.0*f[77]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*(f[78]-1.0*f[77]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*(f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*(f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[6][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[6][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_8 
  vcutSq_i = 0.025*(1.732050807568877*(21.21320343559643*phiWall[18]-21.21320343559643*phi[18]-14.14213562373095*phiWall[3]+14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[15]-2.828427124746191*phi[15]-1.414213562373095*phiWall[12]+1.414213562373095*phi[12])+7.071067811865476*phiWall[14]-7.071067811865476*phi[14]-5.656854249492382*phiWall[13]+5.656854249492382*phi[13])+2.23606797749979*(14.14213562373095*phiWall[9]-14.14213562373095*phi[9]-7.071067811865476*phiWall[8]+7.071067811865476*phi[8]+5.656854249492382*phiWall[7]-5.656854249492382*phi[7]+1.732050807568877*(8.485281374238571*phi[5]-8.485281374238571*phiWall[5])+8.485281374238571*phiWall[1]-8.485281374238571*phi[1])+14.14213562373095*phiWall[0]-14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[7][0] = 0.0; 
  fReflXYQuad[7][1] = 0.0; 
  fReflXYQuad[7][2] = 0.0; 
  fReflXYQuad[7][3] = 0.0; 
  fReflXYQuad[7][4] = 0.0; 
  fReflXYQuad[7][5] = 0.0; 
  fReflXYQuad[7][6] = 0.0; 
  fReflXYQuad[7][7] = 0.0; 
  fReflXYQuad[7][8] = 0.0; 
  fReflXYQuad[7][9] = 0.0; 
  fReflXYQuad[7][10] = 0.0; 
  fReflXYQuad[7][11] = 0.0; 
  fReflXYQuad[7][12] = 0.0; 
  fReflXYQuad[7][13] = 0.0; 
  fReflXYQuad[7][14] = 0.0; 
  fReflXYQuad[7][15] = 0.0; 
  fReflXYQuad[7][16] = 0.0; 
  fReflXYQuad[7][17] = 0.0; 
  fReflXYQuad[7][18] = 0.0; 
  fReflXYQuad[7][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[7][0] = -0.05*(2.23606797749979*(6.708203932499369*f[32]+5.0*f[17]-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]); 
  fReflXYQuad[7][1] = -0.01666666666666667*(45.0*f[57]+6.708203932499369*(5.0*f[34]-4.0*f[33])-6.0*(6.708203932499369*f[7]+5.0*f[3])); 
  fReflXYQuad[7][2] = -0.01666666666666667*(45.0*f[60]+6.708203932499369*(5.0*f[38]-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])); 
  fReflXYQuad[7][3] = -0.01666666666666667*(45.0*f[69]+6.708203932499369*(5.0*f[44]-4.0*f[43])-6.0*(6.708203932499369*f[12]+5.0*f[5])); 
  fReflXYQuad[7][4] = -0.05*(2.23606797749979*(6.708203932499369*f[88]+5.0*f[62]-2.0*(2.0*f[61]+3.0*f[23]))-10.0*f[11]); 
  fReflXYQuad[7][5] = -0.05*(2.23606797749979*(6.708203932499369*f[92]+5.0*f[71]-2.0*(2.0*f[70]+3.0*f[26]))-10.0*f[14]); 
  fReflXYQuad[7][6] = -0.05*(2.23606797749979*(6.708203932499369*f[95]+5.0*f[75]-2.0*(2.0*f[74]+3.0*f[28]))-10.0*f[15]); 
  fReflXYQuad[7][7] = 0.1*(6.708203932499369*f[35]+5.0*f[18]); 
  fReflXYQuad[7][8] = 0.1*(6.708203932499369*f[40]+5.0*f[19]); 
  fReflXYQuad[7][9] = 0.1*(6.708203932499369*f[47]+5.0*f[20]); 
  fReflXYQuad[7][10] = -0.01666666666666667*(45.0*f[108]+6.708203932499369*(5.0*f[97]-4.0*f[96])-6.0*(6.708203932499369*f[54]+5.0*f[30])); 
  fReflXYQuad[7][11] = 0.1*(6.708203932499369*f[63]+5.0*f[39]); 
  fReflXYQuad[7][12] = 0.1*(6.708203932499369*f[66]+5.0*f[42]); 
  fReflXYQuad[7][13] = 0.1*(6.708203932499369*f[72]+5.0*f[45]); 
  fReflXYQuad[7][14] = 0.1*(6.708203932499369*f[77]+5.0*f[46]); 
  fReflXYQuad[7][15] = 0.1*(6.708203932499369*f[81]+5.0*f[49]); 
  fReflXYQuad[7][16] = 0.1*(6.708203932499369*f[83]+5.0*f[50]); 
  fReflXYQuad[7][17] = 0.1*(6.708203932499369*f[98]+5.0*f[76]); 
  fReflXYQuad[7][18] = 0.1*(6.708203932499369*f[101]+5.0*f[79]); 
  fReflXYQuad[7][19] = 0.1*(6.708203932499369*f[105]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[42]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[42]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*(f[79]-1.0*f[77]+f[66])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*(f[79]-1.0*f[77]+f[66])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[46]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[46]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[40]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[40]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[79]+f[77])-1.0*(9.0*f[66]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[79]+f[77])-1.0*(9.0*f[66]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))/(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))+5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))+5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[7][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[7][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_9 
  vcutSq_i = 0.01*(2.23606797749979*(63.63961030678928*phiWall[19]-63.63961030678928*phi[19]+35.35533905932738*phiWall[9]-35.35533905932738*phi[9]+14.14213562373095*phiWall[8]-14.14213562373095*phi[8]+14.14213562373095*phiWall[7]-14.14213562373095*phi[7]+1.732050807568877*((-21.21320343559643*phiWall[6])+21.21320343559643*phi[6]-21.21320343559643*phiWall[5]+21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*phi[2]+21.21320343559643*phiWall[1]-21.21320343559643*phi[1])+1.732050807568877*((-42.42640687119286*phiWall[18])+42.42640687119286*phi[18]-42.42640687119286*phiWall[17]+42.42640687119286*phi[17]-63.63961030678928*phiWall[10]+63.63961030678928*phi[10]-35.35533905932738*phiWall[3]+35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*phi[16]+7.071067811865476*phiWall[15]-7.071067811865476*phi[15]+2.828427124746191*phiWall[12]-2.828427124746191*phi[12]+2.828427124746191*phiWall[11]-2.828427124746191*phi[11])-14.14213562373095*phiWall[14]+14.14213562373095*phi[14]-14.14213562373095*phiWall[13]+14.14213562373095*phi[13])+63.63961030678928*phiWall[4]-63.63961030678928*phi[4]+35.35533905932738*phiWall[0]-35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[8][0] = 0.0; 
  fReflXYQuad[8][1] = 0.0; 
  fReflXYQuad[8][2] = 0.0; 
  fReflXYQuad[8][3] = 0.0; 
  fReflXYQuad[8][4] = 0.0; 
  fReflXYQuad[8][5] = 0.0; 
  fReflXYQuad[8][6] = 0.0; 
  fReflXYQuad[8][7] = 0.0; 
  fReflXYQuad[8][8] = 0.0; 
  fReflXYQuad[8][9] = 0.0; 
  fReflXYQuad[8][10] = 0.0; 
  fReflXYQuad[8][11] = 0.0; 
  fReflXYQuad[8][12] = 0.0; 
  fReflXYQuad[8][13] = 0.0; 
  fReflXYQuad[8][14] = 0.0; 
  fReflXYQuad[8][15] = 0.0; 
  fReflXYQuad[8][16] = 0.0; 
  fReflXYQuad[8][17] = 0.0; 
  fReflXYQuad[8][18] = 0.0; 
  fReflXYQuad[8][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[8][0] = 0.02*(2.23606797749979*(13.41640786499874*(f[32]+f[31])+5.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+5.0*(9.0*f[6]+5.0*f[0])); 
  fReflXYQuad[8][1] = 0.03333333333333333*(2.0*(9.0*(f[57]+f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(3.0*f[21]+2.23606797749979*(f[8]+f[7]))+5.0*f[3])); 
  fReflXYQuad[8][2] = 0.03333333333333333*(2.0*(9.0*(f[60]+f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4])); 
  fReflXYQuad[8][3] = 0.03333333333333333*(2.0*(9.0*(f[69]+f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(3.0*f[25]+2.23606797749979*(f[13]+f[12]))+5.0*f[5])); 
  fReflXYQuad[8][4] = 0.02*(2.23606797749979*(13.41640786499874*(f[88]+f[87])+5.0*(2.0*(f[62]+f[61])+3.0*(f[24]+f[23])))+5.0*(9.0*f[51]+5.0*f[11])); 
  fReflXYQuad[8][5] = 0.02*(2.23606797749979*(13.41640786499874*(f[92]+f[91])+5.0*(2.0*(f[71]+f[70])+3.0*(f[27]+f[26])))+5.0*(9.0*f[52]+5.0*f[14])); 
  fReflXYQuad[8][6] = 0.02*(2.23606797749979*(13.41640786499874*(f[95]+f[94])+5.0*(2.0*(f[75]+f[74])+3.0*(f[29]+f[28])))+5.0*(9.0*f[53]+5.0*f[15])); 
  fReflXYQuad[8][7] = 0.1*(9.0*f[58]+6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[8][8] = 0.1*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[8][9] = 0.1*(9.0*f[80]+6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[8][10] = 0.03333333333333333*(2.0*(9.0*(f[108]+f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(3.0*f[86]+2.23606797749979*(f[55]+f[54]))+5.0*f[30])); 
  fReflXYQuad[8][11] = 0.1*(9.0*f[89]+6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[8][12] = 0.1*(9.0*f[90]+6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[8][13] = 0.1*(9.0*f[93]+6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[8][14] = 0.1*(9.0*f[100]+6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[8][15] = 0.1*(9.0*f[103]+6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[8][16] = 0.1*(9.0*f[104]+6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[8][17] = 0.1*(9.0*f[109]+6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[8][18] = 0.1*(9.0*f[110]+6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[8][19] = 0.1*(9.0*f[111]+6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[67]+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[67]+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[78]+f[77])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[78]+f[77])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))/(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))+5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))+5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))+5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))+5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[8][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[8][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
  fRefl[0] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0])+8.0*(5.0*fReflXYQuad[5][0]+8.0*fReflXYQuad[4][0])+5.0*(8.0*fReflXYQuad[3][0]+5.0*fReflXYQuad[2][0]+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0])); 
  fRefl[1] = 0.0414086662499961*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0]-1.0*(5.0*fReflXYQuad[2][0]+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0])); 
  fRefl[2] = 0.0414086662499961*(5.0*(fReflXYQuad[8][0]-1.0*fReflXYQuad[6][0])+8.0*(fReflXYQuad[5][0]-1.0*fReflXYQuad[3][0])+5.0*(fReflXYQuad[2][0]-1.0*fReflXYQuad[0][0])); 
  fRefl[3] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1])+8.0*(5.0*fReflXYQuad[5][1]+8.0*fReflXYQuad[4][1])+5.0*(8.0*fReflXYQuad[3][1]+5.0*fReflXYQuad[2][1]+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1])); 
  fRefl[4] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2])+8.0*(5.0*fReflXYQuad[5][2]+8.0*fReflXYQuad[4][2])+5.0*(8.0*fReflXYQuad[3][2]+5.0*fReflXYQuad[2][2]+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2])); 
  fRefl[5] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3])+8.0*(5.0*fReflXYQuad[5][3]+8.0*fReflXYQuad[4][3])+5.0*(8.0*fReflXYQuad[3][3]+5.0*fReflXYQuad[2][3]+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3])); 
  fRefl[6] = 0.2777777777777778*(fReflXYQuad[8][0]-1.0*(fReflXYQuad[6][0]+fReflXYQuad[2][0])+fReflXYQuad[0][0]); 
  fRefl[7] = 0.0414086662499961*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1]-1.0*(5.0*fReflXYQuad[2][1]+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1])); 
  fRefl[8] = 0.0414086662499961*(5.0*(fReflXYQuad[8][1]-1.0*fReflXYQuad[6][1])+8.0*(fReflXYQuad[5][1]-1.0*fReflXYQuad[3][1])+5.0*(fReflXYQuad[2][1]-1.0*fReflXYQuad[0][1])); 
  fRefl[9] = 0.0414086662499961*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2]-1.0*(5.0*fReflXYQuad[2][2]+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2])); 
  fRefl[10] = 0.0414086662499961*(5.0*(fReflXYQuad[8][2]-1.0*fReflXYQuad[6][2])+8.0*(fReflXYQuad[5][2]-1.0*fReflXYQuad[3][2])+5.0*(fReflXYQuad[2][2]-1.0*fReflXYQuad[0][2])); 
  fRefl[11] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4])+8.0*(5.0*fReflXYQuad[5][4]+8.0*fReflXYQuad[4][4])+5.0*(8.0*fReflXYQuad[3][4]+5.0*fReflXYQuad[2][4]+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4])); 
  fRefl[12] = 0.0414086662499961*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3]-1.0*(5.0*fReflXYQuad[2][3]+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3])); 
  fRefl[13] = 0.0414086662499961*(5.0*(fReflXYQuad[8][3]-1.0*fReflXYQuad[6][3])+8.0*(fReflXYQuad[5][3]-1.0*fReflXYQuad[3][3])+5.0*(fReflXYQuad[2][3]-1.0*fReflXYQuad[0][3])); 
  fRefl[14] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5])+8.0*(5.0*fReflXYQuad[5][5]+8.0*fReflXYQuad[4][5])+5.0*(8.0*fReflXYQuad[3][5]+5.0*fReflXYQuad[2][5]+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5])); 
  fRefl[15] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6])+8.0*(5.0*fReflXYQuad[5][6]+8.0*fReflXYQuad[4][6])+5.0*(8.0*fReflXYQuad[3][6]+5.0*fReflXYQuad[2][6]+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6])); 
  fRefl[16] = 0.0276057774999974*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0]-2.0*(5.0*fReflXYQuad[5][0]+8.0*fReflXYQuad[4][0])+5.0*(fReflXYQuad[2][0]-2.0*fReflXYQuad[3][0])+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0]); 
  fRefl[17] = 0.0276057774999974*(5.0*(fReflXYQuad[8][0]-2.0*fReflXYQuad[7][0]+fReflXYQuad[6][0])+8.0*(fReflXYQuad[5][0]-2.0*fReflXYQuad[4][0]+fReflXYQuad[3][0])+5.0*(fReflXYQuad[2][0]-2.0*fReflXYQuad[1][0]+fReflXYQuad[0][0])); 
  fRefl[18] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][7]+8.0*fReflXYQuad[7][7]+5.0*fReflXYQuad[6][7])+8.0*(5.0*fReflXYQuad[5][7]+8.0*fReflXYQuad[4][7])+5.0*(8.0*fReflXYQuad[3][7]+5.0*fReflXYQuad[2][7]+8.0*fReflXYQuad[1][7]+5.0*fReflXYQuad[0][7])); 
  fRefl[19] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][8]+8.0*fReflXYQuad[7][8]+5.0*fReflXYQuad[6][8])+8.0*(5.0*fReflXYQuad[5][8]+8.0*fReflXYQuad[4][8])+5.0*(8.0*fReflXYQuad[3][8]+5.0*fReflXYQuad[2][8]+8.0*fReflXYQuad[1][8]+5.0*fReflXYQuad[0][8])); 
  fRefl[20] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][9]+8.0*fReflXYQuad[7][9]+5.0*fReflXYQuad[6][9])+8.0*(5.0*fReflXYQuad[5][9]+8.0*fReflXYQuad[4][9])+5.0*(8.0*fReflXYQuad[3][9]+5.0*fReflXYQuad[2][9]+8.0*fReflXYQuad[1][9]+5.0*fReflXYQuad[0][9])); 
  fRefl[21] = 0.2777777777777778*(fReflXYQuad[8][1]-1.0*(fReflXYQuad[6][1]+fReflXYQuad[2][1])+fReflXYQuad[0][1]); 
  fRefl[22] = 0.2777777777777778*(fReflXYQuad[8][2]-1.0*(fReflXYQuad[6][2]+fReflXYQuad[2][2])+fReflXYQuad[0][2]); 
  fRefl[23] = 0.0414086662499961*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4]-1.0*(5.0*fReflXYQuad[2][4]+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4])); 
  fRefl[24] = 0.0414086662499961*(5.0*(fReflXYQuad[8][4]-1.0*fReflXYQuad[6][4])+8.0*(fReflXYQuad[5][4]-1.0*fReflXYQuad[3][4])+5.0*(fReflXYQuad[2][4]-1.0*fReflXYQuad[0][4])); 
  fRefl[25] = 0.2777777777777778*(fReflXYQuad[8][3]-1.0*(fReflXYQuad[6][3]+fReflXYQuad[2][3])+fReflXYQuad[0][3]); 
  fRefl[26] = 0.0414086662499961*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5]-1.0*(5.0*fReflXYQuad[2][5]+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5])); 
  fRefl[27] = 0.0414086662499961*(5.0*(fReflXYQuad[8][5]-1.0*fReflXYQuad[6][5])+8.0*(fReflXYQuad[5][5]-1.0*fReflXYQuad[3][5])+5.0*(fReflXYQuad[2][5]-1.0*fReflXYQuad[0][5])); 
  fRefl[28] = 0.0414086662499961*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6]-1.0*(5.0*fReflXYQuad[2][6]+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6])); 
  fRefl[29] = 0.0414086662499961*(5.0*(fReflXYQuad[8][6]-1.0*fReflXYQuad[6][6])+8.0*(fReflXYQuad[5][6]-1.0*fReflXYQuad[3][6])+5.0*(fReflXYQuad[2][6]-1.0*fReflXYQuad[0][6])); 
  fRefl[30] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10])+8.0*(5.0*fReflXYQuad[5][10]+8.0*fReflXYQuad[4][10])+5.0*(8.0*fReflXYQuad[3][10]+5.0*fReflXYQuad[2][10]+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10])); 
  fRefl[31] = 0.1851851851851853*(fReflXYQuad[8][0]-1.0*fReflXYQuad[6][0]+2.0*(fReflXYQuad[3][0]-1.0*fReflXYQuad[5][0])+fReflXYQuad[2][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[32] = 0.1851851851851853*(fReflXYQuad[8][0]-2.0*fReflXYQuad[7][0]+fReflXYQuad[6][0]-1.0*fReflXYQuad[2][0]+2.0*fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[33] = 0.02760577749999742*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1]-2.0*(5.0*fReflXYQuad[5][1]+8.0*fReflXYQuad[4][1])+5.0*(fReflXYQuad[2][1]-2.0*fReflXYQuad[3][1])+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1]); 
  fRefl[34] = 0.02760577749999742*(5.0*(fReflXYQuad[8][1]-2.0*fReflXYQuad[7][1]+fReflXYQuad[6][1])+8.0*(fReflXYQuad[5][1]-2.0*fReflXYQuad[4][1]+fReflXYQuad[3][1])+5.0*(fReflXYQuad[2][1]-2.0*fReflXYQuad[1][1]+fReflXYQuad[0][1])); 
  fRefl[35] = 0.04140866624999612*(5.0*fReflXYQuad[8][7]+8.0*fReflXYQuad[7][7]+5.0*fReflXYQuad[6][7]-1.0*(5.0*fReflXYQuad[2][7]+8.0*fReflXYQuad[1][7]+5.0*fReflXYQuad[0][7])); 
  fRefl[36] = 0.04140866624999612*(5.0*(fReflXYQuad[8][7]-1.0*fReflXYQuad[6][7])+8.0*(fReflXYQuad[5][7]-1.0*fReflXYQuad[3][7])+5.0*(fReflXYQuad[2][7]-1.0*fReflXYQuad[0][7])); 
  fRefl[37] = 0.02760577749999742*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2]-2.0*(5.0*fReflXYQuad[5][2]+8.0*fReflXYQuad[4][2])+5.0*(fReflXYQuad[2][2]-2.0*fReflXYQuad[3][2])+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2]); 
  fRefl[38] = 0.02760577749999742*(5.0*(fReflXYQuad[8][2]-2.0*fReflXYQuad[7][2]+fReflXYQuad[6][2])+8.0*(fReflXYQuad[5][2]-2.0*fReflXYQuad[4][2]+fReflXYQuad[3][2])+5.0*(fReflXYQuad[2][2]-2.0*fReflXYQuad[1][2]+fReflXYQuad[0][2])); 
  fRefl[39] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][11]+8.0*fReflXYQuad[7][11]+5.0*fReflXYQuad[6][11])+8.0*(5.0*fReflXYQuad[5][11]+8.0*fReflXYQuad[4][11])+5.0*(8.0*fReflXYQuad[3][11]+5.0*fReflXYQuad[2][11]+8.0*fReflXYQuad[1][11]+5.0*fReflXYQuad[0][11])); 
  fRefl[40] = 0.04140866624999612*(5.0*fReflXYQuad[8][8]+8.0*fReflXYQuad[7][8]+5.0*fReflXYQuad[6][8]-1.0*(5.0*fReflXYQuad[2][8]+8.0*fReflXYQuad[1][8]+5.0*fReflXYQuad[0][8])); 
  fRefl[41] = 0.04140866624999612*(5.0*(fReflXYQuad[8][8]-1.0*fReflXYQuad[6][8])+8.0*(fReflXYQuad[5][8]-1.0*fReflXYQuad[3][8])+5.0*(fReflXYQuad[2][8]-1.0*fReflXYQuad[0][8])); 
  fRefl[42] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][12]+8.0*fReflXYQuad[7][12]+5.0*fReflXYQuad[6][12])+8.0*(5.0*fReflXYQuad[5][12]+8.0*fReflXYQuad[4][12])+5.0*(8.0*fReflXYQuad[3][12]+5.0*fReflXYQuad[2][12]+8.0*fReflXYQuad[1][12]+5.0*fReflXYQuad[0][12])); 
  fRefl[43] = 0.02760577749999742*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3]-2.0*(5.0*fReflXYQuad[5][3]+8.0*fReflXYQuad[4][3])+5.0*(fReflXYQuad[2][3]-2.0*fReflXYQuad[3][3])+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3]); 
  fRefl[44] = 0.02760577749999742*(5.0*(fReflXYQuad[8][3]-2.0*fReflXYQuad[7][3]+fReflXYQuad[6][3])+8.0*(fReflXYQuad[5][3]-2.0*fReflXYQuad[4][3]+fReflXYQuad[3][3])+5.0*(fReflXYQuad[2][3]-2.0*fReflXYQuad[1][3]+fReflXYQuad[0][3])); 
  fRefl[45] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][13]+8.0*fReflXYQuad[7][13]+5.0*fReflXYQuad[6][13])+8.0*(5.0*fReflXYQuad[5][13]+8.0*fReflXYQuad[4][13])+5.0*(8.0*fReflXYQuad[3][13]+5.0*fReflXYQuad[2][13]+8.0*fReflXYQuad[1][13]+5.0*fReflXYQuad[0][13])); 
  fRefl[46] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][14]+8.0*fReflXYQuad[7][14]+5.0*fReflXYQuad[6][14])+8.0*(5.0*fReflXYQuad[5][14]+8.0*fReflXYQuad[4][14])+5.0*(8.0*fReflXYQuad[3][14]+5.0*fReflXYQuad[2][14]+8.0*fReflXYQuad[1][14]+5.0*fReflXYQuad[0][14])); 
  fRefl[47] = 0.04140866624999612*(5.0*fReflXYQuad[8][9]+8.0*fReflXYQuad[7][9]+5.0*fReflXYQuad[6][9]-1.0*(5.0*fReflXYQuad[2][9]+8.0*fReflXYQuad[1][9]+5.0*fReflXYQuad[0][9])); 
  fRefl[48] = 0.04140866624999612*(5.0*(fReflXYQuad[8][9]-1.0*fReflXYQuad[6][9])+8.0*(fReflXYQuad[5][9]-1.0*fReflXYQuad[3][9])+5.0*(fReflXYQuad[2][9]-1.0*fReflXYQuad[0][9])); 
  fRefl[49] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][15]+8.0*fReflXYQuad[7][15]+5.0*fReflXYQuad[6][15])+8.0*(5.0*fReflXYQuad[5][15]+8.0*fReflXYQuad[4][15])+5.0*(8.0*fReflXYQuad[3][15]+5.0*fReflXYQuad[2][15]+8.0*fReflXYQuad[1][15]+5.0*fReflXYQuad[0][15])); 
  fRefl[50] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][16]+8.0*fReflXYQuad[7][16]+5.0*fReflXYQuad[6][16])+8.0*(5.0*fReflXYQuad[5][16]+8.0*fReflXYQuad[4][16])+5.0*(8.0*fReflXYQuad[3][16]+5.0*fReflXYQuad[2][16]+8.0*fReflXYQuad[1][16]+5.0*fReflXYQuad[0][16])); 
  fRefl[51] = 0.2777777777777778*(fReflXYQuad[8][4]-1.0*(fReflXYQuad[6][4]+fReflXYQuad[2][4])+fReflXYQuad[0][4]); 
  fRefl[52] = 0.2777777777777778*(fReflXYQuad[8][5]-1.0*(fReflXYQuad[6][5]+fReflXYQuad[2][5])+fReflXYQuad[0][5]); 
  fRefl[53] = 0.2777777777777778*(fReflXYQuad[8][6]-1.0*(fReflXYQuad[6][6]+fReflXYQuad[2][6])+fReflXYQuad[0][6]); 
  fRefl[54] = 0.0414086662499961*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10]-1.0*(5.0*fReflXYQuad[2][10]+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10])); 
  fRefl[55] = 0.0414086662499961*(5.0*(fReflXYQuad[8][10]-1.0*fReflXYQuad[6][10])+8.0*(fReflXYQuad[5][10]-1.0*fReflXYQuad[3][10])+5.0*(fReflXYQuad[2][10]-1.0*fReflXYQuad[0][10])); 
  fRefl[56] = 0.1851851851851852*(fReflXYQuad[8][1]-1.0*fReflXYQuad[6][1]+2.0*(fReflXYQuad[3][1]-1.0*fReflXYQuad[5][1])+fReflXYQuad[2][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[57] = 0.1851851851851852*(fReflXYQuad[8][1]-2.0*fReflXYQuad[7][1]+fReflXYQuad[6][1]-1.0*fReflXYQuad[2][1]+2.0*fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[58] = 0.2777777777777778*(fReflXYQuad[8][7]-1.0*(fReflXYQuad[6][7]+fReflXYQuad[2][7])+fReflXYQuad[0][7]); 
  fRefl[59] = 0.1851851851851852*(fReflXYQuad[8][2]-1.0*fReflXYQuad[6][2]+2.0*(fReflXYQuad[3][2]-1.0*fReflXYQuad[5][2])+fReflXYQuad[2][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[60] = 0.1851851851851852*(fReflXYQuad[8][2]-2.0*fReflXYQuad[7][2]+fReflXYQuad[6][2]-1.0*fReflXYQuad[2][2]+2.0*fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[61] = 0.0276057774999974*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4]-2.0*(5.0*fReflXYQuad[5][4]+8.0*fReflXYQuad[4][4])+5.0*(fReflXYQuad[2][4]-2.0*fReflXYQuad[3][4])+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4]); 
  fRefl[62] = 0.0276057774999974*(5.0*(fReflXYQuad[8][4]-2.0*fReflXYQuad[7][4]+fReflXYQuad[6][4])+8.0*(fReflXYQuad[5][4]-2.0*fReflXYQuad[4][4]+fReflXYQuad[3][4])+5.0*(fReflXYQuad[2][4]-2.0*fReflXYQuad[1][4]+fReflXYQuad[0][4])); 
  fRefl[63] = 0.04140866624999612*(5.0*fReflXYQuad[8][11]+8.0*fReflXYQuad[7][11]+5.0*fReflXYQuad[6][11]-1.0*(5.0*fReflXYQuad[2][11]+8.0*fReflXYQuad[1][11]+5.0*fReflXYQuad[0][11])); 
  fRefl[64] = 0.04140866624999612*(5.0*(fReflXYQuad[8][11]-1.0*fReflXYQuad[6][11])+8.0*(fReflXYQuad[5][11]-1.0*fReflXYQuad[3][11])+5.0*(fReflXYQuad[2][11]-1.0*fReflXYQuad[0][11])); 
  fRefl[65] = 0.2777777777777778*(fReflXYQuad[8][8]-1.0*(fReflXYQuad[6][8]+fReflXYQuad[2][8])+fReflXYQuad[0][8]); 
  fRefl[66] = 0.04140866624999612*(5.0*fReflXYQuad[8][12]+8.0*fReflXYQuad[7][12]+5.0*fReflXYQuad[6][12]-1.0*(5.0*fReflXYQuad[2][12]+8.0*fReflXYQuad[1][12]+5.0*fReflXYQuad[0][12])); 
  fRefl[67] = 0.04140866624999612*(5.0*(fReflXYQuad[8][12]-1.0*fReflXYQuad[6][12])+8.0*(fReflXYQuad[5][12]-1.0*fReflXYQuad[3][12])+5.0*(fReflXYQuad[2][12]-1.0*fReflXYQuad[0][12])); 
  fRefl[68] = 0.1851851851851852*(fReflXYQuad[8][3]-1.0*fReflXYQuad[6][3]+2.0*(fReflXYQuad[3][3]-1.0*fReflXYQuad[5][3])+fReflXYQuad[2][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[69] = 0.1851851851851852*(fReflXYQuad[8][3]-2.0*fReflXYQuad[7][3]+fReflXYQuad[6][3]-1.0*fReflXYQuad[2][3]+2.0*fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[70] = 0.0276057774999974*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5]-2.0*(5.0*fReflXYQuad[5][5]+8.0*fReflXYQuad[4][5])+5.0*(fReflXYQuad[2][5]-2.0*fReflXYQuad[3][5])+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5]); 
  fRefl[71] = 0.0276057774999974*(5.0*(fReflXYQuad[8][5]-2.0*fReflXYQuad[7][5]+fReflXYQuad[6][5])+8.0*(fReflXYQuad[5][5]-2.0*fReflXYQuad[4][5]+fReflXYQuad[3][5])+5.0*(fReflXYQuad[2][5]-2.0*fReflXYQuad[1][5]+fReflXYQuad[0][5])); 
  fRefl[72] = 0.04140866624999612*(5.0*fReflXYQuad[8][13]+8.0*fReflXYQuad[7][13]+5.0*fReflXYQuad[6][13]-1.0*(5.0*fReflXYQuad[2][13]+8.0*fReflXYQuad[1][13]+5.0*fReflXYQuad[0][13])); 
  fRefl[73] = 0.04140866624999612*(5.0*(fReflXYQuad[8][13]-1.0*fReflXYQuad[6][13])+8.0*(fReflXYQuad[5][13]-1.0*fReflXYQuad[3][13])+5.0*(fReflXYQuad[2][13]-1.0*fReflXYQuad[0][13])); 
  fRefl[74] = 0.0276057774999974*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6]-2.0*(5.0*fReflXYQuad[5][6]+8.0*fReflXYQuad[4][6])+5.0*(fReflXYQuad[2][6]-2.0*fReflXYQuad[3][6])+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6]); 
  fRefl[75] = 0.0276057774999974*(5.0*(fReflXYQuad[8][6]-2.0*fReflXYQuad[7][6]+fReflXYQuad[6][6])+8.0*(fReflXYQuad[5][6]-2.0*fReflXYQuad[4][6]+fReflXYQuad[3][6])+5.0*(fReflXYQuad[2][6]-2.0*fReflXYQuad[1][6]+fReflXYQuad[0][6])); 
  fRefl[76] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][17]+8.0*fReflXYQuad[7][17]+5.0*fReflXYQuad[6][17])+8.0*(5.0*fReflXYQuad[5][17]+8.0*fReflXYQuad[4][17])+5.0*(8.0*fReflXYQuad[3][17]+5.0*fReflXYQuad[2][17]+8.0*fReflXYQuad[1][17]+5.0*fReflXYQuad[0][17])); 
  fRefl[77] = 0.04140866624999612*(5.0*fReflXYQuad[8][14]+8.0*fReflXYQuad[7][14]+5.0*fReflXYQuad[6][14]-1.0*(5.0*fReflXYQuad[2][14]+8.0*fReflXYQuad[1][14]+5.0*fReflXYQuad[0][14])); 
  fRefl[78] = 0.04140866624999612*(5.0*(fReflXYQuad[8][14]-1.0*fReflXYQuad[6][14])+8.0*(fReflXYQuad[5][14]-1.0*fReflXYQuad[3][14])+5.0*(fReflXYQuad[2][14]-1.0*fReflXYQuad[0][14])); 
  fRefl[79] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][18]+8.0*fReflXYQuad[7][18]+5.0*fReflXYQuad[6][18])+8.0*(5.0*fReflXYQuad[5][18]+8.0*fReflXYQuad[4][18])+5.0*(8.0*fReflXYQuad[3][18]+5.0*fReflXYQuad[2][18]+8.0*fReflXYQuad[1][18]+5.0*fReflXYQuad[0][18])); 
  fRefl[80] = 0.2777777777777778*(fReflXYQuad[8][9]-1.0*(fReflXYQuad[6][9]+fReflXYQuad[2][9])+fReflXYQuad[0][9]); 
  fRefl[81] = 0.04140866624999612*(5.0*fReflXYQuad[8][15]+8.0*fReflXYQuad[7][15]+5.0*fReflXYQuad[6][15]-1.0*(5.0*fReflXYQuad[2][15]+8.0*fReflXYQuad[1][15]+5.0*fReflXYQuad[0][15])); 
  fRefl[82] = 0.04140866624999612*(5.0*(fReflXYQuad[8][15]-1.0*fReflXYQuad[6][15])+8.0*(fReflXYQuad[5][15]-1.0*fReflXYQuad[3][15])+5.0*(fReflXYQuad[2][15]-1.0*fReflXYQuad[0][15])); 
  fRefl[83] = 0.04140866624999612*(5.0*fReflXYQuad[8][16]+8.0*fReflXYQuad[7][16]+5.0*fReflXYQuad[6][16]-1.0*(5.0*fReflXYQuad[2][16]+8.0*fReflXYQuad[1][16]+5.0*fReflXYQuad[0][16])); 
  fRefl[84] = 0.04140866624999612*(5.0*(fReflXYQuad[8][16]-1.0*fReflXYQuad[6][16])+8.0*(fReflXYQuad[5][16]-1.0*fReflXYQuad[3][16])+5.0*(fReflXYQuad[2][16]-1.0*fReflXYQuad[0][16])); 
  fRefl[85] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][19]+8.0*fReflXYQuad[7][19]+5.0*fReflXYQuad[6][19])+8.0*(5.0*fReflXYQuad[5][19]+8.0*fReflXYQuad[4][19])+5.0*(8.0*fReflXYQuad[3][19]+5.0*fReflXYQuad[2][19]+8.0*fReflXYQuad[1][19]+5.0*fReflXYQuad[0][19])); 
  fRefl[86] = 0.2777777777777778*(fReflXYQuad[8][10]-1.0*(fReflXYQuad[6][10]+fReflXYQuad[2][10])+fReflXYQuad[0][10]); 
  fRefl[87] = 0.1851851851851853*(fReflXYQuad[8][4]-1.0*fReflXYQuad[6][4]+2.0*(fReflXYQuad[3][4]-1.0*fReflXYQuad[5][4])+fReflXYQuad[2][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[88] = 0.1851851851851853*(fReflXYQuad[8][4]-2.0*fReflXYQuad[7][4]+fReflXYQuad[6][4]-1.0*fReflXYQuad[2][4]+2.0*fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[89] = 0.2777777777777778*(fReflXYQuad[8][11]-1.0*(fReflXYQuad[6][11]+fReflXYQuad[2][11])+fReflXYQuad[0][11]); 
  fRefl[90] = 0.2777777777777778*(fReflXYQuad[8][12]-1.0*(fReflXYQuad[6][12]+fReflXYQuad[2][12])+fReflXYQuad[0][12]); 
  fRefl[91] = 0.1851851851851853*(fReflXYQuad[8][5]-1.0*fReflXYQuad[6][5]+2.0*(fReflXYQuad[3][5]-1.0*fReflXYQuad[5][5])+fReflXYQuad[2][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[92] = 0.1851851851851853*(fReflXYQuad[8][5]-2.0*fReflXYQuad[7][5]+fReflXYQuad[6][5]-1.0*fReflXYQuad[2][5]+2.0*fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[93] = 0.2777777777777778*(fReflXYQuad[8][13]-1.0*(fReflXYQuad[6][13]+fReflXYQuad[2][13])+fReflXYQuad[0][13]); 
  fRefl[94] = 0.1851851851851853*(fReflXYQuad[8][6]-1.0*fReflXYQuad[6][6]+2.0*(fReflXYQuad[3][6]-1.0*fReflXYQuad[5][6])+fReflXYQuad[2][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[95] = 0.1851851851851853*(fReflXYQuad[8][6]-2.0*fReflXYQuad[7][6]+fReflXYQuad[6][6]-1.0*fReflXYQuad[2][6]+2.0*fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[96] = 0.02760577749999742*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10]-2.0*(5.0*fReflXYQuad[5][10]+8.0*fReflXYQuad[4][10])+5.0*(fReflXYQuad[2][10]-2.0*fReflXYQuad[3][10])+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10]); 
  fRefl[97] = 0.02760577749999742*(5.0*(fReflXYQuad[8][10]-2.0*fReflXYQuad[7][10]+fReflXYQuad[6][10])+8.0*(fReflXYQuad[5][10]-2.0*fReflXYQuad[4][10]+fReflXYQuad[3][10])+5.0*(fReflXYQuad[2][10]-2.0*fReflXYQuad[1][10]+fReflXYQuad[0][10])); 
  fRefl[98] = 0.04140866624999612*(5.0*fReflXYQuad[8][17]+8.0*fReflXYQuad[7][17]+5.0*fReflXYQuad[6][17]-1.0*(5.0*fReflXYQuad[2][17]+8.0*fReflXYQuad[1][17]+5.0*fReflXYQuad[0][17])); 
  fRefl[99] = 0.04140866624999612*(5.0*(fReflXYQuad[8][17]-1.0*fReflXYQuad[6][17])+8.0*(fReflXYQuad[5][17]-1.0*fReflXYQuad[3][17])+5.0*(fReflXYQuad[2][17]-1.0*fReflXYQuad[0][17])); 
  fRefl[100] = 0.2777777777777778*(fReflXYQuad[8][14]-1.0*(fReflXYQuad[6][14]+fReflXYQuad[2][14])+fReflXYQuad[0][14]); 
  fRefl[101] = 0.04140866624999612*(5.0*fReflXYQuad[8][18]+8.0*fReflXYQuad[7][18]+5.0*fReflXYQuad[6][18]-1.0*(5.0*fReflXYQuad[2][18]+8.0*fReflXYQuad[1][18]+5.0*fReflXYQuad[0][18])); 
  fRefl[102] = 0.04140866624999612*(5.0*(fReflXYQuad[8][18]-1.0*fReflXYQuad[6][18])+8.0*(fReflXYQuad[5][18]-1.0*fReflXYQuad[3][18])+5.0*(fReflXYQuad[2][18]-1.0*fReflXYQuad[0][18])); 
  fRefl[103] = 0.2777777777777778*(fReflXYQuad[8][15]-1.0*(fReflXYQuad[6][15]+fReflXYQuad[2][15])+fReflXYQuad[0][15]); 
  fRefl[104] = 0.2777777777777778*(fReflXYQuad[8][16]-1.0*(fReflXYQuad[6][16]+fReflXYQuad[2][16])+fReflXYQuad[0][16]); 
  fRefl[105] = 0.04140866624999612*(5.0*fReflXYQuad[8][19]+8.0*fReflXYQuad[7][19]+5.0*fReflXYQuad[6][19]-1.0*(5.0*fReflXYQuad[2][19]+8.0*fReflXYQuad[1][19]+5.0*fReflXYQuad[0][19])); 
  fRefl[106] = 0.04140866624999612*(5.0*(fReflXYQuad[8][19]-1.0*fReflXYQuad[6][19])+8.0*(fReflXYQuad[5][19]-1.0*fReflXYQuad[3][19])+5.0*(fReflXYQuad[2][19]-1.0*fReflXYQuad[0][19])); 
  fRefl[107] = 0.1851851851851852*(fReflXYQuad[8][10]-1.0*fReflXYQuad[6][10]+2.0*(fReflXYQuad[3][10]-1.0*fReflXYQuad[5][10])+fReflXYQuad[2][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[108] = 0.1851851851851852*(fReflXYQuad[8][10]-2.0*fReflXYQuad[7][10]+fReflXYQuad[6][10]-1.0*fReflXYQuad[2][10]+2.0*fReflXYQuad[1][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[109] = 0.2777777777777778*(fReflXYQuad[8][17]-1.0*(fReflXYQuad[6][17]+fReflXYQuad[2][17])+fReflXYQuad[0][17]); 
  fRefl[110] = 0.2777777777777778*(fReflXYQuad[8][18]-1.0*(fReflXYQuad[6][18]+fReflXYQuad[2][18])+fReflXYQuad[0][18]); 
  fRefl[111] = 0.2777777777777778*(fReflXYQuad[8][19]-1.0*(fReflXYQuad[6][19]+fReflXYQuad[2][19])+fReflXYQuad[0][19]); 
}

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; double xc, b, xbarVal, fac; 
  double fReflXYQuad[9][20]; 
  double fReflXYZMuQuad[9][8]; 
  

// node (x,y)_1 
  vcutSq_i = -0.01*(2.23606797749979*((-63.63961030678928*phiWall[19])+63.63961030678928*phi[19]-35.35533905932738*phiWall[9]+35.35533905932738*phi[9]-14.14213562373095*phiWall[8]+14.14213562373095*phi[8]-14.14213562373095*phiWall[7]+14.14213562373095*phi[7]+1.732050807568877*(21.21320343559643*phiWall[6]-21.21320343559643*phi[6]+21.21320343559643*phiWall[5]-21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*phi[2]+21.21320343559643*phiWall[1]-21.21320343559643*phi[1])+1.732050807568877*(42.42640687119286*phiWall[18]-42.42640687119286*phi[18]+42.42640687119286*phiWall[17]-42.42640687119286*phi[17]-63.63961030678928*phiWall[10]+63.63961030678928*phi[10]-35.35533905932738*phiWall[3]+35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*phi[16]+7.071067811865476*phiWall[15]-7.071067811865476*phi[15]+2.828427124746191*phiWall[12]-2.828427124746191*phi[12]+2.828427124746191*phiWall[11]-2.828427124746191*phi[11])-14.14213562373095*phiWall[14]+14.14213562373095*phi[14]-14.14213562373095*phiWall[13]+14.14213562373095*phi[13])-63.63961030678928*phiWall[4]+63.63961030678928*phi[4]-35.35533905932738*phiWall[0]+35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[0][0] = 0.0; 
  fReflXYQuad[0][1] = 0.0; 
  fReflXYQuad[0][2] = 0.0; 
  fReflXYQuad[0][3] = 0.0; 
  fReflXYQuad[0][4] = 0.0; 
  fReflXYQuad[0][5] = 0.0; 
  fReflXYQuad[0][6] = 0.0; 
  fReflXYQuad[0][7] = 0.0; 
  fReflXYQuad[0][8] = 0.0; 
  fReflXYQuad[0][9] = 0.0; 
  fReflXYQuad[0][10] = 0.0; 
  fReflXYQuad[0][11] = 0.0; 
  fReflXYQuad[0][12] = 0.0; 
  fReflXYQuad[0][13] = 0.0; 
  fReflXYQuad[0][14] = 0.0; 
  fReflXYQuad[0][15] = 0.0; 
  fReflXYQuad[0][16] = 0.0; 
  fReflXYQuad[0][17] = 0.0; 
  fReflXYQuad[0][18] = 0.0; 
  fReflXYQuad[0][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = -0.02*(30.0*(f[32]+f[31])-1.0*(11.18033988749895*(2.0*(f[17]+f[16])-3.0*(f[2]+f[1]))+5.0*(9.0*f[6]+5.0*f[0]))); 
  fReflXYQuad[0][1] = -0.03333333333333333*(2.0*(9.0*(f[57]+f[56])-6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[8]+f[7])-3.0*f[21])-5.0*f[3])); 
  fReflXYQuad[0][2] = -0.03333333333333333*(2.0*(9.0*(f[60]+f[59])-6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[10]+f[9])-3.0*f[22])-5.0*f[4])); 
  fReflXYQuad[0][3] = -0.03333333333333333*(2.0*(9.0*(f[69]+f[68])-6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[13]+f[12])-3.0*f[25])-5.0*f[5])); 
  fReflXYQuad[0][4] = -0.02*(30.0*(f[88]+f[87])-1.0*(11.18033988749895*(2.0*(f[62]+f[61])-3.0*(f[24]+f[23]))+5.0*(9.0*f[51]+5.0*f[11]))); 
  fReflXYQuad[0][5] = -0.02*(30.0*(f[92]+f[91])-1.0*(11.18033988749895*(2.0*(f[71]+f[70])-3.0*(f[27]+f[26]))+5.0*(9.0*f[52]+5.0*f[14]))); 
  fReflXYQuad[0][6] = -0.02*(30.0*(f[95]+f[94])-1.0*(11.18033988749895*(2.0*(f[75]+f[74])-3.0*(f[29]+f[28]))+5.0*(9.0*f[53]+5.0*f[15]))); 
  fReflXYQuad[0][7] = 0.1*(9.0*f[58]-6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[0][8] = 0.1*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[0][9] = 0.1*(9.0*f[80]-6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[0][10] = -0.03333333333333333*(2.0*(9.0*(f[108]+f[107])-6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[55]+f[54])-3.0*f[86])-5.0*f[30])); 
  fReflXYQuad[0][11] = 0.1*(9.0*f[89]-6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[0][12] = 0.1*(9.0*f[90]-6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[0][13] = 0.1*(9.0*f[93]-6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[0][14] = 0.1*(9.0*f[100]-6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[0][15] = 0.1*(9.0*f[103]-6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[0][16] = 0.1*(9.0*f[104]-6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[0][17] = 0.1*(9.0*f[109]-6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[0][18] = 0.1*(9.0*f[110]-6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[0][19] = 0.1*(9.0*f[111]-6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))/(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(30.0*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(11.18033988749895*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])-3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]+f[1])))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))-5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))-5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]+f[81]+f[80])+2.0*(3.0*f[21]-2.0*(f[58]+f[57]+f[56])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]+f[1])-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])+2.0*(5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])-3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))-5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))-5.0*(9.0*(2.0*(f[82]+f[81]+f[80]-1.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56])))+3.0*(f[27]+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]+f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12])-1.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*f[111]-1.0*(81.0*(f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*(f[30]+f[29]+f[28])-1.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))-5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])-3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))+2.0*(5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])-3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))-5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))/(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(11.18033988749895*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])+2.0*(3.0*(f[2]+f[1])-2.0*(f[17]+f[16])))-1.0*(15.0*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31]))+10.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(2.0*(3.0*(2.0*(f[60]+f[59])-3.0*f[22]+2.23606797749979*(f[10]+f[9]))-5.0*f[4])-15.0*(f[84]+f[83]+f[64]+f[63]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]-6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]-1.0*(f[2]+f[1]))))))+10.0*(9.0*(f[13]+f[12])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[29]+f[28])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]-1.0*f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])-6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*(f[58]-1.0*(f[57]+f[56]))+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*(f[2]+f[1])))-5.0*f[20])))+10.0*(9.0*(f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109])-1.0*(81.0*(f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_2 
  vcutSq_i = -0.025*(1.732050807568877*((-21.21320343559643*phiWall[18])+21.21320343559643*phi[18]-14.14213562373095*phiWall[3]+14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[15]-2.828427124746191*phi[15]-1.414213562373095*phiWall[12]+1.414213562373095*phi[12])+7.071067811865476*phiWall[14]-7.071067811865476*phi[14]-5.656854249492382*phiWall[13]+5.656854249492382*phi[13])+2.23606797749979*((-14.14213562373095*phiWall[9])+14.14213562373095*phi[9]+7.071067811865476*phiWall[8]-7.071067811865476*phi[8]-5.656854249492382*phiWall[7]+5.656854249492382*phi[7]+1.732050807568877*(8.485281374238571*phiWall[5]-8.485281374238571*phi[5])+8.485281374238571*phiWall[1]-8.485281374238571*phi[1])-14.14213562373095*phiWall[0]+14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[1][0] = 0.0; 
  fReflXYQuad[1][1] = 0.0; 
  fReflXYQuad[1][2] = 0.0; 
  fReflXYQuad[1][3] = 0.0; 
  fReflXYQuad[1][4] = 0.0; 
  fReflXYQuad[1][5] = 0.0; 
  fReflXYQuad[1][6] = 0.0; 
  fReflXYQuad[1][7] = 0.0; 
  fReflXYQuad[1][8] = 0.0; 
  fReflXYQuad[1][9] = 0.0; 
  fReflXYQuad[1][10] = 0.0; 
  fReflXYQuad[1][11] = 0.0; 
  fReflXYQuad[1][12] = 0.0; 
  fReflXYQuad[1][13] = 0.0; 
  fReflXYQuad[1][14] = 0.0; 
  fReflXYQuad[1][15] = 0.0; 
  fReflXYQuad[1][16] = 0.0; 
  fReflXYQuad[1][17] = 0.0; 
  fReflXYQuad[1][18] = 0.0; 
  fReflXYQuad[1][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = 0.05*(2.23606797749979*(6.708203932499369*f[32]-1.0*(5.0*f[17]+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]); 
  fReflXYQuad[1][1] = 0.01666666666666667*(45.0*f[57]-6.708203932499369*(5.0*f[34]-4.0*f[33])+6.0*(5.0*f[3]-6.708203932499369*f[7])); 
  fReflXYQuad[1][2] = 0.01666666666666667*(45.0*f[60]-6.708203932499369*(5.0*f[38]-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])); 
  fReflXYQuad[1][3] = 0.01666666666666667*(45.0*f[69]-6.708203932499369*(5.0*f[44]-4.0*f[43])+6.0*(5.0*f[5]-6.708203932499369*f[12])); 
  fReflXYQuad[1][4] = 0.05*(2.23606797749979*(6.708203932499369*f[88]-1.0*(5.0*f[62]+2.0*(3.0*f[23]-2.0*f[61])))+10.0*f[11]); 
  fReflXYQuad[1][5] = 0.05*(2.23606797749979*(6.708203932499369*f[92]-1.0*(5.0*f[71]+2.0*(3.0*f[26]-2.0*f[70])))+10.0*f[14]); 
  fReflXYQuad[1][6] = 0.05*(2.23606797749979*(6.708203932499369*f[95]-1.0*(5.0*f[75]+2.0*(3.0*f[28]-2.0*f[74])))+10.0*f[15]); 
  fReflXYQuad[1][7] = -0.1*(6.708203932499369*f[35]-5.0*f[18]); 
  fReflXYQuad[1][8] = -0.1*(6.708203932499369*f[40]-5.0*f[19]); 
  fReflXYQuad[1][9] = -0.1*(6.708203932499369*f[47]-5.0*f[20]); 
  fReflXYQuad[1][10] = 0.01666666666666667*(45.0*f[108]-6.708203932499369*(5.0*f[97]-4.0*f[96])+6.0*(5.0*f[30]-6.708203932499369*f[54])); 
  fReflXYQuad[1][11] = -0.1*(6.708203932499369*f[63]-5.0*f[39]); 
  fReflXYQuad[1][12] = -0.1*(6.708203932499369*f[66]-5.0*f[42]); 
  fReflXYQuad[1][13] = -0.1*(6.708203932499369*f[72]-5.0*f[45]); 
  fReflXYQuad[1][14] = -0.1*(6.708203932499369*f[77]-5.0*f[46]); 
  fReflXYQuad[1][15] = -0.1*(6.708203932499369*f[81]-5.0*f[49]); 
  fReflXYQuad[1][16] = -0.1*(6.708203932499369*f[83]-5.0*f[50]); 
  fReflXYQuad[1][17] = -0.1*(6.708203932499369*f[98]-5.0*f[76]); 
  fReflXYQuad[1][18] = -0.1*(6.708203932499369*f[101]-5.0*f[79]); 
  fReflXYQuad[1][19] = -0.1*(6.708203932499369*f[105]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))-5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-1.0*(5.0*(f[69]+f[57])+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]-3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))-5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))/(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]+2.0*(3.0*(f[3]+f[1])-2.0*f[16]))-1.0*(15.0*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32])+10.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])-1.0*(9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]-1.0*(f[3]+f[1]))))))+10.0*(9.0*(f[14]+f[12])-1.0*(9.0*f[7]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[30]+f[28])-1.0*(9.0*f[23]+5.0*f[4]))-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])+2.0*(3.0*(f[5]+f[1])-2.0*f[16]))+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32])))-10.0*(9.0*f[12]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))/(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[40]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])-1.0*(5.0*(f[20]+f[18]+f[17])+2.0*(3.0*f[1]-2.0*f[16])))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[83]+f[63]+f[60])-6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[9])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[40]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[77]-1.0*(6.708203932499369*(f[46]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[5]-1.0*f[1])))-6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[77]-1.0*(6.708203932499369*(f[46]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]-1.0*f[3]+f[1])-2.0*f[16])))+10.0*(9.0*(f[14]-1.0*f[12]+f[7])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*(f[30]-1.0*f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[66]-1.0*(6.708203932499369*(f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[3]-1.0*f[1])))-6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[66]-1.0*(6.708203932499369*(f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])-1.0*(9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])-6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[1]-1.0*(f[5]+f[3]))-2.0*f[16]))))+10.0*(9.0*(f[14]-1.0*(f[12]+f[7]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*(f[30]-1.0*(f[28]+f[23]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_3 
  vcutSq_i = 0.01*(2.23606797749979*((-63.63961030678928*phiWall[19])+63.63961030678928*phi[19]+35.35533905932738*phiWall[9]-35.35533905932738*phi[9]+14.14213562373095*phiWall[8]-14.14213562373095*phi[8]+14.14213562373095*phiWall[7]-14.14213562373095*phi[7]+1.732050807568877*(21.21320343559643*phiWall[6]-21.21320343559643*(phi[6]+phiWall[5])+21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*(phi[2]+phiWall[1])+21.21320343559643*phi[1])+1.732050807568877*((-42.42640687119286*phiWall[18])+42.42640687119286*(phi[18]+phiWall[17])-42.42640687119286*phi[17]-63.63961030678928*phiWall[10]+63.63961030678928*phi[10]+35.35533905932738*phiWall[3]-35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*(phi[16]+phiWall[15])+7.071067811865476*phi[15]-2.828427124746191*phiWall[12]+2.828427124746191*(phi[12]+phiWall[11])-2.828427124746191*phi[11])+14.14213562373095*phiWall[14]-14.14213562373095*phi[14]+14.14213562373095*phiWall[13]-14.14213562373095*phi[13])-63.63961030678928*phiWall[4]+63.63961030678928*phi[4]+35.35533905932738*phiWall[0]-35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[2][0] = 0.0; 
  fReflXYQuad[2][1] = 0.0; 
  fReflXYQuad[2][2] = 0.0; 
  fReflXYQuad[2][3] = 0.0; 
  fReflXYQuad[2][4] = 0.0; 
  fReflXYQuad[2][5] = 0.0; 
  fReflXYQuad[2][6] = 0.0; 
  fReflXYQuad[2][7] = 0.0; 
  fReflXYQuad[2][8] = 0.0; 
  fReflXYQuad[2][9] = 0.0; 
  fReflXYQuad[2][10] = 0.0; 
  fReflXYQuad[2][11] = 0.0; 
  fReflXYQuad[2][12] = 0.0; 
  fReflXYQuad[2][13] = 0.0; 
  fReflXYQuad[2][14] = 0.0; 
  fReflXYQuad[2][15] = 0.0; 
  fReflXYQuad[2][16] = 0.0; 
  fReflXYQuad[2][17] = 0.0; 
  fReflXYQuad[2][18] = 0.0; 
  fReflXYQuad[2][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.02*(2.23606797749979*(13.41640786499874*(f[32]-1.0*f[31])-5.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))+5.0*(9.0*f[6]-5.0*f[0])); 
  fReflXYQuad[2][1] = -0.03333333333333333*(2.0*(9.0*f[57]-1.0*(9.0*f[56]+6.708203932499369*(f[34]+f[33])))+3.0*(9.0*f[21]-1.0*(6.708203932499369*(f[8]-1.0*f[7])+5.0*f[3]))); 
  fReflXYQuad[2][2] = -0.03333333333333333*(2.0*(9.0*f[60]-1.0*(9.0*f[59]+6.708203932499369*(f[38]+f[37])))+3.0*(9.0*f[22]-1.0*(6.708203932499369*(f[10]-1.0*f[9])+5.0*f[4]))); 
  fReflXYQuad[2][3] = -0.03333333333333333*(2.0*(9.0*f[69]-1.0*(9.0*f[68]+6.708203932499369*(f[44]+f[43])))+3.0*(9.0*f[25]-1.0*(6.708203932499369*(f[13]-1.0*f[12])+5.0*f[5]))); 
  fReflXYQuad[2][4] = -0.02*(2.23606797749979*(13.41640786499874*(f[88]-1.0*f[87])-5.0*(2.0*(f[62]+f[61])+3.0*(f[24]-1.0*f[23])))+5.0*(9.0*f[51]-5.0*f[11])); 
  fReflXYQuad[2][5] = -0.02*(2.23606797749979*(13.41640786499874*(f[92]-1.0*f[91])-5.0*(2.0*(f[71]+f[70])+3.0*(f[27]-1.0*f[26])))+5.0*(9.0*f[52]-5.0*f[14])); 
  fReflXYQuad[2][6] = -0.02*(2.23606797749979*(13.41640786499874*(f[95]-1.0*f[94])-5.0*(2.0*(f[75]+f[74])+3.0*(f[29]-1.0*f[28])))+5.0*(9.0*f[53]-5.0*f[15])); 
  fReflXYQuad[2][7] = -0.1*(9.0*f[58]-1.0*(6.708203932499369*(f[36]-1.0*f[35])+5.0*f[18])); 
  fReflXYQuad[2][8] = -0.1*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])); 
  fReflXYQuad[2][9] = -0.1*(9.0*f[80]-1.0*(6.708203932499369*(f[48]-1.0*f[47])+5.0*f[20])); 
  fReflXYQuad[2][10] = -0.03333333333333333*(2.0*(9.0*f[108]-1.0*(9.0*f[107]+6.708203932499369*(f[97]+f[96])))+3.0*(9.0*f[86]-1.0*(6.708203932499369*(f[55]-1.0*f[54])+5.0*f[30]))); 
  fReflXYQuad[2][11] = -0.1*(9.0*f[89]-1.0*(6.708203932499369*(f[64]-1.0*f[63])+5.0*f[39])); 
  fReflXYQuad[2][12] = -0.1*(9.0*f[90]-1.0*(6.708203932499369*(f[67]-1.0*f[66])+5.0*f[42])); 
  fReflXYQuad[2][13] = -0.1*(9.0*f[93]-1.0*(6.708203932499369*(f[73]-1.0*f[72])+5.0*f[45])); 
  fReflXYQuad[2][14] = -0.1*(9.0*f[100]-1.0*(6.708203932499369*(f[78]-1.0*f[77])+5.0*f[46])); 
  fReflXYQuad[2][15] = -0.1*(9.0*f[103]-1.0*(6.708203932499369*(f[82]-1.0*f[81])+5.0*f[49])); 
  fReflXYQuad[2][16] = -0.1*(9.0*f[104]-1.0*(6.708203932499369*(f[84]-1.0*f[83])+5.0*f[50])); 
  fReflXYQuad[2][17] = -0.1*(9.0*f[109]-1.0*(6.708203932499369*(f[99]-1.0*f[98])+5.0*f[76])); 
  fReflXYQuad[2][18] = -0.1*(9.0*f[110]-1.0*(6.708203932499369*(f[102]-1.0*f[101])+5.0*f[79])); 
  fReflXYQuad[2][19] = -0.1*(9.0*f[111]-1.0*(6.708203932499369*(f[106]-1.0*f[105])+5.0*f[85])); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]-1.0*(f[25]+f[21])))+5.0*(3.0*(f[5]+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])-6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*((-1.0*f[58])+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]-1.0*f[1]))-5.0*f[20])))+10.0*(9.0*(f[8]-1.0*f[7]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[24]-1.0*f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*(f[27]-1.0*f[26]+f[25]-1.0*f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[2]+f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))+9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-1.0*(5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]-1.0*f[2]+f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*(f[13]-1.0*f[12]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*(f[29]-1.0*f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[2]-1.0*f[1])))-6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*((-1.0*f[84])+f[83]-1.0*f[64]+f[63])+2.0*(3.0*(2.0*(f[59]-1.0*f[60])-3.0*f[22]+2.23606797749979*(f[10]-1.0*f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]-1.0*(6.708203932499369*(f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))-1.0*(9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]-6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[1]-1.0*(f[5]+f[2]))-2.0*(f[17]+f[16])))))+10.0*(9.0*(f[13]-1.0*(f[12]+f[6]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*(f[29]-1.0*(f[28]+f[22]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]-1.0*f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))-5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))-6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57]-1.0*f[56])+3.0*f[21]))+5.0*(5.0*f[20]+2.0*(3.0*(f[1]-1.0*(f[3]+f[2]))-2.0*(f[18]+f[17]+f[16])))))+10.0*(9.0*(f[8]-1.0*(f[7]+f[6]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]-6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*(f[24]-1.0*(f[23]+f[22]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))-5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]+f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108])-1.0*(81.0*f[107]+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37]))))+3.0*(3.0*(27.0*f[86]+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]-6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_4 
  vcutSq_i = -0.025*(1.732050807568877*((-21.21320343559643*phiWall[17])+21.21320343559643*phi[17]-14.14213562373095*phiWall[3]+14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[16]-2.828427124746191*phi[16]-1.414213562373095*phiWall[11]+1.414213562373095*phi[11])-5.656854249492382*phiWall[14]+5.656854249492382*phi[14]+7.071067811865476*phiWall[13]-7.071067811865476*phi[13])+2.23606797749979*((-14.14213562373095*phiWall[9])+14.14213562373095*phi[9]-5.656854249492382*phiWall[8]+5.656854249492382*phi[8]+7.071067811865476*phiWall[7]-7.071067811865476*phi[7]+1.732050807568877*(8.485281374238571*phiWall[6]-8.485281374238571*phi[6])+8.485281374238571*phiWall[2]-8.485281374238571*phi[2])-14.14213562373095*phiWall[0]+14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[3][0] = 0.0; 
  fReflXYQuad[3][1] = 0.0; 
  fReflXYQuad[3][2] = 0.0; 
  fReflXYQuad[3][3] = 0.0; 
  fReflXYQuad[3][4] = 0.0; 
  fReflXYQuad[3][5] = 0.0; 
  fReflXYQuad[3][6] = 0.0; 
  fReflXYQuad[3][7] = 0.0; 
  fReflXYQuad[3][8] = 0.0; 
  fReflXYQuad[3][9] = 0.0; 
  fReflXYQuad[3][10] = 0.0; 
  fReflXYQuad[3][11] = 0.0; 
  fReflXYQuad[3][12] = 0.0; 
  fReflXYQuad[3][13] = 0.0; 
  fReflXYQuad[3][14] = 0.0; 
  fReflXYQuad[3][15] = 0.0; 
  fReflXYQuad[3][16] = 0.0; 
  fReflXYQuad[3][17] = 0.0; 
  fReflXYQuad[3][18] = 0.0; 
  fReflXYQuad[3][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.05*(2.23606797749979*(6.708203932499369*f[31]+4.0*f[17]-1.0*(5.0*f[16]+6.0*f[2]))+10.0*f[0]); 
  fReflXYQuad[3][1] = 0.01666666666666667*(45.0*f[56]+6.708203932499369*(4.0*f[34]-5.0*f[33])+6.0*(5.0*f[3]-6.708203932499369*f[8])); 
  fReflXYQuad[3][2] = 0.01666666666666667*(45.0*f[59]+6.708203932499369*(4.0*f[38]-5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])); 
  fReflXYQuad[3][3] = 0.01666666666666667*(45.0*f[68]+6.708203932499369*(4.0*f[44]-5.0*f[43])+6.0*(5.0*f[5]-6.708203932499369*f[13])); 
  fReflXYQuad[3][4] = 0.05*(2.23606797749979*(6.708203932499369*f[87]+4.0*f[62]-1.0*(5.0*f[61]+6.0*f[24]))+10.0*f[11]); 
  fReflXYQuad[3][5] = 0.05*(2.23606797749979*(6.708203932499369*f[91]+4.0*f[71]-1.0*(5.0*f[70]+6.0*f[27]))+10.0*f[14]); 
  fReflXYQuad[3][6] = 0.05*(2.23606797749979*(6.708203932499369*f[94]+4.0*f[75]-1.0*(5.0*f[74]+6.0*f[29]))+10.0*f[15]); 
  fReflXYQuad[3][7] = -0.1*(6.708203932499369*f[36]-5.0*f[18]); 
  fReflXYQuad[3][8] = -0.1*(6.708203932499369*f[41]-5.0*f[19]); 
  fReflXYQuad[3][9] = -0.1*(6.708203932499369*f[48]-5.0*f[20]); 
  fReflXYQuad[3][10] = 0.01666666666666667*(45.0*f[107]+6.708203932499369*(4.0*f[97]-5.0*f[96])+6.0*(5.0*f[30]-6.708203932499369*f[55])); 
  fReflXYQuad[3][11] = -0.1*(6.708203932499369*f[64]-5.0*f[39]); 
  fReflXYQuad[3][12] = -0.1*(6.708203932499369*f[67]-5.0*f[42]); 
  fReflXYQuad[3][13] = -0.1*(6.708203932499369*f[73]-5.0*f[45]); 
  fReflXYQuad[3][14] = -0.1*(6.708203932499369*f[78]-5.0*f[46]); 
  fReflXYQuad[3][15] = -0.1*(6.708203932499369*f[82]-5.0*f[49]); 
  fReflXYQuad[3][16] = -0.1*(6.708203932499369*f[84]-5.0*f[50]); 
  fReflXYQuad[3][17] = -0.1*(6.708203932499369*f[99]-5.0*f[76]); 
  fReflXYQuad[3][18] = -0.1*(6.708203932499369*f[102]-5.0*f[79]); 
  fReflXYQuad[3][19] = -0.1*(6.708203932499369*f[106]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))-5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]+f[73]+f[71])-1.0*(5.0*(f[70]+f[68]+f[56])+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-1.0*(5.0*f[16]+6.0*(f[5]+f[3]+f[2]))))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))-5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))/(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(11.18033988749895*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]+6.0*(f[3]+f[2]))-1.0*(15.0*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))+10.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]-1.0*f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]-1.0*(f[3]+f[2])))))+10.0*(9.0*(f[14]+f[13])-1.0*(9.0*f[8]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[30]+f[29])-1.0*(9.0*f[24]+5.0*f[4]))-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]+6.0*(f[5]+f[2]))+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31])))-10.0*(9.0*f[13]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))/(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[41]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])-1.0*(5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]+6.0*f[2]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(45.0*(f[84]+f[64]+f[59])-6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(6.708203932499369*f[41]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[78]-1.0*(6.708203932499369*(f[46]-1.0*f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[5]-1.0*f[2]))-6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.05*(9.0*f[78]-1.0*(6.708203932499369*(f[46]-1.0*f[41])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))+9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*(f[68]-1.0*f[70])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]-1.0*f[3]+f[2])))+10.0*(9.0*(f[14]-1.0*f[13]+f[8])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*(f[30]-1.0*f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[67]-1.0*(6.708203932499369*(f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[3]-1.0*f[2]))-6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.05*(9.0*f[67]-1.0*(6.708203932499369*(f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))-1.0*(9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56]))+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[2]-1.0*(f[5]+f[3])))))+10.0*(9.0*(f[14]-1.0*(f[13]+f[8]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*(f[30]-1.0*(f[29]+f[24]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_5 
  vcutSq_i = -0.125*(3.872983346207417*(1.414213562373095*phiWall[14]-1.414213562373095*phi[14]+1.414213562373095*phiWall[13]-1.414213562373095*phi[13])+2.23606797749979*((-2.828427124746191*phiWall[9])+2.828427124746191*phi[9]+1.414213562373095*phiWall[8]-1.414213562373095*phi[8]+1.414213562373095*phiWall[7]-1.414213562373095*phi[7])+1.732050807568877*(2.828427124746191*phi[3]-2.828427124746191*phiWall[3])-2.828427124746191*phiWall[0]+2.828427124746191*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[4][0] = 0.0; 
  fReflXYQuad[4][1] = 0.0; 
  fReflXYQuad[4][2] = 0.0; 
  fReflXYQuad[4][3] = 0.0; 
  fReflXYQuad[4][4] = 0.0; 
  fReflXYQuad[4][5] = 0.0; 
  fReflXYQuad[4][6] = 0.0; 
  fReflXYQuad[4][7] = 0.0; 
  fReflXYQuad[4][8] = 0.0; 
  fReflXYQuad[4][9] = 0.0; 
  fReflXYQuad[4][10] = 0.0; 
  fReflXYQuad[4][11] = 0.0; 
  fReflXYQuad[4][12] = 0.0; 
  fReflXYQuad[4][13] = 0.0; 
  fReflXYQuad[4][14] = 0.0; 
  fReflXYQuad[4][15] = 0.0; 
  fReflXYQuad[4][16] = 0.0; 
  fReflXYQuad[4][17] = 0.0; 
  fReflXYQuad[4][18] = 0.0; 
  fReflXYQuad[4][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[4][0] = -0.25*(2.23606797749979*(f[17]+f[16])-2.0*f[0]); 
  fReflXYQuad[4][1] = -0.08333333333333333*(6.708203932499369*(f[34]+f[33])-6.0*f[3]); 
  fReflXYQuad[4][2] = -0.08333333333333333*(6.708203932499369*(f[38]+f[37])-6.0*f[4]); 
  fReflXYQuad[4][3] = -0.08333333333333333*(6.708203932499369*(f[44]+f[43])-6.0*f[5]); 
  fReflXYQuad[4][4] = -0.25*(2.23606797749979*(f[62]+f[61])-2.0*f[11]); 
  fReflXYQuad[4][5] = -0.25*(2.23606797749979*(f[71]+f[70])-2.0*f[14]); 
  fReflXYQuad[4][6] = -0.25*(2.23606797749979*(f[75]+f[74])-2.0*f[15]); 
  fReflXYQuad[4][7] = 0.5*f[18]; 
  fReflXYQuad[4][8] = 0.5*f[19]; 
  fReflXYQuad[4][9] = 0.5*f[20]; 
  fReflXYQuad[4][10] = -0.08333333333333333*(6.708203932499369*(f[97]+f[96])-6.0*f[30]); 
  fReflXYQuad[4][11] = 0.5*f[39]; 
  fReflXYQuad[4][12] = 0.5*f[42]; 
  fReflXYQuad[4][13] = 0.5*f[45]; 
  fReflXYQuad[4][14] = 0.5*f[46]; 
  fReflXYQuad[4][15] = 0.5*f[49]; 
  fReflXYQuad[4][16] = 0.5*f[50]; 
  fReflXYQuad[4][17] = 0.5*f[76]; 
  fReflXYQuad[4][18] = 0.5*f[79]; 
  fReflXYQuad[4][19] = 0.5*f[85]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.05*(9.0*f[79]-6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])+6.0*(f[5]+f[3]))+6.708203932499369*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33])))-10.0*(9.0*f[14]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]+f[76])-5.0*(f[75]+f[74]+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]+f[11])-3.0*f[30])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.05*(9.0*f[79]-6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))/(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(6.708203932499369*f[42]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])-1.0*(5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])+6.0*f[3]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(45.0*(f[85]+f[62]+f[61])-6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[11])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(6.708203932499369*f[42]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.05*(9.0*f[79]-1.0*(6.708203932499369*(f[46]-1.0*f[42])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[5]-1.0*f[3]))-6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[76]-1.0*f[85])+5.0*((-1.0*(f[75]+f[74]))+f[62]+f[61]))+2.0*(3.0*(2.23606797749979*(f[15]-1.0*f[11])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.05*(9.0*f[79]-1.0*(6.708203932499369*(f[46]-1.0*f[42])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))/(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(6.708203932499369*f[46]-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.025*(2.23606797749979*(6.708203932499369*(f[45]+f[44]+f[43])+4.0*f[20]-1.0*(5.0*(f[18]+f[17]+f[16])+6.0*f[5]))+10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(45.0*(f[76]+f[75]+f[74])+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(5.0*f[4]-6.708203932499369*f[15])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(6.708203932499369*f[46]-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))/(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.04166666666666666*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.25*f[19])*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.125*(2.23606797749979*(f[20]+f[18]+f[17]+f[16])-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.04166666666666666*(6.708203932499369*(f[50]+f[39]+f[38]+f[37])-6.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.25*f[19])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))/(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(6.708203932499369*f[46]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.025*(15.0*(f[45]+f[44]+f[43])-1.0*(2.23606797749979*(4.0*f[20]-5.0*(f[18]+f[17]+f[16])+6.0*f[5])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(45.0*(f[76]+f[75]+f[74])-1.0*(6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]+f[37]))+6.0*(6.708203932499369*f[15]+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(6.708203932499369*f[46]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[79]+6.708203932499369*(f[46]-1.0*f[42])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[71]+f[70])+4.0*(f[20]+f[18])-5.0*(f[17]+f[16])+6.0*(f[3]-1.0*f[5]))+6.708203932499369*(4.0*(f[49]-1.0*f[45])+5.0*(f[44]+f[43]-1.0*(f[34]+f[33]))))+10.0*(5.0*f[0]-9.0*f[14])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])+4.0*(f[50]+f[39])-5.0*(f[38]+f[37]))+3.0*(3.0*(4.0*(f[85]-1.0*f[76])+5.0*(f[75]+f[74]-1.0*(f[62]+f[61])))+2.0*(3.0*(2.23606797749979*(f[11]-1.0*f[15])-3.0*f[30])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[79]+6.708203932499369*(f[46]-1.0*f[42])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(6.708203932499369*f[42]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[49]+f[34]+f[33])+5.0*f[20]-4.0*f[18]+5.0*(f[17]+f[16])-6.0*f[3])-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(45.0*(f[85]+f[62]+f[61])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*(f[38]+f[37]))-6.0*(6.708203932499369*f[11]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(6.708203932499369*f[42]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.05*(9.0*f[79]+6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(11.18033988749895*(9.0*(f[71]+f[70])-4.0*(f[20]+f[18])+5.0*(f[17]+f[16])-6.0*(f[5]+f[3]))-1.0*(15.0*(4.0*(f[49]+f[45])-5.0*(f[44]+f[43]+f[34]+f[33]))+10.0*(9.0*f[14]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[97]+f[96])-4.0*(f[50]+f[39])+5.0*(f[38]+f[37]))+3.0*(3.0*(5.0*(f[75]+f[74]+f[62]+f[61])-4.0*(f[85]+f[76]))-2.0*(3.0*(3.0*f[30]+2.23606797749979*(f[15]+f[11]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.05*(9.0*f[79]+6.708203932499369*(f[46]+f[42])+5.0*f[19]))*fac; 
   } 
  } 
  fReflXYQuad[4][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[4][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[4][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_6 
  vcutSq_i = 0.025*(1.732050807568877*((-21.21320343559643*phiWall[17])+21.21320343559643*phi[17]+14.14213562373095*phiWall[3]-14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[16]-2.828427124746191*phi[16]-1.414213562373095*phiWall[11]+1.414213562373095*phi[11])+5.656854249492382*phiWall[14]-5.656854249492382*phi[14]-7.071067811865476*phiWall[13]+7.071067811865476*phi[13])+2.23606797749979*(14.14213562373095*phiWall[9]-14.14213562373095*phi[9]+5.656854249492382*phiWall[8]-5.656854249492382*phi[8]-7.071067811865476*phiWall[7]+7.071067811865476*phi[7]+1.732050807568877*(8.485281374238571*phiWall[6]-8.485281374238571*phi[6])+8.485281374238571*phiWall[2]-8.485281374238571*phi[2])+14.14213562373095*phiWall[0]-14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[5][0] = 0.0; 
  fReflXYQuad[5][1] = 0.0; 
  fReflXYQuad[5][2] = 0.0; 
  fReflXYQuad[5][3] = 0.0; 
  fReflXYQuad[5][4] = 0.0; 
  fReflXYQuad[5][5] = 0.0; 
  fReflXYQuad[5][6] = 0.0; 
  fReflXYQuad[5][7] = 0.0; 
  fReflXYQuad[5][8] = 0.0; 
  fReflXYQuad[5][9] = 0.0; 
  fReflXYQuad[5][10] = 0.0; 
  fReflXYQuad[5][11] = 0.0; 
  fReflXYQuad[5][12] = 0.0; 
  fReflXYQuad[5][13] = 0.0; 
  fReflXYQuad[5][14] = 0.0; 
  fReflXYQuad[5][15] = 0.0; 
  fReflXYQuad[5][16] = 0.0; 
  fReflXYQuad[5][17] = 0.0; 
  fReflXYQuad[5][18] = 0.0; 
  fReflXYQuad[5][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[5][0] = -0.05*(15.0*f[31]-1.0*(2.23606797749979*(4.0*f[17]-5.0*f[16]+6.0*f[2])+10.0*f[0])); 
  fReflXYQuad[5][1] = -0.01666666666666667*(45.0*f[56]-1.0*(6.708203932499369*(4.0*f[34]-5.0*f[33])+6.0*(6.708203932499369*f[8]+5.0*f[3]))); 
  fReflXYQuad[5][2] = -0.01666666666666667*(45.0*f[59]-1.0*(6.708203932499369*(4.0*f[38]-5.0*f[37])+6.0*(6.708203932499369*f[10]+5.0*f[4]))); 
  fReflXYQuad[5][3] = -0.01666666666666667*(45.0*f[68]-1.0*(6.708203932499369*(4.0*f[44]-5.0*f[43])+6.0*(6.708203932499369*f[13]+5.0*f[5]))); 
  fReflXYQuad[5][4] = -0.05*(15.0*f[87]-1.0*(2.23606797749979*(4.0*f[62]-5.0*f[61]+6.0*f[24])+10.0*f[11])); 
  fReflXYQuad[5][5] = -0.05*(15.0*f[91]-1.0*(2.23606797749979*(4.0*f[71]-5.0*f[70]+6.0*f[27])+10.0*f[14])); 
  fReflXYQuad[5][6] = -0.05*(15.0*f[94]-1.0*(2.23606797749979*(4.0*f[75]-5.0*f[74]+6.0*f[29])+10.0*f[15])); 
  fReflXYQuad[5][7] = 0.1*(6.708203932499369*f[36]+5.0*f[18]); 
  fReflXYQuad[5][8] = 0.1*(6.708203932499369*f[41]+5.0*f[19]); 
  fReflXYQuad[5][9] = 0.1*(6.708203932499369*f[48]+5.0*f[20]); 
  fReflXYQuad[5][10] = -0.01666666666666667*(45.0*f[107]-1.0*(6.708203932499369*(4.0*f[97]-5.0*f[96])+6.0*(6.708203932499369*f[55]+5.0*f[30]))); 
  fReflXYQuad[5][11] = 0.1*(6.708203932499369*f[64]+5.0*f[39]); 
  fReflXYQuad[5][12] = 0.1*(6.708203932499369*f[67]+5.0*f[42]); 
  fReflXYQuad[5][13] = 0.1*(6.708203932499369*f[73]+5.0*f[45]); 
  fReflXYQuad[5][14] = 0.1*(6.708203932499369*f[78]+5.0*f[46]); 
  fReflXYQuad[5][15] = 0.1*(6.708203932499369*f[82]+5.0*f[49]); 
  fReflXYQuad[5][16] = 0.1*(6.708203932499369*f[84]+5.0*f[50]); 
  fReflXYQuad[5][17] = 0.1*(6.708203932499369*f[99]+5.0*f[76]); 
  fReflXYQuad[5][18] = 0.1*(6.708203932499369*f[102]+5.0*f[79]); 
  fReflXYQuad[5][19] = 0.1*(6.708203932499369*f[106]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])+5.0*(f[31]-1.0*f[33]))+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*f[70]-1.0*(5.0*(f[68]+f[56])+6.0*f[27]))+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*(f[5]+f[3]-1.0*f[2])))+10.0*(9.0*((-1.0*f[14])+f[13]+f[8])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]+f[99]-1.0*f[97])+5.0*(9.0*f[96]-1.0*(9.0*(f[94]+f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]+f[24])-5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[10]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[41]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[42]-1.0*f[41])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]+6.0*(f[2]-1.0*f[3]))+6.708203932499369*(5.0*(f[49]-1.0*f[48])+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[8])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])-5.0*f[50]+4.0*(f[39]+f[38])-5.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[42]-1.0*f[41])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*(f[79]-1.0*f[78]+f[67])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*(f[33]-1.0*f[31]))-1.0*(9.0*(4.0*(f[82]-1.0*f[73]+f[71])+5.0*((-1.0*f[70])+f[68]-1.0*f[56])+6.0*f[27])+5.0*((-4.0*(f[20]+f[18]+f[17]))+5.0*f[16]+6.0*((-1.0*f[5])+f[3]-1.0*f[2]))))+10.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[8])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]-1.0*f[99]+f[97])+5.0*(9.0*(f[94]-1.0*f[96])-1.0*(9.0*f[87]+4.0*(f[50]+f[39]+f[38])-5.0*f[37])))+3.0*(15.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[24])+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*((-1.0*f[15])+f[11]-1.0*f[10]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*((-1.0*f[46])+f[42]-1.0*f[41]))+5.0*(9.0*(f[79]-1.0*f[78]+f[67])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[46]-1.0*f[41])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[73]+f[68])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]+6.0*(f[2]-1.0*f[5]))+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44]+5.0*(f[43]-1.0*(f[36]+f[31]))))+10.0*(5.0*f[0]-9.0*f[13])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])+4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+3.0*(3.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+2.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[46]-1.0*f[41])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[41]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[48]+f[36]+f[31])+5.0*(f[20]+f[18])-4.0*f[17]+5.0*f[16]-6.0*f[2])-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[84]+f[64]+f[59])+6.708203932499369*(5.0*(f[50]+f[39])-4.0*f[38]+5.0*f[37])-6.0*(6.708203932499369*f[10]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[41]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[73]+f[68])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]-6.0*(f[5]+f[2]))-1.0*(15.0*(4.0*f[48]-5.0*f[45]+4.0*f[44]-5.0*(f[43]+f[36]+f[31]))+10.0*(9.0*f[13]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[99]+f[94])-4.0*f[50]+5.0*f[39]-4.0*f[38]+5.0*f[37])+3.0*(3.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))-2.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[79]+f[78])-1.0*(9.0*f[67]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[91]+4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34])-5.0*(f[33]+f[31]))+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68])-1.0*(5.0*f[56]+6.0*f[27]))+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*((-1.0*f[5])+f[3]+f[2])))+10.0*(9.0*(f[8]-1.0*(f[14]+f[13]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[107]+6.708203932499369*(36.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(9.0*(f[96]+f[94]-1.0*f[87])+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+2.0*(5.0*(9.0*(f[24]-1.0*(f[30]+f[29]))+5.0*f[4])-6.708203932499369*(9.0*f[55]+5.0*(f[15]-1.0*(f[11]+f[10])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[79]+f[78])-1.0*(9.0*f[67]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[82]+f[56])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]-6.0*(f[3]+f[2]))+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34])+5.0*(f[33]+f[31])))-10.0*(9.0*f[8]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[106]+f[87])+5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))/(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))+5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[91]-4.0*(f[49]+f[48]+f[45]+f[44])+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56])+6.0*f[27])+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]+6.0*(f[5]+f[3]+f[2])))+10.0*(9.0*(f[14]+f[13]+f[8])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[107]-6.708203932499369*(36.0*(f[106]+f[99]+f[97])+5.0*((-9.0*(f[96]+f[94]+f[87]))+4.0*(f[50]+f[39]+f[38])-5.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-2.0*(6.708203932499369*(9.0*f[55]+5.0*(f[15]+f[11]+f[10]))+5.0*(9.0*(f[30]+f[29]+f[24])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[102]+5.0*(f[46]+f[42]+f[41]))+5.0*(9.0*(f[79]+f[78]+f[67])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[5][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[5][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[5][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_7 
  vcutSq_i = -0.01*(2.23606797749979*(63.63961030678928*phiWall[19]-63.63961030678928*phi[19]-35.35533905932738*phiWall[9]+35.35533905932738*phi[9]-14.14213562373095*phiWall[8]+14.14213562373095*phi[8]-14.14213562373095*phiWall[7]+14.14213562373095*phi[7]+1.732050807568877*(21.21320343559643*phiWall[6]-21.21320343559643*(phi[6]+phiWall[5])+21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*(phi[2]+phiWall[1])+21.21320343559643*phi[1])+1.732050807568877*((-42.42640687119286*phiWall[18])+42.42640687119286*(phi[18]+phiWall[17])-42.42640687119286*phi[17]+63.63961030678928*phiWall[10]-63.63961030678928*phi[10]-35.35533905932738*phiWall[3]+35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*(phi[16]+phiWall[15])+7.071067811865476*phi[15]-2.828427124746191*phiWall[12]+2.828427124746191*(phi[12]+phiWall[11])-2.828427124746191*phi[11])-14.14213562373095*phiWall[14]+14.14213562373095*phi[14]-14.14213562373095*phiWall[13]+14.14213562373095*phi[13])+63.63961030678928*phiWall[4]-63.63961030678928*phi[4]-35.35533905932738*phiWall[0]+35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[6][0] = 0.0; 
  fReflXYQuad[6][1] = 0.0; 
  fReflXYQuad[6][2] = 0.0; 
  fReflXYQuad[6][3] = 0.0; 
  fReflXYQuad[6][4] = 0.0; 
  fReflXYQuad[6][5] = 0.0; 
  fReflXYQuad[6][6] = 0.0; 
  fReflXYQuad[6][7] = 0.0; 
  fReflXYQuad[6][8] = 0.0; 
  fReflXYQuad[6][9] = 0.0; 
  fReflXYQuad[6][10] = 0.0; 
  fReflXYQuad[6][11] = 0.0; 
  fReflXYQuad[6][12] = 0.0; 
  fReflXYQuad[6][13] = 0.0; 
  fReflXYQuad[6][14] = 0.0; 
  fReflXYQuad[6][15] = 0.0; 
  fReflXYQuad[6][16] = 0.0; 
  fReflXYQuad[6][17] = 0.0; 
  fReflXYQuad[6][18] = 0.0; 
  fReflXYQuad[6][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[6][0] = 0.02*(2.23606797749979*(13.41640786499874*(f[32]-1.0*f[31])+5.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+5.0*(5.0*f[0]-9.0*f[6])); 
  fReflXYQuad[6][1] = 0.03333333333333333*(2.0*(9.0*(f[57]-1.0*f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[7]-1.0*f[8])-3.0*f[21])+5.0*f[3])); 
  fReflXYQuad[6][2] = 0.03333333333333333*(2.0*(9.0*(f[60]-1.0*f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[10])-3.0*f[22])+5.0*f[4])); 
  fReflXYQuad[6][3] = 0.03333333333333333*(2.0*(9.0*(f[69]-1.0*f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[12]-1.0*f[13])-3.0*f[25])+5.0*f[5])); 
  fReflXYQuad[6][4] = 0.02*(2.23606797749979*(13.41640786499874*(f[88]-1.0*f[87])+5.0*(2.0*(f[62]+f[61])+3.0*(f[23]-1.0*f[24])))+5.0*(5.0*f[11]-9.0*f[51])); 
  fReflXYQuad[6][5] = 0.02*(2.23606797749979*(13.41640786499874*(f[92]-1.0*f[91])+5.0*(2.0*(f[71]+f[70])+3.0*(f[26]-1.0*f[27])))+5.0*(5.0*f[14]-9.0*f[52])); 
  fReflXYQuad[6][6] = 0.02*(2.23606797749979*(13.41640786499874*(f[95]-1.0*f[94])+5.0*(2.0*(f[75]+f[74])+3.0*(f[28]-1.0*f[29])))+5.0*(5.0*f[15]-9.0*f[53])); 
  fReflXYQuad[6][7] = -0.1*(9.0*f[58]+6.708203932499369*(f[36]-1.0*f[35])-5.0*f[18]); 
  fReflXYQuad[6][8] = -0.1*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]); 
  fReflXYQuad[6][9] = -0.1*(9.0*f[80]+6.708203932499369*(f[48]-1.0*f[47])-5.0*f[20]); 
  fReflXYQuad[6][10] = 0.03333333333333333*(2.0*(9.0*(f[108]-1.0*f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[54]-1.0*f[55])-3.0*f[86])+5.0*f[30])); 
  fReflXYQuad[6][11] = -0.1*(9.0*f[89]+6.708203932499369*(f[64]-1.0*f[63])-5.0*f[39]); 
  fReflXYQuad[6][12] = -0.1*(9.0*f[90]+6.708203932499369*(f[67]-1.0*f[66])-5.0*f[42]); 
  fReflXYQuad[6][13] = -0.1*(9.0*f[93]+6.708203932499369*(f[73]-1.0*f[72])-5.0*f[45]); 
  fReflXYQuad[6][14] = -0.1*(9.0*f[100]+6.708203932499369*(f[78]-1.0*f[77])-5.0*f[46]); 
  fReflXYQuad[6][15] = -0.1*(9.0*f[103]+6.708203932499369*(f[82]-1.0*f[81])-5.0*f[49]); 
  fReflXYQuad[6][16] = -0.1*(9.0*f[104]+6.708203932499369*(f[84]-1.0*f[83])-5.0*f[50]); 
  fReflXYQuad[6][17] = -0.1*(9.0*f[109]+6.708203932499369*(f[99]-1.0*f[98])-5.0*f[76]); 
  fReflXYQuad[6][18] = -0.1*(9.0*f[110]+6.708203932499369*(f[102]-1.0*f[101])-5.0*f[79]); 
  fReflXYQuad[6][19] = -0.1*(9.0*f[111]+6.708203932499369*(f[106]-1.0*f[105])-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]-1.0*f[91])+5.0*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))+3.0*((-1.0*f[27])+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])-3.0*(f[5]+f[3]+f[2]-1.0*f[1]))))+5.0*(5.0*(9.0*(f[14]+f[13]-1.0*f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[83]-1.0*(f[85]+f[84]))-1.0*(10.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]+f[51]))+5.0*(f[15]+f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*(f[30]+f[29]-1.0*f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(f[46]+f[42]+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[47]-1.0*(f[49]+f[48]))+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+2.0*(2.0*(f[58]+f[57])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(5.0*f[20]+2.0*(3.0*(f[3]+f[2]-1.0*f[1])-2.0*(f[18]+f[17]+f[16]))))+10.0*(9.0*((-1.0*f[8])+f[7]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(3.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-2.23606797749979*(9.0*f[51]+5.0*(f[9]-1.0*(f[11]+f[10]))))+5.0*(9.0*((-1.0*f[24])+f[23]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[40]-1.0*(f[42]+f[41])))+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+5.0*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+6.708203932499369*(9.0*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+2.23606797749979*(9.0*(f[55]-1.0*(f[54]+f[53]-1.0*f[51]))+5.0*(f[15]-1.0*(f[11]+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(f[46]-1.0*(f[42]+f[41]-1.0*f[40])))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45])+4.0*(f[44]+f[43])+5.0*(f[35]-1.0*f[36])+4.0*(f[31]-1.0*f[32]))+9.0*(4.0*f[80]+5.0*(f[73]-1.0*f[72])+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]+6.0*f[25]))+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*(f[5]+f[2]-1.0*f[1])-2.0*(f[17]+f[16]))))+10.0*(9.0*((-1.0*f[13])+f[12]+f[6])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(3.0*(10.0*(f[59]-1.0*f[60])-2.23606797749979*(9.0*f[53]+5.0*(f[9]-1.0*(f[15]+f[10]))))+5.0*(9.0*((-1.0*f[29])+f[28]+f[22])-5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[40]-1.0*(f[46]+f[41])))+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])-5.0*(f[20]+f[18])+2.0*(2.0*(f[17]+f[16])+3.0*(f[1]-1.0*f[2])))+6.708203932499369*(5.0*(f[48]-1.0*f[47]+f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31])))+10.0*(5.0*f[0]-9.0*f[6])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])-5.0*(f[50]+f[39])+4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]-1.0*f[83]+f[64]-1.0*f[63])+2.0*(3.0*(2.0*f[60]-1.0*(2.0*f[59]+3.0*f[22]-2.23606797749979*(f[9]-1.0*f[10])))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[65]+6.708203932499369*(f[41]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*(f[78]-1.0*f[77]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*(f[36]-1.0*f[35])+4.0*(f[32]-1.0*f[31]))-1.0*(9.0*(4.0*f[80]+5.0*(f[72]-1.0*f[73])+4.0*(f[68]-1.0*f[69])-5.0*f[58]+6.0*f[25])+5.0*((-4.0*f[20])+5.0*f[18]+2.0*(3.0*((-1.0*f[5])+f[2]-1.0*f[1])-2.0*(f[17]+f[16])))))+10.0*(9.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95]))+5.0*((-1.0*(9.0*f[89]+4.0*f[50]))+5.0*f[39]-4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(3.0*(10.0*f[60]-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[53]+5.0*((-1.0*f[15])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[29])+f[28]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*((-1.0*f[46])+f[41]-1.0*f[40]))+5.0*(9.0*(f[78]-1.0*f[77]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+5.0*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])+3.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*(f[5]-1.0*f[3]+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64])-1.0*(10.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59])+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]-1.0*f[51])+5.0*((-1.0*f[15])+f[11]-1.0*f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*((-1.0*f[46])+f[42]-1.0*f[41]+f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*(f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*((-1.0*f[49])+f[48]-1.0*f[47])+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+9.0*(5.0*(f[82]-1.0*f[81]+f[80])+2.0*(2.0*(f[57]-1.0*f[58])-1.0*(2.0*f[56]+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]-1.0*f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*((-1.0*f[8])+f[7]-1.0*f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(3.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+2.23606797749979*(9.0*f[51]+5.0*((-1.0*f[11])+f[10]-1.0*f[9]))))+5.0*(9.0*((-1.0*f[24])+f[23]-1.0*f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*((-1.0*f[42])+f[41]-1.0*f[40]))+5.0*(9.0*(f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*f[92]+f[91])+5.0*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31])))+5.0*(9.0*(2.0*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+3.0*(f[27]-1.0*f[26]+f[25]+f[21]))+5.0*(3.0*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1])-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])-5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+6.708203932499369*(9.0*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87]))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+2.23606797749979*(9.0*(f[55]-1.0*f[54]+f[53]+f[51])+5.0*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9])))+5.0*(9.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])-5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*((-1.0*(f[46]+f[42]))+f[41]-1.0*f[40]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[6][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[6][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[6][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_8 
  vcutSq_i = 0.025*(1.732050807568877*((-21.21320343559643*phiWall[18])+21.21320343559643*phi[18]+14.14213562373095*phiWall[3]-14.14213562373095*phi[3])+3.872983346207417*(3.872983346207417*(2.828427124746191*phiWall[15]-2.828427124746191*phi[15]-1.414213562373095*phiWall[12]+1.414213562373095*phi[12])-7.071067811865476*phiWall[14]+7.071067811865476*phi[14]+5.656854249492382*phiWall[13]-5.656854249492382*phi[13])+2.23606797749979*(14.14213562373095*phiWall[9]-14.14213562373095*phi[9]-7.071067811865476*phiWall[8]+7.071067811865476*phi[8]+5.656854249492382*phiWall[7]-5.656854249492382*phi[7]+1.732050807568877*(8.485281374238571*phiWall[5]-8.485281374238571*phi[5])+8.485281374238571*phiWall[1]-8.485281374238571*phi[1])+14.14213562373095*phiWall[0]-14.14213562373095*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[7][0] = 0.0; 
  fReflXYQuad[7][1] = 0.0; 
  fReflXYQuad[7][2] = 0.0; 
  fReflXYQuad[7][3] = 0.0; 
  fReflXYQuad[7][4] = 0.0; 
  fReflXYQuad[7][5] = 0.0; 
  fReflXYQuad[7][6] = 0.0; 
  fReflXYQuad[7][7] = 0.0; 
  fReflXYQuad[7][8] = 0.0; 
  fReflXYQuad[7][9] = 0.0; 
  fReflXYQuad[7][10] = 0.0; 
  fReflXYQuad[7][11] = 0.0; 
  fReflXYQuad[7][12] = 0.0; 
  fReflXYQuad[7][13] = 0.0; 
  fReflXYQuad[7][14] = 0.0; 
  fReflXYQuad[7][15] = 0.0; 
  fReflXYQuad[7][16] = 0.0; 
  fReflXYQuad[7][17] = 0.0; 
  fReflXYQuad[7][18] = 0.0; 
  fReflXYQuad[7][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[7][0] = -0.05*(2.23606797749979*(6.708203932499369*f[32]+5.0*f[17]-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]); 
  fReflXYQuad[7][1] = -0.01666666666666667*(45.0*f[57]+6.708203932499369*(5.0*f[34]-4.0*f[33])-6.0*(6.708203932499369*f[7]+5.0*f[3])); 
  fReflXYQuad[7][2] = -0.01666666666666667*(45.0*f[60]+6.708203932499369*(5.0*f[38]-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])); 
  fReflXYQuad[7][3] = -0.01666666666666667*(45.0*f[69]+6.708203932499369*(5.0*f[44]-4.0*f[43])-6.0*(6.708203932499369*f[12]+5.0*f[5])); 
  fReflXYQuad[7][4] = -0.05*(2.23606797749979*(6.708203932499369*f[88]+5.0*f[62]-2.0*(2.0*f[61]+3.0*f[23]))-10.0*f[11]); 
  fReflXYQuad[7][5] = -0.05*(2.23606797749979*(6.708203932499369*f[92]+5.0*f[71]-2.0*(2.0*f[70]+3.0*f[26]))-10.0*f[14]); 
  fReflXYQuad[7][6] = -0.05*(2.23606797749979*(6.708203932499369*f[95]+5.0*f[75]-2.0*(2.0*f[74]+3.0*f[28]))-10.0*f[15]); 
  fReflXYQuad[7][7] = 0.1*(6.708203932499369*f[35]+5.0*f[18]); 
  fReflXYQuad[7][8] = 0.1*(6.708203932499369*f[40]+5.0*f[19]); 
  fReflXYQuad[7][9] = 0.1*(6.708203932499369*f[47]+5.0*f[20]); 
  fReflXYQuad[7][10] = -0.01666666666666667*(45.0*f[108]+6.708203932499369*(5.0*f[97]-4.0*f[96])-6.0*(6.708203932499369*f[54]+5.0*f[30])); 
  fReflXYQuad[7][11] = 0.1*(6.708203932499369*f[63]+5.0*f[39]); 
  fReflXYQuad[7][12] = 0.1*(6.708203932499369*f[66]+5.0*f[42]); 
  fReflXYQuad[7][13] = 0.1*(6.708203932499369*f[72]+5.0*f[45]); 
  fReflXYQuad[7][14] = 0.1*(6.708203932499369*f[77]+5.0*f[46]); 
  fReflXYQuad[7][15] = 0.1*(6.708203932499369*f[81]+5.0*f[49]); 
  fReflXYQuad[7][16] = 0.1*(6.708203932499369*f[83]+5.0*f[50]); 
  fReflXYQuad[7][17] = 0.1*(6.708203932499369*f[98]+5.0*f[76]); 
  fReflXYQuad[7][18] = 0.1*(6.708203932499369*f[101]+5.0*f[79]); 
  fReflXYQuad[7][19] = 0.1*(6.708203932499369*f[105]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34])+4.0*f[33]+5.0*f[32])+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])+6.0*f[26]))+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*(f[5]+f[3]-1.0*f[1])-2.0*f[16])))+10.0*(9.0*((-1.0*f[14])+f[12]+f[7])-5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96])+5.0*((-1.0*(9.0*(f[95]+f[88])+4.0*(f[50]+f[39])-5.0*f[38]))-4.0*f[37]))+3.0*(15.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]+f[23])-5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[9]-1.0*(f[15]+f[11])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[40]-1.0*(f[46]+f[42])))+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[42]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])-5.0*f[20]+4.0*f[18]-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[3])))+6.708203932499369*(5.0*(f[49]-1.0*f[47])+4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))+10.0*(5.0*f[0]-9.0*f[7])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])-5.0*f[50]+4.0*f[39]-5.0*f[38]+4.0*f[37])+3.0*(3.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[42]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*(f[79]-1.0*f[77]+f[66])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32]))-1.0*(9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*(f[69]-1.0*f[57])+6.0*f[26])+5.0*((-4.0*(f[20]+f[18]))+5.0*f[17]+2.0*(3.0*((-1.0*f[5])+f[3]-1.0*f[1])-2.0*f[16]))))+10.0*(9.0*((-1.0*f[14])+f[12]-1.0*f[7])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96]))+5.0*(9.0*f[95]-1.0*(9.0*f[88]+4.0*(f[50]+f[39])-5.0*f[38])-4.0*f[37]))+3.0*(15.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+2.0*(5.0*(9.0*((-1.0*f[30])+f[28]-1.0*f[23])+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*((-1.0*f[15])+f[11]-1.0*f[9]))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*((-1.0*f[46])+f[42]-1.0*f[40]))+5.0*(9.0*(f[79]-1.0*f[77]+f[66])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[46]-1.0*f[40])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(5.0*(9.0*(f[72]+f[69])+4.0*f[20]-5.0*(f[18]+f[17])+2.0*(2.0*f[16]+3.0*(f[1]-1.0*f[5])))+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-1.0*(4.0*f[43]+5.0*(f[35]+f[32]))))+10.0*(5.0*f[0]-9.0*f[12])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])+4.0*f[50]-5.0*(f[39]+f[38])+4.0*f[37])+3.0*(3.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+2.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[46]-1.0*f[40])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))/(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[40]+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.025*(2.23606797749979*(6.708203932499369*(f[47]+f[35]+f[32])+5.0*(f[20]+f[18]+f[17])-2.0*(2.0*f[16]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(45.0*(f[83]+f[63]+f[60])+6.708203932499369*(5.0*(f[50]+f[39]+f[38])-4.0*f[37])-6.0*(6.708203932499369*f[9]+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(6.708203932499369*f[40]+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))/(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(11.18033988749895*(9.0*(f[72]+f[69])-4.0*f[20]+5.0*(f[18]+f[17])-2.0*(2.0*f[16]+3.0*(f[5]+f[1])))-1.0*(15.0*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43]-5.0*(f[35]+f[32]))+10.0*(9.0*f[12]+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[98]+f[95])-4.0*f[50]+5.0*(f[39]+f[38])-4.0*f[37])+3.0*(3.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))-2.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))/(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[79]+f[77])-1.0*(9.0*f[66]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[92]+4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33]-5.0*f[32])+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69]-1.0*(5.0*f[57]+6.0*f[26]))+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*((-1.0*f[5])+f[3]+f[1]))))+10.0*(9.0*(f[7]-1.0*(f[14]+f[12]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[108]+6.708203932499369*(9.0*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96])+5.0*(9.0*(f[95]-1.0*f[88])+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61]-5.0*f[60])+2.0*(5.0*(9.0*(f[23]-1.0*(f[30]+f[28]))+5.0*f[4])-6.708203932499369*(9.0*f[54]+5.0*(f[15]-1.0*(f[11]+f[9])))))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]-1.0*(f[42]+f[40])))+5.0*(9.0*(f[79]+f[77])-1.0*(9.0*f[66]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[81]+f[57])+5.0*f[20]-4.0*f[18]+5.0*f[17]-2.0*(2.0*f[16]+3.0*(f[3]+f[1])))+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33]+5.0*f[32]))-10.0*(9.0*f[7]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[105]+f[88])+5.0*f[50]-4.0*f[39]+5.0*f[38]-4.0*f[37])+3.0*(3.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))/(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))+5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (-0.005*(15.0*(9.0*f[92]-4.0*(f[49]+f[47]+f[45])+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-1.0*(2.23606797749979*(9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57])+6.0*f[26])+5.0*(4.0*(f[20]+f[18])-5.0*f[17]+2.0*(2.0*f[16]+3.0*(f[5]+f[3]+f[1]))))+10.0*(9.0*(f[14]+f[12]+f[7])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (-0.001666666666666667*(405.0*f[108]-6.708203932499369*(9.0*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96])+5.0*((-9.0*(f[95]+f[88]))+4.0*(f[50]+f[39])-5.0*f[38]+4.0*f[37]))+3.0*(15.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-2.0*(6.708203932499369*(9.0*f[54]+5.0*(f[15]+f[11]+f[9]))+5.0*(9.0*(f[30]+f[28]+f[23])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(6.708203932499369*(9.0*f[101]+5.0*(f[46]+f[42]+f[40]))+5.0*(9.0*(f[79]+f[77]+f[66])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[7][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[7][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[7][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
// node (x,y)_9 
  vcutSq_i = 0.01*(2.23606797749979*(63.63961030678928*phiWall[19]-63.63961030678928*phi[19]+35.35533905932738*phiWall[9]-35.35533905932738*phi[9]+14.14213562373095*phiWall[8]-14.14213562373095*phi[8]+14.14213562373095*phiWall[7]-14.14213562373095*phi[7]+1.732050807568877*(21.21320343559643*phiWall[6]-21.21320343559643*phi[6]+21.21320343559643*phiWall[5]-21.21320343559643*phi[5])+21.21320343559643*phiWall[2]-21.21320343559643*phi[2]+21.21320343559643*phiWall[1]-21.21320343559643*phi[1])+1.732050807568877*(42.42640687119286*phiWall[18]-42.42640687119286*phi[18]+42.42640687119286*phiWall[17]-42.42640687119286*phi[17]+63.63961030678928*phiWall[10]-63.63961030678928*phi[10]+35.35533905932738*phiWall[3]-35.35533905932738*phi[3])+3.872983346207417*(3.872983346207417*(7.071067811865476*phiWall[16]-7.071067811865476*phi[16]+7.071067811865476*phiWall[15]-7.071067811865476*phi[15]+2.828427124746191*phiWall[12]-2.828427124746191*phi[12]+2.828427124746191*phiWall[11]-2.828427124746191*phi[11])+14.14213562373095*phiWall[14]-14.14213562373095*phi[14]+14.14213562373095*phiWall[13]-14.14213562373095*phi[13])+63.63961030678928*phiWall[4]-63.63961030678928*phi[4]+35.35533905932738*phiWall[0]-35.35533905932738*phi[0])*q2Dm; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[8][0] = 0.0; 
  fReflXYQuad[8][1] = 0.0; 
  fReflXYQuad[8][2] = 0.0; 
  fReflXYQuad[8][3] = 0.0; 
  fReflXYQuad[8][4] = 0.0; 
  fReflXYQuad[8][5] = 0.0; 
  fReflXYQuad[8][6] = 0.0; 
  fReflXYQuad[8][7] = 0.0; 
  fReflXYQuad[8][8] = 0.0; 
  fReflXYQuad[8][9] = 0.0; 
  fReflXYQuad[8][10] = 0.0; 
  fReflXYQuad[8][11] = 0.0; 
  fReflXYQuad[8][12] = 0.0; 
  fReflXYQuad[8][13] = 0.0; 
  fReflXYQuad[8][14] = 0.0; 
  fReflXYQuad[8][15] = 0.0; 
  fReflXYQuad[8][16] = 0.0; 
  fReflXYQuad[8][17] = 0.0; 
  fReflXYQuad[8][18] = 0.0; 
  fReflXYQuad[8][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[8][0] = 0.02*(2.23606797749979*(13.41640786499874*(f[32]+f[31])+5.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+5.0*(9.0*f[6]+5.0*f[0])); 
  fReflXYQuad[8][1] = 0.03333333333333333*(2.0*(9.0*(f[57]+f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(3.0*f[21]+2.23606797749979*(f[8]+f[7]))+5.0*f[3])); 
  fReflXYQuad[8][2] = 0.03333333333333333*(2.0*(9.0*(f[60]+f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4])); 
  fReflXYQuad[8][3] = 0.03333333333333333*(2.0*(9.0*(f[69]+f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(3.0*f[25]+2.23606797749979*(f[13]+f[12]))+5.0*f[5])); 
  fReflXYQuad[8][4] = 0.02*(2.23606797749979*(13.41640786499874*(f[88]+f[87])+5.0*(2.0*(f[62]+f[61])+3.0*(f[24]+f[23])))+5.0*(9.0*f[51]+5.0*f[11])); 
  fReflXYQuad[8][5] = 0.02*(2.23606797749979*(13.41640786499874*(f[92]+f[91])+5.0*(2.0*(f[71]+f[70])+3.0*(f[27]+f[26])))+5.0*(9.0*f[52]+5.0*f[14])); 
  fReflXYQuad[8][6] = 0.02*(2.23606797749979*(13.41640786499874*(f[95]+f[94])+5.0*(2.0*(f[75]+f[74])+3.0*(f[29]+f[28])))+5.0*(9.0*f[53]+5.0*f[15])); 
  fReflXYQuad[8][7] = 0.1*(9.0*f[58]+6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[8][8] = 0.1*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[8][9] = 0.1*(9.0*f[80]+6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[8][10] = 0.03333333333333333*(2.0*(9.0*(f[108]+f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(3.0*f[86]+2.23606797749979*(f[55]+f[54]))+5.0*f[30])); 
  fReflXYQuad[8][11] = 0.1*(9.0*f[89]+6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[8][12] = 0.1*(9.0*f[90]+6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[8][13] = 0.1*(9.0*f[93]+6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[8][14] = 0.1*(9.0*f[100]+6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[8][15] = 0.1*(9.0*f[103]+6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[8][16] = 0.1*(9.0*f[104]+6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[8][17] = 0.1*(9.0*f[109]+6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[8][18] = 0.1*(9.0*f[110]+6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[8][19] = 0.1*(9.0*f[111]+6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]-1.0*(f[92]+f[91]))+5.0*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31])))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))+3.0*((-1.0*(f[27]+f[26]))+f[25]+f[21]))+5.0*(3.0*(f[5]+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(5.0*(9.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7])-1.0*(9.0*f[6]+5.0*f[0]))-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])))-1.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]-1.0*(f[53]+f[51]))+5.0*((-1.0*(f[15]+f[11]))+f[10]+f[9]))))+5.0*(9.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*((-1.0*(f[46]+f[42]))+f[41]+f[40]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[67]+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+9.0*(5.0*(f[82]+f[81]-1.0*f[80])+2.0*(2.0*f[58]-1.0*(2.0*(f[57]+f[56])+3.0*f[21])))+5.0*(2.0*(2.0*(f[18]+f[17]+f[16])+3.0*((-1.0*f[3])+f[2]+f[1]))-5.0*f[20]))+10.0*(9.0*(f[6]-1.0*(f[8]+f[7]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+5.0*(4.0*(f[39]+f[38]+f[37])-5.0*f[50]))+3.0*(75.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(3.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-2.23606797749979*(9.0*f[51]+5.0*(f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[24]+f[23]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[67]+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*f[93]+f[92]+f[91])+5.0*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))))+5.0*(9.0*(2.0*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+3.0*(f[27]+f[26]-1.0*f[25]+f[21]))+5.0*(3.0*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1]))-2.0*(f[20]+f[18]+f[17]+f[16]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7])-1.0*(9.0*f[6]+5.0*f[0])))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87])))-5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+2.23606797749979*(9.0*(f[55]+f[54]-1.0*f[53]+f[51])+5.0*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23])-1.0*(9.0*f[22]+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*((-1.0*f[46])+f[42]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[78]+f[77])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(6.708203932499369*(9.0*f[93]+4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))+4.0*(f[32]+f[31]))+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]+6.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*((-1.0*f[5])+f[2]+f[1]))))+10.0*(9.0*(f[6]-1.0*(f[13]+f[12]))+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[109]+6.708203932499369*(9.0*(4.0*f[104]+5.0*(f[99]+f[98])-4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(3.0*(10.0*(f[60]+f[59])-2.23606797749979*(9.0*f[53]+5.0*(f[15]-1.0*(f[10]+f[9]))))+5.0*(9.0*(f[22]-1.0*(f[29]+f[28]))+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]-1.0*(f[41]+f[40])))+5.0*(9.0*(f[78]+f[77])-1.0*(9.0*f[65]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))/(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(5.0*(9.0*(f[80]+f[58])+5.0*(f[20]+f[18])-2.0*(2.0*(f[17]+f[16])+3.0*(f[2]+f[1])))+6.708203932499369*(5.0*(f[48]+f[47]+f[36]+f[35])-4.0*(f[32]+f[31])))-10.0*(9.0*f[6]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(6.708203932499369*(9.0*(f[104]+f[89])+5.0*(f[50]+f[39])-4.0*(f[38]+f[37]))+3.0*(15.0*(f[84]+f[83]+f[64]+f[63])-2.0*(3.0*(2.0*(f[60]+f[59])+3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))/(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))+5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(15.0*(9.0*f[93]-4.0*(f[48]+f[47])+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35])-4.0*(f[32]+f[31]))-1.0*(2.23606797749979*(9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58]+6.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]+2.0*(2.0*(f[17]+f[16])+3.0*(f[5]+f[2]+f[1]))))+10.0*(9.0*(f[13]+f[12]+f[6])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[109]-6.708203932499369*(9.0*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94]))+5.0*((-9.0*f[89])+4.0*f[50]-5.0*f[39]+4.0*(f[38]+f[37])))+3.0*(15.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(3.0*(10.0*(f[60]+f[59])+2.23606797749979*(9.0*f[53]+5.0*(f[15]+f[10]+f[9])))+5.0*(9.0*(f[29]+f[28]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(6.708203932499369*(9.0*f[100]+5.0*(f[46]+f[41]+f[40]))+5.0*(9.0*(f[78]+f[77]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]-1.0*(f[93]+f[92]+f[91]))+5.0*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80])-1.0*(2.0*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+3.0*(f[27]+f[26]+f[25]-1.0*f[21])))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*((-1.0*f[5])+f[3]+f[2]+f[1]))))+5.0*(5.0*(9.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+5.0*f[0])-81.0*f[52])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+6.708203932499369*(9.0*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-2.23606797749979*(9.0*(f[55]+f[54]+f[53]-1.0*f[51])+5.0*(f[15]-1.0*(f[11]+f[10]+f[9]))))+5.0*(9.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(f[46]-1.0*(f[42]+f[41]+f[40])))+5.0*(9.0*(f[79]+f[78]+f[77])-1.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))/(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))+5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(6.708203932499369*(9.0*f[103]+5.0*(f[49]+f[48]+f[47])-4.0*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+9.0*(5.0*(f[82]+f[81]+f[80])-2.0*(2.0*(f[58]+f[57]+f[56])+3.0*f[21]))+5.0*(5.0*f[20]-2.0*(2.0*(f[18]+f[17]+f[16])+3.0*(f[3]+f[2]+f[1]))))-10.0*(9.0*(f[8]+f[7]+f[6])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[111]+6.708203932499369*(9.0*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+5.0*(5.0*f[50]-4.0*(f[39]+f[38]+f[37])))+3.0*(75.0*(f[85]+f[84]+f[83])-2.0*(3.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*f[51]+5.0*(f[11]+f[10]+f[9])))+5.0*(9.0*(f[24]+f[23]+f[22])+5.0*f[4])))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(6.708203932499369*(9.0*f[90]+5.0*(f[42]+f[41]+f[40]))+5.0*(9.0*(f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))/(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))) <= 0. || fabs(xbarVal)>=.95) { 
  fReflXYZMuQuad[8][0] = 0.0; 
  fReflXYZMuQuad[8][1] = 0.0; 
  fReflXYZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))))*fac; 
    fac = (b>500 || fabs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3-(2*(b*b+3*(b+1))*exp(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][0] = (0.002*(2.23606797749979*(13.41640786499874*(9.0*(f[103]+f[93]+f[92]+f[91])+5.0*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31]))+5.0*(9.0*(2.0*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+3.0*(f[27]+f[26]+f[25]+f[21]))+5.0*(2.0*(f[20]+f[18]+f[17]+f[16])+3.0*(f[5]+f[3]+f[2]+f[1]))))+5.0*(81.0*f[52]+5.0*(9.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+5.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabs(b)<1e-8)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+6.708203932499369*(9.0*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87])+5.0*(f[50]+f[39]+f[38]+f[37])))+3.0*(3.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+2.23606797749979*(9.0*(f[55]+f[54]+f[53]+f[51])+5.0*(f[15]+f[11]+f[10]+f[9])))+5.0*(9.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*exp(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %G, b=%G, xbarVal=%G \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[8][2] = (0.01*(81.0*f[110]+6.708203932499369*(9.0*(f[102]+f[101]+f[100]+f[90])+5.0*(f[46]+f[42]+f[41]+f[40]))+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[8][0] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0])+8.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(8.0*fReflXYZMuQuad[3][0]+5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][1] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-1.0*(5.0*fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][2] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1])+8.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(8.0*fReflXYZMuQuad[3][1]+5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][3] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-1.0*fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][4] = 0.0414086662499961*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-1.0*(5.0*fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][5] = 0.2777777777777778*(fReflXYZMuQuad[8][0]-1.0*(fReflXYZMuQuad[6][0]+fReflXYZMuQuad[2][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][6] = 0.0414086662499961*(5.0*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-1.0*fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][7] = 0.0276057774999974*(5.0*fReflXYZMuQuad[8][0]+8.0*fReflXYZMuQuad[7][0]+5.0*fReflXYZMuQuad[6][0]-2.0*(5.0*fReflXYZMuQuad[5][0]+8.0*fReflXYZMuQuad[4][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[3][0])+8.0*fReflXYZMuQuad[1][0]+5.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][8] = 0.006172839506172839*(5.0*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2])+8.0*(5.0*fReflXYZMuQuad[5][2]+8.0*fReflXYZMuQuad[4][2])+5.0*(8.0*fReflXYZMuQuad[3][2]+5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][9] = 0.0276057774999974*(5.0*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0])+8.0*(fReflXYZMuQuad[5][0]-2.0*fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+5.0*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[8][10] = 0.2777777777777778*(fReflXYZMuQuad[8][1]-1.0*(fReflXYZMuQuad[6][1]+fReflXYZMuQuad[2][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][11] = 0.02760577749999742*(5.0*fReflXYZMuQuad[8][1]+8.0*fReflXYZMuQuad[7][1]+5.0*fReflXYZMuQuad[6][1]-2.0*(5.0*fReflXYZMuQuad[5][1]+8.0*fReflXYZMuQuad[4][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[3][1])+8.0*fReflXYZMuQuad[1][1]+5.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][12] = 0.04140866624999612*(5.0*fReflXYZMuQuad[8][2]+8.0*fReflXYZMuQuad[7][2]+5.0*fReflXYZMuQuad[6][2]-1.0*(5.0*fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+5.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][13] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-1.0*fReflXYZMuQuad[6][0]+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[5][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][14] = 0.04140866624999612*(5.0*(fReflXYZMuQuad[8][2]-1.0*fReflXYZMuQuad[6][2])+8.0*(fReflXYZMuQuad[5][2]-1.0*fReflXYZMuQuad[3][2])+5.0*(fReflXYZMuQuad[2][2]-1.0*fReflXYZMuQuad[0][2])); 
  fReflXYQuad[8][15] = 0.1851851851851853*(fReflXYZMuQuad[8][0]-2.0*fReflXYZMuQuad[7][0]+fReflXYZMuQuad[6][0]-1.0*fReflXYZMuQuad[2][0]+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[8][16] = 0.02760577749999742*(5.0*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1])+8.0*(fReflXYZMuQuad[5][1]-2.0*fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+5.0*(fReflXYZMuQuad[2][1]-2.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[8][17] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-1.0*fReflXYZMuQuad[6][1]+2.0*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[5][1])+fReflXYZMuQuad[2][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[8][18] = 0.2777777777777778*(fReflXYZMuQuad[8][2]-1.0*(fReflXYZMuQuad[6][2]+fReflXYZMuQuad[2][2])+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[8][19] = 0.1851851851851852*(fReflXYZMuQuad[8][1]-2.0*fReflXYZMuQuad[7][1]+fReflXYZMuQuad[6][1]-1.0*fReflXYZMuQuad[2][1]+2.0*fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  } 

 
  fRefl[0] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0])+8.0*(5.0*fReflXYQuad[5][0]+8.0*fReflXYQuad[4][0])+5.0*(8.0*fReflXYQuad[3][0]+5.0*fReflXYQuad[2][0]+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0])); 
  fRefl[1] = 0.0414086662499961*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0]-1.0*(5.0*fReflXYQuad[2][0]+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0])); 
  fRefl[2] = 0.0414086662499961*(5.0*(fReflXYQuad[8][0]-1.0*fReflXYQuad[6][0])+8.0*(fReflXYQuad[5][0]-1.0*fReflXYQuad[3][0])+5.0*(fReflXYQuad[2][0]-1.0*fReflXYQuad[0][0])); 
  fRefl[3] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1])+8.0*(5.0*fReflXYQuad[5][1]+8.0*fReflXYQuad[4][1])+5.0*(8.0*fReflXYQuad[3][1]+5.0*fReflXYQuad[2][1]+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1])); 
  fRefl[4] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2])+8.0*(5.0*fReflXYQuad[5][2]+8.0*fReflXYQuad[4][2])+5.0*(8.0*fReflXYQuad[3][2]+5.0*fReflXYQuad[2][2]+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2])); 
  fRefl[5] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3])+8.0*(5.0*fReflXYQuad[5][3]+8.0*fReflXYQuad[4][3])+5.0*(8.0*fReflXYQuad[3][3]+5.0*fReflXYQuad[2][3]+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3])); 
  fRefl[6] = 0.2777777777777778*(fReflXYQuad[8][0]-1.0*(fReflXYQuad[6][0]+fReflXYQuad[2][0])+fReflXYQuad[0][0]); 
  fRefl[7] = 0.0414086662499961*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1]-1.0*(5.0*fReflXYQuad[2][1]+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1])); 
  fRefl[8] = 0.0414086662499961*(5.0*(fReflXYQuad[8][1]-1.0*fReflXYQuad[6][1])+8.0*(fReflXYQuad[5][1]-1.0*fReflXYQuad[3][1])+5.0*(fReflXYQuad[2][1]-1.0*fReflXYQuad[0][1])); 
  fRefl[9] = 0.0414086662499961*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2]-1.0*(5.0*fReflXYQuad[2][2]+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2])); 
  fRefl[10] = 0.0414086662499961*(5.0*(fReflXYQuad[8][2]-1.0*fReflXYQuad[6][2])+8.0*(fReflXYQuad[5][2]-1.0*fReflXYQuad[3][2])+5.0*(fReflXYQuad[2][2]-1.0*fReflXYQuad[0][2])); 
  fRefl[11] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4])+8.0*(5.0*fReflXYQuad[5][4]+8.0*fReflXYQuad[4][4])+5.0*(8.0*fReflXYQuad[3][4]+5.0*fReflXYQuad[2][4]+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4])); 
  fRefl[12] = 0.0414086662499961*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3]-1.0*(5.0*fReflXYQuad[2][3]+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3])); 
  fRefl[13] = 0.0414086662499961*(5.0*(fReflXYQuad[8][3]-1.0*fReflXYQuad[6][3])+8.0*(fReflXYQuad[5][3]-1.0*fReflXYQuad[3][3])+5.0*(fReflXYQuad[2][3]-1.0*fReflXYQuad[0][3])); 
  fRefl[14] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5])+8.0*(5.0*fReflXYQuad[5][5]+8.0*fReflXYQuad[4][5])+5.0*(8.0*fReflXYQuad[3][5]+5.0*fReflXYQuad[2][5]+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5])); 
  fRefl[15] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6])+8.0*(5.0*fReflXYQuad[5][6]+8.0*fReflXYQuad[4][6])+5.0*(8.0*fReflXYQuad[3][6]+5.0*fReflXYQuad[2][6]+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6])); 
  fRefl[16] = 0.0276057774999974*(5.0*fReflXYQuad[8][0]+8.0*fReflXYQuad[7][0]+5.0*fReflXYQuad[6][0]-2.0*(5.0*fReflXYQuad[5][0]+8.0*fReflXYQuad[4][0])+5.0*(fReflXYQuad[2][0]-2.0*fReflXYQuad[3][0])+8.0*fReflXYQuad[1][0]+5.0*fReflXYQuad[0][0]); 
  fRefl[17] = 0.0276057774999974*(5.0*(fReflXYQuad[8][0]-2.0*fReflXYQuad[7][0]+fReflXYQuad[6][0])+8.0*(fReflXYQuad[5][0]-2.0*fReflXYQuad[4][0]+fReflXYQuad[3][0])+5.0*(fReflXYQuad[2][0]-2.0*fReflXYQuad[1][0]+fReflXYQuad[0][0])); 
  fRefl[18] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][7]+8.0*fReflXYQuad[7][7]+5.0*fReflXYQuad[6][7])+8.0*(5.0*fReflXYQuad[5][7]+8.0*fReflXYQuad[4][7])+5.0*(8.0*fReflXYQuad[3][7]+5.0*fReflXYQuad[2][7]+8.0*fReflXYQuad[1][7]+5.0*fReflXYQuad[0][7])); 
  fRefl[19] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][8]+8.0*fReflXYQuad[7][8]+5.0*fReflXYQuad[6][8])+8.0*(5.0*fReflXYQuad[5][8]+8.0*fReflXYQuad[4][8])+5.0*(8.0*fReflXYQuad[3][8]+5.0*fReflXYQuad[2][8]+8.0*fReflXYQuad[1][8]+5.0*fReflXYQuad[0][8])); 
  fRefl[20] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][9]+8.0*fReflXYQuad[7][9]+5.0*fReflXYQuad[6][9])+8.0*(5.0*fReflXYQuad[5][9]+8.0*fReflXYQuad[4][9])+5.0*(8.0*fReflXYQuad[3][9]+5.0*fReflXYQuad[2][9]+8.0*fReflXYQuad[1][9]+5.0*fReflXYQuad[0][9])); 
  fRefl[21] = 0.2777777777777778*(fReflXYQuad[8][1]-1.0*(fReflXYQuad[6][1]+fReflXYQuad[2][1])+fReflXYQuad[0][1]); 
  fRefl[22] = 0.2777777777777778*(fReflXYQuad[8][2]-1.0*(fReflXYQuad[6][2]+fReflXYQuad[2][2])+fReflXYQuad[0][2]); 
  fRefl[23] = 0.0414086662499961*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4]-1.0*(5.0*fReflXYQuad[2][4]+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4])); 
  fRefl[24] = 0.0414086662499961*(5.0*(fReflXYQuad[8][4]-1.0*fReflXYQuad[6][4])+8.0*(fReflXYQuad[5][4]-1.0*fReflXYQuad[3][4])+5.0*(fReflXYQuad[2][4]-1.0*fReflXYQuad[0][4])); 
  fRefl[25] = 0.2777777777777778*(fReflXYQuad[8][3]-1.0*(fReflXYQuad[6][3]+fReflXYQuad[2][3])+fReflXYQuad[0][3]); 
  fRefl[26] = 0.0414086662499961*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5]-1.0*(5.0*fReflXYQuad[2][5]+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5])); 
  fRefl[27] = 0.0414086662499961*(5.0*(fReflXYQuad[8][5]-1.0*fReflXYQuad[6][5])+8.0*(fReflXYQuad[5][5]-1.0*fReflXYQuad[3][5])+5.0*(fReflXYQuad[2][5]-1.0*fReflXYQuad[0][5])); 
  fRefl[28] = 0.0414086662499961*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6]-1.0*(5.0*fReflXYQuad[2][6]+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6])); 
  fRefl[29] = 0.0414086662499961*(5.0*(fReflXYQuad[8][6]-1.0*fReflXYQuad[6][6])+8.0*(fReflXYQuad[5][6]-1.0*fReflXYQuad[3][6])+5.0*(fReflXYQuad[2][6]-1.0*fReflXYQuad[0][6])); 
  fRefl[30] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10])+8.0*(5.0*fReflXYQuad[5][10]+8.0*fReflXYQuad[4][10])+5.0*(8.0*fReflXYQuad[3][10]+5.0*fReflXYQuad[2][10]+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10])); 
  fRefl[31] = 0.1851851851851853*(fReflXYQuad[8][0]-1.0*fReflXYQuad[6][0]+2.0*(fReflXYQuad[3][0]-1.0*fReflXYQuad[5][0])+fReflXYQuad[2][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[32] = 0.1851851851851853*(fReflXYQuad[8][0]-2.0*fReflXYQuad[7][0]+fReflXYQuad[6][0]-1.0*fReflXYQuad[2][0]+2.0*fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[33] = 0.02760577749999742*(5.0*fReflXYQuad[8][1]+8.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[6][1]-2.0*(5.0*fReflXYQuad[5][1]+8.0*fReflXYQuad[4][1])+5.0*(fReflXYQuad[2][1]-2.0*fReflXYQuad[3][1])+8.0*fReflXYQuad[1][1]+5.0*fReflXYQuad[0][1]); 
  fRefl[34] = 0.02760577749999742*(5.0*(fReflXYQuad[8][1]-2.0*fReflXYQuad[7][1]+fReflXYQuad[6][1])+8.0*(fReflXYQuad[5][1]-2.0*fReflXYQuad[4][1]+fReflXYQuad[3][1])+5.0*(fReflXYQuad[2][1]-2.0*fReflXYQuad[1][1]+fReflXYQuad[0][1])); 
  fRefl[35] = 0.04140866624999612*(5.0*fReflXYQuad[8][7]+8.0*fReflXYQuad[7][7]+5.0*fReflXYQuad[6][7]-1.0*(5.0*fReflXYQuad[2][7]+8.0*fReflXYQuad[1][7]+5.0*fReflXYQuad[0][7])); 
  fRefl[36] = 0.04140866624999612*(5.0*(fReflXYQuad[8][7]-1.0*fReflXYQuad[6][7])+8.0*(fReflXYQuad[5][7]-1.0*fReflXYQuad[3][7])+5.0*(fReflXYQuad[2][7]-1.0*fReflXYQuad[0][7])); 
  fRefl[37] = 0.02760577749999742*(5.0*fReflXYQuad[8][2]+8.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[6][2]-2.0*(5.0*fReflXYQuad[5][2]+8.0*fReflXYQuad[4][2])+5.0*(fReflXYQuad[2][2]-2.0*fReflXYQuad[3][2])+8.0*fReflXYQuad[1][2]+5.0*fReflXYQuad[0][2]); 
  fRefl[38] = 0.02760577749999742*(5.0*(fReflXYQuad[8][2]-2.0*fReflXYQuad[7][2]+fReflXYQuad[6][2])+8.0*(fReflXYQuad[5][2]-2.0*fReflXYQuad[4][2]+fReflXYQuad[3][2])+5.0*(fReflXYQuad[2][2]-2.0*fReflXYQuad[1][2]+fReflXYQuad[0][2])); 
  fRefl[39] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][11]+8.0*fReflXYQuad[7][11]+5.0*fReflXYQuad[6][11])+8.0*(5.0*fReflXYQuad[5][11]+8.0*fReflXYQuad[4][11])+5.0*(8.0*fReflXYQuad[3][11]+5.0*fReflXYQuad[2][11]+8.0*fReflXYQuad[1][11]+5.0*fReflXYQuad[0][11])); 
  fRefl[40] = 0.04140866624999612*(5.0*fReflXYQuad[8][8]+8.0*fReflXYQuad[7][8]+5.0*fReflXYQuad[6][8]-1.0*(5.0*fReflXYQuad[2][8]+8.0*fReflXYQuad[1][8]+5.0*fReflXYQuad[0][8])); 
  fRefl[41] = 0.04140866624999612*(5.0*(fReflXYQuad[8][8]-1.0*fReflXYQuad[6][8])+8.0*(fReflXYQuad[5][8]-1.0*fReflXYQuad[3][8])+5.0*(fReflXYQuad[2][8]-1.0*fReflXYQuad[0][8])); 
  fRefl[42] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][12]+8.0*fReflXYQuad[7][12]+5.0*fReflXYQuad[6][12])+8.0*(5.0*fReflXYQuad[5][12]+8.0*fReflXYQuad[4][12])+5.0*(8.0*fReflXYQuad[3][12]+5.0*fReflXYQuad[2][12]+8.0*fReflXYQuad[1][12]+5.0*fReflXYQuad[0][12])); 
  fRefl[43] = 0.02760577749999742*(5.0*fReflXYQuad[8][3]+8.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[6][3]-2.0*(5.0*fReflXYQuad[5][3]+8.0*fReflXYQuad[4][3])+5.0*(fReflXYQuad[2][3]-2.0*fReflXYQuad[3][3])+8.0*fReflXYQuad[1][3]+5.0*fReflXYQuad[0][3]); 
  fRefl[44] = 0.02760577749999742*(5.0*(fReflXYQuad[8][3]-2.0*fReflXYQuad[7][3]+fReflXYQuad[6][3])+8.0*(fReflXYQuad[5][3]-2.0*fReflXYQuad[4][3]+fReflXYQuad[3][3])+5.0*(fReflXYQuad[2][3]-2.0*fReflXYQuad[1][3]+fReflXYQuad[0][3])); 
  fRefl[45] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][13]+8.0*fReflXYQuad[7][13]+5.0*fReflXYQuad[6][13])+8.0*(5.0*fReflXYQuad[5][13]+8.0*fReflXYQuad[4][13])+5.0*(8.0*fReflXYQuad[3][13]+5.0*fReflXYQuad[2][13]+8.0*fReflXYQuad[1][13]+5.0*fReflXYQuad[0][13])); 
  fRefl[46] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][14]+8.0*fReflXYQuad[7][14]+5.0*fReflXYQuad[6][14])+8.0*(5.0*fReflXYQuad[5][14]+8.0*fReflXYQuad[4][14])+5.0*(8.0*fReflXYQuad[3][14]+5.0*fReflXYQuad[2][14]+8.0*fReflXYQuad[1][14]+5.0*fReflXYQuad[0][14])); 
  fRefl[47] = 0.04140866624999612*(5.0*fReflXYQuad[8][9]+8.0*fReflXYQuad[7][9]+5.0*fReflXYQuad[6][9]-1.0*(5.0*fReflXYQuad[2][9]+8.0*fReflXYQuad[1][9]+5.0*fReflXYQuad[0][9])); 
  fRefl[48] = 0.04140866624999612*(5.0*(fReflXYQuad[8][9]-1.0*fReflXYQuad[6][9])+8.0*(fReflXYQuad[5][9]-1.0*fReflXYQuad[3][9])+5.0*(fReflXYQuad[2][9]-1.0*fReflXYQuad[0][9])); 
  fRefl[49] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][15]+8.0*fReflXYQuad[7][15]+5.0*fReflXYQuad[6][15])+8.0*(5.0*fReflXYQuad[5][15]+8.0*fReflXYQuad[4][15])+5.0*(8.0*fReflXYQuad[3][15]+5.0*fReflXYQuad[2][15]+8.0*fReflXYQuad[1][15]+5.0*fReflXYQuad[0][15])); 
  fRefl[50] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][16]+8.0*fReflXYQuad[7][16]+5.0*fReflXYQuad[6][16])+8.0*(5.0*fReflXYQuad[5][16]+8.0*fReflXYQuad[4][16])+5.0*(8.0*fReflXYQuad[3][16]+5.0*fReflXYQuad[2][16]+8.0*fReflXYQuad[1][16]+5.0*fReflXYQuad[0][16])); 
  fRefl[51] = 0.2777777777777778*(fReflXYQuad[8][4]-1.0*(fReflXYQuad[6][4]+fReflXYQuad[2][4])+fReflXYQuad[0][4]); 
  fRefl[52] = 0.2777777777777778*(fReflXYQuad[8][5]-1.0*(fReflXYQuad[6][5]+fReflXYQuad[2][5])+fReflXYQuad[0][5]); 
  fRefl[53] = 0.2777777777777778*(fReflXYQuad[8][6]-1.0*(fReflXYQuad[6][6]+fReflXYQuad[2][6])+fReflXYQuad[0][6]); 
  fRefl[54] = 0.0414086662499961*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10]-1.0*(5.0*fReflXYQuad[2][10]+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10])); 
  fRefl[55] = 0.0414086662499961*(5.0*(fReflXYQuad[8][10]-1.0*fReflXYQuad[6][10])+8.0*(fReflXYQuad[5][10]-1.0*fReflXYQuad[3][10])+5.0*(fReflXYQuad[2][10]-1.0*fReflXYQuad[0][10])); 
  fRefl[56] = 0.1851851851851852*(fReflXYQuad[8][1]-1.0*fReflXYQuad[6][1]+2.0*(fReflXYQuad[3][1]-1.0*fReflXYQuad[5][1])+fReflXYQuad[2][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[57] = 0.1851851851851852*(fReflXYQuad[8][1]-2.0*fReflXYQuad[7][1]+fReflXYQuad[6][1]-1.0*fReflXYQuad[2][1]+2.0*fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[58] = 0.2777777777777778*(fReflXYQuad[8][7]-1.0*(fReflXYQuad[6][7]+fReflXYQuad[2][7])+fReflXYQuad[0][7]); 
  fRefl[59] = 0.1851851851851852*(fReflXYQuad[8][2]-1.0*fReflXYQuad[6][2]+2.0*(fReflXYQuad[3][2]-1.0*fReflXYQuad[5][2])+fReflXYQuad[2][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[60] = 0.1851851851851852*(fReflXYQuad[8][2]-2.0*fReflXYQuad[7][2]+fReflXYQuad[6][2]-1.0*fReflXYQuad[2][2]+2.0*fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[61] = 0.0276057774999974*(5.0*fReflXYQuad[8][4]+8.0*fReflXYQuad[7][4]+5.0*fReflXYQuad[6][4]-2.0*(5.0*fReflXYQuad[5][4]+8.0*fReflXYQuad[4][4])+5.0*(fReflXYQuad[2][4]-2.0*fReflXYQuad[3][4])+8.0*fReflXYQuad[1][4]+5.0*fReflXYQuad[0][4]); 
  fRefl[62] = 0.0276057774999974*(5.0*(fReflXYQuad[8][4]-2.0*fReflXYQuad[7][4]+fReflXYQuad[6][4])+8.0*(fReflXYQuad[5][4]-2.0*fReflXYQuad[4][4]+fReflXYQuad[3][4])+5.0*(fReflXYQuad[2][4]-2.0*fReflXYQuad[1][4]+fReflXYQuad[0][4])); 
  fRefl[63] = 0.04140866624999612*(5.0*fReflXYQuad[8][11]+8.0*fReflXYQuad[7][11]+5.0*fReflXYQuad[6][11]-1.0*(5.0*fReflXYQuad[2][11]+8.0*fReflXYQuad[1][11]+5.0*fReflXYQuad[0][11])); 
  fRefl[64] = 0.04140866624999612*(5.0*(fReflXYQuad[8][11]-1.0*fReflXYQuad[6][11])+8.0*(fReflXYQuad[5][11]-1.0*fReflXYQuad[3][11])+5.0*(fReflXYQuad[2][11]-1.0*fReflXYQuad[0][11])); 
  fRefl[65] = 0.2777777777777778*(fReflXYQuad[8][8]-1.0*(fReflXYQuad[6][8]+fReflXYQuad[2][8])+fReflXYQuad[0][8]); 
  fRefl[66] = 0.04140866624999612*(5.0*fReflXYQuad[8][12]+8.0*fReflXYQuad[7][12]+5.0*fReflXYQuad[6][12]-1.0*(5.0*fReflXYQuad[2][12]+8.0*fReflXYQuad[1][12]+5.0*fReflXYQuad[0][12])); 
  fRefl[67] = 0.04140866624999612*(5.0*(fReflXYQuad[8][12]-1.0*fReflXYQuad[6][12])+8.0*(fReflXYQuad[5][12]-1.0*fReflXYQuad[3][12])+5.0*(fReflXYQuad[2][12]-1.0*fReflXYQuad[0][12])); 
  fRefl[68] = 0.1851851851851852*(fReflXYQuad[8][3]-1.0*fReflXYQuad[6][3]+2.0*(fReflXYQuad[3][3]-1.0*fReflXYQuad[5][3])+fReflXYQuad[2][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[69] = 0.1851851851851852*(fReflXYQuad[8][3]-2.0*fReflXYQuad[7][3]+fReflXYQuad[6][3]-1.0*fReflXYQuad[2][3]+2.0*fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[70] = 0.0276057774999974*(5.0*fReflXYQuad[8][5]+8.0*fReflXYQuad[7][5]+5.0*fReflXYQuad[6][5]-2.0*(5.0*fReflXYQuad[5][5]+8.0*fReflXYQuad[4][5])+5.0*(fReflXYQuad[2][5]-2.0*fReflXYQuad[3][5])+8.0*fReflXYQuad[1][5]+5.0*fReflXYQuad[0][5]); 
  fRefl[71] = 0.0276057774999974*(5.0*(fReflXYQuad[8][5]-2.0*fReflXYQuad[7][5]+fReflXYQuad[6][5])+8.0*(fReflXYQuad[5][5]-2.0*fReflXYQuad[4][5]+fReflXYQuad[3][5])+5.0*(fReflXYQuad[2][5]-2.0*fReflXYQuad[1][5]+fReflXYQuad[0][5])); 
  fRefl[72] = 0.04140866624999612*(5.0*fReflXYQuad[8][13]+8.0*fReflXYQuad[7][13]+5.0*fReflXYQuad[6][13]-1.0*(5.0*fReflXYQuad[2][13]+8.0*fReflXYQuad[1][13]+5.0*fReflXYQuad[0][13])); 
  fRefl[73] = 0.04140866624999612*(5.0*(fReflXYQuad[8][13]-1.0*fReflXYQuad[6][13])+8.0*(fReflXYQuad[5][13]-1.0*fReflXYQuad[3][13])+5.0*(fReflXYQuad[2][13]-1.0*fReflXYQuad[0][13])); 
  fRefl[74] = 0.0276057774999974*(5.0*fReflXYQuad[8][6]+8.0*fReflXYQuad[7][6]+5.0*fReflXYQuad[6][6]-2.0*(5.0*fReflXYQuad[5][6]+8.0*fReflXYQuad[4][6])+5.0*(fReflXYQuad[2][6]-2.0*fReflXYQuad[3][6])+8.0*fReflXYQuad[1][6]+5.0*fReflXYQuad[0][6]); 
  fRefl[75] = 0.0276057774999974*(5.0*(fReflXYQuad[8][6]-2.0*fReflXYQuad[7][6]+fReflXYQuad[6][6])+8.0*(fReflXYQuad[5][6]-2.0*fReflXYQuad[4][6]+fReflXYQuad[3][6])+5.0*(fReflXYQuad[2][6]-2.0*fReflXYQuad[1][6]+fReflXYQuad[0][6])); 
  fRefl[76] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][17]+8.0*fReflXYQuad[7][17]+5.0*fReflXYQuad[6][17])+8.0*(5.0*fReflXYQuad[5][17]+8.0*fReflXYQuad[4][17])+5.0*(8.0*fReflXYQuad[3][17]+5.0*fReflXYQuad[2][17]+8.0*fReflXYQuad[1][17]+5.0*fReflXYQuad[0][17])); 
  fRefl[77] = 0.04140866624999612*(5.0*fReflXYQuad[8][14]+8.0*fReflXYQuad[7][14]+5.0*fReflXYQuad[6][14]-1.0*(5.0*fReflXYQuad[2][14]+8.0*fReflXYQuad[1][14]+5.0*fReflXYQuad[0][14])); 
  fRefl[78] = 0.04140866624999612*(5.0*(fReflXYQuad[8][14]-1.0*fReflXYQuad[6][14])+8.0*(fReflXYQuad[5][14]-1.0*fReflXYQuad[3][14])+5.0*(fReflXYQuad[2][14]-1.0*fReflXYQuad[0][14])); 
  fRefl[79] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][18]+8.0*fReflXYQuad[7][18]+5.0*fReflXYQuad[6][18])+8.0*(5.0*fReflXYQuad[5][18]+8.0*fReflXYQuad[4][18])+5.0*(8.0*fReflXYQuad[3][18]+5.0*fReflXYQuad[2][18]+8.0*fReflXYQuad[1][18]+5.0*fReflXYQuad[0][18])); 
  fRefl[80] = 0.2777777777777778*(fReflXYQuad[8][9]-1.0*(fReflXYQuad[6][9]+fReflXYQuad[2][9])+fReflXYQuad[0][9]); 
  fRefl[81] = 0.04140866624999612*(5.0*fReflXYQuad[8][15]+8.0*fReflXYQuad[7][15]+5.0*fReflXYQuad[6][15]-1.0*(5.0*fReflXYQuad[2][15]+8.0*fReflXYQuad[1][15]+5.0*fReflXYQuad[0][15])); 
  fRefl[82] = 0.04140866624999612*(5.0*(fReflXYQuad[8][15]-1.0*fReflXYQuad[6][15])+8.0*(fReflXYQuad[5][15]-1.0*fReflXYQuad[3][15])+5.0*(fReflXYQuad[2][15]-1.0*fReflXYQuad[0][15])); 
  fRefl[83] = 0.04140866624999612*(5.0*fReflXYQuad[8][16]+8.0*fReflXYQuad[7][16]+5.0*fReflXYQuad[6][16]-1.0*(5.0*fReflXYQuad[2][16]+8.0*fReflXYQuad[1][16]+5.0*fReflXYQuad[0][16])); 
  fRefl[84] = 0.04140866624999612*(5.0*(fReflXYQuad[8][16]-1.0*fReflXYQuad[6][16])+8.0*(fReflXYQuad[5][16]-1.0*fReflXYQuad[3][16])+5.0*(fReflXYQuad[2][16]-1.0*fReflXYQuad[0][16])); 
  fRefl[85] = 0.006172839506172839*(5.0*(5.0*fReflXYQuad[8][19]+8.0*fReflXYQuad[7][19]+5.0*fReflXYQuad[6][19])+8.0*(5.0*fReflXYQuad[5][19]+8.0*fReflXYQuad[4][19])+5.0*(8.0*fReflXYQuad[3][19]+5.0*fReflXYQuad[2][19]+8.0*fReflXYQuad[1][19]+5.0*fReflXYQuad[0][19])); 
  fRefl[86] = 0.2777777777777778*(fReflXYQuad[8][10]-1.0*(fReflXYQuad[6][10]+fReflXYQuad[2][10])+fReflXYQuad[0][10]); 
  fRefl[87] = 0.1851851851851853*(fReflXYQuad[8][4]-1.0*fReflXYQuad[6][4]+2.0*(fReflXYQuad[3][4]-1.0*fReflXYQuad[5][4])+fReflXYQuad[2][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[88] = 0.1851851851851853*(fReflXYQuad[8][4]-2.0*fReflXYQuad[7][4]+fReflXYQuad[6][4]-1.0*fReflXYQuad[2][4]+2.0*fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[89] = 0.2777777777777778*(fReflXYQuad[8][11]-1.0*(fReflXYQuad[6][11]+fReflXYQuad[2][11])+fReflXYQuad[0][11]); 
  fRefl[90] = 0.2777777777777778*(fReflXYQuad[8][12]-1.0*(fReflXYQuad[6][12]+fReflXYQuad[2][12])+fReflXYQuad[0][12]); 
  fRefl[91] = 0.1851851851851853*(fReflXYQuad[8][5]-1.0*fReflXYQuad[6][5]+2.0*(fReflXYQuad[3][5]-1.0*fReflXYQuad[5][5])+fReflXYQuad[2][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[92] = 0.1851851851851853*(fReflXYQuad[8][5]-2.0*fReflXYQuad[7][5]+fReflXYQuad[6][5]-1.0*fReflXYQuad[2][5]+2.0*fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[93] = 0.2777777777777778*(fReflXYQuad[8][13]-1.0*(fReflXYQuad[6][13]+fReflXYQuad[2][13])+fReflXYQuad[0][13]); 
  fRefl[94] = 0.1851851851851853*(fReflXYQuad[8][6]-1.0*fReflXYQuad[6][6]+2.0*(fReflXYQuad[3][6]-1.0*fReflXYQuad[5][6])+fReflXYQuad[2][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[95] = 0.1851851851851853*(fReflXYQuad[8][6]-2.0*fReflXYQuad[7][6]+fReflXYQuad[6][6]-1.0*fReflXYQuad[2][6]+2.0*fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[96] = 0.02760577749999742*(5.0*fReflXYQuad[8][10]+8.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[6][10]-2.0*(5.0*fReflXYQuad[5][10]+8.0*fReflXYQuad[4][10])+5.0*(fReflXYQuad[2][10]-2.0*fReflXYQuad[3][10])+8.0*fReflXYQuad[1][10]+5.0*fReflXYQuad[0][10]); 
  fRefl[97] = 0.02760577749999742*(5.0*(fReflXYQuad[8][10]-2.0*fReflXYQuad[7][10]+fReflXYQuad[6][10])+8.0*(fReflXYQuad[5][10]-2.0*fReflXYQuad[4][10]+fReflXYQuad[3][10])+5.0*(fReflXYQuad[2][10]-2.0*fReflXYQuad[1][10]+fReflXYQuad[0][10])); 
  fRefl[98] = 0.04140866624999612*(5.0*fReflXYQuad[8][17]+8.0*fReflXYQuad[7][17]+5.0*fReflXYQuad[6][17]-1.0*(5.0*fReflXYQuad[2][17]+8.0*fReflXYQuad[1][17]+5.0*fReflXYQuad[0][17])); 
  fRefl[99] = 0.04140866624999612*(5.0*(fReflXYQuad[8][17]-1.0*fReflXYQuad[6][17])+8.0*(fReflXYQuad[5][17]-1.0*fReflXYQuad[3][17])+5.0*(fReflXYQuad[2][17]-1.0*fReflXYQuad[0][17])); 
  fRefl[100] = 0.2777777777777778*(fReflXYQuad[8][14]-1.0*(fReflXYQuad[6][14]+fReflXYQuad[2][14])+fReflXYQuad[0][14]); 
  fRefl[101] = 0.04140866624999612*(5.0*fReflXYQuad[8][18]+8.0*fReflXYQuad[7][18]+5.0*fReflXYQuad[6][18]-1.0*(5.0*fReflXYQuad[2][18]+8.0*fReflXYQuad[1][18]+5.0*fReflXYQuad[0][18])); 
  fRefl[102] = 0.04140866624999612*(5.0*(fReflXYQuad[8][18]-1.0*fReflXYQuad[6][18])+8.0*(fReflXYQuad[5][18]-1.0*fReflXYQuad[3][18])+5.0*(fReflXYQuad[2][18]-1.0*fReflXYQuad[0][18])); 
  fRefl[103] = 0.2777777777777778*(fReflXYQuad[8][15]-1.0*(fReflXYQuad[6][15]+fReflXYQuad[2][15])+fReflXYQuad[0][15]); 
  fRefl[104] = 0.2777777777777778*(fReflXYQuad[8][16]-1.0*(fReflXYQuad[6][16]+fReflXYQuad[2][16])+fReflXYQuad[0][16]); 
  fRefl[105] = 0.04140866624999612*(5.0*fReflXYQuad[8][19]+8.0*fReflXYQuad[7][19]+5.0*fReflXYQuad[6][19]-1.0*(5.0*fReflXYQuad[2][19]+8.0*fReflXYQuad[1][19]+5.0*fReflXYQuad[0][19])); 
  fRefl[106] = 0.04140866624999612*(5.0*(fReflXYQuad[8][19]-1.0*fReflXYQuad[6][19])+8.0*(fReflXYQuad[5][19]-1.0*fReflXYQuad[3][19])+5.0*(fReflXYQuad[2][19]-1.0*fReflXYQuad[0][19])); 
  fRefl[107] = 0.1851851851851852*(fReflXYQuad[8][10]-1.0*fReflXYQuad[6][10]+2.0*(fReflXYQuad[3][10]-1.0*fReflXYQuad[5][10])+fReflXYQuad[2][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[108] = 0.1851851851851852*(fReflXYQuad[8][10]-2.0*fReflXYQuad[7][10]+fReflXYQuad[6][10]-1.0*fReflXYQuad[2][10]+2.0*fReflXYQuad[1][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[109] = 0.2777777777777778*(fReflXYQuad[8][17]-1.0*(fReflXYQuad[6][17]+fReflXYQuad[2][17])+fReflXYQuad[0][17]); 
  fRefl[110] = 0.2777777777777778*(fReflXYQuad[8][18]-1.0*(fReflXYQuad[6][18]+fReflXYQuad[2][18])+fReflXYQuad[0][18]); 
  fRefl[111] = 0.2777777777777778*(fReflXYQuad[8][19]-1.0*(fReflXYQuad[6][19]+fReflXYQuad[2][19])+fReflXYQuad[0][19]); 
}
