#include <gkyl_bc_sheath_gyrokinetic.h> 


void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
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
  if(-0.7071067811865475*(f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[3]+f[2]))/(f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.7071067811865475*(f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
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

void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
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
  if(-0.7071067811865475*(f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (-0.7071067811865475*(f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (-0.7071067811865475*(f[3]-1.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.04714045207910316*(15.0*f[5]-15.0*f[4]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[3]+f[2]))/(f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.7071067811865475*(f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.04714045207910316*(15.0*f[5]+15.0*f[4]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (0.7071067811865475*(f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (0.7071067811865475*(f[3]+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
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

void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
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
  if (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))/(f[5]-1.0*f[3]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]-1.0*(f[4]+f[2])))/(f[5]+f[3]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]+f[4]+f[2]))/(f[5]+f[3]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]+f[3]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
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

void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
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
  if (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.5*(f[5]-1.0*(f[3]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.5*(f[7]-1.0*(f[6]+f[4])+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(15.0*f[11]-15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))/(f[5]-1.0*f[3]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.5*(f[5]-1.0*f[3]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.5*(f[7]-1.0*f[6]+f[4]-1.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(15.0*f[11]+15.0*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]-1.0*(f[4]+f[2])))/(f[5]+f[3]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.5*(f[5]+f[3]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.5*(f[7]+f[6]-1.0*(f[4]+f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(15.0*f[11]+15.0*f[10]-1.0*(15.0*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+f[6]+f[4]+f[2]))/(f[5]+f[3]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.5*(f[5]+f[3]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(15.0*f[11]+15.0*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.5*(f[5]+f[3]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.5*(f[7]+f[6]+f[4]+f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
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

void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; long double xc, b, xbarVal, fac; 
  double fReflXYQuad[4][12]; 
  double fReflXYZMuQuad[4][6]; 
  

// node (x,y)_1 
  vcutSq_i = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*phi[6]-2.449489742783178*phiWall[5]+2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*phi[2]+1.414213562373095*phiWall[1]-1.414213562373095*(phi[1]+phiWall[0])+1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = 0.5*(f[6]-1.0*(f[2]+f[1])+f[0]); 
  fReflXYQuad[0][1] = 0.5*(f[16]-1.0*(f[8]+f[7])+f[3]); 
  fReflXYQuad[0][2] = 0.5*(f[17]-1.0*(f[10]+f[9])+f[4]); 
  fReflXYQuad[0][3] = 0.5*(f[20]-1.0*(f[13]+f[12])+f[5]); 
  fReflXYQuad[0][4] = 0.5*(f[26]-1.0*(f[19]+f[18])+f[11]); 
  fReflXYQuad[0][5] = 0.5*(f[27]-1.0*(f[22]+f[21])+f[14]); 
  fReflXYQuad[0][6] = 0.5*(f[28]-1.0*(f[24]+f[23])+f[15]); 
  fReflXYQuad[0][7] = 0.5*(f[31]-1.0*(f[30]+f[29])+f[25]); 
  fReflXYQuad[0][8] = 0.03333333333333333*(15.0*f[37]-15.0*(f[34]+f[33])+15.0*f[32]); 
  fReflXYQuad[0][9] = 0.03333333333333333*(15.0*f[43]-15.0*(f[39]+f[38])+15.0*f[35]); 
  fReflXYQuad[0][10] = 0.03333333333333333*(15.0*f[44]-15.0*(f[41]+f[40])+15.0*f[36]); 
  fReflXYQuad[0][11] = 0.03333333333333333*(15.0*f[47]-15.0*(f[46]+f[45])+15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))/(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]-1.0*f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[36]-1.0*f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]-1.0*f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[36]-1.0*f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))/(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*(f[36]+f[35]-1.0*(f[34]+f[33]))+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*(f[36]+f[35]-1.0*(f[34]+f[33]))+15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_2 
  vcutSq_i = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*(phi[6]+phiWall[5])-2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*(phi[2]+phiWall[1])+1.414213562373095*(phi[1]+phiWall[0])-1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = -0.5*(f[6]-1.0*f[2]+f[1]-1.0*f[0]); 
  fReflXYQuad[1][1] = -0.5*(f[16]-1.0*f[8]+f[7]-1.0*f[3]); 
  fReflXYQuad[1][2] = -0.5*(f[17]-1.0*f[10]+f[9]-1.0*f[4]); 
  fReflXYQuad[1][3] = -0.5*(f[20]-1.0*f[13]+f[12]-1.0*f[5]); 
  fReflXYQuad[1][4] = -0.5*(f[26]-1.0*f[19]+f[18]-1.0*f[11]); 
  fReflXYQuad[1][5] = -0.5*(f[27]-1.0*f[22]+f[21]-1.0*f[14]); 
  fReflXYQuad[1][6] = -0.5*(f[28]-1.0*f[24]+f[23]-1.0*f[15]); 
  fReflXYQuad[1][7] = -0.5*(f[31]-1.0*f[30]+f[29]-1.0*f[25]); 
  fReflXYQuad[1][8] = -0.03333333333333333*(15.0*f[37]+15.0*(f[33]-1.0*f[34])-15.0*f[32]); 
  fReflXYQuad[1][9] = -0.03333333333333333*(15.0*f[43]+15.0*(f[38]-1.0*f[39])-15.0*f[35]); 
  fReflXYQuad[1][10] = -0.03333333333333333*(15.0*f[44]+15.0*(f[40]-1.0*f[41])-15.0*f[36]); 
  fReflXYQuad[1][11] = -0.03333333333333333*(15.0*f[47]+15.0*(f[45]-1.0*f[46])-15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*(f[44]+f[43]))+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[36]+f[35]-1.0*f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*(f[44]+f[43]))+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[36]+f[35]-1.0*f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))/(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*(f[36]-1.0*f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*(f[36]-1.0*f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))/(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]-1.0*f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]-1.0*f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))/(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]+f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]-1.0*f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]+f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]-1.0*f[33])-15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_3 
  vcutSq_i = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*(phi[6]+phiWall[5])+2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*(phi[2]+phiWall[1])-1.414213562373095*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.5*(f[6]+f[2]-1.0*(f[1]+f[0])); 
  fReflXYQuad[2][1] = -0.5*(f[16]+f[8]-1.0*(f[7]+f[3])); 
  fReflXYQuad[2][2] = -0.5*(f[17]+f[10]-1.0*(f[9]+f[4])); 
  fReflXYQuad[2][3] = -0.5*(f[20]+f[13]-1.0*(f[12]+f[5])); 
  fReflXYQuad[2][4] = -0.5*(f[26]+f[19]-1.0*(f[18]+f[11])); 
  fReflXYQuad[2][5] = -0.5*(f[27]+f[22]-1.0*(f[21]+f[14])); 
  fReflXYQuad[2][6] = -0.5*(f[28]+f[24]-1.0*(f[23]+f[15])); 
  fReflXYQuad[2][7] = -0.5*(f[31]+f[30]-1.0*(f[29]+f[25])); 
  fReflXYQuad[2][8] = -0.03333333333333333*(15.0*f[37]+15.0*f[34]-1.0*(15.0*f[33]+15.0*f[32])); 
  fReflXYQuad[2][9] = -0.03333333333333333*(15.0*f[43]+15.0*f[39]-1.0*(15.0*f[38]+15.0*f[35])); 
  fReflXYQuad[2][10] = -0.03333333333333333*(15.0*f[44]+15.0*f[41]-1.0*(15.0*f[40]+15.0*f[36])); 
  fReflXYQuad[2][11] = -0.03333333333333333*(15.0*f[47]+15.0*f[46]-1.0*(15.0*f[45]+15.0*f[42])); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))/(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]+f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34])-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]+f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34])-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))/(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]-1.0*f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*(f[36]-1.0*(f[35]+f[34]-1.0*f[33]))+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]-1.0*f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*(f[36]-1.0*(f[35]+f[34]-1.0*f[33]))+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))/(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]-1.0*f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]-1.0*f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]-1.0*f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]-1.0*f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))/(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[34]-1.0*(f[36]+f[35]))-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[34]-1.0*(f[36]+f[35]))-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_4 
  vcutSq_i = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*phi[6]+2.449489742783178*phiWall[5]-2.449489742783178*phi[5]-1.414213562373095*phiWall[4]+1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*phi[2]-1.414213562373095*phiWall[1]+1.414213562373095*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.5*(f[6]+f[2]+f[1]+f[0]); 
  fReflXYQuad[3][1] = 0.5*(f[16]+f[8]+f[7]+f[3]); 
  fReflXYQuad[3][2] = 0.5*(f[17]+f[10]+f[9]+f[4]); 
  fReflXYQuad[3][3] = 0.5*(f[20]+f[13]+f[12]+f[5]); 
  fReflXYQuad[3][4] = 0.5*(f[26]+f[19]+f[18]+f[11]); 
  fReflXYQuad[3][5] = 0.5*(f[27]+f[22]+f[21]+f[14]); 
  fReflXYQuad[3][6] = 0.5*(f[28]+f[24]+f[23]+f[15]); 
  fReflXYQuad[3][7] = 0.5*(f[31]+f[30]+f[29]+f[25]); 
  fReflXYQuad[3][8] = 0.03333333333333333*(15.0*f[37]+15.0*(f[34]+f[33])+15.0*f[32]); 
  fReflXYQuad[3][9] = 0.03333333333333333*(15.0*f[43]+15.0*(f[39]+f[38])+15.0*f[35]); 
  fReflXYQuad[3][10] = 0.03333333333333333*(15.0*f[44]+15.0*(f[41]+f[40])+15.0*f[36]); 
  fReflXYQuad[3][11] = 0.03333333333333333*(15.0*f[47]+15.0*(f[46]+f[45])+15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))/(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*(f[44]+f[43]))+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*((-1.0*(f[36]+f[35]))+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*(f[44]+f[43]))+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*((-1.0*(f[36]+f[35]))+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))/(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[35]-1.0*f[36])-1.0*(15.0*(f[34]+f[33])+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[35]-1.0*f[36])-1.0*(15.0*(f[34]+f[33])+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))/(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*f[36]-1.0*(15.0*(f[35]+f[34]+f[33])+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*f[36]-1.0*(15.0*(f[35]+f[34]+f[33])+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))/(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
  fRefl[0] = 0.5*(fReflXYQuad[3][0]+fReflXYQuad[2][0]+fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[1] = 0.5*(fReflXYQuad[3][0]+fReflXYQuad[2][0]-1.0*(fReflXYQuad[1][0]+fReflXYQuad[0][0])); 
  fRefl[2] = 0.5*(fReflXYQuad[3][0]-1.0*fReflXYQuad[2][0]+fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[3] = 0.5*(fReflXYQuad[3][1]+fReflXYQuad[2][1]+fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[4] = 0.5*(fReflXYQuad[3][2]+fReflXYQuad[2][2]+fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[5] = 0.5*(fReflXYQuad[3][3]+fReflXYQuad[2][3]+fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[6] = 0.5*(fReflXYQuad[3][0]-1.0*(fReflXYQuad[2][0]+fReflXYQuad[1][0])+fReflXYQuad[0][0]); 
  fRefl[7] = 0.5*(fReflXYQuad[3][1]+fReflXYQuad[2][1]-1.0*(fReflXYQuad[1][1]+fReflXYQuad[0][1])); 
  fRefl[8] = 0.5*(fReflXYQuad[3][1]-1.0*fReflXYQuad[2][1]+fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[9] = 0.5*(fReflXYQuad[3][2]+fReflXYQuad[2][2]-1.0*(fReflXYQuad[1][2]+fReflXYQuad[0][2])); 
  fRefl[10] = 0.5*(fReflXYQuad[3][2]-1.0*fReflXYQuad[2][2]+fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[11] = 0.5*(fReflXYQuad[3][4]+fReflXYQuad[2][4]+fReflXYQuad[1][4]+fReflXYQuad[0][4]); 
  fRefl[12] = 0.5*(fReflXYQuad[3][3]+fReflXYQuad[2][3]-1.0*(fReflXYQuad[1][3]+fReflXYQuad[0][3])); 
  fRefl[13] = 0.5*(fReflXYQuad[3][3]-1.0*fReflXYQuad[2][3]+fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[14] = 0.5*(fReflXYQuad[3][5]+fReflXYQuad[2][5]+fReflXYQuad[1][5]+fReflXYQuad[0][5]); 
  fRefl[15] = 0.5*(fReflXYQuad[3][6]+fReflXYQuad[2][6]+fReflXYQuad[1][6]+fReflXYQuad[0][6]); 
  fRefl[16] = 0.5*(fReflXYQuad[3][1]-1.0*(fReflXYQuad[2][1]+fReflXYQuad[1][1])+fReflXYQuad[0][1]); 
  fRefl[17] = 0.5*(fReflXYQuad[3][2]-1.0*(fReflXYQuad[2][2]+fReflXYQuad[1][2])+fReflXYQuad[0][2]); 
  fRefl[18] = 0.5*(fReflXYQuad[3][4]+fReflXYQuad[2][4]-1.0*(fReflXYQuad[1][4]+fReflXYQuad[0][4])); 
  fRefl[19] = 0.5*(fReflXYQuad[3][4]-1.0*fReflXYQuad[2][4]+fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[20] = 0.5*(fReflXYQuad[3][3]-1.0*(fReflXYQuad[2][3]+fReflXYQuad[1][3])+fReflXYQuad[0][3]); 
  fRefl[21] = 0.5*(fReflXYQuad[3][5]+fReflXYQuad[2][5]-1.0*(fReflXYQuad[1][5]+fReflXYQuad[0][5])); 
  fRefl[22] = 0.5*(fReflXYQuad[3][5]-1.0*fReflXYQuad[2][5]+fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[23] = 0.5*(fReflXYQuad[3][6]+fReflXYQuad[2][6]-1.0*(fReflXYQuad[1][6]+fReflXYQuad[0][6])); 
  fRefl[24] = 0.5*(fReflXYQuad[3][6]-1.0*fReflXYQuad[2][6]+fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[25] = 0.5*(fReflXYQuad[3][7]+fReflXYQuad[2][7]+fReflXYQuad[1][7]+fReflXYQuad[0][7]); 
  fRefl[26] = 0.5*(fReflXYQuad[3][4]-1.0*(fReflXYQuad[2][4]+fReflXYQuad[1][4])+fReflXYQuad[0][4]); 
  fRefl[27] = 0.5*(fReflXYQuad[3][5]-1.0*(fReflXYQuad[2][5]+fReflXYQuad[1][5])+fReflXYQuad[0][5]); 
  fRefl[28] = 0.5*(fReflXYQuad[3][6]-1.0*(fReflXYQuad[2][6]+fReflXYQuad[1][6])+fReflXYQuad[0][6]); 
  fRefl[29] = 0.5*(fReflXYQuad[3][7]+fReflXYQuad[2][7]-1.0*(fReflXYQuad[1][7]+fReflXYQuad[0][7])); 
  fRefl[30] = 0.5*(fReflXYQuad[3][7]-1.0*fReflXYQuad[2][7]+fReflXYQuad[1][7]-1.0*fReflXYQuad[0][7]); 
  fRefl[31] = 0.5*(fReflXYQuad[3][7]-1.0*(fReflXYQuad[2][7]+fReflXYQuad[1][7])+fReflXYQuad[0][7]); 
  fRefl[32] = 0.5*(fReflXYQuad[3][8]+fReflXYQuad[2][8]+fReflXYQuad[1][8]+fReflXYQuad[0][8]); 
  fRefl[33] = 0.5000000000000001*(fReflXYQuad[3][8]+fReflXYQuad[2][8]-1.0*(fReflXYQuad[1][8]+fReflXYQuad[0][8])); 
  fRefl[34] = 0.5000000000000001*(fReflXYQuad[3][8]-1.0*fReflXYQuad[2][8]+fReflXYQuad[1][8]-1.0*fReflXYQuad[0][8]); 
  fRefl[35] = 0.5*(fReflXYQuad[3][9]+fReflXYQuad[2][9]+fReflXYQuad[1][9]+fReflXYQuad[0][9]); 
  fRefl[36] = 0.5*(fReflXYQuad[3][10]+fReflXYQuad[2][10]+fReflXYQuad[1][10]+fReflXYQuad[0][10]); 
  fRefl[37] = 0.5*(fReflXYQuad[3][8]-1.0*(fReflXYQuad[2][8]+fReflXYQuad[1][8])+fReflXYQuad[0][8]); 
  fRefl[38] = 0.5000000000000001*(fReflXYQuad[3][9]+fReflXYQuad[2][9]-1.0*(fReflXYQuad[1][9]+fReflXYQuad[0][9])); 
  fRefl[39] = 0.5000000000000001*(fReflXYQuad[3][9]-1.0*fReflXYQuad[2][9]+fReflXYQuad[1][9]-1.0*fReflXYQuad[0][9]); 
  fRefl[40] = 0.5000000000000001*(fReflXYQuad[3][10]+fReflXYQuad[2][10]-1.0*(fReflXYQuad[1][10]+fReflXYQuad[0][10])); 
  fRefl[41] = 0.5000000000000001*(fReflXYQuad[3][10]-1.0*fReflXYQuad[2][10]+fReflXYQuad[1][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[42] = 0.5*(fReflXYQuad[3][11]+fReflXYQuad[2][11]+fReflXYQuad[1][11]+fReflXYQuad[0][11]); 
  fRefl[43] = 0.5*(fReflXYQuad[3][9]-1.0*(fReflXYQuad[2][9]+fReflXYQuad[1][9])+fReflXYQuad[0][9]); 
  fRefl[44] = 0.5*(fReflXYQuad[3][10]-1.0*(fReflXYQuad[2][10]+fReflXYQuad[1][10])+fReflXYQuad[0][10]); 
  fRefl[45] = 0.5000000000000001*(fReflXYQuad[3][11]+fReflXYQuad[2][11]-1.0*(fReflXYQuad[1][11]+fReflXYQuad[0][11])); 
  fRefl[46] = 0.5000000000000001*(fReflXYQuad[3][11]-1.0*fReflXYQuad[2][11]+fReflXYQuad[1][11]-1.0*fReflXYQuad[0][11]); 
  fRefl[47] = 0.5*(fReflXYQuad[3][11]-1.0*(fReflXYQuad[2][11]+fReflXYQuad[1][11])+fReflXYQuad[0][11]); 
}

void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; long double xc, b, xbarVal, fac; 
  double fReflXYQuad[4][12]; 
  double fReflXYZMuQuad[4][6]; 
  

// node (x,y)_1 
  vcutSq_i = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*phi[6]-2.449489742783178*phiWall[5]+2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*phi[2]-1.414213562373095*phiWall[1]+1.414213562373095*(phi[1]+phiWall[0])-1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = 0.5*(f[6]-1.0*(f[2]+f[1])+f[0]); 
  fReflXYQuad[0][1] = 0.5*(f[16]-1.0*(f[8]+f[7])+f[3]); 
  fReflXYQuad[0][2] = 0.5*(f[17]-1.0*(f[10]+f[9])+f[4]); 
  fReflXYQuad[0][3] = 0.5*(f[20]-1.0*(f[13]+f[12])+f[5]); 
  fReflXYQuad[0][4] = 0.5*(f[26]-1.0*(f[19]+f[18])+f[11]); 
  fReflXYQuad[0][5] = 0.5*(f[27]-1.0*(f[22]+f[21])+f[14]); 
  fReflXYQuad[0][6] = 0.5*(f[28]-1.0*(f[24]+f[23])+f[15]); 
  fReflXYQuad[0][7] = 0.5*(f[31]-1.0*(f[30]+f[29])+f[25]); 
  fReflXYQuad[0][8] = 0.03333333333333333*(15.0*f[37]-15.0*(f[34]+f[33])+15.0*f[32]); 
  fReflXYQuad[0][9] = 0.03333333333333333*(15.0*f[43]-15.0*(f[39]+f[38])+15.0*f[35]); 
  fReflXYQuad[0][10] = 0.03333333333333333*(15.0*f[44]-15.0*(f[41]+f[40])+15.0*f[36]); 
  fReflXYQuad[0][11] = 0.03333333333333333*(15.0*f[47]-15.0*(f[46]+f[45])+15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))/(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]-1.0*(f[22]+f[21]+f[20]+f[16])+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]-1.0*(f[30]+f[29]+f[28]+f[26])+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10]+f[9])+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]-1.0*(f[22]+f[21]+f[20])+f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]+f[5])+f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]-1.0*(f[30]+f[29]+f[28])+f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17]+f[15])+f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]-15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]-1.0*f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[36]-1.0*f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]-1.0*(f[22]+f[21])+f[20]-1.0*f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]-1.0*(f[30]+f[29])+f[28]-1.0*f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]-1.0*f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[36]-1.0*f[35]+f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))/(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*(f[36]+f[35]-1.0*(f[34]+f[33]))+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]-1.0*(f[22]+f[21])+f[20]+f[16]+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]+f[5]+f[3]-1.0*(f[2]+f[1])+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]-1.0*(f[30]+f[29])+f[28]+f[26]+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]+f[15]+f[11]-1.0*(f[10]+f[9])+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*(f[46]+f[45]))+f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*(f[36]+f[35]-1.0*(f[34]+f[33]))+15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[0][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_2 
  vcutSq_i = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*(phi[6]+phiWall[5])-2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]-1.414213562373095*phiWall[2]+1.414213562373095*(phi[2]+phiWall[1])-1.414213562373095*(phi[1]+phiWall[0])+1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = -0.5*(f[6]-1.0*f[2]+f[1]-1.0*f[0]); 
  fReflXYQuad[1][1] = -0.5*(f[16]-1.0*f[8]+f[7]-1.0*f[3]); 
  fReflXYQuad[1][2] = -0.5*(f[17]-1.0*f[10]+f[9]-1.0*f[4]); 
  fReflXYQuad[1][3] = -0.5*(f[20]-1.0*f[13]+f[12]-1.0*f[5]); 
  fReflXYQuad[1][4] = -0.5*(f[26]-1.0*f[19]+f[18]-1.0*f[11]); 
  fReflXYQuad[1][5] = -0.5*(f[27]-1.0*f[22]+f[21]-1.0*f[14]); 
  fReflXYQuad[1][6] = -0.5*(f[28]-1.0*f[24]+f[23]-1.0*f[15]); 
  fReflXYQuad[1][7] = -0.5*(f[31]-1.0*f[30]+f[29]-1.0*f[25]); 
  fReflXYQuad[1][8] = -0.03333333333333333*(15.0*f[37]+15.0*(f[33]-1.0*f[34])-15.0*f[32]); 
  fReflXYQuad[1][9] = -0.03333333333333333*(15.0*f[43]+15.0*(f[38]-1.0*f[39])-15.0*f[35]); 
  fReflXYQuad[1][10] = -0.03333333333333333*(15.0*f[44]+15.0*(f[40]-1.0*f[41])-15.0*f[36]); 
  fReflXYQuad[1][11] = -0.03333333333333333*(15.0*f[47]+15.0*(f[45]-1.0*f[46])-15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))/(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*(f[44]+f[43]))+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[36]+f[35]-1.0*f[34]+f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]-1.0*f[22]+f[21]-1.0*(f[20]+f[16]+f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]+f[5]+f[3]-1.0*f[2]+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]-1.0*f[30]+f[29]-1.0*(f[28]+f[26]+f[25])+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]+f[15]+f[11]-1.0*f[10]+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*(f[44]+f[43]))+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[36]+f[35]-1.0*f[34]+f[33])-15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))/(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*(f[36]-1.0*f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]-1.0*f[22]+f[21]-1.0*f[20]+f[16]-1.0*f[14]+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*f[6]+f[5]-1.0*f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]-1.0*f[30]+f[29]-1.0*f[28]+f[26]-1.0*f[25]+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*f[17]+f[15]-1.0*f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]-1.0*f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*(f[36]-1.0*f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))/(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]-1.0*f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]-1.0*f[22]+f[21]+f[20]-1.0*(f[16]+f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6]+f[5])+f[3]+f[2]-1.0*f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]-1.0*f[30]+f[29]+f[28]-1.0*(f[26]+f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17]+f[15])+f[11]+f[10]-1.0*f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]-1.0*f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]+f[34]-1.0*f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))/(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]+f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]-1.0*f[33])-15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]-1.0*f[22]+f[21]+f[20]+f[16]-1.0*(f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]-1.0*(f[5]+f[3]+f[2])+f[1]-1.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]-1.0*f[30]+f[29]+f[28]+f[26]-1.0*(f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]-1.0*(f[15]+f[11]+f[10])+f[9]-1.0*f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*((-1.0*f[46])+f[45]+f[44]+f[43])+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])-15.0*(f[36]+f[35]+f[34]-1.0*f[33])-15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[1][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_3 
  vcutSq_i = -0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*(phi[6]+phiWall[5])+2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*(phi[2]+phiWall[1])+1.414213562373095*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.5*(f[6]+f[2]-1.0*(f[1]+f[0])); 
  fReflXYQuad[2][1] = -0.5*(f[16]+f[8]-1.0*(f[7]+f[3])); 
  fReflXYQuad[2][2] = -0.5*(f[17]+f[10]-1.0*(f[9]+f[4])); 
  fReflXYQuad[2][3] = -0.5*(f[20]+f[13]-1.0*(f[12]+f[5])); 
  fReflXYQuad[2][4] = -0.5*(f[26]+f[19]-1.0*(f[18]+f[11])); 
  fReflXYQuad[2][5] = -0.5*(f[27]+f[22]-1.0*(f[21]+f[14])); 
  fReflXYQuad[2][6] = -0.5*(f[28]+f[24]-1.0*(f[23]+f[15])); 
  fReflXYQuad[2][7] = -0.5*(f[31]+f[30]-1.0*(f[29]+f[25])); 
  fReflXYQuad[2][8] = -0.03333333333333333*(15.0*f[37]+15.0*f[34]-1.0*(15.0*f[33]+15.0*f[32])); 
  fReflXYQuad[2][9] = -0.03333333333333333*(15.0*f[43]+15.0*f[39]-1.0*(15.0*f[38]+15.0*f[35])); 
  fReflXYQuad[2][10] = -0.03333333333333333*(15.0*f[44]+15.0*f[41]-1.0*(15.0*f[40]+15.0*f[36])); 
  fReflXYQuad[2][11] = -0.03333333333333333*(15.0*f[47]+15.0*f[46]-1.0*(15.0*f[45]+15.0*f[42])); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))/(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]+f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34])-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.25*(f[27]+f[22]-1.0*(f[21]+f[20]+f[16]+f[14]+f[13])+f[12]-1.0*f[8]+f[7]+f[6]+f[5]+f[3]+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.25*(f[31]+f[30]-1.0*(f[29]+f[28]+f[26]+f[25]+f[24])+f[23]-1.0*f[19]+f[18]+f[17]+f[15]+f[11]+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]+f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]-1.0*f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34])-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))/(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]-1.0*f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*(f[36]-1.0*(f[35]+f[34]-1.0*f[33]))+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.25*(f[27]+f[22]-1.0*(f[21]+f[20])+f[16]-1.0*(f[14]+f[13])+f[12]+f[8]-1.0*(f[7]+f[6])+f[5]-1.0*(f[3]+f[2])+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.25*(f[31]+f[30]-1.0*(f[29]+f[28])+f[26]-1.0*(f[25]+f[24])+f[23]+f[19]-1.0*(f[18]+f[17])+f[15]-1.0*(f[11]+f[10])+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*(f[45]+f[44]-1.0*f[43]))+15.0*((-1.0*(f[42]+f[41]))+f[40]+f[39]-1.0*(f[38]+f[37]))+15.0*(f[36]-1.0*(f[35]+f[34]-1.0*f[33]))+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))/(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]-1.0*f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]-1.0*f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.25*(f[27]+f[22]-1.0*f[21]+f[20]-1.0*(f[16]+f[14])+f[13]-1.0*(f[12]+f[8])+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.25*(f[31]+f[30]-1.0*f[29]+f[28]-1.0*(f[26]+f[25])+f[24]-1.0*(f[23]+f[19])+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]-1.0*f[43])+15.0*((-1.0*f[42])+f[41]-1.0*(f[40]+f[39]-1.0*f[38]+f[37]))+15.0*((-1.0*f[36])+f[35]-1.0*f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))/(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[34]-1.0*(f[36]+f[35]))-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.25*(f[27]+f[22]-1.0*f[21]+f[20]+f[16]-1.0*f[14]+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6]-1.0*(f[5]+f[3])+f[2]-1.0*(f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.25*(f[31]+f[30]-1.0*f[29]+f[28]+f[26]-1.0*f[25]+f[24]-1.0*f[23]+f[19]-1.0*f[18]+f[17]-1.0*(f[15]+f[11])+f[10]-1.0*(f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]-1.0*f[45]+f[44]+f[43])+15.0*((-1.0*f[42])+f[41]-1.0*f[40]+f[39]-1.0*f[38]+f[37])+15.0*(f[34]-1.0*(f[36]+f[35]))-1.0*(15.0*f[33]+15.0*f[32])))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[2][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
// node (x,y)_4 
  vcutSq_i = 0.25*(2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*phi[6]+2.449489742783178*phiWall[5]-2.449489742783178*phi[5]+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3]+1.414213562373095*phiWall[2]-1.414213562373095*phi[2]+1.414213562373095*phiWall[1]-1.414213562373095*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
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
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.5*(f[6]+f[2]+f[1]+f[0]); 
  fReflXYQuad[3][1] = 0.5*(f[16]+f[8]+f[7]+f[3]); 
  fReflXYQuad[3][2] = 0.5*(f[17]+f[10]+f[9]+f[4]); 
  fReflXYQuad[3][3] = 0.5*(f[20]+f[13]+f[12]+f[5]); 
  fReflXYQuad[3][4] = 0.5*(f[26]+f[19]+f[18]+f[11]); 
  fReflXYQuad[3][5] = 0.5*(f[27]+f[22]+f[21]+f[14]); 
  fReflXYQuad[3][6] = 0.5*(f[28]+f[24]+f[23]+f[15]); 
  fReflXYQuad[3][7] = 0.5*(f[31]+f[30]+f[29]+f[25]); 
  fReflXYQuad[3][8] = 0.03333333333333333*(15.0*f[37]+15.0*(f[34]+f[33])+15.0*f[32]); 
  fReflXYQuad[3][9] = 0.03333333333333333*(15.0*f[43]+15.0*(f[39]+f[38])+15.0*f[35]); 
  fReflXYQuad[3][10] = 0.03333333333333333*(15.0*f[44]+15.0*(f[41]+f[40])+15.0*f[36]); 
  fReflXYQuad[3][11] = 0.03333333333333333*(15.0*f[47]+15.0*(f[46]+f[45])+15.0*f[42]); 
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))/(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*(f[44]+f[43]))+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*((-1.0*(f[36]+f[35]))+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.25*(f[27]+f[22]+f[21]-1.0*(f[20]+f[16])+f[14]-1.0*(f[13]+f[12]+f[8]+f[7])+f[6]-1.0*(f[5]+f[3])+f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.25*(f[31]+f[30]+f[29]-1.0*(f[28]+f[26])+f[25]-1.0*(f[24]+f[23]+f[19]+f[18])+f[17]-1.0*(f[15]+f[11])+f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*(f[44]+f[43]))+15.0*(f[42]-1.0*(f[41]+f[40]+f[39]+f[38]-1.0*f[37]))+15.0*((-1.0*(f[36]+f[35]))+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))/(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[35]-1.0*f[36])-1.0*(15.0*(f[34]+f[33])+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.25*(f[27]+f[22]+f[21]-1.0*f[20]+f[16]+f[14]-1.0*(f[13]+f[12])+f[8]+f[7]-1.0*(f[6]+f[5])+f[3]-1.0*(f[2]+f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.25*(f[31]+f[30]+f[29]-1.0*f[28]+f[26]+f[25]-1.0*(f[24]+f[23])+f[19]+f[18]-1.0*(f[17]+f[15])+f[11]-1.0*(f[10]+f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]-1.0*f[44]+f[43])+15.0*(f[42]-1.0*(f[41]+f[40]-1.0*f[39])+f[38]-1.0*f[37])+15.0*(f[35]-1.0*f[36])-1.0*(15.0*(f[34]+f[33])+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))/(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*f[36]-1.0*(15.0*(f[35]+f[34]+f[33])+15.0*f[32])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.25*(f[27]+f[22]+f[21]+f[20]-1.0*f[16]+f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6])+f[5]-1.0*(f[3]+f[2]+f[1]+f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.25*(f[31]+f[30]+f[29]+f[28]-1.0*f[26]+f[25]+f[24]+f[23]-1.0*(f[19]+f[18]+f[17])+f[15]-1.0*(f[11]+f[10]+f[9]+f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]-1.0*f[43])+15.0*(f[42]+f[41]+f[40]-1.0*(f[39]+f[38]+f[37]))+15.0*f[36]-1.0*(15.0*(f[35]+f[34]+f[33])+15.0*f[32])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))/(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.25*(f[27]+f[22]+f[21]+f[20]+f[16]+f[14]+f[13]+f[12]+f[8]+f[7]+f[6]+f[5]+f[3]+f[2]+f[1]+f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.25*(f[31]+f[30]+f[29]+f[28]+f[26]+f[25]+f[24]+f[23]+f[19]+f[18]+f[17]+f[15]+f[11]+f[10]+f[9]+f[4]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01666666666666667*(15.0*f[47]+15.0*(f[46]+f[45]+f[44]+f[43])+15.0*(f[42]+f[41]+f[40]+f[39]+f[38]+f[37])+15.0*(f[36]+f[35]+f[34]+f[33])+15.0*f[32]))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][1] = 0.5*(fReflXYZMuQuad[3][0]+fReflXYZMuQuad[2][0]-1.0*(fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][2] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][3] = 0.5*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][4] = 0.5*(fReflXYZMuQuad[3][1]+fReflXYZMuQuad[2][1]-1.0*(fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][5] = 0.5*(fReflXYZMuQuad[3][0]-1.0*(fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][6] = 0.5*(fReflXYZMuQuad[3][1]-1.0*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1]-1.0*fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][7] = 0.5*(fReflXYZMuQuad[3][1]-1.0*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][8] = 0.5*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][9] = 0.5000000000000001*(fReflXYZMuQuad[3][2]+fReflXYZMuQuad[2][2]-1.0*(fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2])); 
  fReflXYQuad[3][10] = 0.5000000000000001*(fReflXYZMuQuad[3][2]-1.0*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2]-1.0*fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][11] = 0.5*(fReflXYZMuQuad[3][2]-1.0*(fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])+fReflXYZMuQuad[0][2]); 
  } 

 
  fRefl[0] = 0.5*(fReflXYQuad[3][0]+fReflXYQuad[2][0]+fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[1] = 0.5*(fReflXYQuad[3][0]+fReflXYQuad[2][0]-1.0*(fReflXYQuad[1][0]+fReflXYQuad[0][0])); 
  fRefl[2] = 0.5*(fReflXYQuad[3][0]-1.0*fReflXYQuad[2][0]+fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[3] = 0.5*(fReflXYQuad[3][1]+fReflXYQuad[2][1]+fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[4] = 0.5*(fReflXYQuad[3][2]+fReflXYQuad[2][2]+fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[5] = 0.5*(fReflXYQuad[3][3]+fReflXYQuad[2][3]+fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[6] = 0.5*(fReflXYQuad[3][0]-1.0*(fReflXYQuad[2][0]+fReflXYQuad[1][0])+fReflXYQuad[0][0]); 
  fRefl[7] = 0.5*(fReflXYQuad[3][1]+fReflXYQuad[2][1]-1.0*(fReflXYQuad[1][1]+fReflXYQuad[0][1])); 
  fRefl[8] = 0.5*(fReflXYQuad[3][1]-1.0*fReflXYQuad[2][1]+fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[9] = 0.5*(fReflXYQuad[3][2]+fReflXYQuad[2][2]-1.0*(fReflXYQuad[1][2]+fReflXYQuad[0][2])); 
  fRefl[10] = 0.5*(fReflXYQuad[3][2]-1.0*fReflXYQuad[2][2]+fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[11] = 0.5*(fReflXYQuad[3][4]+fReflXYQuad[2][4]+fReflXYQuad[1][4]+fReflXYQuad[0][4]); 
  fRefl[12] = 0.5*(fReflXYQuad[3][3]+fReflXYQuad[2][3]-1.0*(fReflXYQuad[1][3]+fReflXYQuad[0][3])); 
  fRefl[13] = 0.5*(fReflXYQuad[3][3]-1.0*fReflXYQuad[2][3]+fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[14] = 0.5*(fReflXYQuad[3][5]+fReflXYQuad[2][5]+fReflXYQuad[1][5]+fReflXYQuad[0][5]); 
  fRefl[15] = 0.5*(fReflXYQuad[3][6]+fReflXYQuad[2][6]+fReflXYQuad[1][6]+fReflXYQuad[0][6]); 
  fRefl[16] = 0.5*(fReflXYQuad[3][1]-1.0*(fReflXYQuad[2][1]+fReflXYQuad[1][1])+fReflXYQuad[0][1]); 
  fRefl[17] = 0.5*(fReflXYQuad[3][2]-1.0*(fReflXYQuad[2][2]+fReflXYQuad[1][2])+fReflXYQuad[0][2]); 
  fRefl[18] = 0.5*(fReflXYQuad[3][4]+fReflXYQuad[2][4]-1.0*(fReflXYQuad[1][4]+fReflXYQuad[0][4])); 
  fRefl[19] = 0.5*(fReflXYQuad[3][4]-1.0*fReflXYQuad[2][4]+fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[20] = 0.5*(fReflXYQuad[3][3]-1.0*(fReflXYQuad[2][3]+fReflXYQuad[1][3])+fReflXYQuad[0][3]); 
  fRefl[21] = 0.5*(fReflXYQuad[3][5]+fReflXYQuad[2][5]-1.0*(fReflXYQuad[1][5]+fReflXYQuad[0][5])); 
  fRefl[22] = 0.5*(fReflXYQuad[3][5]-1.0*fReflXYQuad[2][5]+fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[23] = 0.5*(fReflXYQuad[3][6]+fReflXYQuad[2][6]-1.0*(fReflXYQuad[1][6]+fReflXYQuad[0][6])); 
  fRefl[24] = 0.5*(fReflXYQuad[3][6]-1.0*fReflXYQuad[2][6]+fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[25] = 0.5*(fReflXYQuad[3][7]+fReflXYQuad[2][7]+fReflXYQuad[1][7]+fReflXYQuad[0][7]); 
  fRefl[26] = 0.5*(fReflXYQuad[3][4]-1.0*(fReflXYQuad[2][4]+fReflXYQuad[1][4])+fReflXYQuad[0][4]); 
  fRefl[27] = 0.5*(fReflXYQuad[3][5]-1.0*(fReflXYQuad[2][5]+fReflXYQuad[1][5])+fReflXYQuad[0][5]); 
  fRefl[28] = 0.5*(fReflXYQuad[3][6]-1.0*(fReflXYQuad[2][6]+fReflXYQuad[1][6])+fReflXYQuad[0][6]); 
  fRefl[29] = 0.5*(fReflXYQuad[3][7]+fReflXYQuad[2][7]-1.0*(fReflXYQuad[1][7]+fReflXYQuad[0][7])); 
  fRefl[30] = 0.5*(fReflXYQuad[3][7]-1.0*fReflXYQuad[2][7]+fReflXYQuad[1][7]-1.0*fReflXYQuad[0][7]); 
  fRefl[31] = 0.5*(fReflXYQuad[3][7]-1.0*(fReflXYQuad[2][7]+fReflXYQuad[1][7])+fReflXYQuad[0][7]); 
  fRefl[32] = 0.5*(fReflXYQuad[3][8]+fReflXYQuad[2][8]+fReflXYQuad[1][8]+fReflXYQuad[0][8]); 
  fRefl[33] = 0.5000000000000001*(fReflXYQuad[3][8]+fReflXYQuad[2][8]-1.0*(fReflXYQuad[1][8]+fReflXYQuad[0][8])); 
  fRefl[34] = 0.5000000000000001*(fReflXYQuad[3][8]-1.0*fReflXYQuad[2][8]+fReflXYQuad[1][8]-1.0*fReflXYQuad[0][8]); 
  fRefl[35] = 0.5*(fReflXYQuad[3][9]+fReflXYQuad[2][9]+fReflXYQuad[1][9]+fReflXYQuad[0][9]); 
  fRefl[36] = 0.5*(fReflXYQuad[3][10]+fReflXYQuad[2][10]+fReflXYQuad[1][10]+fReflXYQuad[0][10]); 
  fRefl[37] = 0.5*(fReflXYQuad[3][8]-1.0*(fReflXYQuad[2][8]+fReflXYQuad[1][8])+fReflXYQuad[0][8]); 
  fRefl[38] = 0.5000000000000001*(fReflXYQuad[3][9]+fReflXYQuad[2][9]-1.0*(fReflXYQuad[1][9]+fReflXYQuad[0][9])); 
  fRefl[39] = 0.5000000000000001*(fReflXYQuad[3][9]-1.0*fReflXYQuad[2][9]+fReflXYQuad[1][9]-1.0*fReflXYQuad[0][9]); 
  fRefl[40] = 0.5000000000000001*(fReflXYQuad[3][10]+fReflXYQuad[2][10]-1.0*(fReflXYQuad[1][10]+fReflXYQuad[0][10])); 
  fRefl[41] = 0.5000000000000001*(fReflXYQuad[3][10]-1.0*fReflXYQuad[2][10]+fReflXYQuad[1][10]-1.0*fReflXYQuad[0][10]); 
  fRefl[42] = 0.5*(fReflXYQuad[3][11]+fReflXYQuad[2][11]+fReflXYQuad[1][11]+fReflXYQuad[0][11]); 
  fRefl[43] = 0.5*(fReflXYQuad[3][9]-1.0*(fReflXYQuad[2][9]+fReflXYQuad[1][9])+fReflXYQuad[0][9]); 
  fRefl[44] = 0.5*(fReflXYQuad[3][10]-1.0*(fReflXYQuad[2][10]+fReflXYQuad[1][10])+fReflXYQuad[0][10]); 
  fRefl[45] = 0.5000000000000001*(fReflXYQuad[3][11]+fReflXYQuad[2][11]-1.0*(fReflXYQuad[1][11]+fReflXYQuad[0][11])); 
  fRefl[46] = 0.5000000000000001*(fReflXYQuad[3][11]-1.0*fReflXYQuad[2][11]+fReflXYQuad[1][11]-1.0*fReflXYQuad[0][11]); 
  fRefl[47] = 0.5*(fReflXYQuad[3][11]-1.0*(fReflXYQuad[2][11]+fReflXYQuad[1][11])+fReflXYQuad[0][11]); 
}
