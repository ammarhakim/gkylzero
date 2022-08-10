#include <gkyl_bc_sheath_gyrokinetic.h> 


void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_tensor_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZQuad[3][9]; 
  

  vcutSq = 0.5*(2.23606797749979*(1.414213562373095*phiWall[2]-1.414213562373095*phi[2])+1.732050807568877*(1.414213562373095*phi[1]-1.414213562373095*phiWall[1])+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
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
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))/(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (0.1414213562373095*(4.47213595499958*f[8]-6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (0.1414213562373095*(4.47213595499958*f[8]-6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*f[6]-6.0*f[2]))/(2.23606797749979*f[4]-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (-0.3535533905932737*(2.23606797749979*f[8]-2.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (-0.3535533905932737*(2.23606797749979*f[8]-2.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))/(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[2][0] = 0.0; 
  fReflZQuad[2][1] = 0.0; 
  fReflZQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(4.47213595499958*f[8]+6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(4.47213595499958*f[8]+6.708203932499369*f[7]+5.0*f[5]))*fac; 
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
  fRefl[8] = 0.3513641844631533*(fReflZQuad[2][2]-2.0*fReflZQuad[1][2]+fReflZQuad[0][2]); 
  } 

 
}

void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_tensor_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZQuad[3][9]; 
  

  vcutSq = 0.5*(2.23606797749979*(1.414213562373095*phiWall[2]-1.414213562373095*phi[2])+1.732050807568877*(1.414213562373095*phiWall[1]-1.414213562373095*phi[1])+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if (vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
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
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))/(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (0.1414213562373095*(4.47213595499958*f[8]-6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (0.1414213562373095*(4.47213595499958*f[8]-6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*f[6]-6.0*f[2]))/(2.23606797749979*f[4]-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (-0.3535533905932737*(2.23606797749979*f[8]-2.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (-0.3535533905932737*(2.23606797749979*f[8]-2.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))/(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZQuad[2][0] = 0.0; 
  fReflZQuad[2][1] = 0.0; 
  fReflZQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if (wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(4.47213595499958*f[8]+6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(4.47213595499958*f[8]+6.708203932499369*f[7]+5.0*f[5]))*fac; 
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
  fRefl[8] = 0.3513641844631533*(fReflZQuad[2][2]-2.0*fReflZQuad[1][2]+fReflZQuad[0][2]); 
  } 

 
}

void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_tensor_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[9][9]; 
  

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
  fRefl[20] = 0.0; 
  fRefl[21] = 0.0; 
  fRefl[22] = 0.0; 
  fRefl[23] = 0.0; 
  fRefl[24] = 0.0; 
  fRefl[25] = 0.0; 
  fRefl[26] = 0.0; 
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
  fRefl[20] = f[20]; 
  fRefl[21] = f[21]; 
  fRefl[22] = f[22]; 
  fRefl[23] = f[23]; 
  fRefl[24] = f[24]; 
  fRefl[25] = f[25]; 
  fRefl[26] = f[26]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(2.0*(2.0*f[26]-3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(2.0*(2.0*f[26]-3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))/(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.01666666666666667*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[25])+2.23606797749979*(5.0*f[22]-4.0*f[20])+2.0*(6.708203932499369*f[12]-5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.01666666666666667*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[25])+2.23606797749979*(5.0*f[22]-4.0*f[20])+2.0*(6.708203932499369*f[12]-5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[23]-1.0*f[25])+2.23606797749979*(f[22]+f[20]))-9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[23]-1.0*f[25])+2.23606797749979*(f[22]+f[20]))-9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))/(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (-0.01666666666666667*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[23])+2.23606797749979*(5.0*f[20]-4.0*f[22])+2.0*(6.708203932499369*f[14]-5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (-0.01666666666666667*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[23])+2.23606797749979*(5.0*f[20]-4.0*f[22])+2.0*(6.708203932499369*f[14]-5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))/(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))); 
  // if f is not realizable, no reflection from this node 
  if (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (0.04166666666666666*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.125*(5.0*f[26]+2.0*(2.0*f[8]-2.23606797749979*(f[22]+f[20]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (0.04166666666666666*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.125*(5.0*f[26]+2.0*(2.0*f[8]-2.23606797749979*(f[22]+f[20]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))/(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[23])-1.0*(2.23606797749979*(4.0*f[22]-5.0*f[20])+2.0*(6.708203932499369*f[14]+5.0*f[8]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[23])-1.0*(2.23606797749979*(4.0*f[22]-5.0*f[20])+2.0*(6.708203932499369*f[14]+5.0*f[8]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]-1.0*f[23])+2.23606797749979*(f[22]+f[20]))-1.0*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12]))+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]-1.0*f[23])+2.23606797749979*(f[22]+f[20]))-1.0*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12]))+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))/(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[25])-1.0*(2.23606797749979*(4.0*f[20]-5.0*f[22])+2.0*(6.708203932499369*f[12]+5.0*f[8]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[25])-1.0*(2.23606797749979*(4.0*f[20]-5.0*f[22])+2.0*(6.708203932499369*f[12]+5.0*f[8]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[8][0] = 0.0; 
  fReflZMuQuad[8][1] = 0.0; 
  fReflZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
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
  fRefl[20] = 0.0276057774999974*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2]-2.0*(5.0*fReflZMuQuad[5][2]+8.0*fReflZMuQuad[4][2])+5.0*(fReflZMuQuad[2][2]-2.0*fReflZMuQuad[3][2])+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2]); 
  fRefl[21] = 0.1234567901234568*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0]+2.0*((-1.0*fReflZMuQuad[5][0])+2.0*fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])+fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[22] = 0.0276057774999974*(5.0*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2])+8.0*(fReflZMuQuad[5][2]-2.0*fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+5.0*(fReflZMuQuad[2][2]-2.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[23] = 0.1851851851851852*(fReflZMuQuad[8][2]-1.0*fReflZMuQuad[6][2]+2.0*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[5][2])+fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[24] = 0.1234567901234568*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1]+2.0*((-1.0*fReflZMuQuad[5][1])+2.0*fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])+fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[25] = 0.1851851851851852*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2]-1.0*fReflZMuQuad[2][2]+2.0*fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[26] = 0.1234567901234568*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2]+2.0*((-1.0*fReflZMuQuad[5][2])+2.0*fReflZMuQuad[4][2]-1.0*fReflZMuQuad[3][2])+fReflZMuQuad[2][2]-2.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  } 

 
}

void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_tensor_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[9][9]; 
  

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
  fRefl[20] = 0.0; 
  fRefl[21] = 0.0; 
  fRefl[22] = 0.0; 
  fRefl[23] = 0.0; 
  fRefl[24] = 0.0; 
  fRefl[25] = 0.0; 
  fRefl[26] = 0.0; 
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
  fRefl[20] = f[20]; 
  fRefl[21] = f[21]; 
  fRefl[22] = f[22]; 
  fRefl[23] = f[23]; 
  fRefl[24] = f[24]; 
  fRefl[25] = f[25]; 
  fRefl[26] = f[26]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(2.0*(2.0*f[26]-3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])-3.0*(f[3]+f[1]))-13.41640786499874*(f[15]+f[13]))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]-3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]-2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(2.0*(2.0*f[26]-3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))/(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.01666666666666667*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[25])+2.23606797749979*(5.0*f[22]-4.0*f[20])+2.0*(6.708203932499369*f[12]-5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.05*(10.0*f[21]+2.23606797749979*((-6.708203932499369*f[15])+5.0*f[9]+2.0*(3.0*f[1]-2.0*f[7]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.01666666666666667*(15.0*(2.0*f[24]-3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])+6.0*(6.708203932499369*f[4]-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[25])+2.23606797749979*(5.0*f[22]-4.0*f[20])+2.0*(6.708203932499369*f[12]-5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[23]-1.0*f[25])+2.23606797749979*(f[22]+f[20]))-9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(20.0*f[21]+2.23606797749979*(5.0*(2.0*(f[9]+f[7])+3.0*(f[3]-1.0*f[1]))-13.41640786499874*(f[15]-1.0*f[13]))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[17]-1.0*f[19]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]-1.0*f[4])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[23]-1.0*f[25])+2.23606797749979*(f[22]+f[20]))-9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))/(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (-0.01666666666666667*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[23])+2.23606797749979*(5.0*f[20]-4.0*f[22])+2.0*(6.708203932499369*f[14]-5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (-0.05*(10.0*f[21]-1.0*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-1.0*(5.0*f[7]+6.0*f[3]))+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (-0.01666666666666667*(30.0*f[24]-1.0*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11]))+6.0*(6.708203932499369*f[6]-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.05*(5.0*(2.0*f[26]-3.0*f[23])+2.23606797749979*(5.0*f[20]-4.0*f[22])+2.0*(6.708203932499369*f[14]-5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))/(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))); 
  // if f is not realizable, no reflection from this node 
  if (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (0.04166666666666666*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.125*(5.0*f[26]+2.0*(2.0*f[8]-2.23606797749979*(f[22]+f[20]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (0.125*(5.0*f[21]+2.0*(2.0*f[0]-2.23606797749979*(f[9]+f[7]))))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (0.04166666666666666*(15.0*f[24]+2.0*(6.0*f[2]-6.708203932499369*(f[16]+f[11]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.125*(5.0*f[26]+2.0*(2.0*f[8]-2.23606797749979*(f[22]+f[20]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))/(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[23])-1.0*(2.23606797749979*(4.0*f[22]-5.0*f[20])+2.0*(6.708203932499369*f[14]+5.0*f[8]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.05*(10.0*f[21]+15.0*f[13]-1.0*(2.23606797749979*(4.0*f[9]-5.0*f[7]+6.0*f[3])+10.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[17])-1.0*(6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(6.708203932499369*f[6]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[23])-1.0*(2.23606797749979*(4.0*f[22]-5.0*f[20])+2.0*(6.708203932499369*f[14]+5.0*f[8]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]-1.0*f[23])+2.23606797749979*(f[22]+f[20]))-1.0*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12]))+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]-1.0*f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[1]-1.0*f[3])))+5.0*(5.0*f[0]-9.0*f[5])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]-1.0*f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]-1.0*f[23])+2.23606797749979*(f[22]+f[20]))-1.0*(9.0*f[18]+6.708203932499369*(f[14]-1.0*f[12]))+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))/(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[25])-1.0*(2.23606797749979*(4.0*f[20]-5.0*f[22])+2.0*(6.708203932499369*f[12]+5.0*f[8]))))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (-0.05*(10.0*f[21]+2.23606797749979*(6.708203932499369*f[15]+5.0*f[9]-2.0*(2.0*f[7]+3.0*f[1]))-10.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (-0.01666666666666667*(15.0*(2.0*f[24]+3.0*f[19])+6.708203932499369*(5.0*f[16]-4.0*f[11])-6.0*(6.708203932499369*f[4]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (-0.05*(5.0*(2.0*f[26]+3.0*f[25])-1.0*(2.23606797749979*(4.0*f[20]-5.0*f[22])+2.0*(6.708203932499369*f[12]+5.0*f[8]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[8][0] = 0.0; 
  fReflZMuQuad[8][1] = 0.0; 
  fReflZMuQuad[8][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][0] = (0.02*(20.0*f[21]+2.23606797749979*(13.41640786499874*(f[15]+f[13])+5.0*(2.0*(f[9]+f[7])+3.0*(f[3]+f[1])))+5.0*(9.0*f[5]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*coshl(b)-sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][1] = (0.03333333333333333*(2.0*(3.0*(2.0*f[24]+3.0*(f[19]+f[17]))+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*coshl(b) + 4/3*(3 + b*b)*sinhl(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[8][2] = (0.1*(2.0*(2.0*f[26]+3.0*(f[25]+f[23])+2.23606797749979*(f[22]+f[20]))+9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
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
  fRefl[20] = 0.0276057774999974*(5.0*fReflZMuQuad[8][2]+8.0*fReflZMuQuad[7][2]+5.0*fReflZMuQuad[6][2]-2.0*(5.0*fReflZMuQuad[5][2]+8.0*fReflZMuQuad[4][2])+5.0*(fReflZMuQuad[2][2]-2.0*fReflZMuQuad[3][2])+8.0*fReflZMuQuad[1][2]+5.0*fReflZMuQuad[0][2]); 
  fRefl[21] = 0.1234567901234568*(fReflZMuQuad[8][0]-2.0*fReflZMuQuad[7][0]+fReflZMuQuad[6][0]+2.0*((-1.0*fReflZMuQuad[5][0])+2.0*fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])+fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[22] = 0.0276057774999974*(5.0*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2])+8.0*(fReflZMuQuad[5][2]-2.0*fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+5.0*(fReflZMuQuad[2][2]-2.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[23] = 0.1851851851851852*(fReflZMuQuad[8][2]-1.0*fReflZMuQuad[6][2]+2.0*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[5][2])+fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[24] = 0.1234567901234568*(fReflZMuQuad[8][1]-2.0*fReflZMuQuad[7][1]+fReflZMuQuad[6][1]+2.0*((-1.0*fReflZMuQuad[5][1])+2.0*fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])+fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[25] = 0.1851851851851852*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2]-1.0*fReflZMuQuad[2][2]+2.0*fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[26] = 0.1234567901234568*(fReflZMuQuad[8][2]-2.0*fReflZMuQuad[7][2]+fReflZMuQuad[6][2]+2.0*((-1.0*fReflZMuQuad[5][2])+2.0*fReflZMuQuad[4][2]-1.0*fReflZMuQuad[3][2])+fReflZMuQuad[2][2]-2.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  } 

 
}
