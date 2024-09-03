#include <gkyl_positivity_shift_gyrokinetic_kernels.h> 

GKYL_CU_DH bool positivity_shift_gyrokinetic_1x1v_ser_p1(double ffloor, double *distf) 
{ 
  // ffloor: Distribution function floor to shift to when f<0.
  // distf: distribution function.

  bool shifted = false;

  double fnod[6];
  fnod[0] = -(0.44721359549995804*distf[5])+0.4472135954999579*distf[4]+0.6708203932499369*distf[3]-0.6708203932499369*distf[2]-0.5*distf[1]+0.5*distf[0]; 
  fnod[1] = 0.5590169943749476*distf[5]-0.5590169943749475*distf[4]-0.5*distf[1]+0.5*distf[0]; 
  fnod[2] = -(0.44721359549995804*distf[5])+0.4472135954999579*distf[4]-0.6708203932499369*distf[3]+0.6708203932499369*distf[2]-0.5*distf[1]+0.5*distf[0]; 
  fnod[3] = 0.44721359549995804*distf[5]+0.4472135954999579*distf[4]-0.6708203932499369*distf[3]-0.6708203932499369*distf[2]+0.5*distf[1]+0.5*distf[0]; 
  fnod[4] = -(0.5590169943749476*distf[5])-0.5590169943749475*distf[4]+0.5*distf[1]+0.5*distf[0]; 
  fnod[5] = 0.44721359549995804*distf[5]+0.4472135954999579*distf[4]+0.6708203932499369*distf[3]+0.6708203932499369*distf[2]+0.5*distf[1]+0.5*distf[0]; 

  // If f < 0. at check nodes, set it to ffloor.
  if (fnod[0] < 0.) {
    fnod[0] = ffloor;
    shifted = true;
  }
  if (fnod[1] < 0.) {
    fnod[1] = ffloor;
    shifted = true;
  }
  if (fnod[2] < 0.) {
    fnod[2] = ffloor;
    shifted = true;
  }
  if (fnod[3] < 0.) {
    fnod[3] = ffloor;
    shifted = true;
  }
  if (fnod[4] < 0.) {
    fnod[4] = ffloor;
    shifted = true;
  }
  if (fnod[5] < 0.) {
    fnod[5] = ffloor;
    shifted = true;
  }

  if (shifted) {
  distf[0] = 0.2777777777777778*fnod[5]+0.4444444444444444*fnod[4]+0.2777777777777778*fnod[3]+0.2777777777777778*fnod[2]+0.4444444444444444*fnod[1]+0.2777777777777778*fnod[0]; 
  distf[1] = 0.2777777777777778*fnod[5]+0.4444444444444444*fnod[4]+0.2777777777777778*fnod[3]-0.2777777777777778*fnod[2]-0.4444444444444444*fnod[1]-0.2777777777777778*fnod[0]; 
  distf[2] = 0.37267799624996495*fnod[5]-0.37267799624996495*fnod[3]+0.37267799624996495*fnod[2]-0.37267799624996495*fnod[0]; 
  distf[3] = 0.37267799624996495*fnod[5]-0.37267799624996495*fnod[3]-0.37267799624996495*fnod[2]+0.37267799624996495*fnod[0]; 
  distf[4] = 0.24845199749997662*fnod[5]-0.49690399499995325*fnod[4]+0.24845199749997662*fnod[3]+0.24845199749997662*fnod[2]-0.49690399499995325*fnod[1]+0.24845199749997662*fnod[0]; 
  distf[5] = 0.24845199749997673*fnod[5]-0.49690399499995347*fnod[4]+0.24845199749997673*fnod[3]-0.24845199749997673*fnod[2]+0.49690399499995347*fnod[1]-0.24845199749997673*fnod[0]; 
  }

  return shifted;

}
