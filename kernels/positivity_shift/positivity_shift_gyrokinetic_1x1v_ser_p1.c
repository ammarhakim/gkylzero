#include <gkyl_positivity_shift_gyrokinetic_kernels.h> 

GKYL_CU_DH bool positivity_shift_gyrokinetic_1x1v_ser_p1(double ffloor, double *distf, double *Deltaf) 
{ 
  // ffloor: Distribution function floor to shift to when f<0.
  // distf: distribution function.
  // Deltaf: Change in the distribution function.

  bool shifted = false;

  Deltaf[0] = -1.0*distf[0]; 
  Deltaf[1] = -1.0*distf[1]; 
  Deltaf[2] = -1.0*distf[2]; 
  Deltaf[3] = -1.0*distf[3]; 
  Deltaf[4] = -1.0*distf[4]; 
  Deltaf[5] = -1.0*distf[5]; 

  double fnod[6];
  fnod[0] = (-0.447213595499958*distf[5])+0.4472135954999579*distf[4]+0.6708203932499369*distf[3]-0.6708203932499369*distf[2]-0.5*distf[1]+0.5*distf[0]; 
  fnod[1] = 0.5590169943749476*distf[5]-0.5590169943749475*distf[4]-0.5*distf[1]+0.5*distf[0]; 
  fnod[2] = (-0.447213595499958*distf[5])+0.4472135954999579*distf[4]-0.6708203932499369*distf[3]+0.6708203932499369*distf[2]-0.5*distf[1]+0.5*distf[0]; 
  fnod[3] = 0.447213595499958*distf[5]+0.4472135954999579*distf[4]-0.6708203932499369*distf[3]-0.6708203932499369*distf[2]+0.5*distf[1]+0.5*distf[0]; 
  fnod[4] = (-0.5590169943749476*distf[5])-0.5590169943749475*distf[4]+0.5*distf[1]+0.5*distf[0]; 
  fnod[5] = 0.447213595499958*distf[5]+0.4472135954999579*distf[4]+0.6708203932499369*distf[3]+0.6708203932499369*distf[2]+0.5*distf[1]+0.5*distf[0]; 

  // If f < 0. at positivity control points, set it to ffloor.
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

  distf[0] = 0.2777777777777778*fnod[5]+0.4444444444444444*fnod[4]+0.2777777777777778*fnod[3]+0.2777777777777778*fnod[2]+0.4444444444444444*fnod[1]+0.2777777777777778*fnod[0]; 
  distf[1] = 0.2777777777777778*fnod[5]+0.4444444444444444*fnod[4]+0.2777777777777778*fnod[3]-0.2777777777777778*fnod[2]-0.4444444444444444*fnod[1]-0.2777777777777778*fnod[0]; 
  distf[2] = 0.3726779962499649*fnod[5]-0.3726779962499649*fnod[3]+0.3726779962499649*fnod[2]-0.3726779962499649*fnod[0]; 
  distf[3] = 0.3726779962499649*fnod[5]-0.3726779962499649*fnod[3]-0.3726779962499649*fnod[2]+0.3726779962499649*fnod[0]; 
  distf[4] = 0.2484519974999766*fnod[5]-0.4969039949999532*fnod[4]+0.2484519974999766*fnod[3]+0.2484519974999766*fnod[2]-0.4969039949999532*fnod[1]+0.2484519974999766*fnod[0]; 
  distf[5] = 0.2484519974999767*fnod[5]-0.4969039949999535*fnod[4]+0.2484519974999767*fnod[3]-0.2484519974999767*fnod[2]+0.4969039949999535*fnod[1]-0.2484519974999767*fnod[0]; 

  if (shifted) {
  Deltaf[0] += distf[0]; 
  Deltaf[1] += distf[1]; 
  Deltaf[2] += distf[2]; 
  Deltaf[3] += distf[3]; 
  Deltaf[4] += distf[4]; 
  Deltaf[5] += distf[5]; 
  }

  return shifted;

}
