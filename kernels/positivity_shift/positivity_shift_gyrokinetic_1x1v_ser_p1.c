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
  fnod[0] = 0.2151657414559676*distf[5]-0.3726779962499649*distf[4]+0.1666666666666667*distf[3]-0.2886751345948129*distf[2]-0.2886751345948129*distf[1]+0.5*distf[0]; 
  fnod[1] = 0.3227486121839514*distf[5]-0.5590169943749475*distf[4]-0.2886751345948129*distf[1]+0.5*distf[0]; 
  fnod[2] = 0.2151657414559676*distf[5]-0.3726779962499649*distf[4]-0.1666666666666667*distf[3]+0.2886751345948129*distf[2]-0.2886751345948129*distf[1]+0.5*distf[0]; 
  fnod[3] = (-0.2151657414559676*distf[5])-0.3726779962499649*distf[4]-0.1666666666666667*distf[3]-0.2886751345948129*distf[2]+0.2886751345948129*distf[1]+0.5*distf[0]; 
  fnod[4] = (-0.3227486121839514*distf[5])-0.5590169943749475*distf[4]+0.2886751345948129*distf[1]+0.5*distf[0]; 
  fnod[5] = (-0.2151657414559676*distf[5])-0.3726779962499649*distf[4]+0.1666666666666667*distf[3]+0.2886751345948129*distf[2]+0.2886751345948129*distf[1]+0.5*distf[0]; 

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

  distf[0] = 1.5*fnod[5]-2.0*fnod[4]+1.5*fnod[3]+1.5*fnod[2]-2.0*fnod[1]+1.5*fnod[0]; 
  distf[1] = 2.598076211353316*fnod[5]-3.464101615137754*fnod[4]+2.598076211353316*fnod[3]-2.598076211353316*fnod[2]+3.464101615137754*fnod[1]-2.598076211353316*fnod[0]; 
  distf[2] = 0.8660254037844386*fnod[5]-0.8660254037844386*fnod[3]+0.8660254037844386*fnod[2]-0.8660254037844386*fnod[0]; 
  distf[3] = 1.5*fnod[5]-1.5*fnod[3]-1.5*fnod[2]+1.5*fnod[0]; 
  distf[4] = 1.341640786499874*fnod[5]-2.683281572999748*fnod[4]+1.341640786499874*fnod[3]+1.341640786499874*fnod[2]-2.683281572999748*fnod[1]+1.341640786499874*fnod[0]; 
  distf[5] = 2.32379000772445*fnod[5]-4.6475800154489*fnod[4]+2.32379000772445*fnod[3]-2.32379000772445*fnod[2]+4.6475800154489*fnod[1]-2.32379000772445*fnod[0]; 

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
