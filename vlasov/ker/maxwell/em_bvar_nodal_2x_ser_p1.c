#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH int em_bvar_nodal_2x_ser_p1(const double *em, double* GKYL_RESTRICT b) 
{ 
  // em: Input electromagnetic fields [Ex, Ey, Ez, Bx, By, Bz, phi, psi]. 
  // b:  Output magnetic field unit vector [bx, by, bz] = B/|B|, normalized at the 
  //     4 Gauss-Legendre nodes and interpolated to the p=1 basis. 
  // Returns the number of nodes where |B| = 0 (b is set to zero at such nodes). 
 
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  double *bx = &b[0]; 
  double *by = &b[4]; 
  double *bz = &b[8]; 
 
  double bx_n[4] = {0.0}; 
  double by_n[4] = {0.0}; 
  double bz_n[4] = {0.0}; 
  int num_null_b = 0; 
 
  { 
    const double Bx_q = 0.5*B_x[3]-0.5*B_x[2]-0.5*B_x[1]+0.5*B_x[0]; 
    const double By_q = 0.5*B_y[3]-0.5*B_y[2]-0.5*B_y[1]+0.5*B_y[0]; 
    const double Bz_q = 0.5*B_z[3]-0.5*B_z[2]-0.5*B_z[1]+0.5*B_z[0]; 
    const double magB_q = sqrt(Bx_q*Bx_q + By_q*By_q + Bz_q*Bz_q); 
    if (magB_q > 0.0) { 
      bx_n[0] = Bx_q/magB_q; 
      by_n[0] = By_q/magB_q; 
      bz_n[0] = Bz_q/magB_q; 
    } else { 
      num_null_b += 1; 
    } 
  } 
  { 
    const double Bx_q = -(0.5*B_x[3])+0.5*B_x[2]-0.5*B_x[1]+0.5*B_x[0]; 
    const double By_q = -(0.5*B_y[3])+0.5*B_y[2]-0.5*B_y[1]+0.5*B_y[0]; 
    const double Bz_q = -(0.5*B_z[3])+0.5*B_z[2]-0.5*B_z[1]+0.5*B_z[0]; 
    const double magB_q = sqrt(Bx_q*Bx_q + By_q*By_q + Bz_q*Bz_q); 
    if (magB_q > 0.0) { 
      bx_n[1] = Bx_q/magB_q; 
      by_n[1] = By_q/magB_q; 
      bz_n[1] = Bz_q/magB_q; 
    } else { 
      num_null_b += 1; 
    } 
  } 
  { 
    const double Bx_q = -(0.5*B_x[3])-0.5*B_x[2]+0.5*B_x[1]+0.5*B_x[0]; 
    const double By_q = -(0.5*B_y[3])-0.5*B_y[2]+0.5*B_y[1]+0.5*B_y[0]; 
    const double Bz_q = -(0.5*B_z[3])-0.5*B_z[2]+0.5*B_z[1]+0.5*B_z[0]; 
    const double magB_q = sqrt(Bx_q*Bx_q + By_q*By_q + Bz_q*Bz_q); 
    if (magB_q > 0.0) { 
      bx_n[2] = Bx_q/magB_q; 
      by_n[2] = By_q/magB_q; 
      bz_n[2] = Bz_q/magB_q; 
    } else { 
      num_null_b += 1; 
    } 
  } 
  { 
    const double Bx_q = 0.5*B_x[3]+0.5*B_x[2]+0.5*B_x[1]+0.5*B_x[0]; 
    const double By_q = 0.5*B_y[3]+0.5*B_y[2]+0.5*B_y[1]+0.5*B_y[0]; 
    const double Bz_q = 0.5*B_z[3]+0.5*B_z[2]+0.5*B_z[1]+0.5*B_z[0]; 
    const double magB_q = sqrt(Bx_q*Bx_q + By_q*By_q + Bz_q*Bz_q); 
    if (magB_q > 0.0) { 
      bx_n[3] = Bx_q/magB_q; 
      by_n[3] = By_q/magB_q; 
      bz_n[3] = Bz_q/magB_q; 
    } else { 
      num_null_b += 1; 
    } 
  } 
 
  bx[0] = 0.5*bx_n[3]+0.5*bx_n[2]+0.5*bx_n[1]+0.5*bx_n[0]; 
  bx[1] = 0.5*bx_n[3]+0.5*bx_n[2]-0.5*bx_n[1]-0.5*bx_n[0]; 
  bx[2] = 0.5*bx_n[3]-0.5*bx_n[2]+0.5*bx_n[1]-0.5*bx_n[0]; 
  bx[3] = 0.5*bx_n[3]-0.5*bx_n[2]-0.5*bx_n[1]+0.5*bx_n[0]; 
  by[0] = 0.5*by_n[3]+0.5*by_n[2]+0.5*by_n[1]+0.5*by_n[0]; 
  by[1] = 0.5*by_n[3]+0.5*by_n[2]-0.5*by_n[1]-0.5*by_n[0]; 
  by[2] = 0.5*by_n[3]-0.5*by_n[2]+0.5*by_n[1]-0.5*by_n[0]; 
  by[3] = 0.5*by_n[3]-0.5*by_n[2]-0.5*by_n[1]+0.5*by_n[0]; 
  bz[0] = 0.5*bz_n[3]+0.5*bz_n[2]+0.5*bz_n[1]+0.5*bz_n[0]; 
  bz[1] = 0.5*bz_n[3]+0.5*bz_n[2]-0.5*bz_n[1]-0.5*bz_n[0]; 
  bz[2] = 0.5*bz_n[3]-0.5*bz_n[2]+0.5*bz_n[1]-0.5*bz_n[0]; 
  bz[3] = 0.5*bz_n[3]-0.5*bz_n[2]-0.5*bz_n[1]+0.5*bz_n[0]; 
 
  return num_null_b; 
} 
