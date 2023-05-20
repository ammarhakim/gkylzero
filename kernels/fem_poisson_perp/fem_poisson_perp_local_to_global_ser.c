#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+6; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+6; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+17; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+17; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+9; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+9; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+19; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+16; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+19; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+16; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+4*idx[1]; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+4*idx[1]+4; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+4*idx[1]+4; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+4*idx[1]+6; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+4*idx[1]+6; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+12*idx[1]; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+12*idx[1]+3; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+12*idx[1]+12; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+12*idx[1]+12; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+7; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+12*idx[1]+17; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+12*idx[1]+17; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+10; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+12*idx[1]+7; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+12*idx[1]+10; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+12*idx[1]+19; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+12*idx[1]+19; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+11; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+6; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+6; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+10; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+8; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+6; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+17; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+17; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+10; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+25; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+22; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+17; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+9; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+29; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+24; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+19; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+16; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

    globalIdxs[1] = 2*idx[1]; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = 2*idx[1]+2; 

  }  else {
    globalIdxs[3] = 2*idx[1]+2; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

    globalIdxs[5] = 2*idx[1]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 2*idx[1]+3; 

  }  else {
    globalIdxs[7] = 2*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

    globalIdxs[2] = 7*idx[1]; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

    globalIdxs[4] = 7*idx[1]+2; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+8; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 7*idx[1]+7; 

  }  else {
    globalIdxs[7] = 7*idx[1]+7; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

    globalIdxs[9] = 7*idx[1]+3; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+10; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = 7*idx[1]+10; 

  }  else {
    globalIdxs[11] = 7*idx[1]+10; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

    globalIdxs[14] = 7*idx[1]+4; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

    globalIdxs[16] = 7*idx[1]+6; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+12; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = 7*idx[1]+11; 

  }  else {
    globalIdxs[19] = 7*idx[1]+11; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

    globalIdxs[1] = 2*idx[1]; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = 2*idx[1]+2; 

  }  else {
    globalIdxs[3] = 2*idx[1]+2; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

    globalIdxs[5] = 2*idx[1]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 2*idx[1]+4; 

  }  else {
    globalIdxs[7] = 2*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

    globalIdxs[2] = 7*idx[1]; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

    globalIdxs[4] = 7*idx[1]+2; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 7*idx[1]+7; 

  }  else {
    globalIdxs[7] = 7*idx[1]+7; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

    globalIdxs[9] = 7*idx[1]+3; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = 7*idx[1]+12; 

  }  else {
    globalIdxs[11] = 7*idx[1]+10; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

    globalIdxs[14] = 7*idx[1]+4; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

    globalIdxs[16] = 7*idx[1]+6; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = 7*idx[1]+14; 

  }  else {
    globalIdxs[19] = 7*idx[1]+11; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+4*idx[1]; 

    globalIdxs[1] = 2*idx[0]*numCells[1]+4*idx[1]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+4*idx[1]+4; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+4*idx[1]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = 2*idx[0]*numCells[1]+4*idx[1]+5; 

  }  else {
    globalIdxs[3] = 2*idx[0]*numCells[1]+4*idx[1]+5; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+4*idx[1]+2; 

    globalIdxs[5] = 2*idx[0]*numCells[1]+4*idx[1]+3; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+4*idx[1]+6; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+4*idx[1]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 2*idx[0]*numCells[1]+4*idx[1]+7; 

  }  else {
    globalIdxs[7] = 2*idx[0]*numCells[1]+4*idx[1]+7; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+12*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+12*idx[1]+1; 

    globalIdxs[2] = 7*idx[0]*numCells[1]+12*idx[1]+2; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+12*idx[1]+3; 

    globalIdxs[4] = 7*idx[0]*numCells[1]+12*idx[1]+4; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+12*idx[1]+12; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+12*idx[1]+12; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+12*idx[1]+13; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+12*idx[1]+13; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 7*idx[0]*numCells[1]+12*idx[1]+14; 

  }  else {
    globalIdxs[7] = 7*idx[0]*numCells[1]+12*idx[1]+14; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+12*idx[1]+5; 

    globalIdxs[9] = 7*idx[0]*numCells[1]+12*idx[1]+6; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+12*idx[1]+17; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+12*idx[1]+17; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = 7*idx[0]*numCells[1]+12*idx[1]+18; 

  }  else {
    globalIdxs[11] = 7*idx[0]*numCells[1]+12*idx[1]+18; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+12*idx[1]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+12*idx[1]+8; 

    globalIdxs[14] = 7*idx[0]*numCells[1]+12*idx[1]+9; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+12*idx[1]+10; 

    globalIdxs[16] = 7*idx[0]*numCells[1]+12*idx[1]+11; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+12*idx[1]+19; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+12*idx[1]+19; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+12*idx[1]+20; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+12*idx[1]+20; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = 7*idx[0]*numCells[1]+12*idx[1]+21; 

  }  else {
    globalIdxs[19] = 7*idx[0]*numCells[1]+12*idx[1]+21; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]; 

    globalIdxs[1] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[2] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[3] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+5; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+2; 

    globalIdxs[5] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+3; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+8; 

  }  else {
    globalIdxs[6] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+9; 

  }  else {
    globalIdxs[7] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+7; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+1; 

    globalIdxs[2] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+2; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+3; 

    globalIdxs[4] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+4; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[5] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+12; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+13; 

  }  else {
    globalIdxs[6] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+13; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+14; 

  }  else {
    globalIdxs[7] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+14; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+5; 

    globalIdxs[9] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+6; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[10] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+20; 

  }  else {
    globalIdxs[10] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+17; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[11] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+21; 

  }  else {
    globalIdxs[11] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+18; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+8; 

    globalIdxs[14] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+9; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+10; 

    globalIdxs[16] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+11; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[17] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+24; 

  }  else {
    globalIdxs[17] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+19; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[18] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+25; 

  }  else {
    globalIdxs[18] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+20; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[19] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+26; 

  }  else {
    globalIdxs[19] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+21; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }

    globalIdxs[2] = 2*idx[0]*numCells[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }

    globalIdxs[6] = 2*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+1; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }

    globalIdxs[5] = 7*idx[0]*numCells[1]; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }

    globalIdxs[10] = 7*idx[0]*numCells[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+3; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+3; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }

    globalIdxs[17] = 7*idx[0]*numCells[1]+4; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+4; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+4; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }

    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+9; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }

    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+4*idx[1]; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]; 

  }

    globalIdxs[2] = 2*idx[0]*numCells[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+1; 

  }

    globalIdxs[6] = 2*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+12*idx[1]; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+12*idx[1]+3; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+2; 

  }

    globalIdxs[5] = 7*idx[0]*numCells[1]; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+3; 

  }

    globalIdxs[10] = 7*idx[0]*numCells[1]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+5; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+3; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+12*idx[1]+7; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+4; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+12*idx[1]+10; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+6; 

  }

    globalIdxs[17] = 7*idx[0]*numCells[1]+4; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+4; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+3; 

  }

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+6; 

  }  else {
    globalIdxs[5] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+4*idx[1]+2*idx[0]+8; 

  }  else {
    globalIdxs[7] = (2*idx[0]+2)*numCells[1]+2*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+5; 

  }

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+7; 

  }

    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+10; 

  }  else {
    globalIdxs[7] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+8; 

  }

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+13; 

  }  else {
    globalIdxs[9] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+10; 

  }

    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+15; 

  }  else {
    globalIdxs[11] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+11; 

  }

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+17; 

  }  else {
    globalIdxs[14] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+12; 

  }

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+9; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+20; 

  }  else {
    globalIdxs[16] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+14; 

  }

    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+12*idx[1]+5*idx[0]+22; 

  }  else {
    globalIdxs[19] = (7*idx[0]+7)*numCells[1]+7*idx[1]+5*idx[0]+15; 

  }


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]; 

    globalIdxs[1] = 2*idx[1]; 

    globalIdxs[2] = 2*idx[0]*numCells[1]; 

    globalIdxs[3] = 0; 

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+1; 

    globalIdxs[5] = 2*idx[1]+1; 

    globalIdxs[6] = 2*idx[0]*numCells[1]+1; 

    globalIdxs[7] = 1; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+1; 

    globalIdxs[2] = 7*idx[1]; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+2; 

    globalIdxs[4] = 7*idx[1]+2; 

    globalIdxs[5] = 7*idx[0]*numCells[1]; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+1; 

    globalIdxs[7] = 0; 

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+3; 

    globalIdxs[9] = 7*idx[1]+3; 

    globalIdxs[10] = 7*idx[0]*numCells[1]+3; 

    globalIdxs[11] = 3; 

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+4; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5; 

    globalIdxs[14] = 7*idx[1]+4; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+6; 

    globalIdxs[16] = 7*idx[1]+6; 

    globalIdxs[17] = 7*idx[0]*numCells[1]+4; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+5; 

    globalIdxs[19] = 4; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]; 

    globalIdxs[1] = 2*idx[1]; 

    globalIdxs[2] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+1; 

    globalIdxs[3] = 2*idx[1]+1; 

    globalIdxs[4] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+2; 

    globalIdxs[5] = 2*idx[1]+2; 

    globalIdxs[6] = 2*idx[0]*numCells[1]+2*idx[1]+2*idx[0]+3; 

    globalIdxs[7] = 2*idx[1]+3; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+1; 

    globalIdxs[2] = 7*idx[1]; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+2; 

    globalIdxs[4] = 7*idx[1]+2; 

    globalIdxs[5] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+3; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+4; 

    globalIdxs[7] = 7*idx[1]+3; 

    globalIdxs[8] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+5; 

    globalIdxs[9] = 7*idx[1]+5; 

    globalIdxs[10] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+6; 

    globalIdxs[11] = 7*idx[1]+6; 

    globalIdxs[12] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+8; 

    globalIdxs[14] = 7*idx[1]+7; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+9; 

    globalIdxs[16] = 7*idx[1]+9; 

    globalIdxs[17] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+10; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+7*idx[1]+5*idx[0]+11; 

    globalIdxs[19] = 7*idx[1]+10; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+4*idx[1]; 

    globalIdxs[1] = 2*idx[0]*numCells[1]+4*idx[1]+1; 

    globalIdxs[2] = 2*idx[0]*numCells[1]; 

    globalIdxs[3] = 2*idx[0]*numCells[1]+1; 

    globalIdxs[4] = 2*idx[0]*numCells[1]+4*idx[1]+2; 

    globalIdxs[5] = 2*idx[0]*numCells[1]+4*idx[1]+3; 

    globalIdxs[6] = 2*idx[0]*numCells[1]+2; 

    globalIdxs[7] = 2*idx[0]*numCells[1]+3; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+12*idx[1]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+12*idx[1]+1; 

    globalIdxs[2] = 7*idx[0]*numCells[1]+12*idx[1]+2; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+12*idx[1]+3; 

    globalIdxs[4] = 7*idx[0]*numCells[1]+12*idx[1]+4; 

    globalIdxs[5] = 7*idx[0]*numCells[1]; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+1; 

    globalIdxs[7] = 7*idx[0]*numCells[1]+2; 

    globalIdxs[8] = 7*idx[0]*numCells[1]+12*idx[1]+5; 

    globalIdxs[9] = 7*idx[0]*numCells[1]+12*idx[1]+6; 

    globalIdxs[10] = 7*idx[0]*numCells[1]+5; 

    globalIdxs[11] = 7*idx[0]*numCells[1]+6; 

    globalIdxs[12] = 7*idx[0]*numCells[1]+12*idx[1]+7; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+12*idx[1]+8; 

    globalIdxs[14] = 7*idx[0]*numCells[1]+12*idx[1]+9; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+12*idx[1]+10; 

    globalIdxs[16] = 7*idx[0]*numCells[1]+12*idx[1]+11; 

    globalIdxs[17] = 7*idx[0]*numCells[1]+7; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+8; 

    globalIdxs[19] = 7*idx[0]*numCells[1]+9; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]; 

    globalIdxs[1] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+1; 

    globalIdxs[2] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+2; 

    globalIdxs[3] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+3; 

    globalIdxs[4] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+4; 

    globalIdxs[5] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+5; 

    globalIdxs[6] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+6; 

    globalIdxs[7] = 2*idx[0]*numCells[1]+4*idx[1]+2*idx[0]+7; 


}

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]; 

    globalIdxs[1] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+1; 

    globalIdxs[2] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+2; 

    globalIdxs[3] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+3; 

    globalIdxs[4] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+4; 

    globalIdxs[5] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+5; 

    globalIdxs[6] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+6; 

    globalIdxs[7] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+7; 

    globalIdxs[8] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+8; 

    globalIdxs[9] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+9; 

    globalIdxs[10] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+10; 

    globalIdxs[11] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+11; 

    globalIdxs[12] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+12; 

    globalIdxs[13] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+13; 

    globalIdxs[14] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+14; 

    globalIdxs[15] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+15; 

    globalIdxs[16] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+16; 

    globalIdxs[17] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+17; 

    globalIdxs[18] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+18; 

    globalIdxs[19] = 7*idx[0]*numCells[1]+12*idx[1]+5*idx[0]+19; 


}

