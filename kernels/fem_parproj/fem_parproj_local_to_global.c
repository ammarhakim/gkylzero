#include <gkyl_fem_parproj_kernels.h> 
 
void fem_parproj_local_to_global_1x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = parIdx; 
  globalIdxs[1] = parIdx+1; 

}

void fem_parproj_local_to_global_1x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 2*parIdx; 
  globalIdxs[1] = 2*parIdx+1; 
  globalIdxs[2] = 2*parIdx+2; 

}

void fem_parproj_local_to_global_1x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 3*parIdx; 
  globalIdxs[1] = 3*parIdx+1; 
  globalIdxs[2] = 3*parIdx+2; 
  globalIdxs[3] = 3*parIdx+3; 

}

void fem_parproj_local_to_global_3x_ser_p1(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 4*parIdx; 
  globalIdxs[1] = 4*parIdx+1; 
  globalIdxs[2] = 4*parIdx+2; 
  globalIdxs[3] = 4*parIdx+3; 
  globalIdxs[4] = 4*parIdx+4; 
  globalIdxs[5] = 4*parIdx+5; 
  globalIdxs[6] = 4*parIdx+6; 
  globalIdxs[7] = 4*parIdx+7; 

}

void fem_parproj_local_to_global_3x_ser_p2(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 12*parIdx; 
  globalIdxs[1] = 12*parIdx+1; 
  globalIdxs[2] = 12*parIdx+2; 
  globalIdxs[3] = 12*parIdx+3; 
  globalIdxs[4] = 12*parIdx+4; 
  globalIdxs[5] = 12*parIdx+5; 
  globalIdxs[6] = 12*parIdx+6; 
  globalIdxs[7] = 12*parIdx+7; 
  globalIdxs[8] = 12*parIdx+8; 
  globalIdxs[9] = 12*parIdx+9; 
  globalIdxs[10] = 12*parIdx+10; 
  globalIdxs[11] = 12*parIdx+11; 
  globalIdxs[12] = 12*parIdx+12; 
  globalIdxs[13] = 12*parIdx+13; 
  globalIdxs[14] = 12*parIdx+14; 
  globalIdxs[15] = 12*parIdx+15; 
  globalIdxs[16] = 12*parIdx+16; 
  globalIdxs[17] = 12*parIdx+17; 
  globalIdxs[18] = 12*parIdx+18; 
  globalIdxs[19] = 12*parIdx+19; 

}

void fem_parproj_local_to_global_3x_ser_p3(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 20*parIdx; 
  globalIdxs[1] = 20*parIdx+1; 
  globalIdxs[2] = 20*parIdx+2; 
  globalIdxs[3] = 20*parIdx+3; 
  globalIdxs[4] = 20*parIdx+4; 
  globalIdxs[5] = 20*parIdx+5; 
  globalIdxs[6] = 20*parIdx+6; 
  globalIdxs[7] = 20*parIdx+7; 
  globalIdxs[8] = 20*parIdx+8; 
  globalIdxs[9] = 20*parIdx+9; 
  globalIdxs[10] = 20*parIdx+10; 
  globalIdxs[11] = 20*parIdx+11; 
  globalIdxs[12] = 20*parIdx+12; 
  globalIdxs[13] = 20*parIdx+13; 
  globalIdxs[14] = 20*parIdx+14; 
  globalIdxs[15] = 20*parIdx+15; 
  globalIdxs[16] = 20*parIdx+16; 
  globalIdxs[17] = 20*parIdx+17; 
  globalIdxs[18] = 20*parIdx+18; 
  globalIdxs[19] = 20*parIdx+19; 
  globalIdxs[20] = 20*parIdx+20; 
  globalIdxs[21] = 20*parIdx+21; 
  globalIdxs[22] = 20*parIdx+22; 
  globalIdxs[23] = 20*parIdx+23; 
  globalIdxs[24] = 20*parIdx+24; 
  globalIdxs[25] = 20*parIdx+25; 
  globalIdxs[26] = 20*parIdx+26; 
  globalIdxs[27] = 20*parIdx+27; 
  globalIdxs[28] = 20*parIdx+28; 
  globalIdxs[29] = 20*parIdx+29; 
  globalIdxs[30] = 20*parIdx+30; 
  globalIdxs[31] = 20*parIdx+31; 

}

void fem_parproj_local_to_global_1x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = parIdx; 
  globalIdxs[1] = parIdx+1; 

}

void fem_parproj_local_to_global_1x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 2*parIdx; 
  globalIdxs[1] = 2*parIdx+1; 
  globalIdxs[2] = 2*parIdx+2; 

}

void fem_parproj_local_to_global_3x_tensor_p1(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 4*parIdx; 
  globalIdxs[1] = 4*parIdx+1; 
  globalIdxs[2] = 4*parIdx+2; 
  globalIdxs[3] = 4*parIdx+3; 
  globalIdxs[4] = 4*parIdx+4; 
  globalIdxs[5] = 4*parIdx+5; 
  globalIdxs[6] = 4*parIdx+6; 
  globalIdxs[7] = 4*parIdx+7; 

}

void fem_parproj_local_to_global_3x_tensor_p2(const int numCellsPar, const int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = 18*parIdx; 
  globalIdxs[1] = 18*parIdx+1; 
  globalIdxs[2] = 18*parIdx+2; 
  globalIdxs[3] = 18*parIdx+3; 
  globalIdxs[4] = 18*parIdx+4; 
  globalIdxs[5] = 18*parIdx+5; 
  globalIdxs[6] = 18*parIdx+6; 
  globalIdxs[7] = 18*parIdx+7; 
  globalIdxs[8] = 18*parIdx+8; 
  globalIdxs[9] = 18*parIdx+9; 
  globalIdxs[10] = 18*parIdx+10; 
  globalIdxs[11] = 18*parIdx+11; 
  globalIdxs[12] = 18*parIdx+12; 
  globalIdxs[13] = 18*parIdx+13; 
  globalIdxs[14] = 18*parIdx+14; 
  globalIdxs[15] = 18*parIdx+15; 
  globalIdxs[16] = 18*parIdx+16; 
  globalIdxs[17] = 18*parIdx+17; 
  globalIdxs[18] = 18*parIdx+18; 
  globalIdxs[19] = 18*parIdx+19; 
  globalIdxs[20] = 18*parIdx+20; 
  globalIdxs[21] = 18*parIdx+21; 
  globalIdxs[22] = 18*parIdx+22; 
  globalIdxs[23] = 18*parIdx+23; 
  globalIdxs[24] = 18*parIdx+24; 
  globalIdxs[25] = 18*parIdx+25; 
  globalIdxs[26] = 18*parIdx+26; 

}

