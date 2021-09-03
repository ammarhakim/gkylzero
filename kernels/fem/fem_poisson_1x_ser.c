#include <gkyl_fem_poisson.h> 
 
int fem_poisson_num_nodes_global_1x_ser_p1_periodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[0];
}

int fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[0]+1;
}

void fem_poisson_local_to_global_1x_ser_p1_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = idx[0]+1; 
  } else {
    globalIdxs[1] = idx[0]+1; 
  }


}

void fem_poisson_local_to_global_1x_ser_p1_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = idx[0]+1; 
  } else {
    globalIdxs[1] = idx[0]+1; 
  }


}

void fem_poisson_local_to_global_1x_ser_p1_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]; 
  globalIdxs[1] = 0; 

}

void fem_poisson_local_to_global_1x_ser_p1_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]; 
  globalIdxs[1] = idx[0]+1; 

}

