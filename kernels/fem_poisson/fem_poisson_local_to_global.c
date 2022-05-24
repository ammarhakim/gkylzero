#include <gkyl_fem_poisson_kernels.h> 
 
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

void fem_poisson_local_to_global_2x_ser_p1_periodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  }

  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  }

  if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  } else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  } else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  } else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  }


}

void fem_poisson_local_to_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  }

  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  }

  if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  } else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  } else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  } else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  }


}

void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  }

  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  }

  if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+2; 
  } else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+2; 
  } else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  } else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 
  }


}

void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+1; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  }

  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  }

  if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+3; 
  } else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+3; 
  } else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  } else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 
  }


}

void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  globalIdxs[1] = idx[1]; 
  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 
  }

  globalIdxs[3] = idx[1]+1; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  globalIdxs[1] = idx[1]; 
  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  }

  globalIdxs[3] = idx[1]+1; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]; 
  globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+1; 
  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+2; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+2; 
  }

  globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+3; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]+idx[0]; 
  globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+idx[0]+1; 
  if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 
  } else {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 
  }

  globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+idx[0]+3; 

}

void fem_poisson_local_to_global_2x_ser_p1_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  }

  globalIdxs[2] = idx[0]*numCells[1]; 
  globalIdxs[3] = (idx[0]+1)*numCells[1]; 

}

void fem_poisson_local_to_global_2x_ser_p1_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  }

  globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

}

void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 
  }

  globalIdxs[2] = idx[0]*numCells[1]; 
  globalIdxs[3] = (idx[0]+1)*numCells[1]; 

}

void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+1; 
  } else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 
  }

  globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 
  globalIdxs[1] = idx[1]; 
  globalIdxs[2] = idx[0]*numCells[1]; 
  globalIdxs[3] = 0; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 
  globalIdxs[1] = idx[1]; 
  globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 
  globalIdxs[3] = idx[1]+1; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]; 
  globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+1; 
  globalIdxs[2] = idx[0]*numCells[1]; 
  globalIdxs[3] = idx[0]*numCells[1]+1; 

}

void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

  globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]+idx[0]; 
  globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+idx[0]+1; 
  globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 
  globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+idx[0]+3; 

}

