#include <gkyl_fem_poisson_kernels.h> 
 
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

int fem_poisson_num_nodes_global_2x_ser_p1_periodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[1]+(numCells[0]-1)*(numCells[1]-1)+numCells[0]-1;
}

int fem_poisson_num_nodes_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[1]+(numCells[0]-1)*(numCells[1]-1)+2*(numCells[0]-1)+1;
}

int fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return (numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+numCells[0]+1;
}

int fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return (numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+2*(numCells[0]-1)+4;
}

