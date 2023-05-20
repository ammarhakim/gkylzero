#include <gkyl_fem_poisson_perp_kernels.h> 
 
long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+2*(numCells[0]-1)+2;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p2_periodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 7*(numCells[0]-1)*(numCells[1]-1)+7*(numCells[1]-1)+7*(numCells[0]-1)+7;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+4*(numCells[0]-1)+4;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p2_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 7*(numCells[0]-1)*(numCells[1]-1)+7*(numCells[1]-1)+12*(numCells[0]-1)+12;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+4*(numCells[1]-1)+2*(numCells[0]-1)+4;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p2_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 7*(numCells[0]-1)*(numCells[1]-1)+12*(numCells[1]-1)+7*(numCells[0]-1)+12;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+4*(numCells[1]-1)+4*(numCells[0]-1)+8;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p2_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 7*(numCells[0]-1)*(numCells[1]-1)+12*(numCells[1]-1)+12*(numCells[0]-1)+20;
}

