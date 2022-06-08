#include <gkyl_fem_parproj_kernels.h> 
 
long fem_parproj_num_nodes_global_1x_ser_p1(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return numCellsPar+1;
}

long fem_parproj_num_nodes_global_1x_ser_p2(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+3;
}

long fem_parproj_num_nodes_global_1x_ser_p3(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 3*(numCellsPar-1)+4;
}

long fem_parproj_num_nodes_global_3x_ser_p1(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 4*(numCellsPar-1)+8;
}

long fem_parproj_num_nodes_global_3x_ser_p2(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 12*(numCellsPar-1)+20;
}

long fem_parproj_num_nodes_global_3x_ser_p3(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 20*(numCellsPar-1)+32;
}

long fem_parproj_num_nodes_global_1x_tensor_p1(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return numCellsPar+1;
}

long fem_parproj_num_nodes_global_1x_tensor_p2(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+3;
}

long fem_parproj_num_nodes_global_3x_tensor_p1(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 4*(numCellsPar-1)+8;
}

long fem_parproj_num_nodes_global_3x_tensor_p2(const int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 18*(numCellsPar-1)+27;
}

