#include <gkyl_fem_poisson_kernels.h> 
 
GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_lox_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_lox_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_lox_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_lox_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_lox_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[0]*bcVals[2])-0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[0]*bcVals[2])-0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_lox_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-1.0*rdx2[0]*bcVals[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-1.0*rdx2[0]*bcVals[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_lox_robinx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[0]*bcVals[2])/bcVals[1])-0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[0]*bcVals[2])/bcVals[1])-0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_lox_robinx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_upx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_upx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_upx_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_upx_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_upx_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]+0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]+0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_upx_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[0]*bcVals[5]+0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[0]*bcVals[5]+0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p1_upx_robinx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]+0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_1x_ser_p2_upx_robinx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[1]; 
  rdx2[0] = 2.0/dx[0]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[0]*bcVals[5])/bcVals[4]+0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[2]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[0]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[3]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[2]));
  #else
  bsrc[globalIdxs[5]] = bcVals[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-1.0*rdx2[0]*bcVals[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-1.333333333333333*rdx2[0]*bcVals[2]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.3333333333333333*rdx2[0]*bcVals[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_lox_robinx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-(1.0*rdx2[0]*bcVals[2])/bcVals[1]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_lox_robinx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-(1.333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-(0.3333333333333333*rdx2[0]*bcVals[2])/bcVals[1]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])+rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])+rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[0]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[1]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[8]));
  #else
  bsrc[globalIdxs[2]] = bcVals[8];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-1.0*rdx2[1]*bcVals[8])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.0*rdx2[1]*bcVals[8])+(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.0*rdx2[1]*bcVals[8])+(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.3333333333333333*rdx2[1]*bcVals[8])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-1.333333333333333*rdx2[1]*bcVals[8])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.3333333333333333*rdx2[1]*bcVals[8])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.0*rdx2[1]*bcVals[8])/bcVals[7])+(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_loy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])-0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-(1.333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-(0.3333333333333333*rdx2[1]*bcVals[8])/bcVals[7])+0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[3]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[2]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[4]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += rdx2[0]*bcVals[5]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+rdx2[0]*bcVals[5]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+1.333333333333333*rdx2[0]*bcVals[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.3333333333333333*rdx2[0]*bcVals[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[2]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[3]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[5]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[6]] = bcVals[11];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(bcVals[11]));
  #else
  bsrc[globalIdxs[7]] = bcVals[11];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += rdx2[1]*bcVals[11]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],rdx2[1]*bcVals[11]+(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += rdx2[1]*bcVals[11]+(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += 0.3333333333333333*rdx2[1]*bcVals[11]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += 1.333333333333333*rdx2[1]*bcVals[11]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += 0.3333333333333333*rdx2[1]*bcVals[11]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p1_upx_robinx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[0]] += 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],(rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[1]] += (rdx2[0]*bcVals[5])/bcVals[4]-0.1666666666666667*rho[3]-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],(rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[2]] += (rdx2[1]*bcVals[11])/bcVals[10]-0.1666666666666667*rho[3]+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],(rdx2[1]*bcVals[11])/bcVals[10]+(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[globalIdxs[3]] += (rdx2[1]*bcVals[11])/bcVals[10]+(rdx2[0]*bcVals[5])/bcVals[4]+0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_src_stencil_vareps_2x_ser_p2_upx_robinx_upy_robiny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[0]],(-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[0]] += (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[1]],0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[1]] += 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[2]],0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[2]] += 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[3]],0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[3]] += 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[4]],(-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[4]] += (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+(1.333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[5]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[5]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.08606629658238707*rho[7]+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[6]],(1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  #else
  bsrc[globalIdxs[6]] += (1.333333333333333*rdx2[1]*bcVals[11])/bcVals[10]-0.1721325931647741*rho[6]-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[globalIdxs[7]],(0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  #else
  bsrc[globalIdxs[7]] += (0.3333333333333333*rdx2[1]*bcVals[11])/bcVals[10]+0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+(0.3333333333333333*rdx2[0]*bcVals[5])/bcVals[4]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0];
  #endif

}

