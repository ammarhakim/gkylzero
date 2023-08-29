#include <gkyl_fem_parproj_kernels.h> 
 
GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1_inx_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2_inx_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1_lox_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2_lox_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1_lox_dirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(1.224744871391589*phiBC[1]+0.7071067811865476*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[0]] = 1.224744871391589*phiBC[1]+0.7071067811865476*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2_lox_dirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(1.58113883008419*phiBC[2]+1.224744871391589*phiBC[1]+0.7071067811865475*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[0]] = 1.58113883008419*phiBC[2]+1.224744871391589*phiBC[1]+0.7071067811865475*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1_upx_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.408248290463863*rho[1]+0.7071067811865476*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.408248290463863*rho[1]+0.7071067811865476*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2_upx_nondirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p1_upx_dirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.7071067811865476*rho[0]-0.408248290463863*rho[1];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(0.7071067811865476*phiBC[0]-1.224744871391589*phiBC[1]));
  #else
  bsrc[nodeOff+globalIdxs[1]] = 0.7071067811865476*phiBC[0]-1.224744871391589*phiBC[1];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_1x_ser_p2_upx_dirichletx(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.9428090415820636*rho[0]-0.4216370213557841*rho[2];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(1.58113883008419*phiBC[2]-1.224744871391589*phiBC[1]+0.7071067811865475*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[2]] = 1.58113883008419*phiBC[2]-1.224744871391589*phiBC[1]+0.7071067811865475*phiBC[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1_inz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2_inz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],(-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += (-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],(-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += (-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += 0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[8]],(-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[8]] += (-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[9]],0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[9]] += 0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[10]],0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[10]] += 0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[11]],(-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[11]] += (-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[12]],0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[12]] += 0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[13]],0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[13]] += 0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[14]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[14]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[15]],0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[15]] += 0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[16]],(-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[16]] += (-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[17]],(-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[17]] += (-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[18]],(-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[18]] += (-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[19]],0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[19]] += 0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1_loz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2_loz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],(-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += (-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],(-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += (-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += 0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[8]],(-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[8]] += (-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[9]],0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[9]] += 0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[10]],0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[10]] += 0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[11]],(-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[11]] += (-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[12]],0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[12]] += 0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[13]],0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[13]] += 0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[14]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[14]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[15]],0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[15]] += 0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[16]],(-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[16]] += (-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[17]],(-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[17]] += (-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[18]],(-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[18]] += (-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[19]],0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[19]] += 0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1_loz_dirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(1.837117307087384*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[0]] = 1.837117307087384*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong((-1.837117307087384*phiBC[7])-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[1]] = (-1.837117307087384*phiBC[7])-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong((-1.837117307087384*phiBC[7])+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[2]] = (-1.837117307087384*phiBC[7])+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(1.837117307087384*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[3]] = 1.837117307087384*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2_loz_dirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(2.371708245126285*phiBC[19]-2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[0]] = 2.371708245126285*phiBC[19]-2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(1.185854122563142*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[14]-0.6846531968814578*phiBC[13]+0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]-1.060660171779821*phiBC[6]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[1]] = 1.185854122563142*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[14]-0.6846531968814578*phiBC[13]+0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]-1.060660171779821*phiBC[6]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong((-2.371708245126285*phiBC[19])+2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[2]] = (-2.371708245126285*phiBC[19])+2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(1.185854122563142*phiBC[18]-1.369306393762916*phiBC[15]-0.6846531968814578*phiBC[14]+1.369306393762916*phiBC[13]+0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[5]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[3]] = 1.185854122563142*phiBC[18]-1.369306393762916*phiBC[15]-0.6846531968814578*phiBC[14]+1.369306393762916*phiBC[13]+0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[5]+0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong((-1.185854122563142*phiBC[18])+1.369306393762916*phiBC[15]-0.6846531968814578*phiBC[14]+1.369306393762916*phiBC[13]-0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[5]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[4]] = (-1.185854122563142*phiBC[18])+1.369306393762916*phiBC[15]-0.6846531968814578*phiBC[14]+1.369306393762916*phiBC[13]-0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[5]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong((-2.371708245126285*phiBC[19])-2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[5]] = (-2.371708245126285*phiBC[19])-2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[6]],__double_as_longlong((-1.185854122563142*phiBC[17])+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[14]-0.6846531968814578*phiBC[13]-0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]+1.060660171779821*phiBC[6]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[6]] = (-1.185854122563142*phiBC[17])+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[14]-0.6846531968814578*phiBC[13]-0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]+1.060660171779821*phiBC[6]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(2.371708245126285*phiBC[19]+2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[7]] = 2.371708245126285*phiBC[19]+2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]+1.369306393762916*phiBC[14]+1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]+0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[8]],(-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[8]] += (-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[9]],0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[9]] += 0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[10]],0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[10]] += 0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[11]],(-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[11]] += (-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[12]],0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[12]] += 0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[13]],0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[13]] += 0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[14]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[14]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[15]],0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[15]] += 0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[16]],(-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[16]] += (-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[17]],(-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[17]] += (-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[18]],(-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[18]] += (-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[19]],0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[19]] += 0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1_upz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2_upz_nondirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],(-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += (-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],(-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += (-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += 0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[8]],(-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[8]] += (-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[9]],0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[9]] += 0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[10]],0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[10]] += 0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[11]],(-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[11]] += (-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[12]],0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[12]] += 0.03513641844631534*rho[19]-0.03513641844631534*rho[18]-0.03513641844631534*rho[17]-0.06085806194501847*rho[16]-0.06085806194501847*rho[15]+0.06085806194501847*rho[14]+0.06085806194501847*rho[13]-0.06085806194501847*rho[12]-0.06085806194501847*rho[11]+0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.03928371006591933*rho[6]-0.03928371006591933*rho[5]+0.03928371006591933*rho[4]-0.06804138174397718*rho[3]+0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[13]],0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[13]] += 0.07027283689263068*rho[17]-0.1217161238900369*rho[13]+0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]+0.2721655269759087*rho[3]-0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[14]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[14]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]-0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]-0.06085806194501844*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]-0.06804138174397717*rho[3]+0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[15]],0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[15]] += 0.07027283689263068*rho[18]-0.121716123890037*rho[14]+0.121716123890037*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]+0.2721655269759088*rho[3]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[16]],(-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[16]] += (-0.07027283689263068*rho[18])-0.121716123890037*rho[14]-0.121716123890037*rho[12]-0.210818510677892*rho[8]+0.1571348402636773*rho[5]+0.2721655269759088*rho[3]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[17]],(-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[17]] += (-0.03513641844631533*rho[19])-0.03513641844631533*rho[18]+0.03513641844631533*rho[17]+0.06085806194501849*rho[16]-0.06085806194501849*rho[15]+0.06085806194501849*rho[14]+0.06085806194501849*rho[13]-0.06085806194501849*rho[12]+0.06085806194501849*rho[11]-0.06804138174397718*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]-0.03928371006591933*rho[5]-0.03928371006591933*rho[4]-0.06804138174397718*rho[3]-0.06804138174397718*rho[2]+0.06804138174397718*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[18]],(-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[18]] += (-0.07027283689263068*rho[17])-0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]+0.1571348402636773*rho[6]+0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[19]],0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[19]] += 0.03513641844631532*rho[19]+0.03513641844631532*rho[18]+0.03513641844631532*rho[17]+0.06085806194501844*rho[16]+0.06085806194501844*rho[15]+0.06085806194501844*rho[14]+0.06085806194501844*rho[13]+0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.0392837100659193*rho[6]+0.0392837100659193*rho[5]+0.0392837100659193*rho[4]-0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p1_upz_dirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong((-1.837117307087384*phiBC[7])+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[4]] = (-1.837117307087384*phiBC[7])+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(1.837117307087384*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[5]] = 1.837117307087384*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[6]],__double_as_longlong(1.837117307087384*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[6]] = 1.837117307087384*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong((-1.837117307087384*phiBC[7])-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[7]] = (-1.837117307087384*phiBC[7])-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif

}

GKYL_CU_DH void fem_parproj_src_stencil_3x_ser_p2_upz_dirichletz(const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // rho: right side source.
  // phiBC: Dirichlet boundary potential, given as a DG expansion in the ghost cell (volume).
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[0]],0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[0]] += 0.03513641844631534*rho[19]+0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501849*rho[16]-0.06085806194501849*rho[15]-0.06085806194501849*rho[14]-0.06085806194501849*rho[13]-0.06085806194501849*rho[12]-0.06085806194501849*rho[11]-0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591933*rho[6]+0.03928371006591933*rho[5]+0.03928371006591933*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]+0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[1]],(-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[1]] += (-0.07027283689263068*rho[17])+0.121716123890037*rho[13]+0.121716123890037*rho[11]-0.2108185106778921*rho[7]+0.1571348402636773*rho[6]-0.2721655269759088*rho[3]-0.2721655269759088*rho[2]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[2]],(-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[2]] += (-0.03513641844631534*rho[19])-0.03513641844631534*rho[18]+0.03513641844631534*rho[17]-0.06085806194501848*rho[16]+0.06085806194501848*rho[15]-0.06085806194501848*rho[14]-0.06085806194501848*rho[13]+0.06085806194501848*rho[12]-0.06085806194501848*rho[11]+0.0680413817439772*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]+0.03928371006591932*rho[6]-0.03928371006591932*rho[5]-0.03928371006591932*rho[4]+0.0680413817439772*rho[3]+0.0680413817439772*rho[2]-0.0680413817439772*rho[1]-0.3535533905932739*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[3]],(-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[3]] += (-0.07027283689263064*rho[18])+0.1217161238900369*rho[14]+0.1217161238900369*rho[12]-0.2108185106778919*rho[8]+0.1571348402636772*rho[5]-0.2721655269759087*rho[3]-0.2721655269759087*rho[1]+0.4714045207910316*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[4]],0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[4]] += 0.07027283689263068*rho[18]+0.1217161238900369*rho[14]-0.1217161238900369*rho[12]-0.210818510677892*rho[8]-0.1571348402636773*rho[5]-0.2721655269759087*rho[3]+0.2721655269759087*rho[1]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[5]],(-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[5]] += (-0.03513641844631532*rho[19])+0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501844*rho[16]-0.06085806194501844*rho[15]-0.06085806194501844*rho[14]-0.06085806194501844*rho[13]-0.06085806194501844*rho[12]+0.06085806194501844*rho[11]+0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]+0.0392837100659193*rho[5]-0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]+0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[6]],0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[6]] += 0.07027283689263068*rho[17]+0.1217161238900369*rho[13]-0.1217161238900369*rho[11]-0.210818510677892*rho[7]-0.1571348402636773*rho[6]-0.2721655269759087*rho[3]+0.2721655269759087*rho[2]+0.4714045207910318*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[7]],0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[7]] += 0.03513641844631532*rho[19]-0.03513641844631532*rho[18]-0.03513641844631532*rho[17]+0.06085806194501846*rho[16]+0.06085806194501846*rho[15]-0.06085806194501846*rho[14]-0.06085806194501846*rho[13]+0.06085806194501846*rho[12]+0.06085806194501846*rho[11]-0.06804138174397717*rho[10]+0.105409255338946*rho[9]+0.105409255338946*rho[8]+0.105409255338946*rho[7]-0.0392837100659193*rho[6]-0.0392837100659193*rho[5]+0.0392837100659193*rho[4]+0.06804138174397717*rho[3]-0.06804138174397717*rho[2]-0.06804138174397717*rho[1]-0.3535533905932737*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[8]],(-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[8]] += (-0.07027283689263064*rho[19])+0.1217161238900369*rho[16]+0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]-0.2721655269759087*rho[2]-0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[9]],0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[9]] += 0.07027283689263068*rho[19]+0.121716123890037*rho[16]-0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]-0.2721655269759088*rho[2]+0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[10]],0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[10]] += 0.07027283689263068*rho[19]-0.121716123890037*rho[16]+0.121716123890037*rho[15]-0.2108185106778921*rho[9]-0.1571348402636773*rho[4]+0.2721655269759088*rho[2]-0.2721655269759088*rho[1]+0.4714045207910319*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[nodeOff+globalIdxs[11]],(-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0]);
  #else
  bsrc[nodeOff+globalIdxs[11]] += (-0.07027283689263064*rho[19])-0.1217161238900369*rho[16]-0.1217161238900369*rho[15]-0.2108185106778919*rho[9]+0.1571348402636772*rho[4]+0.2721655269759087*rho[2]+0.2721655269759087*rho[1]+0.4714045207910317*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[12]],__double_as_longlong(2.371708245126285*phiBC[19]+2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[12]] = 2.371708245126285*phiBC[19]+2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[13]],__double_as_longlong((-1.185854122563142*phiBC[17])-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[14]+0.6846531968814578*phiBC[13]+0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]+1.060660171779821*phiBC[6]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[13]] = (-1.185854122563142*phiBC[17])-1.369306393762916*phiBC[16]-1.369306393762916*phiBC[14]+0.6846531968814578*phiBC[13]+0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]+1.060660171779821*phiBC[6]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[14]],__double_as_longlong((-2.371708245126285*phiBC[19])-2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[14]] = (-2.371708245126285*phiBC[19])-2.371708245126285*phiBC[18]+2.371708245126285*phiBC[17]-1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]-1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[15]],__double_as_longlong((-1.185854122563142*phiBC[18])-1.369306393762916*phiBC[15]+0.6846531968814578*phiBC[14]-1.369306393762916*phiBC[13]+0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[5]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[15]] = (-1.185854122563142*phiBC[18])-1.369306393762916*phiBC[15]+0.6846531968814578*phiBC[14]-1.369306393762916*phiBC[13]+0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]+1.060660171779821*phiBC[5]-0.6123724356957946*phiBC[3]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[16]],__double_as_longlong(1.185854122563142*phiBC[18]+1.369306393762916*phiBC[15]+0.6846531968814578*phiBC[14]-1.369306393762916*phiBC[13]-0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[5]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[16]] = 1.185854122563142*phiBC[18]+1.369306393762916*phiBC[15]+0.6846531968814578*phiBC[14]-1.369306393762916*phiBC[13]-0.6846531968814578*phiBC[12]+0.7905694150420949*phiBC[9]-0.3952847075210474*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[5]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[17]],__double_as_longlong((-2.371708245126285*phiBC[19])+2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[17]] = (-2.371708245126285*phiBC[19])+2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]-1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]+1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]+1.060660171779821*phiBC[5]-1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]-0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[18]],__double_as_longlong(1.185854122563142*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[14]+0.6846531968814578*phiBC[13]-0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]-1.060660171779821*phiBC[6]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[18]] = 1.185854122563142*phiBC[17]+1.369306393762916*phiBC[16]-1.369306393762916*phiBC[14]+0.6846531968814578*phiBC[13]-0.6846531968814578*phiBC[11]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]-0.3952847075210474*phiBC[7]-1.060660171779821*phiBC[6]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.3535533905932738*phiBC[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[19]],__double_as_longlong(2.371708245126285*phiBC[19]-2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0]));
  #else
  bsrc[nodeOff+globalIdxs[19]] = 2.371708245126285*phiBC[19]-2.371708245126285*phiBC[18]-2.371708245126285*phiBC[17]+1.369306393762916*phiBC[16]+1.369306393762916*phiBC[15]-1.369306393762916*phiBC[14]-1.369306393762916*phiBC[13]+1.369306393762916*phiBC[12]+1.369306393762916*phiBC[11]-1.837117307087384*phiBC[10]+0.7905694150420949*phiBC[9]+0.7905694150420949*phiBC[8]+0.7905694150420949*phiBC[7]-1.060660171779821*phiBC[6]-1.060660171779821*phiBC[5]+1.060660171779821*phiBC[4]-0.6123724356957946*phiBC[3]+0.6123724356957946*phiBC[2]+0.6123724356957946*phiBC[1]+0.3535533905932738*phiBC[0];
  #endif

}

