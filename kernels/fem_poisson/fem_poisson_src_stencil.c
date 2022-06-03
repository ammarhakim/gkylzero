#include <gkyl_fem_poisson_kernels.h> 
 
void fem_poisson_src_stencil_1x_ser_p1_inx_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_inx_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_lox_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_lox_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_lox_dirichletx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_lox_dirichletx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_lox_neumannx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_lox_neumannx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_lox_robinx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_lox_robinx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_upx_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_upx_periodicx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_upx_dirichletx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);

}
void fem_poisson_src_stencil_1x_ser_p2_upx_dirichletx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);

}
void fem_poisson_src_stencil_1x_ser_p1_upx_neumannx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_upx_neumannx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p1_upx_robinx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.7071067811865476*rho[0]-0.408248290463863*rho[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.408248290463863*rho[1]+0.7071067811865476*rho[0]);

}
void fem_poisson_src_stencil_1x_ser_p2_upx_robinx(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.210818510677892*rho[2]-0.408248290463863*rho[1]+0.2357022603955158*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.9428090415820636*rho[0]-0.4216370213557841*rho[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.210818510677892*rho[2]+0.408248290463863*rho[1]+0.2357022603955158*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_iny_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_inx_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_inx_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_dirichletx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_dirichletx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[2]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_neumannx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_neumannx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_lox_robinx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_lox_robinx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_loy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_loy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_insert(tri, globalIdxs[0], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[8]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[8]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_loy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_loy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_periodicx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_dirichletx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_dirichletx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[4], 0, bcVals[5]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[5]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_neumannx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_neumannx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_upy_periodicy(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_upy_dirichlety(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[6], 0, bcVals[11]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], 0, bcVals[11]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_upy_neumanny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p1_upx_robinx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, 0.1666666666666667*rho[3]-0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, (-0.1666666666666667*rho[3])-0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, (-0.1666666666666667*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1666666666666667*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);

}
void fem_poisson_src_stencil_2x_ser_p2_upx_robinx_upy_robiny(const double *rho, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // rho: right side source.
  // bcVals[3*2*dim]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  gkyl_mat_triples_accum(tri, globalIdxs[0], 0, (-0.08606629658238704*rho[7])-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]-0.09622504486493762*rho[2]-0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], 0, 0.1721325931647741*rho[6]-0.298142396999972*rho[4]-0.3849001794597506*rho[2]+0.6666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], 0, 0.08606629658238704*rho[7]-0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]-0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], 0, 0.1721325931647741*rho[7]-0.2981423969999719*rho[5]-0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], 0, (-0.1721325931647741*rho[7])-0.2981423969999719*rho[5]+0.3849001794597505*rho[1]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], 0, (-0.08606629658238707*rho[7])+0.08606629658238707*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]-0.1666666666666667*rho[3]+0.09622504486493766*rho[2]-0.09622504486493766*rho[1]-0.1666666666666667*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], 0, (-0.1721325931647741*rho[6])-0.2981423969999719*rho[4]+0.3849001794597505*rho[2]+0.6666666666666665*rho[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], 0, 0.08606629658238704*rho[7]+0.08606629658238704*rho[6]+0.149071198499986*rho[5]+0.149071198499986*rho[4]+0.1666666666666667*rho[3]+0.09622504486493762*rho[2]+0.09622504486493762*rho[1]-0.1666666666666667*rho[0]);

}
