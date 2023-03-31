#include <gkyl_fem_poisson_kernels.h> 
 
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_dirichletx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_dirichletx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_neumannx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_neumannx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_robinx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.5*(epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_robinx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.1666666666666667*(7.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_periodicx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_dirichletx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_dirichletx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_neumannx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.5*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_neumannx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 1.166666666666667*epsilon[0]*rdx2Sq[0]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_robinx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -0.5*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.5*(epsilon[0]*rdx2Sq[0]*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_robinx(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 1.166666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 2.666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.1666666666666667*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -1.333333333333333*epsilon[0]*rdx2Sq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.1666666666666667*(7.0*epsilon[0]*rdx2Sq[0]*bcVals[4]+6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_iny_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]*bcVals[6]))/(bcVals[1]*bcVals[7]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*bcVals[1]*rdx2Sq[1]*bcVals[6]))/(bcVals[1]*bcVals[7]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_robinx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*(epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -(0.1666666666666667*(2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+2.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*bcVals[1]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]*bcVals[1]-2.0*bcVals[0]*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*bcVals[1]*rdx2Sq[1]*bcVals[9]))/(bcVals[1]*bcVals[10]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_robinx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*(26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.02222222222222222*(80.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-48.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], (0.01111111111111111*(28.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -(0.02222222222222222*(40.0*epsilon[0]*bcVals[1]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]*bcVals[1]+6.0*bcVals[0]*epsilon[0]*rdx2Sq[0]))/bcVals[1]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*bcVals[1]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]*bcVals[1]-12.0*bcVals[0]*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*bcVals[1]*rdx2Sq[1]*bcVals[9]))/(bcVals[1]*bcVals[10]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_loy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_loy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_loy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*(((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3])*bcVals[7]-2.0*epsilon[0]*rdx2Sq[1]*bcVals[4]*bcVals[6]))/(bcVals[4]*bcVals[7]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_loy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-48.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[7]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[7]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[6]))/bcVals[7]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*(((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3])*bcVals[7]-12.0*epsilon[0]*rdx2Sq[1]*bcVals[4]*bcVals[6]))/(bcVals[4]*bcVals[7]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_upy_periodicy(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_upy_dirichlety(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_upy_neumanny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_robinx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.3333333333333333*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1666666666666667*(epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], -0.1666666666666667*(2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.3333333333333333*((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.1666666666666667*(epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], -(0.1666666666666667*((2.0*epsilon[0]*rdx2Sq[1]-1.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-2.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], (0.1666666666666667*((epsilon[0]*rdx2Sq[1]-2.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], (0.3333333333333333*(((epsilon[0]*rdx2Sq[1]+epsilon[0]*rdx2Sq[0])*bcVals[4]+2.0*epsilon[0]*rdx2Sq[0]*bcVals[3])*bcVals[10]+2.0*epsilon[0]*rdx2Sq[1]*bcVals[4]*bcVals[9]))/(bcVals[4]*bcVals[10]));

}
void fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_robinx_upy_robiny(const double *epsilon, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity.
  // dx: cell length in each direction.
  // bcVals[3]: values to impose as BCs, i.e. bcVals[0]*phi+bcVals[1]*d(phi)/dx=bcVals[2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], 0.02222222222222222*(26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[5], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[7], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], 0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[5], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[6], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[7], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.01111111111111111*(17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], 0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[5], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[6], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[7], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.02222222222222222*(80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[4], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[5], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[7], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[0], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[2], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[3], 0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-24.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[4], (0.02222222222222222*((80.0*epsilon[0]*rdx2Sq[1]+24.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+48.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[5], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[4], globalIdxs[7], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[0], 0.01111111111111111*(28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[2], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[3], -0.02222222222222222*(40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[4], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[5], (0.02222222222222222*((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[5], globalIdxs[7], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[0], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[1], -0.02222222222222222*(24.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[2], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[5], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[6], (0.02222222222222222*((24.0*epsilon[0]*rdx2Sq[1]+80.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+48.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[6], globalIdxs[7], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[0], 0.01111111111111111*(23.0*epsilon[0]*rdx2Sq[1]+23.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[1], -0.02222222222222222*(3.0*epsilon[0]*rdx2Sq[1]+20.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[2], (0.01111111111111111*((28.0*epsilon[0]*rdx2Sq[1]+17.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[3], -0.02222222222222222*(20.0*epsilon[0]*rdx2Sq[1]+3.0*epsilon[0]*rdx2Sq[0]));
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[4], -(0.02222222222222222*((40.0*epsilon[0]*rdx2Sq[1]-3.0*epsilon[0]*rdx2Sq[0])*bcVals[4]-6.0*epsilon[0]*rdx2Sq[0]*bcVals[3]))/bcVals[4]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[5], (0.01111111111111111*((17.0*epsilon[0]*rdx2Sq[1]+28.0*epsilon[0]*rdx2Sq[0])*bcVals[10]-6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[6], (0.02222222222222222*((3.0*epsilon[0]*rdx2Sq[1]-40.0*epsilon[0]*rdx2Sq[0])*bcVals[10]+6.0*epsilon[0]*rdx2Sq[1]*bcVals[9]))/bcVals[10]);
  gkyl_mat_triples_accum(tri, globalIdxs[7], globalIdxs[7], (0.02222222222222222*(((26.0*epsilon[0]*rdx2Sq[1]+26.0*epsilon[0]*rdx2Sq[0])*bcVals[4]+12.0*epsilon[0]*rdx2Sq[0]*bcVals[3])*bcVals[10]+12.0*epsilon[0]*rdx2Sq[1]*bcVals[4]*bcVals[9]))/(bcVals[4]*bcVals[10]));

}
