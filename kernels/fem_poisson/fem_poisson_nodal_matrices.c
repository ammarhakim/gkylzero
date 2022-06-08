#include <gkyl_fem_poisson_kernels.h> 
 
GKYL_CU_DH void fem_poisson_stiff_1x_ser_p1(const double *dx, struct gkyl_mat *matout) 
{ 
  // dx: cell length in each direction.
  // matout: local stiffness matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_set(matout,0,0,0.5*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,1,-0.5*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,0,-0.5*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,1,0.5*rdx2Sq[0]); 
}

GKYL_CU_DH void fem_poisson_stiff_1x_ser_p2(const double *dx, struct gkyl_mat *matout) 
{ 
  // dx: cell length in each direction.
  // matout: local stiffness matrix.

  double rdx2Sq[1];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);

  gkyl_mat_set(matout,0,0,1.166666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,1,-1.333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,2,0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,0,-1.333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,1,2.666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,2,-1.333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,0,0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,1,-1.333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,2,1.166666666666667*rdx2Sq[0]); 
}

GKYL_CU_DH void fem_poisson_mass_times_modtonod_1x_ser_p1(struct gkyl_mat *matout) 
{ 
  // matout: mass matrix times modal to nodal matrix.

  gkyl_mat_set(matout,0,0,0.7071067811865475); 
  gkyl_mat_set(matout,0,1,-0.4082482904638631); 
  gkyl_mat_set(matout,1,0,0.7071067811865475); 
  gkyl_mat_set(matout,1,1,0.4082482904638631); 
}

GKYL_CU_DH void fem_poisson_mass_times_modtonod_1x_ser_p2(struct gkyl_mat *matout) 
{ 
  // matout: mass matrix times modal to nodal matrix.

  gkyl_mat_set(matout,0,0,0.2357022603955159); 
  gkyl_mat_set(matout,0,1,-0.4082482904638631); 
  gkyl_mat_set(matout,0,2,0.2108185106778921); 
  gkyl_mat_set(matout,1,0,0.9428090415820638); 
  gkyl_mat_set(matout,1,1,0.0); 
  gkyl_mat_set(matout,1,2,-0.4216370213557841); 
  gkyl_mat_set(matout,2,0,0.2357022603955159); 
  gkyl_mat_set(matout,2,1,0.4082482904638631); 
  gkyl_mat_set(matout,2,2,0.2108185106778921); 
}

GKYL_CU_DH void fem_poisson_nodtomod_1x_ser_p1(struct gkyl_mat *matout) 
{ 
  // matout: nodal to modal matrix.

  gkyl_mat_set(matout,0,0,0.7071067811865475); 
  gkyl_mat_set(matout,0,1,0.7071067811865475); 
  gkyl_mat_set(matout,1,0,-0.408248290463863); 
  gkyl_mat_set(matout,1,1,0.408248290463863); 
}

GKYL_CU_DH void fem_poisson_nodtomod_1x_ser_p2(struct gkyl_mat *matout) 
{ 
  // matout: nodal to modal matrix.

  gkyl_mat_set(matout,0,0,0.2357022603955159); 
  gkyl_mat_set(matout,0,1,0.9428090415820636); 
  gkyl_mat_set(matout,0,2,0.2357022603955159); 
  gkyl_mat_set(matout,1,0,-0.4082482904638631); 
  gkyl_mat_set(matout,1,1,0.0); 
  gkyl_mat_set(matout,1,2,0.4082482904638631); 
  gkyl_mat_set(matout,2,0,0.210818510677892); 
  gkyl_mat_set(matout,2,1,-0.421637021355784); 
  gkyl_mat_set(matout,2,2,0.210818510677892); 
}

GKYL_CU_DH void fem_poisson_stiff_2x_ser_p1(const double *dx, struct gkyl_mat *matout) 
{ 
  // dx: cell length in each direction.
  // matout: local stiffness matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_set(matout,0,0,0.3333333333333333*rdx2Sq[1]+0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,1,0.1666666666666667*rdx2Sq[1]-0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,2,0.1666666666666667*rdx2Sq[0]-0.3333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,0,3,(-0.1666666666666667*rdx2Sq[1])-0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,0,0.1666666666666667*rdx2Sq[1]-0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,1,0.3333333333333333*rdx2Sq[1]+0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,2,(-0.1666666666666667*rdx2Sq[1])-0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,3,0.1666666666666667*rdx2Sq[0]-0.3333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,2,0,0.1666666666666667*rdx2Sq[0]-0.3333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,2,1,(-0.1666666666666667*rdx2Sq[1])-0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,2,0.3333333333333333*rdx2Sq[1]+0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,3,0.1666666666666667*rdx2Sq[1]-0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,0,(-0.1666666666666667*rdx2Sq[1])-0.1666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,1,0.1666666666666667*rdx2Sq[0]-0.3333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,3,2,0.1666666666666667*rdx2Sq[1]-0.3333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,3,0.3333333333333333*rdx2Sq[1]+0.3333333333333333*rdx2Sq[0]); 
}

GKYL_CU_DH void fem_poisson_stiff_2x_ser_p2(const double *dx, struct gkyl_mat *matout) 
{ 
  // dx: cell length in each direction.
  // matout: local stiffness matrix.

  double rdx2Sq[2];
  rdx2Sq[0] = 4.0/(dx[0]*dx[0]);
  rdx2Sq[1] = 4.0/(dx[1]*dx[1]);

  gkyl_mat_set(matout,0,0,0.5777777777777777*rdx2Sq[1]+0.5777777777777777*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,1,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,2,0.1888888888888889*rdx2Sq[1]+0.3111111111111111*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,3,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,0,4,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,5,0.3111111111111111*rdx2Sq[1]+0.1888888888888889*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,6,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,0,7,0.2555555555555555*rdx2Sq[1]+0.2555555555555555*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,0,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,1,0.5333333333333333*rdx2Sq[1]+1.777777777777778*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,2,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,3,0.0); 
  gkyl_mat_set(matout,1,4,0.0); 
  gkyl_mat_set(matout,1,5,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,1,6,0.8888888888888888*rdx2Sq[0]-0.5333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,1,7,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,0,0.1888888888888889*rdx2Sq[1]+0.3111111111111111*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,1,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,2,0.5777777777777777*rdx2Sq[1]+0.5777777777777777*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,3,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,4,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,2,5,0.2555555555555555*rdx2Sq[1]+0.2555555555555555*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,6,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,2,7,0.3111111111111111*rdx2Sq[1]+0.1888888888888889*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,0,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,3,1,0.0); 
  gkyl_mat_set(matout,3,2,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,3,1.777777777777778*rdx2Sq[1]+0.5333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,4,0.8888888888888888*rdx2Sq[1]-0.5333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,3,5,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,3,6,0.0); 
  gkyl_mat_set(matout,3,7,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,4,0,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,4,1,0.0); 
  gkyl_mat_set(matout,4,2,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,4,3,0.8888888888888888*rdx2Sq[1]-0.5333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,4,4,1.777777777777778*rdx2Sq[1]+0.5333333333333333*rdx2Sq[0]); 
  gkyl_mat_set(matout,4,5,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,4,6,0.0); 
  gkyl_mat_set(matout,4,7,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,5,0,0.3111111111111111*rdx2Sq[1]+0.1888888888888889*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,1,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,2,0.2555555555555555*rdx2Sq[1]+0.2555555555555555*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,3,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,5,4,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,5,0.5777777777777777*rdx2Sq[1]+0.5777777777777777*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,6,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,5,7,0.1888888888888889*rdx2Sq[1]+0.3111111111111111*rdx2Sq[0]); 
  gkyl_mat_set(matout,6,0,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,6,1,0.8888888888888888*rdx2Sq[0]-0.5333333333333333*rdx2Sq[1]); 
  gkyl_mat_set(matout,6,2,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,6,3,0.0); 
  gkyl_mat_set(matout,6,4,0.0); 
  gkyl_mat_set(matout,6,5,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,6,6,0.5333333333333333*rdx2Sq[1]+1.777777777777778*rdx2Sq[0]); 
  gkyl_mat_set(matout,6,7,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,0,0.2555555555555555*rdx2Sq[1]+0.2555555555555555*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,1,(-0.06666666666666667*rdx2Sq[1])-0.4444444444444444*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,2,0.3111111111111111*rdx2Sq[1]+0.1888888888888889*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,3,(-0.4444444444444444*rdx2Sq[1])-0.06666666666666667*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,4,0.06666666666666667*rdx2Sq[0]-0.8888888888888888*rdx2Sq[1]); 
  gkyl_mat_set(matout,7,5,0.1888888888888889*rdx2Sq[1]+0.3111111111111111*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,6,0.06666666666666667*rdx2Sq[1]-0.8888888888888888*rdx2Sq[0]); 
  gkyl_mat_set(matout,7,7,0.5777777777777777*rdx2Sq[1]+0.5777777777777777*rdx2Sq[0]); 
}

GKYL_CU_DH void fem_poisson_mass_times_modtonod_2x_ser_p1(struct gkyl_mat *matout) 
{ 
  // matout: mass matrix times modal to nodal matrix.

  gkyl_mat_set(matout,0,0,0.5); 
  gkyl_mat_set(matout,0,1,-0.2886751345948129); 
  gkyl_mat_set(matout,0,2,-0.2886751345948129); 
  gkyl_mat_set(matout,0,3,0.1666666666666667); 
  gkyl_mat_set(matout,1,0,0.5); 
  gkyl_mat_set(matout,1,1,0.2886751345948129); 
  gkyl_mat_set(matout,1,2,-0.2886751345948129); 
  gkyl_mat_set(matout,1,3,-0.1666666666666667); 
  gkyl_mat_set(matout,2,0,0.5); 
  gkyl_mat_set(matout,2,1,-0.2886751345948129); 
  gkyl_mat_set(matout,2,2,0.2886751345948129); 
  gkyl_mat_set(matout,2,3,-0.1666666666666667); 
  gkyl_mat_set(matout,3,0,0.5); 
  gkyl_mat_set(matout,3,1,0.2886751345948129); 
  gkyl_mat_set(matout,3,2,0.2886751345948129); 
  gkyl_mat_set(matout,3,3,0.1666666666666667); 
}

GKYL_CU_DH void fem_poisson_mass_times_modtonod_2x_ser_p2(struct gkyl_mat *matout) 
{ 
  // matout: mass matrix times modal to nodal matrix.

  gkyl_mat_set(matout,0,0,-0.1666666666666667); 
  gkyl_mat_set(matout,0,1,-0.09622504486493766); 
  gkyl_mat_set(matout,0,2,-0.09622504486493766); 
  gkyl_mat_set(matout,0,3,0.1666666666666667); 
  gkyl_mat_set(matout,0,4,0.149071198499986); 
  gkyl_mat_set(matout,0,5,0.149071198499986); 
  gkyl_mat_set(matout,0,6,-0.08606629658238703); 
  gkyl_mat_set(matout,0,7,-0.08606629658238703); 
  gkyl_mat_set(matout,1,0,0.6666666666666666); 
  gkyl_mat_set(matout,1,1,0.0); 
  gkyl_mat_set(matout,1,2,-0.3849001794597506); 
  gkyl_mat_set(matout,1,3,0.0); 
  gkyl_mat_set(matout,1,4,-0.298142396999972); 
  gkyl_mat_set(matout,1,5,0.0); 
  gkyl_mat_set(matout,1,6,0.1721325931647741); 
  gkyl_mat_set(matout,1,7,0.0); 
  gkyl_mat_set(matout,2,0,-0.1666666666666667); 
  gkyl_mat_set(matout,2,1,0.09622504486493766); 
  gkyl_mat_set(matout,2,2,-0.09622504486493766); 
  gkyl_mat_set(matout,2,3,-0.1666666666666667); 
  gkyl_mat_set(matout,2,4,0.149071198499986); 
  gkyl_mat_set(matout,2,5,0.149071198499986); 
  gkyl_mat_set(matout,2,6,-0.08606629658238703); 
  gkyl_mat_set(matout,2,7,0.08606629658238703); 
  gkyl_mat_set(matout,3,0,0.6666666666666666); 
  gkyl_mat_set(matout,3,1,-0.3849001794597505); 
  gkyl_mat_set(matout,3,2,5.551115123125783e-17); 
  gkyl_mat_set(matout,3,3,0.0); 
  gkyl_mat_set(matout,3,4,0.0); 
  gkyl_mat_set(matout,3,5,-0.298142396999972); 
  gkyl_mat_set(matout,3,6,0.0); 
  gkyl_mat_set(matout,3,7,0.1721325931647741); 
  gkyl_mat_set(matout,4,0,0.6666666666666666); 
  gkyl_mat_set(matout,4,1,0.3849001794597505); 
  gkyl_mat_set(matout,4,2,0.0); 
  gkyl_mat_set(matout,4,3,0.0); 
  gkyl_mat_set(matout,4,4,0.0); 
  gkyl_mat_set(matout,4,5,-0.298142396999972); 
  gkyl_mat_set(matout,4,6,0.0); 
  gkyl_mat_set(matout,4,7,-0.1721325931647741); 
  gkyl_mat_set(matout,5,0,-0.1666666666666667); 
  gkyl_mat_set(matout,5,1,-0.09622504486493766); 
  gkyl_mat_set(matout,5,2,0.09622504486493762); 
  gkyl_mat_set(matout,5,3,-0.1666666666666667); 
  gkyl_mat_set(matout,5,4,0.149071198499986); 
  gkyl_mat_set(matout,5,5,0.149071198499986); 
  gkyl_mat_set(matout,5,6,0.08606629658238703); 
  gkyl_mat_set(matout,5,7,-0.08606629658238703); 
  gkyl_mat_set(matout,6,0,0.6666666666666666); 
  gkyl_mat_set(matout,6,1,0.0); 
  gkyl_mat_set(matout,6,2,0.3849001794597506); 
  gkyl_mat_set(matout,6,3,0.0); 
  gkyl_mat_set(matout,6,4,-0.298142396999972); 
  gkyl_mat_set(matout,6,5,0.0); 
  gkyl_mat_set(matout,6,6,-0.1721325931647741); 
  gkyl_mat_set(matout,6,7,0.0); 
  gkyl_mat_set(matout,7,0,-0.1666666666666667); 
  gkyl_mat_set(matout,7,1,0.09622504486493766); 
  gkyl_mat_set(matout,7,2,0.09622504486493766); 
  gkyl_mat_set(matout,7,3,0.1666666666666667); 
  gkyl_mat_set(matout,7,4,0.149071198499986); 
  gkyl_mat_set(matout,7,5,0.149071198499986); 
  gkyl_mat_set(matout,7,6,0.08606629658238703); 
  gkyl_mat_set(matout,7,7,0.08606629658238703); 
}

GKYL_CU_DH void fem_poisson_nodtomod_2x_ser_p1(struct gkyl_mat *matout) 
{ 
  // matout: nodal to modal matrix.

  gkyl_mat_set(matout,0,0,0.5); 
  gkyl_mat_set(matout,0,1,0.5); 
  gkyl_mat_set(matout,0,2,0.5); 
  gkyl_mat_set(matout,0,3,0.5); 
  gkyl_mat_set(matout,1,0,-0.2886751345948129); 
  gkyl_mat_set(matout,1,1,0.2886751345948129); 
  gkyl_mat_set(matout,1,2,-0.2886751345948129); 
  gkyl_mat_set(matout,1,3,0.2886751345948129); 
  gkyl_mat_set(matout,2,0,-0.2886751345948129); 
  gkyl_mat_set(matout,2,1,-0.2886751345948129); 
  gkyl_mat_set(matout,2,2,0.2886751345948129); 
  gkyl_mat_set(matout,2,3,0.2886751345948129); 
  gkyl_mat_set(matout,3,0,0.1666666666666667); 
  gkyl_mat_set(matout,3,1,-0.1666666666666667); 
  gkyl_mat_set(matout,3,2,-0.1666666666666667); 
  gkyl_mat_set(matout,3,3,0.1666666666666667); 
}

GKYL_CU_DH void fem_poisson_nodtomod_2x_ser_p2(struct gkyl_mat *matout) 
{ 
  // matout: nodal to modal matrix.

  gkyl_mat_set(matout,0,0,-0.1666666666666667); 
  gkyl_mat_set(matout,0,1,0.6666666666666666); 
  gkyl_mat_set(matout,0,2,-0.1666666666666667); 
  gkyl_mat_set(matout,0,3,0.6666666666666666); 
  gkyl_mat_set(matout,0,4,0.6666666666666666); 
  gkyl_mat_set(matout,0,5,-0.1666666666666667); 
  gkyl_mat_set(matout,0,6,0.6666666666666666); 
  gkyl_mat_set(matout,0,7,-0.1666666666666667); 
  gkyl_mat_set(matout,1,0,-0.09622504486493762); 
  gkyl_mat_set(matout,1,1,0.0); 
  gkyl_mat_set(matout,1,2,0.09622504486493762); 
  gkyl_mat_set(matout,1,3,-0.3849001794597504); 
  gkyl_mat_set(matout,1,4,0.3849001794597505); 
  gkyl_mat_set(matout,1,5,-0.09622504486493762); 
  gkyl_mat_set(matout,1,6,-4.158640611673197e-18); 
  gkyl_mat_set(matout,1,7,0.09622504486493762); 
  gkyl_mat_set(matout,2,0,-0.09622504486493762); 
  gkyl_mat_set(matout,2,1,-0.3849001794597505); 
  gkyl_mat_set(matout,2,2,-0.09622504486493762); 
  gkyl_mat_set(matout,2,3,4.782436703424177e-17); 
  gkyl_mat_set(matout,2,4,3.326912489338558e-17); 
  gkyl_mat_set(matout,2,5,0.09622504486493759); 
  gkyl_mat_set(matout,2,6,0.3849001794597505); 
  gkyl_mat_set(matout,2,7,0.0962250448649376); 
  gkyl_mat_set(matout,3,0,0.1666666666666667); 
  gkyl_mat_set(matout,3,1,0.0); 
  gkyl_mat_set(matout,3,2,-0.1666666666666667); 
  gkyl_mat_set(matout,3,3,8.317281223346394e-18); 
  gkyl_mat_set(matout,3,4,2.495184367003918e-17); 
  gkyl_mat_set(matout,3,5,-0.1666666666666667); 
  gkyl_mat_set(matout,3,6,0.0); 
  gkyl_mat_set(matout,3,7,0.1666666666666667); 
  gkyl_mat_set(matout,4,0,0.149071198499986); 
  gkyl_mat_set(matout,4,1,-0.2981423969999719); 
  gkyl_mat_set(matout,4,2,0.149071198499986); 
  gkyl_mat_set(matout,4,3,1.663456244669279e-17); 
  gkyl_mat_set(matout,4,4,-1.143626168210129e-17); 
  gkyl_mat_set(matout,4,5,0.149071198499986); 
  gkyl_mat_set(matout,4,6,-0.298142396999972); 
  gkyl_mat_set(matout,4,7,0.149071198499986); 
  gkyl_mat_set(matout,5,0,0.149071198499986); 
  gkyl_mat_set(matout,5,1,2.911048428171238e-17); 
  gkyl_mat_set(matout,5,2,0.149071198499986); 
  gkyl_mat_set(matout,5,3,-0.2981423969999719); 
  gkyl_mat_set(matout,5,4,-0.298142396999972); 
  gkyl_mat_set(matout,5,5,0.149071198499986); 
  gkyl_mat_set(matout,5,6,-3.425906660311779e-17); 
  gkyl_mat_set(matout,5,7,0.149071198499986); 
  gkyl_mat_set(matout,6,0,-0.08606629658238704); 
  gkyl_mat_set(matout,6,1,0.1721325931647741); 
  gkyl_mat_set(matout,6,2,-0.08606629658238706); 
  gkyl_mat_set(matout,6,3,-2.079320305836599e-17); 
  gkyl_mat_set(matout,6,4,9.356941376264694e-18); 
  gkyl_mat_set(matout,6,5,0.08606629658238706); 
  gkyl_mat_set(matout,6,6,-0.1721325931647741); 
  gkyl_mat_set(matout,6,7,0.08606629658238706); 
  gkyl_mat_set(matout,7,0,-0.08606629658238704); 
  gkyl_mat_set(matout,7,1,-2.079320305836599e-17); 
  gkyl_mat_set(matout,7,2,0.08606629658238706); 
  gkyl_mat_set(matout,7,3,0.1721325931647741); 
  gkyl_mat_set(matout,7,4,-0.1721325931647741); 
  gkyl_mat_set(matout,7,5,-0.08606629658238704); 
  gkyl_mat_set(matout,7,6,0.0); 
  gkyl_mat_set(matout,7,7,0.08606629658238706); 
}

