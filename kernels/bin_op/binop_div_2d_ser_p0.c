#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
void binop_div_2d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
  // fdivg:   output field. 
 
  double lhs[1]; 
  lhs[0] = g[0]; 
  gkyl_mat_set(rhs,1,0,f[0]); 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.5*lhs[0]); 

  // Solve the system of equations. 
  long ipiv[1]; 
  gkyl_mat_linsolve_lu(A,rhs,ipiv); 
  for(size_t i=0; i<1; i++) 
  { 
    fdivg[i] = gkyl_mat_get(rhs,i,0); 
  } 
} 
