GKYL_CU_DH static inline void 
ser_2x_p1_sqrt(const double *A, double *ASqrt) 
{ 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A)). 
 
  double AOrd[4] = {0.0}; 

  AOrd[0] = sqrt(0.5*A[3]-0.5*(A[2]+A[1])+0.5*A[0]); 
  AOrd[1] = sqrt((-0.5*A[3])+0.5*A[2]-0.5*A[1]+0.5*A[0]); 
  AOrd[2] = sqrt(0.5*(A[1]+A[0])-0.5*(A[3]+A[2])); 
  AOrd[3] = sqrt(0.5*(A[3]+A[2]+A[1]+A[0])); 
  ASqrt[0] = 0.5*AOrd[3]+0.5*AOrd[2]+0.5*AOrd[1]+0.5*AOrd[0]; 
  ASqrt[1] = 0.5*AOrd[3]+0.5*AOrd[2]-0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[2] = 0.5*AOrd[3]-0.5*AOrd[2]+0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[3] = 0.5*AOrd[3]-0.5*AOrd[2]-0.5*AOrd[1]+0.5*AOrd[0]; 

} 
 
