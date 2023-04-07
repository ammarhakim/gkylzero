GKYL_CU_DH static inline void 
ser_2x_p1_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[9] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = 0.9*A[3]-0.6708203932499369*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = 0.9*ASign[3]-0.6708203932499369*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.5*A[0]-0.6708203932499369*A[1]; 
  temp_sign = 0.5*ASign[0]-0.6708203932499369*ASign[1]; 
  if (temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = (-0.9*A[3])+0.6708203932499369*A[2]-0.6708203932499369*A[1]+0.5*A[0]; 
  temp_sign = (-0.9*ASign[3])+0.6708203932499369*ASign[2]-0.6708203932499369*ASign[1]+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = 0.5*A[0]-0.6708203932499369*A[2]; 
  temp_sign = 0.5*ASign[0]-0.6708203932499369*ASign[2]; 
  if (temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  temp = 0.5*A[0]; 
  temp_sign = 0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[4] = -sqrt(temp); 
  } else { 
  AOrd[4] = sqrt(temp); 
  } 
  temp = 0.6708203932499369*A[2]+0.5*A[0]; 
  temp_sign = 0.6708203932499369*ASign[2]+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[5] = -sqrt(temp); 
  } else { 
  AOrd[5] = sqrt(temp); 
  } 
  temp = (-0.9*A[3])-0.6708203932499369*A[2]+0.6708203932499369*A[1]+0.5*A[0]; 
  temp_sign = (-0.9*ASign[3])-0.6708203932499369*ASign[2]+0.6708203932499369*ASign[1]+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[6] = -sqrt(temp); 
  } else { 
  AOrd[6] = sqrt(temp); 
  } 
  temp = 0.6708203932499369*A[1]+0.5*A[0]; 
  temp_sign = 0.6708203932499369*ASign[1]+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[7] = -sqrt(temp); 
  } else { 
  AOrd[7] = sqrt(temp); 
  } 
  temp = 0.9*A[3]+0.6708203932499369*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = 0.9*ASign[3]+0.6708203932499369*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[8] = -sqrt(temp); 
  } else { 
  AOrd[8] = sqrt(temp); 
  } 
  ASqrt[0] = 0.154320987654321*AOrd[8]+0.2469135802469136*AOrd[7]+0.154320987654321*AOrd[6]+0.2469135802469136*AOrd[5]+0.3950617283950617*AOrd[4]+0.2469135802469136*AOrd[3]+0.154320987654321*AOrd[2]+0.2469135802469136*AOrd[1]+0.154320987654321*AOrd[0]; 
  ASqrt[1] = 0.2070433312499806*AOrd[8]+0.3312693299999688*AOrd[7]+0.2070433312499806*AOrd[6]-0.2070433312499806*AOrd[2]-0.3312693299999688*AOrd[1]-0.2070433312499806*AOrd[0]; 
  ASqrt[2] = 0.2070433312499806*AOrd[8]-0.2070433312499806*AOrd[6]+0.3312693299999688*AOrd[5]-0.3312693299999688*AOrd[3]+0.2070433312499806*AOrd[2]-0.2070433312499806*AOrd[0]; 
  ASqrt[3] = 0.2777777777777778*AOrd[8]-0.2777777777777778*AOrd[6]-0.2777777777777778*AOrd[2]+0.2777777777777778*AOrd[0]; 

} 
 
