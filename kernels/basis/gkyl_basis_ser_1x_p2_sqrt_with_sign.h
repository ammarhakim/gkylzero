GKYL_CU_DH static inline void 
ser_1x_p2_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[5] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = 1.156987065043442*A[2]-1.109839118871799*A[1]+0.7071067811865475*A[0]; 
  temp_sign = 1.156987065043442*ASign[2]-1.109839118871799*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = (-0.1028945116539821*A[2])-0.6594875259537025*A[1]+0.7071067811865475*A[0]; 
  temp_sign = (-0.1028945116539821*ASign[2])-0.6594875259537025*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.7071067811865475*A[0]-0.7905694150420947*A[2]; 
  temp_sign = 0.7071067811865475*ASign[0]-0.7905694150420947*ASign[2]; 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = (-0.1028945116539821*A[2])+0.6594875259537025*A[1]+0.7071067811865475*A[0]; 
  temp_sign = (-0.1028945116539821*ASign[2])+0.6594875259537025*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  temp = 1.156987065043442*A[2]+1.109839118871799*A[1]+0.7071067811865475*A[0]; 
  temp_sign = 1.156987065043442*ASign[2]+1.109839118871799*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[4] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[4] = -sqrt(temp); 
  } else { 
  AOrd[4] = sqrt(temp); 
  } 
  ASqrt[0] = 0.1675326070686369*AOrd[4]+0.3384415785804036*AOrd[3]+0.4022651910750141*AOrd[2]+0.3384415785804036*AOrd[1]+0.1675326070686369*AOrd[0]; 
  ASqrt[1] = 0.262950725347801*AOrd[4]+0.315649637758137*AOrd[3]-0.315649637758137*AOrd[1]-0.262950725347801*AOrd[0]; 
  ASqrt[2] = 0.2741213413710451*AOrd[4]-0.04924826331462695*AOrd[3]-0.4497461561128364*AOrd[2]-0.04924826331462695*AOrd[1]+0.2741213413710451*AOrd[0]; 

} 
 
