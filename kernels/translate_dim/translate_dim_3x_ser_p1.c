#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[2]; 
  ftar[2] = 0.7071067811865475*fdo[3]; 
  ftar[3] = 0.7071067811865475*fdo[6]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[1]; 
  ftar[2] = 0.7071067811865475*fdo[3]; 
  ftar[3] = 0.7071067811865475*fdo[5]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[1]; 
  ftar[2] = 0.7071067811865475*fdo[2]; 
  ftar[3] = 0.7071067811865475*fdo[4]; 
}

