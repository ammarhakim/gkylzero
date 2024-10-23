#include <gkyl_skin_surf_from_ghost_kernels.h> 

GKYL_CU_DH void skin_surf_from_ghost_lower_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 1.224744871391589*fghost[3]+0.7071067811865475*fghost[0]; 
  fghostSurf[1] = 1.224744871391589*fghost[5]+0.7071067811865475*fghost[1]; 
  fghostSurf[2] = 1.224744871391589*fghost[6]+0.7071067811865475*fghost[2]; 
  fghostSurf[3] = 1.224744871391589*fghost[7]+0.7071067811865475*fghost[4]; 

  fskin[0] = 0.8660254037844386*fskin[3]+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskin[1] = 0.8660254037844386*fskin[5]+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskin[2] = 0.8660254037844386*fskin[6]+0.5*fskin[2]+0.7071067811865475*fghostSurf[2]; 
  fskin[3] = 0.5*fskin[3]+0.2886751345948129*fskin[0]-0.408248290463863*fghostSurf[0]; 
  fskin[4] = 0.8660254037844386*fskin[7]+0.5*fskin[4]+0.7071067811865475*fghostSurf[3]; 
  fskin[5] = 0.5*fskin[5]+0.2886751345948129*fskin[1]-0.408248290463863*fghostSurf[1]; 
  fskin[6] = 0.5*fskin[6]+0.2886751345948129*fskin[2]-0.408248290463863*fghostSurf[2]; 
  fskin[7] = 0.5*fskin[7]+0.2886751345948129*fskin[4]-0.408248290463863*fghostSurf[3]; 

}

GKYL_CU_DH void skin_surf_from_ghost_upper_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 0.7071067811865475*fghost[0]-1.224744871391589*fghost[3]; 
  fghostSurf[1] = 0.7071067811865475*fghost[1]-1.224744871391589*fghost[5]; 
  fghostSurf[2] = 0.7071067811865475*fghost[2]-1.224744871391589*fghost[6]; 
  fghostSurf[3] = 0.7071067811865475*fghost[4]-1.224744871391589*fghost[7]; 

  fskin[0] = -(0.8660254037844386*fskin[3])+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskin[1] = -(0.8660254037844386*fskin[5])+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskin[2] = -(0.8660254037844386*fskin[6])+0.5*fskin[2]+0.7071067811865475*fghostSurf[2]; 
  fskin[3] = 0.5*fskin[3]-0.2886751345948129*fskin[0]+0.408248290463863*fghostSurf[0]; 
  fskin[4] = -(0.8660254037844386*fskin[7])+0.5*fskin[4]+0.7071067811865475*fghostSurf[3]; 
  fskin[5] = 0.5*fskin[5]-0.2886751345948129*fskin[1]+0.408248290463863*fghostSurf[1]; 
  fskin[6] = 0.5*fskin[6]-0.2886751345948129*fskin[2]+0.408248290463863*fghostSurf[2]; 
  fskin[7] = 0.5*fskin[7]-0.2886751345948129*fskin[4]+0.408248290463863*fghostSurf[3]; 

}

