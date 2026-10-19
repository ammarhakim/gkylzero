#include <gkyl_dg_fv_interpolate_kernels.h> 

void dg_fv_interpolate_fv_to_modal_2x(const double *fv_data[8], double *dg_out) { 
  // fv_data: List of finite volume data in adjacent cells. Same order as nodes. 
  // out: DG data on 2x coarsened mesh. 
 
  double data[4]; 
  data[0] = fv_data[0][0]; 
  data[1] = fv_data[1][0]; 
  data[2] = fv_data[2][0]; 
  data[3] = fv_data[3][0]; 
  dg_out[0] = 0.5*data[3]+0.5*data[2]+0.5*data[1]+0.5*data[0]; 
  dg_out[1] = 0.5*data[3]+0.5*data[2]-0.5*data[1]-0.5*data[0]; 
  dg_out[2] = 0.5*data[3]-0.5*data[2]+0.5*data[1]-0.5*data[0]; 
  dg_out[3] = 0.5*data[3]-0.5*data[2]-0.5*data[1]+0.5*data[0]; 
}
