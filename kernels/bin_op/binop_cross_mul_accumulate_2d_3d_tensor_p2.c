#include <gkyl_binop_cross_mul_tensor.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_2d_3d_tensor_p2(double a, const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  fg[0] += a*(0.5*f[8]*g[20]+0.5*f[7]*g[12]+0.5*f[6]*g[11]+0.5*f[5]*g[8]+0.5*f[4]*g[7]+0.5*f[3]*g[4]+0.5*f[2]*g[2]+0.5*f[1]*g[1]+0.5*f[0]*g[0]); 
  fg[1] += a*(0.447213595499958*f[7]*g[20]+0.447213595499958*f[8]*g[12]+0.5000000000000001*f[5]*g[12]+0.447213595499958*f[3]*g[11]+0.5000000000000001*f[7]*g[8]+0.4472135954999579*f[1]*g[7]+0.447213595499958*g[4]*f[6]+0.5*f[2]*g[4]+0.4472135954999579*g[1]*f[4]+0.5*g[2]*f[3]+0.5*f[0]*g[1]+0.5*g[0]*f[1]); 
  fg[2] += a*(0.447213595499958*f[6]*g[20]+0.447213595499958*f[3]*g[12]+0.447213595499958*f[8]*g[11]+0.5000000000000001*f[4]*g[11]+0.4472135954999579*f[2]*g[8]+0.5000000000000001*f[6]*g[7]+0.447213595499958*g[4]*f[7]+0.4472135954999579*g[2]*f[5]+0.5*f[1]*g[4]+0.5*g[1]*f[3]+0.5*f[0]*g[2]+0.5*g[0]*f[2]); 
  fg[3] += a*(0.5*f[8]*g[23]+0.5000000000000001*f[7]*g[18]+0.5000000000000001*f[6]*g[17]+0.5000000000000001*f[5]*g[14]+0.5000000000000001*f[4]*g[13]+0.5*f[3]*g[10]+0.5*f[2]*g[6]+0.5*f[1]*g[5]+0.5*f[0]*g[3]); 
  fg[4] += a*(0.4*f[3]*g[20]+0.4*f[6]*g[12]+0.447213595499958*f[2]*g[12]+0.4*f[7]*g[11]+0.447213595499958*f[1]*g[11]+0.4472135954999579*f[3]*g[8]+0.4*g[4]*f[8]+0.4472135954999579*f[3]*g[7]+0.447213595499958*g[2]*f[7]+0.447213595499958*g[1]*f[6]+0.4472135954999579*g[4]*f[5]+0.4472135954999579*f[4]*g[4]+0.5*f[0]*g[4]+0.5*g[0]*f[3]+0.5*f[1]*g[2]+0.5*g[1]*f[2]); 
  fg[5] += a*(0.447213595499958*f[7]*g[23]+0.4472135954999579*f[8]*g[18]+0.5*f[5]*g[18]+0.4472135954999579*f[3]*g[17]+0.5*f[7]*g[14]+0.447213595499958*f[1]*g[13]+0.447213595499958*f[6]*g[10]+0.5*f[2]*g[10]+0.5*f[3]*g[6]+0.4472135954999579*f[4]*g[5]+0.5*f[0]*g[5]+0.5*f[1]*g[3]); 
  fg[6] += a*(0.447213595499958*f[6]*g[23]+0.4472135954999579*f[3]*g[18]+0.4472135954999579*f[8]*g[17]+0.5*f[4]*g[17]+0.447213595499958*f[2]*g[14]+0.5*f[6]*g[13]+0.447213595499958*f[7]*g[10]+0.5*f[1]*g[10]+0.4472135954999579*f[5]*g[6]+0.5*f[0]*g[6]+0.5*f[3]*g[5]+0.5*f[2]*g[3]); 
  fg[7] += a*(0.31943828249997*f[8]*g[20]+0.5*f[5]*g[20]+0.4472135954999579*f[7]*g[12]+0.31943828249997*f[6]*g[11]+0.5000000000000001*f[2]*g[11]+0.5*f[8]*g[8]+0.31943828249997*f[4]*g[7]+0.5*f[0]*g[7]+0.5000000000000001*g[2]*f[6]+0.4472135954999579*f[3]*g[4]+0.5*g[0]*f[4]+0.4472135954999579*f[1]*g[1]); 
  fg[8] += a*(0.31943828249997*f[8]*g[20]+0.5*f[4]*g[20]+0.31943828249997*f[7]*g[12]+0.5000000000000001*f[1]*g[12]+0.4472135954999579*f[6]*g[11]+0.31943828249997*f[5]*g[8]+0.5*f[0]*g[8]+0.5*g[7]*f[8]+0.5000000000000001*g[1]*f[7]+0.5*g[0]*f[5]+0.4472135954999579*f[3]*g[4]+0.4472135954999579*f[2]*g[2]); 
  fg[9] += a*(0.5*f[8]*g[26]+0.5000000000000001*f[7]*g[25]+0.5000000000000001*f[6]*g[24]+0.5*f[5]*g[22]+0.5*f[4]*g[21]+0.5*f[3]*g[19]+0.5000000000000001*f[2]*g[16]+0.5000000000000001*f[1]*g[15]+0.5*f[0]*g[9]); 
  fg[10] += a*(0.4*f[3]*g[23]+0.4*f[6]*g[18]+0.4472135954999579*f[2]*g[18]+0.4*f[7]*g[17]+0.4472135954999579*f[1]*g[17]+0.447213595499958*f[3]*g[14]+0.447213595499958*f[3]*g[13]+0.4*f[8]*g[10]+0.4472135954999579*f[5]*g[10]+0.4472135954999579*f[4]*g[10]+0.5*f[0]*g[10]+0.447213595499958*g[6]*f[7]+0.5*f[1]*g[6]+0.447213595499958*g[5]*f[6]+0.5*f[2]*g[5]+0.5*f[3]*g[3]); 
  fg[11] += a*(0.2857142857142857*f[6]*g[20]+0.447213595499958*f[2]*g[20]+0.4*f[3]*g[12]+0.2857142857142857*f[8]*g[11]+0.4472135954999579*f[5]*g[11]+0.31943828249997*f[4]*g[11]+0.5*f[0]*g[11]+0.4472135954999579*f[6]*g[8]+0.447213595499958*g[2]*f[8]+0.31943828249997*f[6]*g[7]+0.5000000000000001*f[2]*g[7]+0.4*g[4]*f[7]+0.5*g[0]*f[6]+0.447213595499958*f[1]*g[4]+0.5000000000000001*g[2]*f[4]+0.447213595499958*g[1]*f[3]); 
  fg[12] += a*(0.2857142857142857*f[7]*g[20]+0.447213595499958*f[1]*g[20]+0.2857142857142857*f[8]*g[12]+0.31943828249997*f[5]*g[12]+0.4472135954999579*f[4]*g[12]+0.5*f[0]*g[12]+0.4*f[3]*g[11]+0.31943828249997*f[7]*g[8]+0.5000000000000001*f[1]*g[8]+0.447213595499958*g[1]*f[8]+0.4472135954999579*f[7]*g[7]+0.5*g[0]*f[7]+0.4*g[4]*f[6]+0.5000000000000001*g[1]*f[5]+0.447213595499958*f[2]*g[4]+0.447213595499958*g[2]*f[3]); 
  fg[13] += a*(0.31943828249997*f[8]*g[23]+0.5000000000000001*f[5]*g[23]+0.4472135954999579*f[7]*g[18]+0.31943828249997*f[6]*g[17]+0.5000000000000001*f[2]*g[17]+0.5*f[8]*g[14]+0.31943828249997*f[4]*g[13]+0.5*f[0]*g[13]+0.447213595499958*f[3]*g[10]+0.5*f[6]*g[6]+0.447213595499958*f[1]*g[5]+0.5000000000000001*g[3]*f[4]); 
  fg[14] += a*(0.31943828249997*f[8]*g[23]+0.5000000000000001*f[4]*g[23]+0.31943828249997*f[7]*g[18]+0.5000000000000001*f[1]*g[18]+0.4472135954999579*f[6]*g[17]+0.31943828249997*f[5]*g[14]+0.5*f[0]*g[14]+0.5*f[8]*g[13]+0.447213595499958*f[3]*g[10]+0.5*g[5]*f[7]+0.447213595499958*f[2]*g[6]+0.5000000000000001*g[3]*f[5]); 
  fg[15] += a*(0.4472135954999579*f[7]*g[26]+0.447213595499958*f[8]*g[25]+0.5000000000000001*f[5]*g[25]+0.447213595499958*f[3]*g[24]+0.5*f[7]*g[22]+0.447213595499958*f[1]*g[21]+0.4472135954999579*f[6]*g[19]+0.5000000000000001*f[2]*g[19]+0.5*f[3]*g[16]+0.4472135954999579*f[4]*g[15]+0.5*f[0]*g[15]+0.5000000000000001*f[1]*g[9]); 
  fg[16] += a*(0.4472135954999579*f[6]*g[26]+0.447213595499958*f[3]*g[25]+0.447213595499958*f[8]*g[24]+0.5000000000000001*f[4]*g[24]+0.447213595499958*f[2]*g[22]+0.5*f[6]*g[21]+0.4472135954999579*f[7]*g[19]+0.5000000000000001*f[1]*g[19]+0.4472135954999579*f[5]*g[16]+0.5*f[0]*g[16]+0.5*f[3]*g[15]+0.5000000000000001*f[2]*g[9]); 
  fg[17] += a*(0.2857142857142858*f[6]*g[23]+0.4472135954999579*f[2]*g[23]+0.4*f[3]*g[18]+0.2857142857142857*f[8]*g[17]+0.4472135954999579*f[5]*g[17]+0.31943828249997*f[4]*g[17]+0.5*f[0]*g[17]+0.4472135954999579*f[6]*g[14]+0.31943828249997*f[6]*g[13]+0.5000000000000001*f[2]*g[13]+0.4*f[7]*g[10]+0.4472135954999579*f[1]*g[10]+0.4472135954999579*g[6]*f[8]+0.5*f[4]*g[6]+0.5000000000000001*g[3]*f[6]+0.4472135954999579*f[3]*g[5]); 
  fg[18] += a*(0.2857142857142858*f[7]*g[23]+0.4472135954999579*f[1]*g[23]+0.2857142857142857*f[8]*g[18]+0.31943828249997*f[5]*g[18]+0.4472135954999579*f[4]*g[18]+0.5*f[0]*g[18]+0.4*f[3]*g[17]+0.31943828249997*f[7]*g[14]+0.5000000000000001*f[1]*g[14]+0.4472135954999579*f[7]*g[13]+0.4*f[6]*g[10]+0.4472135954999579*f[2]*g[10]+0.4472135954999579*g[5]*f[8]+0.5000000000000001*g[3]*f[7]+0.4472135954999579*f[3]*g[6]+0.5*f[5]*g[5]); 
  fg[19] += a*(0.4*f[3]*g[26]+0.4*f[6]*g[25]+0.4472135954999579*f[2]*g[25]+0.4*f[7]*g[24]+0.4472135954999579*f[1]*g[24]+0.4472135954999579*f[3]*g[22]+0.4472135954999579*f[3]*g[21]+0.4*f[8]*g[19]+0.4472135954999579*f[5]*g[19]+0.4472135954999579*f[4]*g[19]+0.5*f[0]*g[19]+0.4472135954999579*f[7]*g[16]+0.5000000000000001*f[1]*g[16]+0.4472135954999579*f[6]*g[15]+0.5000000000000001*f[2]*g[15]+0.5*f[3]*g[9]); 
  fg[20] += a*(0.2040816326530612*f[8]*g[20]+0.31943828249997*f[5]*g[20]+0.31943828249997*f[4]*g[20]+0.5*f[0]*g[20]+0.2857142857142857*f[7]*g[12]+0.447213595499958*f[1]*g[12]+0.2857142857142857*f[6]*g[11]+0.447213595499958*f[2]*g[11]+0.31943828249997*f[8]*g[8]+0.5*f[4]*g[8]+0.31943828249997*g[7]*f[8]+0.5*g[0]*f[8]+0.5*f[5]*g[7]+0.447213595499958*g[1]*f[7]+0.447213595499958*g[2]*f[6]+0.4*f[3]*g[4]); 
  fg[21] += a*(0.31943828249997*f[8]*g[26]+0.5*f[5]*g[26]+0.447213595499958*f[7]*g[25]+0.31943828249997*f[6]*g[24]+0.5*f[2]*g[24]+0.5*f[8]*g[22]+0.31943828249997*f[4]*g[21]+0.5*f[0]*g[21]+0.4472135954999579*f[3]*g[19]+0.5*f[6]*g[16]+0.447213595499958*f[1]*g[15]+0.5*f[4]*g[9]); 
  fg[22] += a*(0.31943828249997*f[8]*g[26]+0.5*f[4]*g[26]+0.31943828249997*f[7]*g[25]+0.5*f[1]*g[25]+0.447213595499958*f[6]*g[24]+0.31943828249997*f[5]*g[22]+0.5*f[0]*g[22]+0.5*f[8]*g[21]+0.4472135954999579*f[3]*g[19]+0.447213595499958*f[2]*g[16]+0.5*f[7]*g[15]+0.5*f[5]*g[9]); 
  fg[23] += a*(0.2040816326530612*f[8]*g[23]+0.31943828249997*f[5]*g[23]+0.31943828249997*f[4]*g[23]+0.5*f[0]*g[23]+0.2857142857142858*f[7]*g[18]+0.4472135954999579*f[1]*g[18]+0.2857142857142858*f[6]*g[17]+0.4472135954999579*f[2]*g[17]+0.31943828249997*f[8]*g[14]+0.5000000000000001*f[4]*g[14]+0.31943828249997*f[8]*g[13]+0.5000000000000001*f[5]*g[13]+0.4*f[3]*g[10]+0.5*g[3]*f[8]+0.447213595499958*g[5]*f[7]+0.447213595499958*f[6]*g[6]); 
  fg[24] += a*(0.2857142857142858*f[6]*g[26]+0.4472135954999579*f[2]*g[26]+0.4*f[3]*g[25]+0.2857142857142857*f[8]*g[24]+0.4472135954999579*f[5]*g[24]+0.31943828249997*f[4]*g[24]+0.5*f[0]*g[24]+0.447213595499958*f[6]*g[22]+0.31943828249997*f[6]*g[21]+0.5*f[2]*g[21]+0.4*f[7]*g[19]+0.4472135954999579*f[1]*g[19]+0.447213595499958*f[8]*g[16]+0.5000000000000001*f[4]*g[16]+0.447213595499958*f[3]*g[15]+0.5000000000000001*f[6]*g[9]); 
  fg[25] += a*(0.2857142857142858*f[7]*g[26]+0.4472135954999579*f[1]*g[26]+0.2857142857142857*f[8]*g[25]+0.31943828249997*f[5]*g[25]+0.4472135954999579*f[4]*g[25]+0.5*f[0]*g[25]+0.4*f[3]*g[24]+0.31943828249997*f[7]*g[22]+0.5*f[1]*g[22]+0.447213595499958*f[7]*g[21]+0.4*f[6]*g[19]+0.4472135954999579*f[2]*g[19]+0.447213595499958*f[3]*g[16]+0.447213595499958*f[8]*g[15]+0.5000000000000001*f[5]*g[15]+0.5000000000000001*f[7]*g[9]); 
  fg[26] += a*(0.2040816326530612*f[8]*g[26]+0.31943828249997*f[5]*g[26]+0.31943828249997*f[4]*g[26]+0.5*f[0]*g[26]+0.2857142857142858*f[7]*g[25]+0.4472135954999579*f[1]*g[25]+0.2857142857142858*f[6]*g[24]+0.4472135954999579*f[2]*g[24]+0.31943828249997*f[8]*g[22]+0.5*f[4]*g[22]+0.31943828249997*f[8]*g[21]+0.5*f[5]*g[21]+0.4*f[3]*g[19]+0.4472135954999579*f[6]*g[16]+0.4472135954999579*f[7]*g[15]+0.5*f[8]*g[9]); 
} 