#pragma once
#include <gkyl_real_type.h>
#include <math.h> 

void vlasov_M0_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

void vlasov_M0_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M1i_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M2_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_FiveMoments_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2); 
void vlasov_M2ij_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3i_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_M3ijk_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_int_mom_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out); 

