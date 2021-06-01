#pragma once 
#include <math.h>
#include <gkyl_util.h>

void vlasov_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

void vlasov_M0_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M1i_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M2_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_FiveMoments_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict outM0, double* restrict outM1i, double* restrict outM2); 
void vlasov_M2ij_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3i_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_M3ijk_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 
void vlasov_int_mom_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out); 

