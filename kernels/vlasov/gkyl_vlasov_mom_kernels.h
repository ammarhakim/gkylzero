#pragma once

void vlasov_mom_1x1v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_1x1v_m0_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m2_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x1v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_1x2v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_1x2v_m0_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m2_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x2v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_1x3v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_1x3v_m0_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m2_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_1x3v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_2x2v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_2x2v_m0_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m2_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x2v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_2x3v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_2x3v_m0_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m2_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_2x3v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);

void vlasov_mom_3x3v_m0_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_3x3v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_3x3v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_3x3v_m2_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
void vlasov_mom_3x3v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *f, double* restrict out);
