#pragma once 
#include <math.h> 

double vlasov_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x1v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x1v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x1v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x1v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_1x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x2v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_2x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_2x3v_ser_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_3x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x1v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x1v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_1x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x2v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_stream_vol_2x3v_tensor_p2(const double *w, const double *dxv, const double *f, double* restrict out); 
double vlasov_vol_2x3v_tensor_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* restrict out); 
void vlasov_surfx_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfy_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surfz_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvx_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvy_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surfvz_2x3v_tensor_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

