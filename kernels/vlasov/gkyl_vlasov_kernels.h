#pragma once
#include <gkyl_real_type.h>
#include <math.h> 

gkyl_real vlasov_stream_vol_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x1v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x1v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x1v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x1v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_1x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_1x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x2v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x2v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x3v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_2x3v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_3x3v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_3x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x1v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x1v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x1v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_1x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_1x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x2v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x2v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real vlasov_stream_vol_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out); 
gkyl_real vlasov_vol_2x3v_tensor_p2(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out); 
void vlasov_surfx_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfy_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
void vlasov_surfz_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvx_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvy_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real vlasov_surfvz_2x3v_tensor_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

