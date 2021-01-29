#pragma once

#include <gkyl_real_type.h>

#include <math.h> 
typedef struct { gkyl_real c, chi, gamma; } gkyl_maxwell_inp; 

gkyl_real maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_1x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_2x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real maxwell_surfy_2x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_2x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_2x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real maxwell_surfy_2x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_3x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_3x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real maxwell_surfy_3x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real maxwell_surfz_3x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_1x_tensor_p2(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_1x_tensor_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

gkyl_real maxwell_vol_2x_tensor_p2(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out); 
gkyl_real maxwell_surfx_2x_tensor_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 
gkyl_real maxwell_surfy_2x_tensor_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr); 

