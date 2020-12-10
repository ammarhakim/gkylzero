#pragma once 
#include <math.h> 
typedef struct { double c, chi, gamma; } MaxwellEq_t; 
 
double maxwell_vol_1x_p1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_1x_x_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 

double maxwell_vol_1x_p2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_1x_x_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 


double maxwell_vol_2x_p1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_2x_x_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_2x_y_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 

double maxwell_vol_2x_p2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_2x_x_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_2x_y_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 


double maxwell_vol_3x_p1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_3x_x_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_3x_y_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_3x_z_p1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 

double maxwell_vol_3x_p2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double* restrict out); 
double maxwell_surf_3x_x_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_3x_y_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 
double maxwell_surf_3x_z_p2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr); 


