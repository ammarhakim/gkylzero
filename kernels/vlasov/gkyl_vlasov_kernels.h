#pragma once 
#include <math.h> 
double vlasov_vol_1x1v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x1v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x1v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_vol_1x1v_p2(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x1v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x1v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 


double vlasov_vol_1x2v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x2v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x2v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x2v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_vol_1x2v_p2(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x2v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x2v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x2v_vy_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 


double vlasov_vol_1x3v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x3v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vz_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_vol_1x3v_p2(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_1x3v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vy_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_1x3v_vz_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 


double vlasov_vol_2x2v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_2x2v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_2x2v_y_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x2v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x2v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_vol_2x2v_p2(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_2x2v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_2x2v_y_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x2v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x2v_vy_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 


double vlasov_vol_2x3v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_2x3v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_2x3v_y_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vz_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 

double vlasov_vol_2x3v_p2(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_2x3v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_2x3v_y_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vy_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_2x3v_vz_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 


double vlasov_vol_3x3v_p1(const double *w, const double *dxv, const double *EM, const double *f, double* restrict out); 
void vlasov_surf_3x3v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_3x3v_y_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
void vlasov_surf_3x3v_z_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_3x3v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_3x3v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr); 
double vlasov_surf_3x3v_vz_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr);


