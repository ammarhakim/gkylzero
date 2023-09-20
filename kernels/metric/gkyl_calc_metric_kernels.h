#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void gij_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_lo2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_up2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_lo2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_up2_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo2_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up2_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_lo2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_lo2_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_up2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_lo1_up2_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_lo2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_lo2_up3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_up2_lo3_3x_Ser_p1( const double **xyz, double *gij); 
GKYL_CU_DH void gij_up1_up2_up3_3x_Ser_p1( const double **xyz, double *gij); 
EXTERN_C_END 
