#pragma once 
#include <math.h> 
#include <gkyl_util.h> 

// These kernels will need to change void --> double to return cfl.

GKYL_CU_DH static inline double 
sigma_cx_1x1v_ser_p1(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[2]:     neutral particle density. 
  // u_ion[2]:       ion fluid velocity. 
  // u_neut[2]:      neutral fluid velocity. 
  // vt_sq_ion[2]:   ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[2]:  neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:     cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0];
  
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_ion_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
  
  double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b;

  return 0.1666666666666667*m0_neut[0]*v_sigma_cx[0];
  }  
}

GKYL_CU_DH static inline double 
sigma_cx_1x1v_ser_p2(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[3]:         neutral particle density. 
  // u_ion[3]:        ion fluid velocity. 
  // u_neut[3]:       neutral fluid velocity. 
  // vt_sq_ion[3]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0;
  } else {
    double v_in_sq_av = 0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
    
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b;

    return 0.1*m0_neut[0]*v_sigma_cx[0];
  }
} 

GKYL_CU_DH static inline double
sigma_cx_1x2v_ser_p1(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[2]:         neutral particle density. 
  // u_ion[4]:        ion fluid velocity. 
  // u_neut[4]:       neutral fluid velocity. 
  // vt_sq_ion[2]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[2]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.5*pow(u_neut[2],2)-1.0*u_ion[2]*u_neut[2]+0.5*pow(u_ion[2],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b;

    return 0.1666666666666667*m0_neut[0]*v_sigma_cx[0];  
  }
}

GKYL_CU_DH static inline double
sigma_cx_1x2v_ser_p2(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[3]:         neutral particle density. 
  // u_ion[6]:        ion fluid velocity. 
  // u_neut[6]:       neutral fluid velocity. 
  // vt_sq_ion[3]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.5*pow(u_neut[3],2)-1.0*u_ion[3]*u_neut[3]+0.5*pow(u_ion[3],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b;

    return 0.1*m0_neut[0]*v_sigma_cx[0]; 
  }
}

GKYL_CU_DH static inline double
sigma_cx_1x3v_ser_p1(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[2]:         neutral particle density. 
  // u_ion[6]:        ion fluid velocity. 
  // u_neut[6]:       neutral fluid velocity. 
  // vt_sq_ion[2]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[2]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.5*pow(u_neut[4],2)-1.0*u_ion[4]*u_neut[4]+0.5*pow(u_ion[4],2)+0.5*pow(u_neut[2],2)-1.0*u_ion[2]*u_neut[2]+0.5*pow(u_ion[2],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b;

    return 0.1666666666666667*m0_neut[0]*v_sigma_cx[0]; 
  }
}

GKYL_CU_DH static inline double
sigma_cx_1x3v_ser_p2(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[3]:         neutral particle density. 
  // u_ion[9]:        ion fluid velocity. 
  // u_neut[9]:       neutral fluid velocity. 
  // vt_sq_ion[3]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.7071067811865476*m0_neut[0]; 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.5*pow(u_neut[6],2)-1.0*u_ion[6]*u_neut[6]+0.5*pow(u_ion[6],2)+0.5*pow(u_neut[3],2)-1.0*u_ion[3]*u_neut[3]+0.5*pow(u_ion[3],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 1.414213562373095*v_cx*a-1.414213562373095*v_cx*log(v_cx)*b; 

    return 0.1*m0_neut[0]*v_sigma_cx[0]; 
  }
} 

GKYL_CU_DH static inline double
sigma_cx_2x3v_ser_p1(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[4]:         neutral particle density. 
  // u_ion[12]:        ion fluid velocity. 
  // u_neut[12]:       neutral fluid velocity. 
  // vt_sq_ion[4]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[4]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.5*m0_neut[0]; 
  double vt_sq_ion_av = 0.5*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.5*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.25*pow(u_neut[8],2)-0.5*u_ion[8]*u_neut[8]+0.25*pow(u_ion[8],2)+0.25*pow(u_neut[4],2)-0.5*u_ion[4]*u_neut[4]+0.25*pow(u_ion[4],2)+0.25*pow(u_neut[0],2)-0.5*u_ion[0]*u_neut[0]+0.25*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 2.0*v_cx*a-2.0*v_cx*log(v_cx)*b;

    return 0.08333333333333333*m0_neut[0]*v_sigma_cx[0]; 
  }
}

GKYL_CU_DH static inline double
sigma_cx_2x3v_ser_p2(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[8]:         neutral particle density. 
  // u_ion[24]:        ion fluid velocity. 
  // u_neut[24]:       neutral fluid velocity. 
  // vt_sq_ion[8]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[8]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.5*m0_neut[0]; 
  double vt_sq_ion_av = 0.5*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.5*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.25*pow(u_neut[16],2)-0.5*u_ion[16]*u_neut[16]+0.25*pow(u_ion[16],2)+0.25*pow(u_neut[8],2)-0.5*u_ion[8]*u_neut[8]+0.25*pow(u_ion[8],2)+0.25*pow(u_neut[0],2)-0.5*u_ion[0]*u_neut[0]+0.25*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 2.0*v_cx*a-2.0*v_cx*log(v_cx)*b; 

    return 0.05*m0_neut[0]*v_sigma_cx[0]; 
  }
} 

GKYL_CU_DH static inline double
sigma_cx_3x3v_ser_p1(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[8]:         neutral particle density. 
  // u_ion[24]:        ion fluid velocity. 
  // u_neut[24]:       neutral fluid velocity. 
  // vt_sq_ion[8]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[8]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.3535533905932738*m0_neut[0]; 
  double vt_sq_ion_av = 0.3535533905932738*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.3535533905932738*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.125*pow(u_neut[16],2)-0.25*u_ion[16]*u_neut[16]+0.125*pow(u_ion[16],2)+0.125*pow(u_neut[8],2)-0.25*u_ion[8]*u_neut[8]+0.125*pow(u_ion[8],2)+0.125*pow(u_neut[0],2)-0.25*u_ion[0]*u_neut[0]+0.125*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 2.828427124746191*v_cx*a-2.828427124746191*v_cx*log(v_cx)*b;

    return 0.04166666666666666*m0_neut[0]*v_sigma_cx[0]; 
  }
}

GKYL_CU_DH static inline double
sigma_cx_3x3v_ser_p2(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0_neut[20]:         neutral particle density. 
  // u_ion[60]:        ion fluid velocity. 
  // u_neut[60]:       neutral fluid velocity. 
  // vt_sq_ion[20]:     ion squared thermal speed, sqrt(T/m). 
  // vt_sq_neut[20]:    neutral squared thermal speed, sqrt(T/m). 
  // v_sigma_cx:          cell ave cross section fitting eqn. 
 
  double m0_neut_av = 0.3535533905932738*m0_neut[0]; 
  double vt_sq_ion_av = 0.3535533905932738*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.3535533905932738*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
    double v_in_sq_av = 0.125*pow(u_neut[40],2)-0.25*u_ion[40]*u_neut[40]+0.125*pow(u_ion[40],2)+0.125*pow(u_neut[20],2)-0.25*u_ion[20]*u_neut[20]+0.125*pow(u_ion[20],2)+0.125*pow(u_neut[0],2)-0.25*u_ion[0]*u_neut[0]+0.125*pow(u_ion[0],2); 
 
    double v_cx = sqrt(1.273239544735163*vt_sq_neut_av+1.273239544735163*vt_sq_ion_av+v_in_sq_av);
    v_sigma_cx[0] = 2.828427124746191*v_cx*a-2.828427124746191*v_cx*log(v_cx)*b;

    return 0.025*m0_neut[0]*v_sigma_cx[0]; 
  }
}
