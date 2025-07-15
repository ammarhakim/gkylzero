#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 9.0*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 9.0*coeff[1]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffx(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 4.5*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 4.5*coeff[9]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffx(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_fluid_order2_vol_2x_tensor_p2_varcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

