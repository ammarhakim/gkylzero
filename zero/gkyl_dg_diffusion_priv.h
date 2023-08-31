#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*diffusion_surf_t)(const double *w, const double *dx,
  double D, const double *ql, const double *qc, const double *qr,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_diffusion_vol_kern_list;
typedef struct { diffusion_surf_t kernels[3]; } gkyl_dg_diffusion_surf_kern_list;

struct dg_diffusion {
  struct gkyl_dg_eqn eqn;
  diffusion_surf_t surf[3];
  double D; // Constant diffusion coefficient
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_1x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_1x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_2x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_2x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_3x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion_vol_3x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_1x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_1x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_2x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_2x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_3x_ser_p1(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion4_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion4_vol_3x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion6_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion6_vol_1x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion6_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion6_vol_2x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion6_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return dg_diffusion6_vol_3x_ser_p2(xc, dx, diffusion->D, qIn, qRhsOut);
}

// Volume kernel grad^2 diffusion (only returns CFL, same for all basis types and equation systems)
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_dg_diffusion_vol_1x_ser_p1, kernel_dg_diffusion_vol_1x_ser_p2 },
  { NULL, kernel_dg_diffusion_vol_2x_ser_p1, kernel_dg_diffusion_vol_2x_ser_p2 },
  { NULL, kernel_dg_diffusion_vol_3x_ser_p1, kernel_dg_diffusion_vol_3x_ser_p2 },
};

// Volume kernel grad^4 diffusion (only returns CFL, same for all basis types and equation systems)
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol4_kernels[] = {
  { NULL, kernel_dg_diffusion4_vol_1x_ser_p1, kernel_dg_diffusion4_vol_1x_ser_p2 },
  { NULL, kernel_dg_diffusion4_vol_2x_ser_p1, kernel_dg_diffusion4_vol_2x_ser_p2 },
  { NULL, kernel_dg_diffusion4_vol_3x_ser_p1, kernel_dg_diffusion4_vol_3x_ser_p2 },
};

// Volume kernel grad^6 diffusion (only returns CFL, same for all basis types and equation systems)
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol6_kernels[] = {
  { NULL, NULL, kernel_dg_diffusion6_vol_1x_ser_p2 },
  { NULL, NULL, kernel_dg_diffusion6_vol_2x_ser_p2 },
  { NULL, NULL, kernel_dg_diffusion6_vol_3x_ser_p2 },
};

//
// Scalar isotropic diffusion
//

// grad^2 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, dg_diffusion_surfx_1x_ser_p1, dg_diffusion_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_surfx_2x_ser_p1, dg_diffusion_surfx_2x_ser_p2 },
  { NULL, dg_diffusion_surfx_3x_ser_p1, dg_diffusion_surfx_3x_ser_p2 },
};

// grad^2 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_surfy_2x_ser_p1, dg_diffusion_surfy_2x_ser_p2 },
  { NULL, dg_diffusion_surfy_3x_ser_p1, dg_diffusion_surfy_3x_ser_p2 },
};

// grad^2 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_surfz_3x_ser_p1, dg_diffusion_surfz_3x_ser_p2 },
};

// grad^2 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, dg_diffusion_surfx_1x_ser_p1, dg_diffusion_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_surfx_2x_ser_p1, dg_diffusion_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion_surfx_3x_ser_p1, dg_diffusion_surfx_3x_tensor_p2 },
};

// grad^2 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_surfy_2x_ser_p1, dg_diffusion_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion_surfy_3x_ser_p1, dg_diffusion_surfy_3x_tensor_p2 },
};

// grad^2 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_surfz_3x_ser_p1, dg_diffusion_surfz_3x_tensor_p2 },
};

// grad^4 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_x_kernels[] = {
  { NULL, dg_diffusion4_surfx_1x_ser_p1, dg_diffusion4_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_surfx_2x_ser_p1, dg_diffusion4_surfx_2x_ser_p2 },
  { NULL, dg_diffusion4_surfx_3x_ser_p1, dg_diffusion4_surfx_3x_ser_p2 },
};

// grad^4 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_surfy_2x_ser_p1, dg_diffusion4_surfy_2x_ser_p2 },
  { NULL, dg_diffusion4_surfy_3x_ser_p1, dg_diffusion4_surfy_3x_ser_p2 },
};

// grad^4 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_surfz_3x_ser_p1, dg_diffusion4_surfz_3x_ser_p2 },
};

// grad^4 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_x_kernels[] = {
  { NULL, dg_diffusion4_surfx_1x_ser_p1, dg_diffusion4_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_surfx_2x_ser_p1, dg_diffusion4_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion4_surfx_3x_ser_p1, dg_diffusion4_surfx_3x_tensor_p2 },
};

// grad^4 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_surfy_2x_ser_p1, dg_diffusion4_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion4_surfy_3x_ser_p1, dg_diffusion4_surfy_3x_tensor_p2 },
};

// grad^4 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_surfz_3x_ser_p1, dg_diffusion4_surfz_3x_tensor_p2 },
};

// grad^6 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_x_kernels[] = {
  { NULL, NULL, dg_diffusion6_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_surfx_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_surfx_3x_ser_p2 },
};

// grad^6 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_surfy_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_surfy_3x_ser_p2 },
};

// grad^6 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_surfz_3x_ser_p2 },
};

// grad^6 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_x_kernels[] = {
  { NULL, NULL, dg_diffusion6_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_surfx_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_surfx_3x_tensor_p2 },
};

// grad^6 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_surfy_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_surfy_3x_tensor_p2 },
};

// grad^6 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_surfz_3x_tensor_p2 },
};

//
// PKPM isotropic diffusion
//

// PKPM grad^2 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_x_pkpm_kernels[] = {
  { NULL, dg_diffusion_pkpm_surfx_1x_ser_p1, dg_diffusion_pkpm_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_pkpm_surfx_2x_ser_p1, dg_diffusion_pkpm_surfx_2x_ser_p2 },
  { NULL, dg_diffusion_pkpm_surfx_3x_ser_p1, dg_diffusion_pkpm_surfx_3x_ser_p2 },
};

// PKPM grad^2 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_pkpm_surfy_2x_ser_p1, dg_diffusion_pkpm_surfy_2x_ser_p2 },
  { NULL, dg_diffusion_pkpm_surfy_3x_ser_p1, dg_diffusion_pkpm_surfy_3x_ser_p2 },
};

// PKPM grad^2 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_pkpm_surfz_3x_ser_p1, dg_diffusion_pkpm_surfz_3x_ser_p2 },
};

// PKPM grad^2 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_x_pkpm_kernels[] = {
  { NULL, dg_diffusion_pkpm_surfx_1x_ser_p1, dg_diffusion_pkpm_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_pkpm_surfx_2x_ser_p1, dg_diffusion_pkpm_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion_pkpm_surfx_3x_ser_p1, dg_diffusion_pkpm_surfx_3x_tensor_p2 },
};

// PKPM grad^2 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_pkpm_surfy_2x_ser_p1, dg_diffusion_pkpm_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion_pkpm_surfy_3x_ser_p1, dg_diffusion_pkpm_surfy_3x_tensor_p2 },
};

// PKPM grad^2 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_pkpm_surfz_3x_ser_p1, dg_diffusion_pkpm_surfz_3x_tensor_p2 },
};

// PKPM grad^4 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_x_pkpm_kernels[] = {
  { NULL, dg_diffusion4_pkpm_surfx_1x_ser_p1, dg_diffusion4_pkpm_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_pkpm_surfx_2x_ser_p1, dg_diffusion4_pkpm_surfx_2x_ser_p2 },
  { NULL, dg_diffusion4_pkpm_surfx_3x_ser_p1, dg_diffusion4_pkpm_surfx_3x_ser_p2 },
};

// PKPM grad^4 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_pkpm_surfy_2x_ser_p1, dg_diffusion4_pkpm_surfy_2x_ser_p2 },
  { NULL, dg_diffusion4_pkpm_surfy_3x_ser_p1, dg_diffusion4_pkpm_surfy_3x_ser_p2 },
};

// PKPM grad^4 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_pkpm_surfz_3x_ser_p1, dg_diffusion4_pkpm_surfz_3x_ser_p2 },
};

// PKPM grad^4 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_x_pkpm_kernels[] = {
  { NULL, dg_diffusion4_pkpm_surfx_1x_ser_p1, dg_diffusion4_pkpm_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_pkpm_surfx_2x_ser_p1, dg_diffusion4_pkpm_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion4_pkpm_surfx_3x_ser_p1, dg_diffusion4_pkpm_surfx_3x_tensor_p2 },
};

// PKPM grad^4 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_pkpm_surfy_2x_ser_p1, dg_diffusion4_pkpm_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion4_pkpm_surfy_3x_ser_p1, dg_diffusion4_pkpm_surfy_3x_tensor_p2 },
};

// PKPM grad^4 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_pkpm_surfz_3x_ser_p1, dg_diffusion4_pkpm_surfz_3x_tensor_p2 },
};

// PKPM grad^6 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_x_pkpm_kernels[] = {
  { NULL, NULL, dg_diffusion6_pkpm_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfx_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfx_3x_ser_p2 },
};

// PKPM grad^6 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_pkpm_surfy_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfy_3x_ser_p2 },
};

// PKPM grad^6 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_pkpm_surfz_3x_ser_p2 },
};

// PKPM grad^6 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_x_pkpm_kernels[] = {
  { NULL, NULL, dg_diffusion6_pkpm_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfx_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfx_3x_tensor_p2 },
};

// PKPM grad^6 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_y_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_pkpm_surfy_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_pkpm_surfy_3x_tensor_p2 },
};

// PKPM grad^6 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_z_pkpm_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_pkpm_surfz_3x_tensor_p2 },
};

//
// Isothermal Euler isotropic diffusion
//

// Isothermal Euler grad^2 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_x_iso_euler_kernels[] = {
  { NULL, dg_diffusion_iso_euler_surfx_1x_ser_p1, dg_diffusion_iso_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_iso_euler_surfx_2x_ser_p1, dg_diffusion_iso_euler_surfx_2x_ser_p2 },
  { NULL, dg_diffusion_iso_euler_surfx_3x_ser_p1, dg_diffusion_iso_euler_surfx_3x_ser_p2 },
};

// Isothermal Euler grad^2 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_iso_euler_surfy_2x_ser_p1, dg_diffusion_iso_euler_surfy_2x_ser_p2 },
  { NULL, dg_diffusion_iso_euler_surfy_3x_ser_p1, dg_diffusion_iso_euler_surfy_3x_ser_p2 },
};

// Isothermal Euler grad^2 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_iso_euler_surfz_3x_ser_p1, dg_diffusion_iso_euler_surfz_3x_ser_p2 },
};

// Isothermal Euler grad^2 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_x_iso_euler_kernels[] = {
  { NULL, dg_diffusion_iso_euler_surfx_1x_ser_p1, dg_diffusion_iso_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_iso_euler_surfx_2x_ser_p1, dg_diffusion_iso_euler_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion_iso_euler_surfx_3x_ser_p1, dg_diffusion_iso_euler_surfx_3x_tensor_p2 },
};

// Isothermal Euler grad^2 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_iso_euler_surfy_2x_ser_p1, dg_diffusion_iso_euler_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion_iso_euler_surfy_3x_ser_p1, dg_diffusion_iso_euler_surfy_3x_tensor_p2 },
};

// Isothermal Euler grad^2 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_iso_euler_surfz_3x_ser_p1, dg_diffusion_iso_euler_surfz_3x_tensor_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_x_iso_euler_kernels[] = {
  { NULL, dg_diffusion4_iso_euler_surfx_1x_ser_p1, dg_diffusion4_iso_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_iso_euler_surfx_2x_ser_p1, dg_diffusion4_iso_euler_surfx_2x_ser_p2 },
  { NULL, dg_diffusion4_iso_euler_surfx_3x_ser_p1, dg_diffusion4_iso_euler_surfx_3x_ser_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_iso_euler_surfy_2x_ser_p1, dg_diffusion4_iso_euler_surfy_2x_ser_p2 },
  { NULL, dg_diffusion4_iso_euler_surfy_3x_ser_p1, dg_diffusion4_iso_euler_surfy_3x_ser_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_iso_euler_surfz_3x_ser_p1, dg_diffusion4_iso_euler_surfz_3x_ser_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_x_iso_euler_kernels[] = {
  { NULL, dg_diffusion4_iso_euler_surfx_1x_ser_p1, dg_diffusion4_iso_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_iso_euler_surfx_2x_ser_p1, dg_diffusion4_iso_euler_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion4_iso_euler_surfx_3x_ser_p1, dg_diffusion4_iso_euler_surfx_3x_tensor_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_iso_euler_surfy_2x_ser_p1, dg_diffusion4_iso_euler_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion4_iso_euler_surfy_3x_ser_p1, dg_diffusion4_iso_euler_surfy_3x_tensor_p2 },
};

// Isothermal Euler grad^4 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_iso_euler_surfz_3x_ser_p1, dg_diffusion4_iso_euler_surfz_3x_tensor_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_x_iso_euler_kernels[] = {
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_3x_ser_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_iso_euler_surfy_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfy_3x_ser_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_iso_euler_surfz_3x_ser_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_x_iso_euler_kernels[] = {
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfx_3x_tensor_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_y_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_iso_euler_surfy_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_iso_euler_surfy_3x_tensor_p2 },
};

// Isothermal Euler grad^6 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_z_iso_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_iso_euler_surfz_3x_tensor_p2 },
};

//
// Euler isotropic diffusion
//

// Euler grad^2 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_x_euler_kernels[] = {
  { NULL, dg_diffusion_euler_surfx_1x_ser_p1, dg_diffusion_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_euler_surfx_2x_ser_p1, dg_diffusion_euler_surfx_2x_ser_p2 },
  { NULL, dg_diffusion_euler_surfx_3x_ser_p1, dg_diffusion_euler_surfx_3x_ser_p2 },
};

// Euler grad^2 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_euler_surfy_2x_ser_p1, dg_diffusion_euler_surfy_2x_ser_p2 },
  { NULL, dg_diffusion_euler_surfy_3x_ser_p1, dg_diffusion_euler_surfy_3x_ser_p2 },
};

// Euler grad^2 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_euler_surfz_3x_ser_p1, dg_diffusion_euler_surfz_3x_ser_p2 },
};

// Euler grad^2 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_x_euler_kernels[] = {
  { NULL, dg_diffusion_euler_surfx_1x_ser_p1, dg_diffusion_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_euler_surfx_2x_ser_p1, dg_diffusion_euler_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion_euler_surfx_3x_ser_p1, dg_diffusion_euler_surfx_3x_tensor_p2 },
};

// Euler grad^2 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_euler_surfy_2x_ser_p1, dg_diffusion_euler_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion_euler_surfy_3x_ser_p1, dg_diffusion_euler_surfy_3x_tensor_p2 },
};

// Euler grad^2 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_euler_surfz_3x_ser_p1, dg_diffusion_euler_surfz_3x_tensor_p2 },
};

// Euler grad^4 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_x_euler_kernels[] = {
  { NULL, dg_diffusion4_euler_surfx_1x_ser_p1, dg_diffusion4_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_euler_surfx_2x_ser_p1, dg_diffusion4_euler_surfx_2x_ser_p2 },
  { NULL, dg_diffusion4_euler_surfx_3x_ser_p1, dg_diffusion4_euler_surfx_3x_ser_p2 },
};

// Euler grad^4 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_euler_surfy_2x_ser_p1, dg_diffusion4_euler_surfy_2x_ser_p2 },
  { NULL, dg_diffusion4_euler_surfy_3x_ser_p1, dg_diffusion4_euler_surfy_3x_ser_p2 },
};

// Euler grad^4 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf4_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_euler_surfz_3x_ser_p1, dg_diffusion4_euler_surfz_3x_ser_p2 },
};

// Euler grad^4 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_x_euler_kernels[] = {
  { NULL, dg_diffusion4_euler_surfx_1x_ser_p1, dg_diffusion4_euler_surfx_1x_ser_p2 },
  { NULL, dg_diffusion4_euler_surfx_2x_ser_p1, dg_diffusion4_euler_surfx_2x_tensor_p2 },
  { NULL, dg_diffusion4_euler_surfx_3x_ser_p1, dg_diffusion4_euler_surfx_3x_tensor_p2 },
};

// Euler grad^4 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion4_euler_surfy_2x_ser_p1, dg_diffusion4_euler_surfy_2x_tensor_p2 },
  { NULL, dg_diffusion4_euler_surfy_3x_ser_p1, dg_diffusion4_euler_surfy_3x_tensor_p2 },
};

// Euler grad^4 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf4_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion4_euler_surfz_3x_ser_p1, dg_diffusion4_euler_surfz_3x_tensor_p2 },
};

// Euler grad^6 Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_x_euler_kernels[] = {
  { NULL, NULL, dg_diffusion6_euler_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfx_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfx_3x_ser_p2 },
};

// Euler grad^6 Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_euler_surfy_2x_ser_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfy_3x_ser_p2 },
};

// Euler grad^6 Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf6_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_euler_surfz_3x_ser_p2 },
};

// Euler grad^6 Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_x_euler_kernels[] = {
  { NULL, NULL, dg_diffusion6_euler_surfx_1x_ser_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfx_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfx_3x_tensor_p2 },
};

// Euler grad^6 Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_y_euler_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, dg_diffusion6_euler_surfy_2x_tensor_p2 },
  { NULL, NULL, dg_diffusion6_euler_surfy_3x_tensor_p2 },
};

// Euler grad^6 Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ten_surf6_z_euler_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, dg_diffusion6_euler_surfz_3x_tensor_p2 },
};

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  return diffusion->surf[dir](xcC, dxC,
    diffusion->D, 
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  return 0.;
}
