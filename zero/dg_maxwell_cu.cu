/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
}

// Volume kernel list
__device__ static struct { maxwell_vol_t kernels[3]; } p_vol_kernels[] = {
  { NULL, &maxwell_vol_1x_ser_p1, &maxwell_vol_1x_ser_p2 }, // 0
  { NULL, &maxwell_vol_2x_ser_p1, &maxwell_vol_2x_ser_p2 }, // 1
  { NULL, &maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
__device__ static struct { maxwell_surf_t kernels[3]; } p_surf_x_kernels[] = {
  { NULL, &maxwell_surfx_1x_ser_p1, &maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, &maxwell_surfx_2x_ser_p1, &maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, &maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
__device__ static struct { maxwell_surf_t kernels[3]; } p_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, &maxwell_surfy_2x_ser_p1, &maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, &maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
__device__ static struct { maxwell_surf_t kernels[3]; } p_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, &maxwell_surfz_3x_ser_p1, NULL }, // 2
};

// various pointers to functions on device
__device__ static vol_termf_t p_vol = &vol;
__device__ static surf_termf_t p_surf = &surf;
__device__ static boundary_surf_termf_t p_boundary_surf = &boundary_surf;

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  struct dg_maxwell *maxwell = (struct dg_maxwell*) gkyl_malloc(sizeof(struct dg_maxwell));

  // set basic parameters
  maxwell->eqn.num_equations = 8;
  maxwell->maxwell_data.c = lightSpeed;
  maxwell->maxwell_data.chi = lightSpeed*elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = lightSpeed*mgnErrorSpeedFactor;  

  int cdim = cbasis->ndim;
  int polyOrder = cbasis->polyOrder;

  // copy appropriate device function pointers to host to store in struct
  cudaMemcpyFromSymbol(&maxwell->eqn.vol_term, p_vol, sizeof(vol_termf_t));
  cudaMemcpyFromSymbol(&maxwell->eqn.surf_term, p_surf, sizeof(surf_termf_t));
  cudaMemcpyFromSymbol(&maxwell->eqn.boundary_surf_term, p_boundary_surf, sizeof(boundary_surf_termf_t));

  cudaMemcpyFromSymbol(&maxwell->vol, CK(p_vol_kernels, cdim, polyOrder), sizeof(maxwell_vol_t));

  cudaMemcpyFromSymbol(&maxwell->surf[0], CK(p_surf_x_kernels, cdim, polyOrder), sizeof(maxwell_surf_t));
  if (cdim>1)
    cudaMemcpyFromSymbol(&maxwell->surf[1], CK(p_surf_y_kernels, cdim, polyOrder), sizeof(maxwell_surf_t));
  if (cdim>2)
    cudaMemcpyFromSymbol(&maxwell->surf[2], CK(p_surf_z_kernels, cdim, polyOrder), sizeof(maxwell_surf_t));

  // copy the host struct to device struct
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);

  gkyl_free(maxwell);
  
  return &maxwell_cu->eqn;
}
