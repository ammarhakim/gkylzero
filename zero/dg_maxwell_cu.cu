/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
}

#define CK(lst,cdim,polyOrder) lst[cdim-1].kernels[polyOrder]

// various pointers to functions on device
__device__ static vol_termf_t d_vol = &vol;
__device__ static surf_termf_t d_surf = &surf;
__device__ static boundary_surf_termf_t d_boundary_surf = &boundary_surf;

// Volume kernel list
__device__ static struct { maxwell_vol_t kernels[3]; } d_vol_kernels[] = {
  { NULL, &maxwell_vol_1x_ser_p1, &maxwell_vol_1x_ser_p2 }, // 0
  { NULL, &maxwell_vol_2x_ser_p1, &maxwell_vol_2x_ser_p2 }, // 1
  { NULL, &maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
__device__ static struct { maxwell_surf_t kernels[3]; } d_surf_x_kernels[] = {
  { NULL, &maxwell_surfx_1x_ser_p1, &maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, &maxwell_surfx_2x_ser_p1, &maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, &maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
__device__ static struct { maxwell_surf_t kernels[3]; } d_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, &maxwell_surfy_2x_ser_p1, &maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, &maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
__device__ static struct { maxwell_surf_t kernels[3]; } d_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, &maxwell_surfz_3x_ser_p1, NULL }, // 2
};

__global__
void dg_maxwell_set_cu_dev_ptrs(struct dg_maxwell* maxwell, int cdim, int polyOrder)
{
  maxwell->eqn.vol_term = d_vol;
  maxwell->eqn.surf_term = d_surf;
  maxwell->eqn.boundary_surf_term = d_boundary_surf;

  maxwell->vol =  CK(d_vol_kernels, cdim, polyOrder);

  maxwell->surf[0] = CK(d_surf_x_kernels, cdim, polyOrder);
  if (cdim>1)
    maxwell->surf[1] = CK(d_surf_y_kernels, cdim, polyOrder);
  if (cdim>2)
    maxwell->surf[2] = CK(d_surf_z_kernels, cdim, polyOrder);
}

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

  // copy the host struct to device struct
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);

  dg_maxwell_set_cu_dev_ptrs<<<1,1>>>(maxwell_cu, cbasis->ndim, cbasis->polyOrder);

  gkyl_free(maxwell);
  
  return &maxwell_cu->eqn;
}
