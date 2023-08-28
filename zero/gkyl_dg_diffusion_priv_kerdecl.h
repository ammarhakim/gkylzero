#pragma once

#include <gkyl_dg_diffusion_priv_kerdecl_common.h>

// Private header with declaration of kernels for DG diffusion equation object.

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_surfx_1x_ser_p1_constcoeff, dg_diffusion_order2_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfx_2x_ser_p1_constcoeff, dg_diffusion_order2_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_surfx_3x_ser_p1_constcoeff, dg_diffusion_order2_surfx_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_surfx_1x_ser_p1_constcoeff, dg_diffusion_order4_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfx_2x_ser_p1_constcoeff, dg_diffusion_order4_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_surfx_3x_ser_p1_constcoeff, dg_diffusion_order4_surfx_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfx_3x_ser_p2_constcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfy_2x_ser_p1_constcoeff, dg_diffusion_order2_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_surfy_3x_ser_p1_constcoeff, dg_diffusion_order2_surfy_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfy_2x_ser_p1_constcoeff, dg_diffusion_order4_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_surfy_3x_ser_p1_constcoeff, dg_diffusion_order4_surfy_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_3x_ser_p2_constcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfz_3x_ser_p1_constcoeff, dg_diffusion_order2_surfz_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfz_3x_ser_p1_constcoeff, dg_diffusion_order4_surfz_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfz_3x_ser_p2_constcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_boundary_surfx_1x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfx_2x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfx_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_boundary_surfx_1x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfx_2x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfx_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_boundary_surfx_1x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfx_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfx_3x_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_2x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfy_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_2x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfy_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_2x_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_3x_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfz_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfz_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfz_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfz_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfz_3x_ser_p2_constcoeff },},
  },
};

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_surfx_1x_ser_p1_varcoeff, dg_diffusion_order2_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfx_2x_ser_p1_varcoeff, dg_diffusion_order2_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_surfx_3x_ser_p1_varcoeff, dg_diffusion_order2_surfx_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_surfx_1x_ser_p1_varcoeff, dg_diffusion_order4_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfx_2x_ser_p1_varcoeff, dg_diffusion_order4_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_surfx_3x_ser_p1_varcoeff, dg_diffusion_order4_surfx_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfx_3x_ser_p2_varcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfy_2x_ser_p1_varcoeff, dg_diffusion_order2_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_surfy_3x_ser_p1_varcoeff, dg_diffusion_order2_surfy_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfy_2x_ser_p1_varcoeff, dg_diffusion_order4_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_surfy_3x_ser_p1_varcoeff, dg_diffusion_order4_surfy_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_3x_ser_p2_varcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfz_3x_ser_p1_varcoeff, dg_diffusion_order2_surfz_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfz_3x_ser_p1_varcoeff, dg_diffusion_order4_surfz_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfz_3x_ser_p2_varcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_boundary_surfx_1x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfx_2x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfx_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_boundary_surfx_1x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfx_2x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfx_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_boundary_surfx_1x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfx_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfx_3x_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_2x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfy_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_2x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfy_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_2x_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_3x_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfz_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfz_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfz_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfz_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfz_3x_ser_p2_varcoeff },},
  },
};
