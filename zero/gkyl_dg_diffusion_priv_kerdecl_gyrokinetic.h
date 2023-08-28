#pragma once

#include <gkyl_dg_diffusion_priv_kerdecl_common.h>

// Private header with declaration of gyrokinetic kernels for DG diffusion equation object.

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_gyrokinetic_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order2_gyrokinetic_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_gyrokinetic_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_3x2v_ser_p2_constcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfy_3x2v_ser_p2_constcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfz_3x2v_ser_p2_constcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
};

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_gyrokinetic_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order2_gyrokinetic_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_gyrokinetic_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_1x1v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfx_3x2v_ser_p2_varcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfy_3x2v_ser_p2_varcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_gyrokinetic_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_surfz_3x2v_ser_p2_varcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_1x1v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfx_3x2v_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfy_3x2v_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_gyrokinetic_boundary_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_order2_gyrokinetic_boundary_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_gyrokinetic_boundary_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_order4_gyrokinetic_boundary_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_gyrokinetic_boundary_surfz_3x2v_ser_p2_varcoeff },},
  },
};

