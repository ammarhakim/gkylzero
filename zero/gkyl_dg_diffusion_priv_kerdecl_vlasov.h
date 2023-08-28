#pragma once

#include <gkyl_dg_diffusion_priv_kerdecl_common.h>

// Private header with declaration of vlasov kernels for DG diffusion equation object.

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_vlasov_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_1x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_surfx_1x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_1x3v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_surfx_2x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order2_vlasov_surfx_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_vlasov_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_1x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_surfx_1x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_1x3v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_surfx_2x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order4_vlasov_surfx_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_vlasov_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_1x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_1x3v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_2x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_2x3v_ser_p2_constcoeff },
//      { NULL, dg_diffusion_order6_vlasov_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_vlasov_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfy_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_surfy_2x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfy_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order2_vlasov_surfy_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_vlasov_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfy_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_surfy_2x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfy_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order4_vlasov_surfy_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_vlasov_surfy_2x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfy_2x3v_ser_p2_constcoeff },
//      { NULL, dg_diffusion_order6_vlasov_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order2_vlasov_surfz_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order4_vlasov_surfz_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { NULL, dg_diffusion_order6_vlasov_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_vlasov_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_1x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x3v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_2x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order2_vlasov_boundary_surfx_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_vlasov_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_1x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x3v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_2x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order4_vlasov_boundary_surfx_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x3v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_2x3v_ser_p2_constcoeff },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_vlasov_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfy_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfy_2x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfy_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order2_vlasov_boundary_surfy_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p2_constcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfy_2x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfy_2x3v_ser_p2_constcoeff },
//      { dg_diffusion_order4_vlasov_boundary_surfy_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_2x3v_ser_p2_constcoeff },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order2_vlasov_boundary_surfz_3x3v_ser_p1_constcoeff, dg_diffusion_order2_vlasov_boundary_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order4_vlasov_boundary_surfz_3x3v_ser_p1_constcoeff, dg_diffusion_order4_vlasov_boundary_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfz_3x3v_ser_p2_constcoeff },},
      { NULL, NULL },},
  },
};

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_vlasov_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_1x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_surfx_1x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_1x3v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_surfx_2x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order2_vlasov_surfx_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_vlasov_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_1x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_surfx_1x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_1x3v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_surfx_2x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order4_vlasov_surfx_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_vlasov_surfx_1x1v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_1x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_1x3v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_2x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfx_2x3v_ser_p2_varcoeff },
//      { NULL, dg_diffusion_order6_vlasov_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_vlasov_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfy_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_surfy_2x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfy_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order2_vlasov_surfy_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_vlasov_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfy_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_surfy_2x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfy_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order4_vlasov_surfy_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_vlasov_surfy_2x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_surfy_2x3v_ser_p2_varcoeff },
//      { NULL, dg_diffusion_order6_vlasov_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_vlasov_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order2_vlasov_surfz_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order4_vlasov_surfz_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { NULL, dg_diffusion_order6_vlasov_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_vlasov_boundary_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_1x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_1x3v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfx_2x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order2_vlasov_boundary_surfx_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_vlasov_boundary_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_1x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_1x3v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfx_2x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order4_vlasov_boundary_surfx_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x1v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_1x3v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_2x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_2x3v_ser_p2_varcoeff },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfx_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_vlasov_boundary_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfy_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order2_vlasov_boundary_surfy_2x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfy_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order2_vlasov_boundary_surfy_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p2_varcoeff },
      { dg_diffusion_order4_vlasov_boundary_surfy_2x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfy_2x3v_ser_p2_varcoeff },
//      { dg_diffusion_order4_vlasov_boundary_surfy_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_2x2v_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_2x3v_ser_p2_varcoeff },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfy_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_vlasov_boundary_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order2_vlasov_boundary_surfz_3x3v_ser_p1_varcoeff, dg_diffusion_order2_vlasov_boundary_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { dg_diffusion_order4_vlasov_boundary_surfz_3x3v_ser_p1_varcoeff, dg_diffusion_order4_vlasov_boundary_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
//      { NULL, dg_diffusion_order6_vlasov_boundary_surfz_3x3v_ser_p2_varcoeff },},
      { NULL, NULL },},
  },
};

