#include <gkyl_calc_metric.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_metric_kernels.h>
#include <assert.h>
#include <math.h>

typedef void (*metric_kernel)(const double **xyz, double *gFld);



typedef struct { metric_kernel kernels[27]; } metric_kernel_loc_list;
typedef struct { metric_kernel_loc_list list[3]; } metric_kernel_bcz_list;
typedef struct { metric_kernel_bcz_list list[2]; } metric_kernel_bcy_list;
typedef struct { metric_kernel_bcy_list list[2]; } metric_kernel_bcx_list;

GKYL_CU_DH
static const metric_kernel_bcx_list ser_metric_kernel_loc_list[] = {
// periodic x
{ .list = {
    // periodicy
    { .list = {
      // periodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_upz_periodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      // nonperiodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_upy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_upy_periodicy_upz_nonperiodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      }
    },

    // nonperiodicy
    { .list = {
      // periodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_upz_periodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      // nonperiodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_periodicx_upy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_periodicx_upy_nonperiodicy_upz_nonperiodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      }
    },


    }
  },
// nonperiodic x
{ .list = {
    // periodicy
    { .list = {
      // periodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_upz_periodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      // nonperiodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_periodicy_upz_nonperiodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      }
    },

    // nonperiodicy
    { .list = {
      // periodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_upz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_loz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_upz_periodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      // nonperiodic z
      {.list = {
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               { gij_3x_Ser_p1_inx_periodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_inz_periodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_iny_periodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_inx_periodicx_upy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_lox_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_loy_nonperiodicy_upz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_loz_nonperiodicz, gij_3x_Ser_p1_upx_nonperiodicx_upy_nonperiodicy_upz_nonperiodicz },
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL, NULL, NULL},
               },
      },
      }
    },


    }
  },



};





//// Boundary condition types.
//enum gkyl_metric_bc_type {
//  GKYL_metric_PERIODIC = 0,
//  GKYL_metric_NONPERIODIC, // sets the value.
//};
//
//
//struct gkyl_metric_bc {
//  enum gkyl_metric_bc_type type[GKYL_MAX_CDIM];
//};

struct gkyl_calc_metric {
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  struct gkyl_rect_grid* grid;
  bool use_gpu;
  metric_kernel_loc_list kernels;
  int *num_cells;
  int *bcs;
};



GKYL_CU_DH
static inline int idx_to_inloup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[d] == num_cells[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}


GKYL_CU_DH
//static metric_kernel_list metric_choose_kernel(int dim, int basis_type, int poly_order, int lidx)
//static metric_kernel metric_choose_kernel(int lidx)
static metric_kernel_loc_list metric_choose_kernel(int dim, int poly_order, int *bcs)
{
  return ser_metric_kernel_loc_list[bcs[0]].list[bcs[1]].list[bcs[2]].list[poly_order];
  //switch (basis_type) {
  //  case GKYL_BASIS_MODAL_SERENDIPITY:
  //    //return ser_metric_kernel_list[dim].kernels[poly_order];
  //    return ser_metric_kernel_list[lidx];
  //  default:
  //    assert(false);
  //    break;
  //}
}


