#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void skin_surf_from_ghost_lower_2x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upper_2x_ser_p1(const double *fghost, double *fskin);

GKYL_CU_DH void skin_surf_from_ghost_lower_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upper_3x_ser_p1(const double *fghost, double *fskin);


EXTERN_C_END