// Tue Aug 24 13:05:03 2021
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG
GKYL_CU_DH void eval_2d_tensor_p2(const double *z, double *b);
GKYL_CU_DH void flip_sign_2d_tensor_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void eval_3d_tensor_p2(const double *z, double *b);
GKYL_CU_DH void flip_sign_3d_tensor_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void eval_4d_tensor_p2(const double *z, double *b);
GKYL_CU_DH void flip_sign_4d_tensor_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void eval_5d_tensor_p2(const double *z, double *b);
GKYL_CU_DH void flip_sign_5d_tensor_p2(int dir, const double *f, double *fout );
EXTERN_C_END
