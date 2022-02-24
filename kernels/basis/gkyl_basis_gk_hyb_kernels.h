// Wed Feb 23 22:22:36 2022
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG
GKYL_CU_DH void eval_2d_hyb_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_hyb_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_hyb_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_2d_hyb_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_hyb_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_hyb_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_3d_hyb_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_3d_hyb_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_3d_hyb_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_3d_hyb_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_3d_hyb_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_3d_hyb_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_4d_hyb_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_4d_hyb_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_4d_hyb_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_4d_hyb_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_4d_hyb_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_4d_hyb_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_5d_hyb_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_5d_hyb_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_5d_hyb_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_5d_hyb_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_5d_hyb_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_5d_hyb_p1(const double *fnodal, double *fmodal);
EXTERN_C_END
