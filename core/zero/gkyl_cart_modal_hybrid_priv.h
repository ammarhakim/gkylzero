#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_hyb_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, eval_1x1v_hyb_p1, eval_1x2v_hyb_p1, eval_1x3v_hyb_p1,  },
  { NULL, eval_2x1v_hyb_p1, eval_2x2v_hyb_p1, eval_2x3v_hyb_p1,  },
  { NULL, eval_3x1v_hyb_p1, eval_3x2v_hyb_p1, eval_3x3v_hyb_p1,  },
};

// Expansion eval for each dimension: eve_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { double (*ev[4])(const double *z, const double *f); } eve_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, eval_expand_1x1v_hyb_p1, eval_expand_1x2v_hyb_p1, eval_expand_1x3v_hyb_p1 },
  { NULL, eval_expand_2x1v_hyb_p1, eval_expand_2x2v_hyb_p1, eval_expand_2x3v_hyb_p1 },
  { NULL, eval_expand_3x1v_hyb_p1, eval_expand_3x2v_hyb_p1, eval_expand_3x3v_hyb_p1 },
};

// Expansion eval_grad for each dimension: eveg_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { double (*ev[4])(int dir, const double *z, const double *f); } eveg_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, eval_grad_expand_1x1v_hyb_p1, eval_grad_expand_1x2v_hyb_p1, eval_grad_expand_1x3v_hyb_p1 },
  { NULL, eval_grad_expand_2x1v_hyb_p1, eval_grad_expand_2x2v_hyb_p1, eval_grad_expand_2x3v_hyb_p1 },
  { NULL, eval_grad_expand_3x1v_hyb_p1, eval_grad_expand_3x2v_hyb_p1, eval_grad_expand_3x3v_hyb_p1 },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fos_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, flip_odd_sign_1x1v_hyb_p1 , flip_odd_sign_1x2v_hyb_p1 , flip_odd_sign_1x3v_hyb_p1 },
  { NULL, flip_odd_sign_2x1v_hyb_p1 , flip_odd_sign_2x2v_hyb_p1 , flip_odd_sign_2x3v_hyb_p1 },
  { NULL, flip_odd_sign_3x1v_hyb_p1 , flip_odd_sign_3x2v_hyb_p1 , flip_odd_sign_3x3v_hyb_p1 },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fes_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, flip_even_sign_1x1v_hyb_p1, flip_even_sign_1x2v_hyb_p1, flip_even_sign_1x3v_hyb_p1 },
  { NULL, flip_even_sign_2x1v_hyb_p1, flip_even_sign_2x2v_hyb_p1, flip_even_sign_2x3v_hyb_p1 },
  { NULL, flip_even_sign_3x1v_hyb_p1, flip_even_sign_3x2v_hyb_p1, flip_even_sign_3x3v_hyb_p1 },
};

// Number of basis functions: num_basis_list[ndim].count[poly_order]
GKYL_CU_D
static struct { int count[4]; } num_basis_list[] = {
  { 0,  0,  0,  0 },  // No 0x basis functions.
  { 0,  6, 16, 40 },
  { 0, 12, 32, 80 },
  { 0, 24, 64, 160 },
};

// Node list function: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*nl[4])(double * node_list); } nl_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, node_coords_1x1v_hyb_p1, node_coords_1x2v_hyb_p1, node_coords_1x3v_hyb_p1 },
  { NULL, node_coords_2x1v_hyb_p1, node_coords_2x2v_hyb_p1, node_coords_2x3v_hyb_p1 },
  { NULL, node_coords_3x1v_hyb_p1, node_coords_3x2v_hyb_p1, node_coords_3x3v_hyb_p1 },
};

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, nodal_to_modal_1x1v_hyb_p1, nodal_to_modal_1x2v_hyb_p1, nodal_to_modal_1x3v_hyb_p1 },
  { NULL, nodal_to_modal_2x1v_hyb_p1, nodal_to_modal_2x2v_hyb_p1, nodal_to_modal_2x3v_hyb_p1 },
  { NULL, nodal_to_modal_3x1v_hyb_p1, nodal_to_modal_3x2v_hyb_p1, nodal_to_modal_3x3v_hyb_p1 },
};

// Gauss-Legendre quadrature nodes nodal basis -> modal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fquad, double *fmodal, long linc2); } qn2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, quad_to_modal_1x1v_hyb_p1, quad_to_modal_1x2v_hyb_p1, quad_to_modal_1x3v_hyb_p1 },
  { NULL, quad_to_modal_2x1v_hyb_p1, quad_to_modal_2x2v_hyb_p1, quad_to_modal_2x3v_hyb_p1 },
  { NULL, quad_to_modal_3x1v_hyb_p1, quad_to_modal_3x2v_hyb_p1, quad_to_modal_3x3v_hyb_p1 },
};

// modal basis -> Gauss-Legendre quadrature nodes nodal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fmodal, double *fquad, long linc2); } m2qn_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0x basis functions
  { NULL, modal_to_quad_1x1v_hyb_p1, modal_to_quad_1x2v_hyb_p1, modal_to_quad_1x3v_hyb_p1 },
  { NULL, modal_to_quad_2x1v_hyb_p1, modal_to_quad_2x2v_hyb_p1, modal_to_quad_2x3v_hyb_p1 },
  { NULL, modal_to_quad_3x1v_hyb_p1, modal_to_quad_3x2v_hyb_p1, modal_to_quad_3x3v_hyb_p1 },
};
