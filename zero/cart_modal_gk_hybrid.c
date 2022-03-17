#include <assert.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_gk_hyb_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_2d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_3d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, eval_5d_gk_hyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // No 6D basis functions
};

// Expansion eval for each dimension: eve_list[ndim].ev[poly_order]
static struct { double (*ev[4])(const double *z, const double *f); } eve_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_expand_2d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_expand_3d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_expand_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, eval_expand_5d_gk_hyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // No 6D basis functions
};

// Expansion eval_grad for each dimension: eveg_list[ndim].ev[poly_order]
static struct { double (*ev[4])(int dir, const double *z, const double *f); } eveg_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_grad_expand_2d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_grad_expand_3d_gk_hyb_p1, NULL, NULL },
  { NULL, eval_grad_expand_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, eval_grad_expand_5d_gk_hyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // No 6D basis functions
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fos_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, flip_odd_sign_2d_gk_hyb_p1, NULL, NULL },
  { NULL, flip_odd_sign_3d_gk_hyb_p1, NULL, NULL },
  { NULL, flip_odd_sign_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, flip_odd_sign_5d_gk_hyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // No 6D basis functions  
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fes_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, flip_even_sign_2d_gk_hyb_p1, NULL, NULL },
  { NULL, flip_even_sign_3d_gk_hyb_p1, NULL, NULL },
  { NULL, flip_even_sign_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, flip_even_sign_5d_gk_hyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // No 6D basis functions  
};

// Number of basis functions: num_basis_list[ndim].count[poly_order]
static struct { int count[4]; } num_basis_list[] = {
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 6, 0, 0 },
  { 0, 12, 0, 0 },
  { 0, 24, 0, 0 },
  { 0, 48, 0, 0 },
  { 0, 0, 0, 0 },
};

// Node list function: ev_list[ndim].ev[poly_order]
static struct { void (*nl[4])(double * node_list); } nl_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, node_coords_2d_gk_hyb_p1, NULL, NULL },
  { NULL, node_coords_3d_gk_hyb_p1, NULL, NULL },
  { NULL, node_coords_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, node_coords_5d_gk_hyb_p1, NULL, NULL },  
};

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, nodal_to_modal_2d_gk_hyb_p1, NULL, NULL },
  { NULL, nodal_to_modal_3d_gk_hyb_p1, NULL, NULL },
  { NULL, nodal_to_modal_4d_gk_hyb_p1, NULL, NULL },     
  { NULL, nodal_to_modal_5d_gk_hyb_p1, NULL, NULL },
};

void
gkyl_cart_modal_gk_hybrid(struct gkyl_basis *basis, int ndim)
{
  assert(ndim>1 && ndim<6);
  
  basis->ndim = ndim;
  basis->poly_order = 1;
  basis->num_basis = num_basis_list[ndim].count[1];
  strcpy(basis->id, "gk_hybrid");
  basis->b_type = GKYL_BASIS_MODAL_GK_HYBRID;

  // function pointers
  basis->eval = ev_list[ndim].ev[1];
  basis->eval_expand = eve_list[ndim].ev[1];
  basis->eval_grad_expand = eveg_list[ndim].ev[1];
  basis->flip_odd_sign = fos_list[ndim].fs[1];
  basis->flip_even_sign = fes_list[ndim].fs[1];
  basis->node_list = NULL; // TODO nl_list[ndim].nl[1];
  basis->nodal_to_modal = NULL; //TODO n2m_list[ndim].n2m[1];
}
