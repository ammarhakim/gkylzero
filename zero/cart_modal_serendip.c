#include <assert.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_ser_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_1d_ser_p0, eval_1d_ser_p1, eval_1d_ser_p2, eval_1d_ser_p3 },
  { eval_2d_ser_p0, eval_2d_ser_p1, eval_2d_ser_p2, eval_2d_ser_p3 },
  { eval_3d_ser_p0, eval_3d_ser_p1, eval_3d_ser_p2, eval_3d_ser_p3 },
  { eval_4d_ser_p0, eval_4d_ser_p1, eval_4d_ser_p2, eval_4d_ser_p3 },
  { eval_5d_ser_p0, eval_5d_ser_p1, eval_5d_ser_p2, NULL },
  { eval_6d_ser_p0, eval_6d_ser_p1, NULL, NULL },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fs_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { flip_sign_1d_ser_p0, flip_sign_1d_ser_p1, flip_sign_1d_ser_p2, flip_sign_1d_ser_p3 },
  { flip_sign_2d_ser_p0, flip_sign_2d_ser_p1, flip_sign_2d_ser_p2, flip_sign_2d_ser_p3 },
  { flip_sign_3d_ser_p0, flip_sign_3d_ser_p1, flip_sign_3d_ser_p2, flip_sign_3d_ser_p3 },
  { flip_sign_4d_ser_p0, flip_sign_4d_ser_p1, flip_sign_4d_ser_p2, flip_sign_4d_ser_p3 },
  { flip_sign_5d_ser_p0, flip_sign_5d_ser_p1, flip_sign_5d_ser_p2, NULL },
  { flip_sign_6d_ser_p0, flip_sign_6d_ser_p1, NULL, NULL },
};

// Number of basis functions: num_basis_list[ndim].count[poly_order]
static struct { int count[4]; } num_basis_list[] = {
  { 1, 1, 1, 1 },
  { 1, 2, 3, 4 },
  { 1, 4, 8, 12 },
  { 1, 8, 20, 32 },
  { 1, 16, 48, 80 },
  { 1, 32, 112, 192 },
  { 1, 64, 256, 0 },
};

// Node list function: ev_list[ndim].ev[poly_order]
static struct { void (*nl[4])(double * node_list); } nl_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { node_coords_1d_ser_p0, node_coords_1d_ser_p1, node_coords_1d_ser_p2, NULL },
  { node_coords_2d_ser_p0, node_coords_2d_ser_p1, node_coords_2d_ser_p2, NULL },
  { node_coords_3d_ser_p0, node_coords_3d_ser_p1, node_coords_3d_ser_p2, NULL },
  { node_coords_4d_ser_p0, node_coords_4d_ser_p1, node_coords_4d_ser_p2, NULL },
  { node_coords_5d_ser_p0, node_coords_5d_ser_p1, node_coords_5d_ser_p2, NULL },
  { node_coords_6d_ser_p0, node_coords_6d_ser_p1, NULL, NULL },
};

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, nodal_to_modal_1d_ser_p1, nodal_to_modal_1d_ser_p2, NULL },
  { NULL, nodal_to_modal_2d_ser_p1, nodal_to_modal_2d_ser_p2, NULL },
  { NULL, nodal_to_modal_3d_ser_p1, nodal_to_modal_3d_ser_p2, NULL },
  { NULL, nodal_to_modal_4d_ser_p1, nodal_to_modal_4d_ser_p2, NULL },
  { NULL, nodal_to_modal_5d_ser_p1, nodal_to_modal_5d_ser_p2, NULL },
  { NULL, nodal_to_modal_6d_ser_p1, NULL, NULL },
};

void
gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ndim>0 && ndim<=6);
  assert(ev_list[ndim].ev[poly_order]);
  
  basis->ndim = ndim;
  basis->poly_order = poly_order;
  basis->num_basis = num_basis_list[ndim].count[poly_order];
  strcpy(basis->id, "serendipity");
  basis->b_type = GKYL_BASIS_MODAL_SERENDIPITY;
  basis->eval = ev_list[ndim].ev[poly_order];
  basis->flip_sign = fs_list[ndim].fs[poly_order];
  basis->node_list = nl_list[ndim].nl[poly_order];
  basis->nodal_to_modal = n2m_list[ndim].n2m[poly_order];
}
