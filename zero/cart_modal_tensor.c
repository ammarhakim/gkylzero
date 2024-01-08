#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>
#include <gkyl_basis_ser_kernels.h>
#include <gkyl_basis_tensor_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_1d_ser_p0, eval_1d_ser_p1, eval_1d_ser_p2, eval_1d_ser_p3 },
  { eval_2d_ser_p0, eval_2d_ser_p1, eval_2d_tensor_p2, eval_2d_tensor_p3 },
  { eval_3d_ser_p0, eval_3d_ser_p1, eval_3d_tensor_p2, NULL },
  { eval_4d_ser_p0, eval_4d_ser_p1, eval_4d_tensor_p2, NULL },
  { eval_5d_ser_p0, eval_5d_ser_p1, eval_5d_tensor_p2, NULL },
  { eval_6d_ser_p0, eval_6d_ser_p1, NULL, NULL },
};

// Expansion eval for each dimension: eve_list[ndim].ev[poly_order]
static struct { double (*ev[4])(const double *z, const double *f); } eve_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_expand_1d_ser_p0, eval_expand_1d_ser_p1, eval_expand_1d_ser_p2, eval_expand_1d_ser_p3 },
  { eval_expand_2d_ser_p0, eval_expand_2d_ser_p1, eval_expand_2d_tensor_p2,  eval_expand_2d_tensor_p3 },
  { eval_expand_3d_ser_p0, eval_expand_3d_ser_p1, eval_expand_3d_tensor_p2, NULL },
  { eval_expand_4d_ser_p0, eval_expand_4d_ser_p1, eval_expand_4d_tensor_p2, NULL },
  { eval_expand_5d_ser_p0, eval_expand_5d_ser_p1, eval_expand_5d_tensor_p2, NULL },
  { eval_expand_6d_ser_p0, eval_expand_6d_ser_p1, NULL, NULL },
};

// Expansion eval_grad for each dimension: eveg_list[ndim].ev[poly_order]
static struct { double (*ev[4])(int dir, const double *z, const double *f); } eveg_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_grad_expand_1d_ser_p0, eval_grad_expand_1d_ser_p1, eval_grad_expand_1d_ser_p2, eval_grad_expand_1d_ser_p3 },
  { eval_grad_expand_2d_ser_p0, eval_grad_expand_2d_ser_p1, eval_grad_expand_2d_tensor_p2, eval_grad_expand_2d_tensor_p3 },
  { eval_grad_expand_3d_ser_p0, eval_grad_expand_3d_ser_p1, eval_grad_expand_3d_tensor_p2, NULL },
  { eval_grad_expand_4d_ser_p0, eval_grad_expand_4d_ser_p1, eval_grad_expand_4d_tensor_p2, NULL },
  { eval_grad_expand_5d_ser_p0, eval_grad_expand_5d_ser_p1, eval_grad_expand_5d_tensor_p2, NULL },
  { eval_grad_expand_6d_ser_p0, eval_grad_expand_6d_ser_p1, NULL, NULL },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fos_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { flip_odd_sign_1d_ser_p0, flip_odd_sign_1d_ser_p1, flip_odd_sign_1d_ser_p2, flip_odd_sign_1d_ser_p3 },
  { flip_odd_sign_2d_ser_p0, flip_odd_sign_2d_ser_p1, flip_odd_sign_2d_tensor_p2, flip_odd_sign_2d_tensor_p3 },
  { flip_odd_sign_3d_ser_p0, flip_odd_sign_3d_ser_p1, flip_odd_sign_3d_tensor_p2, NULL },
  { flip_odd_sign_4d_ser_p0, flip_odd_sign_4d_ser_p1, flip_odd_sign_4d_tensor_p2, NULL },
  { flip_odd_sign_5d_ser_p0, flip_odd_sign_5d_ser_p1, flip_odd_sign_5d_tensor_p2, NULL },
  { flip_odd_sign_6d_ser_p0, flip_odd_sign_6d_ser_p1, NULL, NULL },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fes_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { flip_even_sign_1d_ser_p0, flip_even_sign_1d_ser_p1, flip_even_sign_1d_ser_p2, flip_even_sign_1d_ser_p3 },
  { flip_even_sign_2d_ser_p0, flip_even_sign_2d_ser_p1, flip_even_sign_2d_tensor_p2, flip_even_sign_2d_tensor_p3 },
  { flip_even_sign_3d_ser_p0, flip_even_sign_3d_ser_p1, flip_even_sign_3d_tensor_p2, NULL },
  { flip_even_sign_4d_ser_p0, flip_even_sign_4d_ser_p1, flip_even_sign_4d_tensor_p2, NULL },
  { flip_even_sign_5d_ser_p0, flip_even_sign_5d_ser_p1, flip_even_sign_5d_tensor_p2, NULL },
  { flip_even_sign_6d_ser_p0, flip_even_sign_6d_ser_p1, NULL, NULL },
};

// Node list function: ev_list[ndim].ev[poly_order]
static struct { void (*nl[4])(double * node_list); } nl_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { node_coords_1d_ser_p0, node_coords_1d_ser_p1, node_coords_1d_ser_p2, NULL },
  { node_coords_2d_ser_p0, node_coords_2d_ser_p1, node_coords_2d_tensor_p2, NULL },
  { node_coords_3d_ser_p0, node_coords_3d_ser_p1, node_coords_3d_tensor_p2, NULL },
  { node_coords_4d_ser_p0, node_coords_4d_ser_p1, node_coords_4d_tensor_p2, NULL },
  { node_coords_5d_ser_p0, node_coords_5d_ser_p1, node_coords_5d_tensor_p2, NULL },
  { node_coords_6d_ser_p0, node_coords_6d_ser_p1, NULL, NULL },
};

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, nodal_to_modal_1d_ser_p1, nodal_to_modal_1d_ser_p2, NULL },
  { NULL, nodal_to_modal_2d_ser_p1, nodal_to_modal_2d_tensor_p2, NULL },
  { NULL, nodal_to_modal_3d_ser_p1, nodal_to_modal_3d_tensor_p2, NULL },
  { NULL, nodal_to_modal_4d_ser_p1, nodal_to_modal_4d_tensor_p2, NULL },
  { NULL, NULL, NULL, NULL }, // TODO
  { NULL, NULL, NULL, NULL }, // TODO
};

void
gkyl_cart_modal_tensor(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ndim>0 && ndim<=6);
  assert(ev_list[ndim].ev[poly_order]);
  
  basis->ndim = ndim;
  basis->poly_order = poly_order;
  basis->num_basis = pow(poly_order+1, ndim);
  strcpy(basis->id, "tensor");
  basis->b_type = GKYL_BASIS_MODAL_TENSOR;

  // function pointers
  basis->eval = ev_list[ndim].ev[poly_order];
  basis->eval_expand = eve_list[ndim].ev[poly_order];
  basis->eval_grad_expand = eveg_list[ndim].ev[poly_order];
  basis->flip_odd_sign = fos_list[ndim].fs[poly_order];
  basis->flip_even_sign = fes_list[ndim].fs[poly_order];
  basis->node_list = nl_list[ndim].nl[poly_order];
  basis->nodal_to_modal = n2m_list[ndim].n2m[poly_order];
}

struct gkyl_basis *
gkyl_cart_modal_tensor_new(int ndim, int poly_order)
{
  struct gkyl_basis *basis = gkyl_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_tensor(basis, ndim, poly_order);
  return basis;
}

#ifndef GKYL_HAVE_CUDA
void
gkyl_cart_modal_tensor_cu_dev(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(false);
}

struct gkyl_basis *
gkyl_cart_modal_tensor_cu_dev_new(int ndim, int poly_order)
{
  assert(false);
}
#endif
