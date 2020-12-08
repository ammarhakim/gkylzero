#include <assert.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include "cart_modal_serendip_funcs.c"

// Basis function eval for each dimension: ev_list[ndim].ev[polyOrder]
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_1d_p0, eval_1d_p1, eval_1d_p2, eval_1d_p3 },
  { eval_2d_p0, eval_2d_p1, eval_2d_p2, eval_2d_p3 },
  { eval_3d_p0, eval_3d_p1, eval_3d_p2, eval_3d_p3 },
  { eval_4d_p0, eval_4d_p1, eval_4d_p2, eval_4d_p3 },
  { eval_5d_p0, eval_5d_p1, eval_5d_p2, NULL },
  { eval_6d_p0, eval_6d_p1, NULL, NULL },
};

// Number of basis functions: numBasis_list[ndim].count[polyOrder]
static struct { int count[4]; } numBasis_list[] = {
  { 1, 1, 1, 1 },
  { 1, 2, 3, 4 },
  { 1, 4, 8, 12 },
  { 1, 8, 20, 32 },
  { 1, 16, 48, 80 },
  { 1, 32, 112, 192 },
  { 1, 64, 256, 0 },
};

void
gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim, int polyOrder)
{
  assert(ndim>0 && ndim<=6);
  assert(ev_list[ndim].ev[polyOrder]);
  
  basis->ndim = ndim;
  basis->polyOrder = polyOrder;
  basis->numBasis = numBasis_list[ndim].count[polyOrder];
  strcpy(basis->id, "serendipity");
  basis->eval = ev_list[ndim].ev[polyOrder];
}
