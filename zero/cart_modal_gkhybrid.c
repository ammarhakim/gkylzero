#include <assert.h>
#include <string.h>
#include <gkyl_util.h>

#include "cart_modal_gkhybrid_priv.h"

void
gkyl_cart_modal_gkhybrid(struct gkyl_basis *basis, int cdim, int vdim)
{
  int ndim = cdim+vdim;
  assert(ndim>1 && ndim<6);
  assert(cdim<4 && vdim>0 && vdim<3);
  
  basis->ndim = ndim;
  basis->poly_order = 1;
  basis->num_basis = num_basis_list[ndim].count[1];
  strcpy(basis->id, "gkhybrid");
  basis->b_type = GKYL_BASIS_MODAL_GKHYBRID;

  // function pointers
  basis->eval = ev_list[ndim].ev[1];
  basis->eval_expand = eve_list[ndim].ev[1];
  basis->eval_grad_expand = eveg_list[ndim].ev[1];
  basis->flip_odd_sign = fos_list[ndim].fs[1];
  basis->flip_even_sign = fes_list[ndim].fs[1];
  basis->node_list = nl_list[ndim].nl[1];
  basis->nodal_to_modal = n2m_list[ndim].n2m[1];
}

#ifndef GKYL_HAVE_CUDA
void
gkyl_cart_modal_gkhybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim)
{
  assert(false);
}
#endif
