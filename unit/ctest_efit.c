#include <stdio.h>
#include <math.h>
#include <string.h>
#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_efit.h>
#include <gkyl_efit_priv.h>


void test_1(){
  char* filepath = "./efit_data/input.geqdsk";
  // RZ basis function

  int rzpoly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  struct gkyl_rect_grid rzgrid;
  struct gkyl_range rzlocal, rzlocal_ext;

  struct gkyl_efit* efit = gkyl_efit_new(filepath, &rzbasis, false);

  gkyl_rect_grid_init(&rzgrid, 2, efit->rzlower, efit->rzupper, efit->rzcells);
  gkyl_create_grid_ranges(&rzgrid, efit->rzghost, &rzlocal_ext, &rzlocal);


  struct gkyl_array* psizr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array* psibyrzr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array* psibyr2zr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);


  gkyl_efit_advance(efit, &rzgrid, &rzlocal, &rzlocal_ext, psizr, psibyrzr, psibyr2zr);

  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psizr, "efit_psi.gkyl");
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psibyrzr, "efit_psibyr.gkyl");
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psibyr2zr, "efit_psibyr2.gkyl");

  gkyl_array_release(psizr);
  gkyl_array_release(psibyrzr);
  gkyl_array_release(psibyr2zr);




}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
