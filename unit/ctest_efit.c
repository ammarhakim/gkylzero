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
  char* filepath = "/home/akash/test_cio/AH_PI4_P5.geqdsk";
  // RZ basis function

  int rzpoly_order = 1;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  struct gkyl_rect_grid rzgrid;
 
  struct gkyl_range rzlocal, rzlocal_ext;

  struct gkyl_efit* efit = gkyl_efit_new(filepath, &rzbasis, &rzgrid, &rzlocal, &rzlocal_ext, false);

  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, efit->psizr, "efit_psi.gkyl");




}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
