#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


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
#include <gkyl_array_rio.h>
#include <gkyl_calc_cart_bmag.h>
#include <gkyl_calc_cart_bmag_priv.h>


void test_1()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  char* filepath = "./data/eqdsk/wham.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_calc_cart_bmag* bcalculator = gkyl_calc_cart_bmag_new(efit->rzbasis, efit->rzgrid, efit->rzlocal, efit->rzlocal_ext, false);
  gkyl_calc_cart_bmag_advance(bcalculator, efit->psizr, efit->psibyrzr, efit->psibyr2zr);
  gkyl_grid_sub_array_write(bcalculator->grid, &bcalculator->local, NULL,  bcalculator->bmag_rz, "cartbmag.gkyl");
  gkyl_grid_sub_array_write(bcalculator->grid, &bcalculator->local, NULL,  bcalculator->br_rz, "cartbr.gkyl");
  gkyl_grid_sub_array_write(bcalculator->grid, &bcalculator->local, NULL,  bcalculator->bz_rz, "cartbz.gkyl");


  double xcart[3] = {0.01, 0.01732, 1.0};
  double bcart[3] = {0.0};
  gkyl_eval_cart_bmag(bcalculator, xcart, bcart);
  double bmag = sqrt(bcart[0]*bcart[0] +  bcart[1]*bcart[1] + bcart[2]*bcart[2]);
  double R = sqrt(xcart[0]*xcart[0] +  xcart[1]*xcart[1]);
  printf("R = %g, Z = %g, bmag = %g\n", R, xcart[2], bmag);

  gkyl_calc_cart_bmag_release(bcalculator);

}


TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
