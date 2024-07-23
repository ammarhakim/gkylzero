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


void test_solovev(){
  char* filepath = "./data/eqdsk/solovev.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order, GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "solovev_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "solovev_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "solovev_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "solovev_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "solovev_q.gkyl");

  gkyl_efit_release(efit);

}

void test_step(){
  char* filepath = "./data/eqdsk/step.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_TENSOR, fluxpoly_order, true, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "step_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "step_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "step_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "step_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "step_q.gkyl");

  gkyl_efit_release(efit);

}

void test_asdex(){
  char* filepath = "./data/eqdsk/asdex.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "asdex_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "asdex_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "asdex_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "asdex_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "asdex_q.gkyl");

  gkyl_efit_release(efit);

}

void test_cerfon(){
  char* filepath = "./data/eqdsk/cerfon.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "cerfon_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "cerfon_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "cerfon_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "cerfon_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "cerfon_q.gkyl");

  gkyl_efit_release(efit);

}

void test_elliptical(){
  char* filepath = "./data/eqdsk/elliptical.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "elliptical_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "elliptical_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "elliptical_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "elliptical_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "elliptical_q.gkyl");

  gkyl_efit_release(efit);

}

void test_wham(){
  char* filepath = "./data/eqdsk/wham.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "wham_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "wham_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "wham_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "wham_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "wham_q.gkyl");

  gkyl_efit_release(efit);

}


void test_tcv(){
  char* filepath = "./data/eqdsk/tcv.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "tcv_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "tcv_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "tcv_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "tcv_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "tcv_q.gkyl");

  gkyl_efit_release(efit);

}

void test_mast(){
  char* filepath = "./data/eqdsk/mast.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "mast_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "mast_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "mast_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "mast_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "mast_q.gkyl");

  gkyl_efit_release(efit);

}

TEST_LIST = {
  { "test_solovev", test_solovev},
  { "test_step", test_step},
  { "test_asdex", test_asdex},
  { "test_cerfon", test_cerfon},
  { "test_elliptical", test_elliptical},
  { "test_wham", test_wham},
  { "test_tcv", test_tcv},
  { "test_mast", test_mast},
  { NULL, NULL },
};
