#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_calc_bmag.h>
#include <gkyl_calc_bmag_kernels.h>
#include <gkyl_gkgeom.h>




void
psi(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  fout[0] = (R-2)*(R-2) + Z*Z/4;
}

void
psibyr(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  psi(t,xn,fout,ctx);
  fout[0] = fout[0]/R;
}

void
psibyr2(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  psi(t,xn,fout,ctx);
  fout[0] = fout[0]/R/R;
}

void
test_1()
{
  // create RZ grid
  double lower[] = { 0.5, -4.0 }, upper[] = { 6.0, 4.0 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 64, 128 };


  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rzpoly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  // allocate psiRZ array, initialize and write it to file
  struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid,
    &rzbasis, 1, &psi, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  //gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psiRZ, "test_bmag_psi.gkyl");


  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal,

      .quad_param = {  .eps = 1e-12 }
    }
  );

  //gkyl gkgeom *app = gkgeom_app_new(&inp);

  //basic_root_test(&inp, app);

  // Computational grid: theta X psi X alpha (only 2D for now)
  double clower[] = { -M_PI/2, 7.0 };
  double cupper[] = { M_PI/2, 12.0 };
  int ccells[] = { 16, 8 };


  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  //int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&cgrid, nghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);

  struct gkyl_gkgeom_geo_inp ginp = {
    .cgrid = &cgrid,
    .cbasis = &cbasis,
    .ftype = GKYL_SOL_DN,
    .rclose = upper[0],
    .zmin = lower[1],
    .zmax = upper[1],
  
    //.write_node_coord_array = true,
    //.node_file_nm = "test_bmag.gkyl"
  }; 

  //make psi
  gkyl_eval_on_nodes *eval_psi = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psi, 0);
  struct gkyl_array* psidg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psi, 0.0, &rzlocal_ext, psidg); //on ghosts with ext_range
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psidg, "test_bmag_psi.gkyl");

  gkyl_eval_on_nodes *eval_psibyr = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr, 0);
  struct gkyl_array* psibyrdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr, 0.0, &rzlocal_ext, psibyrdg); //on ghosts with ext_range
                                                                       //
  gkyl_eval_on_nodes *eval_psibyr2 = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr2, 0);
  struct gkyl_array* psibyr2dg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr2, 0.0, &rzlocal_ext, psibyr2dg); //on ghosts with ext_range
  //make bmag
  struct gkyl_array* bmagdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);


  struct gkyl_array* bmag_compdg= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);

  gkyl_calc_bmag *calculator = gkyl_calc_bmag_new(&cbasis, &rzbasis, &cgrid, &rzgrid, geo, &ginp, false);
  gkyl_calc_bmag_advance(calculator, &clocal, &clocal_ext, &rzlocal, &rzlocal_ext, psidg, psibyrdg, psibyr2dg, bmag_compdg);

  char fileNm[1024];
  do{
    printf("writing the comp bmag file \n");
    const char *fmt = "%s_compbmag.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "test3");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmag_compdg, fileNm);
  } while (0);





  //gkyl_array_release(RZ);
  gkyl_array_release(psidg);
  gkyl_array_release(psibyrdg);
  gkyl_array_release(psibyr2dg);
  gkyl_array_release(bmag_compdg);
  gkyl_array_release(bmagdg);
  gkyl_gkgeom_release(geo);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
