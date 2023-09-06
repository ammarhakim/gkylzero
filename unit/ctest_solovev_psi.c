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
#include <gkyl_geo_gyrokinetic.h>


#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_kernels.h>
#include <gkyl_geo_gyrokinetic.h>



struct mapc2p_ctx{
   struct gkyl_geo_gyrokinetic* app;
   struct gkyl_geo_gyrokinetic_geo_inp* ginp;
};

struct solovev_ctx {
  double B0, R0, k, q0, Ztop;
};

static inline double sq(double x) { return x*x; }



static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}




void
psi(double t, const double *xn, double *fout, void *ctx)
{
  struct solovev_ctx *s = ctx;
  double B0 = s->B0, R0 = s->R0, k = s->k, q0 = s->q0;
  double R = xn[0], Z = xn[1];
  fout[0] = B0*k/(2*sq(R0)*q0)*(sq(R)*sq(Z)/sq(k) + sq(sq(R) - sq(R0))/4);
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


void mapc2p(double t, const double *xn, double* fout, void *ctx)
{
  struct mapc2p_ctx *gc = (struct mapc2p_ctx*) ctx;
  //double RZ[2];
  gkyl_geo_gyrokinetic_mapc2p(gc->app, gc->ginp, xn, fout);
}

void
test_1()
{
  struct solovev_ctx sctx = {
    .B0 = 0.55, .R0 = 0.85, .k = 2, .q0 = 2, .Ztop = 1.5
  };

  double psi_sep = sctx.B0*sctx.k*sq(sctx.R0)/(8*sctx.q0);
  printf("psi_sep = %lg\n", psi_sep);
  

  // create RZ grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
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
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, &psi, &sctx);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  //gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psiRZ, "test_bmag_psi.gkyl");


  gkyl_geo_gyrokinetic *geo = gkyl_geo_gyrokinetic_new(&(struct gkyl_geo_gyrokinetic_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal,

      .quad_param = {  .eps = 1e-12 }
    }
  );

  //gkyl geo_gyrokinetic *app = geo_gyrokinetic_app_new(&inp);

  //basic_root_test(&inp, app);

  // Computational grid: theta X psi X alpha (only 2D for now)
  //double clower[] = { -M_PI + 1e-10, 0.05 };
  //double cupper[] = { M_PI, 0.1 };
  double clower[] = { -M_PI/4, 0.05 };
  double cupper[] = { M_PI/4, 0.1 };
  int ccells[] = { 16, 10 };


  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  //int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&cgrid, nghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);

  struct gkyl_geo_gyrokinetic_geo_inp ginp = {
    .cgrid = &cgrid,
    .cbasis = &cbasis,
    .ftype = GKYL_SOL_DN,
    .rclose = upper[0],
    .zmin = lower[1],
    .zmax = upper[1],
  
    //.write_node_coord_array = true,
    //.node_file_nm = "test_bmag.gkyl"
  }; 


  //check node coords
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ginp.cgrid->ndim, ginp.cbasis->num_basis);
  ginp.cbasis->node_list(gkyl_array_fetch(nodes, 0));
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  struct gkyl_range_iter iter_c;
  gkyl_range_iter_init(&iter_c, &clocal);
  
  while (gkyl_range_iter_next(&iter_c)) {
    gkyl_rect_grid_cell_center(ginp.cgrid, iter_c.idx, xc);

    for (int i=0; i<ginp.cbasis->num_basis; ++i) {
      comp_to_phys(ginp.cgrid->ndim, gkyl_array_cfetch(nodes, i),
        ginp.cgrid->dx, xc, xmu);
      //printf("xc = %g,%g ; xmu = %g,%g\n",xc[0],xc[1],xmu[0],xmu[1]);
      //printf("xmu = %g,%g\n", xmu[0], xmu[1]);
      //up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }
  }



  //Do Ammar's calcgeom

  struct gkyl_array *mapc2p_arr = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
  gkyl_geo_gyrokinetic_calcgeom(geo, &ginp, mapc2p_arr);

  //make psi
  gkyl_eval_on_nodes *eval_psi = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psi, &sctx);
  struct gkyl_array* psidg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psi, 0.0, &rzlocal_ext, psidg); //on ghosts with ext_range
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psidg, "solovev_psi.gkyl");

  gkyl_eval_on_nodes *eval_psibyr = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr, &sctx);
  struct gkyl_array* psibyrdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr, 0.0, &rzlocal_ext, psibyrdg); //on ghosts with ext_range
                                                                       //
  gkyl_eval_on_nodes *eval_psibyr2 = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr2, &sctx);
  struct gkyl_array* psibyr2dg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr2, 0.0, &rzlocal_ext, psibyr2dg); //on ghosts with ext_range
  printf("made the psi arays\n");

  //make bmag
  struct gkyl_array* bmagdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);


  struct gkyl_array* bmag_compdg= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);

  printf("calculating bmag \n");
  gkyl_calc_bmag *calculator = gkyl_calc_bmag_new(&cbasis, &rzbasis, &cgrid, &rzgrid, geo, &ginp, false);
  gkyl_calc_bmag_advance(calculator, &clocal, &clocal_ext, &rzlocal, &rzlocal_ext, psidg, psibyrdg, psibyr2dg, bmag_compdg);

  char fileNm[1024];
  do{
    printf("writing the comp bmag file \n");
    const char *fmt = "%s_compbmag.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "solovev_psi");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmag_compdg, fileNm);
  } while (0);


  //// write mapc2p and get metrics
  printf("calculating mapc2p \n");
  struct mapc2p_ctx *mctx = gkyl_malloc(sizeof(*mctx));
  mctx->app = geo;
  mctx->ginp = &ginp;
  gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&cgrid, &cbasis, 2, mapc2p, mctx);
  struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, clocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &clocal_ext, XYZ);

  printf("writing rz file\n");
  gkyl_grid_sub_array_write(&cgrid, &clocal, XYZ, "solovev_rzfile.gkyl");
  printf("wrote rz file\n");
  
  printf("calculating metrics \n");
  gkyl_calc_metric *mcalculator = gkyl_calc_metric_new(&cbasis, &cgrid, false);
  gkyl_calc_metric_advance( mcalculator, &clocal, XYZ, gFld);
  do{
    printf("writing the gij file \n");
    const char *fmt = "%s_gij.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "solovev_psi");
    gkyl_grid_sub_array_write(&cgrid, &clocal, gFld, fileNm);
  } while (0);

  // done doing metrics






  //gkyl_array_release(RZ);
  gkyl_array_release(psidg);
  gkyl_array_release(psibyrdg);
  gkyl_array_release(psibyr2dg);
  gkyl_array_release(bmag_compdg);
  gkyl_array_release(bmagdg);
  gkyl_geo_gyrokinetic_release(geo);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
