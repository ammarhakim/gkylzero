#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


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


#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_kernels.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_derived_geo_kernels.h>
#include <gkyl_gkgeom.h>



struct mapc2p_ctx{
   struct gkyl_gkgeom* app;
   struct gkyl_gkgeom_geo_inp* ginp;
};


// Cerfon equilibrium
struct cerfon_ctx {
  double R0, psi_prefactor, B0;
};

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }



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
  struct cerfon_ctx *s = ctx;
  double R0 = s->R0, psi_prefactor = s->psi_prefactor;
  double R = xn[0], Z = xn[1];
  double x = R/R0, y = Z/R0;

  fout[0] = psi_prefactor*(0.00373804283369699*hex(x)*log(x) - 0.00574955335438162*hex(x) - 0.0448565140043639*qad(x)*sq(y)*log(x) + 0.0503044260840946*qad(x)*sq(y) + 0.017623348727471*qad(x)*log(x) + 0.0956643504553683*qad(x) + 0.0299043426695759*sq(x)*qad(y)*log(x) - 0.0160920841654771*sq(x)*qad(y) - 0.0704933949098842*sq(x)*sq(y)*log(x) + 0.0644725519961135*sq(x)*sq(y) - 7.00898484784405e-5*sq(x)*log(x) - 0.303766642191745*sq(x) - 0.00199362284463839*hex(y) + 0.0117488991516474*qad(y) + 7.00898484784405e-5*sq(y) + 0.0145368720253975);
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
bphi_func(double t, const double *xn, double *fout, void *ctx)
{
  struct cerfon_ctx *s = ctx;
  double B0 = s->B0, R0 = s->R0;
  double R = xn[0];
  fout[0] = B0*R0/R;
}


void
test_1()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();


  struct cerfon_ctx sctx = {  .R0 = 2.5, .psi_prefactor = 1.0, .B0 = 0.55 };

  
  // create RZ grid
  double lower[] = { 0.01, -6.0 }, upper[] = { 6.0, 6.0 };
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


  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal,
      .B0 = sctx.B0,
      .R0 = sctx.R0,

      .quad_param = {  .eps = 1e-12 }
    }
  );

  // compute outboard SOL geometry
  int npsi = 32, ntheta = 32;
  double psi_min = 0.1, psi_max = 1.2;
  double dpsi = (psi_max-psi_min)/npsi;
  double dtheta = M_PI/ntheta;
  psi_min += dpsi;

  psi_min = 0.03636363636363636; // This gives ghost node on psisep for 32 cells
  //psi_min = 0.07058823529411765; // This gives ghost node on psisep for 16 cells
  printf("psimin = %g\n", psi_min);
  
  // Computational grid: theta X psi X alpha (only 2D for now)
  double clower[] = { psi_min, -0.01, -2.9 };
  double cupper[] = {psi_max, 0.01, 2.9 };
  int ccells[] = { 32,1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_gkgeom_geo_inp ginp = {
    .cgrid = &cgrid,
    .cbasis = &cbasis,
    .ftype = GKYL_SOL_DN,
    .rclose = upper[0],
    .zmin = lower[1],
    .zmax = upper[1],
  
    .write_node_coord_array = true,
    .node_file_nm = "cerfon3d_nodes.gkyl"
  }; 




  //Do Ammar's calcgeom

  struct gkyl_array *mapc2p_arr = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, clocal_ext.volume);
  gkyl_gkgeom_calcgeom(geo, &ginp, mapc2p_arr, &clocal_ext);


  printf("writing mapc2p file from calcgeom\n");
  gkyl_grid_sub_array_write(&cgrid, &clocal, mapc2p_arr, "cerfon3d_mapc2pfile.gkyl");
  printf("wrote mapc2p file\n");

  //make psi
  gkyl_eval_on_nodes *eval_psi = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psi, &sctx);
  struct gkyl_array* psidg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psi, 0.0, &rzlocal_ext, psidg); //on ghosts with ext_range
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psidg, "cerfon3d_psi.gkyl");

  gkyl_eval_on_nodes *eval_psibyr = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr, &sctx);
  struct gkyl_array* psibyrdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr, 0.0, &rzlocal_ext, psibyrdg); //on ghosts with ext_range
                                                                       //
  gkyl_eval_on_nodes *eval_psibyr2 = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr2, &sctx);
  struct gkyl_array* psibyr2dg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_psibyr2, 0.0, &rzlocal_ext, psibyr2dg); //on ghosts with ext_range
  printf("made the psi arays\n");

  gkyl_eval_on_nodes *eval_bphi= gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, bphi_func, &sctx);
  struct gkyl_array* bphidg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_eval_on_nodes_advance(eval_bphi, 0.0, &rzlocal_ext, bphidg); //on ghosts with ext_range
  printf("made the Bphi array\n");

  //make bmag
  struct gkyl_array* bmagdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);


  struct gkyl_array* bmagFld= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);

  printf("calculating bmag \n");
  gkyl_calc_bmag *calculator = gkyl_calc_bmag_new(&cbasis, &rzbasis, &cgrid, &rzgrid, geo, &ginp, false);
  printf("allocated bmag calculator\n");

  //printf("testing a fetch\n");
  //struct gkyl_range_iter citer;
  //gkyl_range_iter_init(&citer, &clocal_ext);
  //citer.idx[0] = 0;
  //citer.idx[1] = 0;
  //citer.idx[2] = 0;
  //printf("citer = %d,%d,%d\n", citer.idx[0], citer.idx[1], citer.idx[2]);

  //long lidx = gkyl_range_idx(&clocal_ext, citer.idx);
  //const double *mcoeffs = gkyl_array_cfetch(mapc2p_arr, lidx);
  //

  //printf("coeffs = ");
  //for(int i = 0; i < cgrid.ndim; i++){
  //  for(int j = 0; j < cbasis.num_basis; j++){
  //    const double *temp = &mcoeffs[i*cbasis.num_basis];
  //    printf(" %g ", temp[j] );
  //  }
  //}
  //printf("\n");
  //printf("done testing a fetch\n");


  gkyl_calc_bmag_advance(calculator, &clocal, &clocal_ext, &rzlocal, &rzlocal_ext, psidg, psibyrdg, psibyr2dg, bphidg, bmagFld, mapc2p_arr);
  printf("advanced bmag calculator\n");



  char fileNm[1024];
  do{
    printf("writing the comp bmag file \n");
    const char *fmt = "%s_compbmag.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmagFld, fileNm);
  } while (0);

  do{
    printf("writing the bphi file \n");
    const char *fmt = "%s_bphi.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
    gkyl_grid_sub_array_write(&rzgrid, &rzlocal, bphidg, fileNm);
  } while (0);


  //// write mapc2p and get metrics
  //printf("calculating mapc2p \n");
  //struct mapc2p_ctx *mctx = gkyl_malloc(sizeof(*mctx));
  //mctx->app = geo;
  //mctx->ginp = &ginp;
  //gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&cgrid, &cbasis, 3, mapc2p, mctx);
  //struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, clocal_ext.volume);
  //gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &clocal_ext, XYZ);

  //printf("writing rz file\n");
  //gkyl_grid_sub_array_write(&cgrid, &clocal, XYZ, "cerfon3d_xyzfile.gkyl");
  //printf("wrote rz file\n");
  
  printf("calculating metrics \n");
  struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 6*cbasis.num_basis, clocal_ext.volume);
  gkyl_calc_metric *mcalculator = gkyl_calc_metric_new(&cbasis, &cgrid, false);
  //gkyl_calc_metric_advance( mcalculator, &clocal, XYZ, gFld);
  gkyl_calc_metric_advance( mcalculator, &clocal, mapc2p_arr, gFld);
  do{
    printf("writing the gij file \n");
    const char *fmt = "%s_g_ij.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
    gkyl_grid_sub_array_write(&cgrid, &clocal, gFld, fileNm);
  } while (0);

  // done doing metrics
  printf("calculating jacobian, upper metrics, etc\n");
  //Now calculate the jacobian and upper metrics
  //struct gkyl_array *bmagFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *jFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *jinvFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *grFld = gkyl_array_new(GKYL_DOUBLE, 6*cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *biFld = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *cmagFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  printf("created geo fields\n");
  gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&cbasis, &cgrid, false);
  printf("made the j calculator \n");
  gkyl_calc_derived_geo_advance( jcalculator, &clocal, gFld, bmagFld, jFld, jinvFld, grFld, biFld, cmagFld);
  

  do{
    printf("writing the  second comp bmag file \n");
    const char *fmt = "%s_compbmag2.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmagFld, fileNm);
  } while (0);

    do{
      printf("writing the j file \n");
      const char *fmt = "%s_j.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
      gkyl_grid_sub_array_write(&cgrid, &clocal, jFld, fileNm);
    } while (0);

    do{
      printf("writing the cmag file \n");
      const char *fmt = "%s_cmag.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
      gkyl_grid_sub_array_write(&cgrid, &clocal, cmagFld, fileNm);
    } while (0);
 
  
    do{
      printf("writing the g^ij file \n");
      const char *fmt = "%s_gij.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
      gkyl_grid_sub_array_write(&cgrid, &clocal, grFld, fileNm);
    } while (0);


    do{
      printf("writing the b_i file \n");
      const char *fmt = "%s_bi.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "cerfon3d");
      gkyl_grid_sub_array_write(&cgrid, &clocal, biFld, fileNm);
    } while (0);








  //gkyl_array_release(RZ);
  gkyl_array_release(psidg);
  gkyl_array_release(psibyrdg);
  gkyl_array_release(psibyr2dg);
  gkyl_array_release(bmagFld);
  gkyl_array_release(bmagdg);
  gkyl_array_release(cmagFld);
  gkyl_array_release(jFld);
  gkyl_array_release(jinvFld);
  gkyl_array_release(biFld);
  gkyl_gkgeom_release(geo);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
