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

#include <gkyl_calc_bmag.h>
#include <gkyl_calc_bmag_kernels.h>
#include <gkyl_geo_gyrokinetic.h>


#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_kernels.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_derived_geo_kernels.h>
#include <gkyl_geo_gyrokinetic.h>
#include <elliptic_integral.h>

#include <gkyl_efit.h>
#include <gkyl_efit_priv.h>


struct mapc2p_ctx{
   struct gkyl_geo_gyrokinetic* app;
   struct gkyl_geo_gyrokinetic_geo_inp* ginp;
};


// step equilibrium
struct step_ctx {
  double R0, B0;
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




double Alist[24] = {1.098170000000000091e+00,  
5.749690000000000190e+00,  
5.254970000000000141e+00,  
7.063550000000000217e+00,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.199999999999999956e+00,  
1.098170000000000091e+00,  
5.749690000000000190e+00,  
5.254970000000000141e+00,  
7.063550000000000217e+00,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.000000000000000056e-01,  
1.199999999999999956e+00,  
2.339999999999999858e+00,  
2.860000000000000320e+00,  
2.600000000000000089e+00,  
2.600000000000000089e+00,  
2.600000000000000089e+00,  
2.600000000000000089e+00};


 double Zlist[24] = {6.067210000000000214e+00,
 8.385300000000000864e+00,
 5.380029999999999646e+00,
 1.933729999999999949e+00,
 2.000000000000000000e+00,
 0.000000000000000000e+00,
 7.500000000000000000e-01,
 9.000000000000000222e-01,
 8.000000000000000000e+00,
 -6.067210000000000214e+00,
 -8.385300000000000864e+00,
 -5.380029999999999646e+00,
 -1.933729999999999949e+00,
 -2.000000000000000000e+00,
 -0.000000000000000000e+00,
 -7.500000000000000000e-01,
 -9.000000000000000222e-01,
 -8.000000000000000000e+00,
 0.000000000000000000e+00,
 0.000000000000000000e+00,
 2.600000000000000089e-01,
 -2.600000000000000089e-01,
 0.000000000000000000e+00,
 0.000000000000000000e+00};

 double current_list[24] = {2.733587001572373323e+06,
 2.576555680242971517e+06,
 -4.769840006204952952e+05,
 -5.848048038329290226e+06,
 -6.272967875940668583e+07,
 -9.374610126189786196e+07,
 -5.444083807201912999e+07,
 -4.781479151113666594e+07,
 2.885451375306561589e+06,
 2.733587001572373323e+06,
 2.576555680242971517e+06,
 -4.769840006204952952e+05,
 -5.848048038329290226e+06,
 -6.272967875940668583e+07,
 -9.374610126189786196e+07,
 -5.444083807201912999e+07,
 -4.781479151113666594e+07,
 2.885451375306561589e+06,
 -2.886111956023474224e+06,
 2.886111956023474224e+06,
 6.830522162672655284e+07,
 6.830522162672655284e+07,
 -1.366104432534531057e+08,
 1.650000000000000000e+07};





double ellipk(double m){
  return Complete_Elliptic_Integral_First_Kind( 'm', m );
}
double ellipe(double m){
  return Complete_Elliptic_Integral_Second_Kind( 'm', m );
}

// input is spherical coords r,theta,a
double Aphi_f(double r, double theta, double a){
    double mu0 = 1.256e-6;
    double ksq = 4*a*r*sin(theta)/(sq(a)+sq(r) + 2*a*r*sin(theta));
    double fac1 = mu0/(4*M_PI);

    double num2 = 4*a;
    double denom2 = sqrt(sq(a) + sq(r) + 2*a*r*sin(theta));
    double fac2 = num2/denom2;

    double num3 = (2-ksq)*ellipk(ksq) - 2*ellipe(ksq);
    double denom3 = ksq;
    double fac3 = num3/denom3;

    return fac1*fac2*fac3;
}

void
psi(double t, const double *xn, double *fout, void *ctx)
{
  struct step_ctx *s = ctx;
  double R0 = s->R0;
  double R = xn[0], Z = xn[1];

  double psiret = 0.0;

  for(int i = 0; i<24; i++){
    double zi = Z - Zlist[i];
    double rs = sqrt(sq(R) + sq(zi));
    double thetas = acos(zi/rs);
    double Aphii = Aphi_f(rs,thetas,Alist[i]);
    double psii = current_list[i]*Aphii*R;
    psiret += psii;
  }
  fout[0] = psiret;
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
  struct step_ctx *s = ctx;
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



  
  // create RZ grid
  //double lower[] = { 0.01, -7.5 }, upper[] = { 7.0, 7.4 };
  //int cells[] = { 64, 128 };

  //struct gkyl_rect_grid rzgrid;
  //gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  //struct gkyl_range rzlocal, rzlocal_ext;
  //int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  //gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  //// RZ basis function
  //int rzpoly_order = 2;
  //struct gkyl_basis rzbasis;
  //gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  // allocate psiRZ array, initialize and write it to file
  //struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  //
  //gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, &psi, &sctx);
  //gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  //gkyl_eval_on_nodes_release(eon);

  char* filepath = "./efit_data/input.geqdsk";
  // RZ basis function
  int rzpoly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);
  struct gkyl_rect_grid rzgrid;
  struct gkyl_range rzlocal, rzlocal_ext;
  // flux basis function
  int fluxpoly_order = 1;
  struct gkyl_basis fluxbasis;
  gkyl_cart_modal_serendip(&fluxbasis, 1, fluxpoly_order);
  struct gkyl_rect_grid fluxgrid;
  struct gkyl_range fluxlocal, fluxlocal_ext;
  struct gkyl_efit* efit = gkyl_efit_new(filepath, &rzbasis, &fluxbasis, false);

  gkyl_rect_grid_init(&rzgrid, 2, efit->rzlower, efit->rzupper, efit->rzcells);
  gkyl_create_grid_ranges(&rzgrid, efit->rzghost, &rzlocal_ext, &rzlocal);

  gkyl_rect_grid_init(&fluxgrid, 1, efit->fluxlower, efit->fluxupper, efit->fluxcells);
  gkyl_create_grid_ranges(&fluxgrid, efit->fluxghost, &fluxlocal_ext, &fluxlocal);



  struct gkyl_array* psizr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array* psibyrzr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array* psibyr2zr = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array* fpolflux = gkyl_array_new(GKYL_DOUBLE, fluxbasis.num_basis, fluxlocal_ext.volume);
  printf("mads psi arrays\n");
  gkyl_efit_advance(efit, &rzgrid, &fluxgrid, &rzlocal, &rzlocal_ext, psizr, psibyrzr, psibyr2zr, &fluxlocal, &fluxlocal_ext, fpolflux);


  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psizr, "efit_psi.gkyl");
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psibyrzr, "efit_psibyr.gkyl");
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psibyr2zr, "efit_psibyr2.gkyl");
  gkyl_grid_sub_array_write(&fluxgrid, &fluxlocal, fpolflux, "stepoutboard_fpol.gkyl");


  //struct step_ctx sctx = {  .R0 = 2.6,  .B0 = 2.1 };
  struct step_ctx sctx = {  .R0 = efit->rcentr,  .B0 = efit->bcentr };
  printf("rcentr, bcentr = %g, %g",efit->rcentr,efit->bcentr);
  



  gkyl_geo_gyrokinetic *geo = gkyl_geo_gyrokinetic_new(&(struct gkyl_geo_gyrokinetic_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psizr,
      .rzlocal = &rzlocal,
      .B0 = sctx.B0,
      .R0 = sctx.R0,

      .quad_param = {  .eps = 1e-10 }
    }
  );

  // compute outboard SOL geometry
  

  double psiSep = 1.50982;

  //double clower[] = { 0.934, -0.01, -3.14 };
  //double cupper[] = {1.4688, 0.01, 3.14 };

  double clower[] = { psiSep, -0.01, -3.14 };
  double cupper[] = {1.9, 0.01, 3.14 };

  int ccells[] = { 1, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext, conversion_range;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  // BCs, 0 is periodic, 1 is nonperiodic
  int bcs[3] = {1,0,1};

  // calcgeom will go into the ghost y cells based on bc. If bc[1]=1 we use ghosts.
  // Need to pass appropriate conversion to modal range depending on the bcs
  if(bcs[1]==1){
    int sublower[3] = {clocal.lower[0], clocal_ext.lower[1], clocal.lower[2]};
    int subupper[3] = {clocal.upper[0], clocal_ext.upper[1], clocal.upper[2]};
    gkyl_sub_range_init(&conversion_range, &clocal_ext, sublower, subupper);
  }
  else{
    conversion_range = clocal;
  }

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_geo_gyrokinetic_geo_inp ginp = {
    .cgrid = &cgrid,
    .cbasis = &cbasis,
    .ftype = GKYL_CORE,
    .bcs = bcs,
    //.zmin = efit->zmin,
    //.zmax = efit->zmax,

    .rclose = efit->rmax,
    .rleft= efit->rmin,
    .rright= efit->rmax,
    .zmin = -6.2,
    .zmax = 6.2,
    .zmaxis = -4.84614e-07,
  
    .write_node_coord_array = true,
    .node_file_nm = "stepclosed_nodes.gkyl"
  }; 




  //Do Ammar's calcgeom

  struct gkyl_array *mapc2p_arr = gkyl_array_new(GKYL_DOUBLE, 3*cbasis.num_basis, clocal_ext.volume);
  gkyl_geo_gyrokinetic_calcgeom(geo, &ginp, mapc2p_arr, &conversion_range);

  printf("writing mapc2p file from calcgeom\n");
  gkyl_grid_sub_array_write(&cgrid, &clocal, mapc2p_arr, "stepclosed_mapc2pfile.gkyl");
  printf("wrote mapc2p file\n");

  //make psi
  //gkyl_eval_on_nodes *eval_psi = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psi, &sctx);
  //struct gkyl_array* psidg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  //gkyl_eval_on_nodes_advance(eval_psi, 0.0, &rzlocal_ext, psidg); //on ghosts with ext_range
  //gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psidg, "stepclosed_psi.gkyl");

  //gkyl_eval_on_nodes *eval_psibyr = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr, &sctx);
  //struct gkyl_array* psibyrdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  //gkyl_eval_on_nodes_advance(eval_psibyr, 0.0, &rzlocal_ext, psibyrdg); //on ghosts with ext_range
  //                                                                     //
  //gkyl_eval_on_nodes *eval_psibyr2 = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, psibyr2, &sctx);
  //struct gkyl_array* psibyr2dg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  //gkyl_eval_on_nodes_advance(eval_psibyr2, 0.0, &rzlocal_ext, psibyr2dg); //on ghosts with ext_range
  //printf("made the psi arays\n");

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

  gkyl_calc_bmag_advance(calculator, &clocal, &clocal_ext, &rzlocal, &rzlocal_ext, psizr, psibyrzr, psibyr2zr, bphidg, bmagFld, mapc2p_arr);
  printf("advanced bmag calculator\n");



  char fileNm[1024];
  do{
    printf("writing the comp bmag file \n");
    const char *fmt = "%s_compbmag.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmagFld, fileNm);
  } while (0);

  do{
    printf("writing the bphi file \n");
    const char *fmt = "%s_bphi.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
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
  //gkyl_grid_sub_array_write(&cgrid, &clocal, XYZ, "stepclosed_xyzfile.gkyl");
  //printf("wrote rz file\n");
  
  printf("calculating metrics \n");
  struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 6*cbasis.num_basis, clocal_ext.volume);
  printf("creating\n");
  gkyl_calc_metric *mcalculator = gkyl_calc_metric_new(&cbasis, &cgrid, bcs, false);
  printf("advancing\n");
  //gkyl_calc_metric_advance( mcalculator, &clocal, mapc2p_arr, gFld);
  gkyl_calc_metric_advance(mcalculator, geo->nrange, geo->mc2p_nodal_fd, geo->dzc, gFld, &conversion_range);

  do{
    printf("writing the gij file \n");
    const char *fmt = "%s_g_ij.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
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
  struct gkyl_array *jtotFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *jtotinvFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *bmaginvFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *bmaginvsqFld = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *gxxJ= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *gxyJ= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  struct gkyl_array *gyyJ= gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  printf("created geo fields\n");
  gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&cbasis, &cgrid, false);
  printf("made the j calculator \n");
  gkyl_calc_derived_geo_advance( jcalculator, &clocal, gFld, bmagFld, jFld, jinvFld, grFld, biFld, cmagFld, jtotFld, jtotinvFld, bmaginvFld, bmaginvsqFld, gxxJ, gxyJ, gyyJ);
  

  do{
    printf("writing the  second comp bmag file \n");
    const char *fmt = "%s_compbmag2.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
    gkyl_grid_sub_array_write(&cgrid, &clocal, bmagFld, fileNm);
  } while (0);

    do{
      printf("writing the j file \n");
      const char *fmt = "%s_j.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
      gkyl_grid_sub_array_write(&cgrid, &clocal, jFld, fileNm);
    } while (0);

    do{
      printf("writing the jtot file \n");
      const char *fmt = "%s_jtot.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
      gkyl_grid_sub_array_write(&cgrid, &clocal, jtotFld, fileNm);
    } while (0);

    do{
      printf("writing the cmag file \n");
      const char *fmt = "%s_cmag.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
      gkyl_grid_sub_array_write(&cgrid, &clocal, cmagFld, fileNm);
    } while (0);
 
  
    do{
      printf("writing the g^ij file \n");
      const char *fmt = "%s_gij.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
      gkyl_grid_sub_array_write(&cgrid, &clocal, grFld, fileNm);
    } while (0);


    do{
      printf("writing the b_i file \n");
      const char *fmt = "%s_bi.gkyl";
      snprintf(fileNm, sizeof fileNm, fmt, "stepclosed");
      gkyl_grid_sub_array_write(&cgrid, &clocal, biFld, fileNm);
    } while (0);








  //gkyl_array_release(RZ);
  gkyl_efit_release(efit);
  gkyl_array_release(psizr);
  gkyl_array_release(psibyrzr);
  gkyl_array_release(psibyr2zr);
  gkyl_array_release(bmagFld);
  gkyl_array_release(bmagdg);
  gkyl_array_release(cmagFld);
  gkyl_array_release(jFld);
  gkyl_array_release(jinvFld);
  gkyl_array_release(biFld);

  gkyl_array_release(jtotFld);
  gkyl_array_release(jtotinvFld);
  gkyl_array_release(bmaginvFld);
  gkyl_array_release(bmaginvsqFld);
  gkyl_array_release(gxxJ);
  gkyl_array_release(gxyJ);
  gkyl_array_release(gyyJ);
  gkyl_geo_gyrokinetic_release(geo);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
