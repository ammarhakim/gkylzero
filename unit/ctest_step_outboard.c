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
#include <gkyl_calc_bmag.h>
#include <gkyl_tok_geo.h>


#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>






// gives Z(R)
double pfunc_upper(double R){
 double dzdr = -0.33333333333332943;
 double c = 9.990999999999982;

 //double dzdr = -0.06600660066006442;
 //double c = 8.634851485148507;
 double Z = dzdr*R + c;
 return Z;
}

double pfunc_lower(double R){
 double dzdr = 0.33333333333332943;
 double c = -9.990999999999982;
 //double dzdr = 0.06600660066006442;
 //double c = -8.634851485148507;
 double Z = dzdr*R + c;
 return Z;
}




void
test_1()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  char* filepath = "./efit_data/input.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath, rzpoly_order, fluxpoly_order, false);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, efit->psizr, "stepoutboard_psi.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, efit->fpolflux, "stepoutboard_fpol.gkyl");


  struct gkyl_tok_geo_inp inp = {
      // psiRZ and related inputs
      .rzgrid = efit->rzgrid,
      .rzbasis = efit->rzbasis,
      .psiRZ = efit->psizr,
      .psibyrRZ = efit->psibyrzr,
      .psibyr2RZ = efit->psibyr2zr,
      .rzlocal = efit->rzlocal,
      .rzlocal_ext = efit->rzlocal_ext,
      .fgrid = efit->fluxgrid,
      .frange = efit->fluxlocal,
      .fpoldg = efit->fpolflux,
      .qdg = efit->qflux,
      .fbasis = efit->fluxbasis,
      .psisep = efit->sibry,
      .plate_spec = false,
      // can set plate func if you dont want a fixed zmin and zmax
      //.plate_func_lower = pfunc_lower,
      //.plate_func_upper = pfunc_upper,
      //.plate_lower_Rl = 4.77,
      //.plate_lower_Rr = 5.073,
      //.plate_upper_Rl = 4.77,
      //.plate_upper_Rr = 5.073,
      .quad_param = {  .eps = 1e-10 }
    };

  //double clower[] = { 0.934, -0.01, -3.14 };
  //double cupper[] = {1.4688, 0.01, 3.14 };

  double clower[] = { 0.934, -0.01, -3.14 };
  double cupper[] = {1.0, 0.01, 3.14 };

  int ccells[] = { 4, 1, 64 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_tok_geo_geo_inp ginp = {
    .cgrid = &cgrid,
    .cbasis = &cbasis,
    .ftype = GKYL_SOL_DN_OUT,
    //.rclose = efit->rmax,
    //.zmin = efit->zmin,
    //.zmax = efit->zmax,

    .rclose = efit->rmax,
    .zmin = -8.3,
    .zmax = 8.3,
    .zmaxis = efit->zmaxis,
  
    .write_node_coord_array = true,
    .node_file_nm = "stepoutboard_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_new(&cgrid, &clocal, &clocal_ext, &cbasis, NULL, &inp, NULL, &ginp, true, false); 
  gkyl_gk_geometry_write(up);


  gkyl_efit_release(efit);
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
