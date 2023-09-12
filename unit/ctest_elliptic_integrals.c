#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


#include <acutest.h>
#include <elliptic_integral.h>
//#include <gkyl_array.h>
//#include <gkyl_array_rio.h>
//#include <gkyl_eval_on_nodes.h>
//#include <gkyl_range.h>
//#include <gkyl_rect_grid.h>
//#include <gkyl_rect_decomp.h>
//#include <gkyl_util.h>
//#include <gkyl_basis.h>
//
//#include <gkyl_calc_bmag.h>
//#include <gkyl_calc_bmag_kernels.h>
//#include <gkyl_geo_gyrokinetic.h>
//
//
//#include <gkyl_calc_metric.h>
//#include <gkyl_calc_metric_kernels.h>
//#include <gkyl_calc_derived_geo.h>
//#include <gkyl_calc_derived_geo_kernels.h>
//#include <gkyl_geo_gyrokinetic.h>



static inline double sq(double x) { return x*x; }



void
test_1()
{

  double m = 0.1;
  double K = Complete_Elliptic_Integral_First_Kind( 'm', m );
  printf("\n K(m) = %12.6f where m = %1.10f\n",K, m);
  double E = Complete_Elliptic_Integral_Second_Kind( 'm', m );
  printf("\n E(m) = %12.6f where m = %1.10f\n",E, m);


  m = 0.3;
  K = Complete_Elliptic_Integral_First_Kind( 'm', m );
  printf("\n K(m) = %12.6f where m = %1.10f\n",K, m);
  E = Complete_Elliptic_Integral_Second_Kind( 'm', m );
  printf("\n E(m) = %12.6f where m = %1.10f\n",E, m);

  m = 0.5;
  K = Complete_Elliptic_Integral_First_Kind( 'm', m );
  printf("\n K(m) = %12.6f where m = %1.10f\n",K, m);
  E = Complete_Elliptic_Integral_Second_Kind( 'm', m );
  printf("\n E(m) = %12.6f where m = %1.10f\n",E, m);

  m = 0.9;
  K = Complete_Elliptic_Integral_First_Kind( 'm', m );
  printf("\n K(m) = %12.6f where m = %1.10f\n",K, m);
  E = Complete_Elliptic_Integral_Second_Kind( 'm', m );
  printf("\n E(m) = %12.6f where m = %1.10f\n",E, m);


}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
