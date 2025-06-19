#include <acutest.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>

static void
test_g1()
{
  double w[gkyl_gauss_max], x[gkyl_gauss_max];
  
  for (unsigned n=1; n<=gkyl_gauss_max; ++n) {
    // use the generic routine to get the ordinates and weights
    gkyl_gauleg(-1, 1, &x[0], &w[0], n);

    // compare with pre-computed values
    const double *xp = gkyl_gauss_ordinates[n];
    const double *wp = gkyl_gauss_weights[n];
    
    for (unsigned i=0; i<n; ++i) {
      TEST_CHECK( gkyl_compare(x[i], xp[i], 1e-12) );
      TEST_CHECK( gkyl_compare(w[i], wp[i], 1e-12) );
    }
  }
}

static void
test_ndim()
{
  double nquad = 4;
  double *x = gkyl_malloc(sizeof(double)*nquad*nquad*2);
  double *w = gkyl_malloc(sizeof(double)*nquad*nquad);

  gkyl_ndim_ordinates_weights(2, x, w, nquad);

  double area = 0.0;
  for (int d=0; d<2; ++d) area += w[d];

  TEST_CHECK( gkyl_compare_double(2.0, area, 1e-14) );

  gkyl_free(x);
  gkyl_free(w);
}

TEST_LIST = {
  { "basic", test_g1 },
  { "ndim", test_ndim },
  { NULL, NULL },
};
