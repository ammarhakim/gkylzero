extern "C" {
#include <float.h>
#include <gkylzero.h>
#include <gkyl_array.h>
#include <gkyl_wv_maxwell_priv.h>
#include <gkyl_wave_prop.h>
int cu_wave_prop_test(gkyl_wave_prop **up, const int ndim);
}

// FIXME: duplicate of gkyl_compare_double in util.c
GKYL_CU_DH static int compare(double a, double b, double eps)
{
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}

/* check primitive values and basic validity of structs */
__global__
void ker_cu_wave_prop_test(const gkyl_wave_prop *up, const int d, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( up->limiter == GKYL_MONOTONIZED_CENTERED, nfail );
  GKYL_CU_CHECK( up->num_up_dirs == 1, nfail );
  GKYL_CU_CHECK( up->update_dirs[0] == d, nfail );
  GKYL_CU_CHECK( up->cfl == 1.0, nfail);

  const struct gkyl_wv_eqn *eqn = up->equation;
  GKYL_CU_CHECK( eqn->num_equations == 8, nfail );
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  GKYL_CU_CHECK( maxwell->c == 299792458.0, nfail );

  const struct gkyl_wave_geom *wg = up->geom;
  double r_inn = 0.25, r_out = 1.25, phi_max = M_PI / 2.;
  double area = 0.5 * (r_out*r_out - r_inn*r_inn);
  double area_c = (r_out - r_inn) * phi_max;
  int idx[] = {1, 1, 1};
  {
    // const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idx);
    const struct gkyl_wave_cell_geom *cg = (const struct gkyl_wave_cell_geom *)
      gkyl_array_cfetch(wg->geom, gkyl_range_idx(&wg->range, idx));

    GKYL_CU_CHECK(compare(cg->kappa, area / area_c, 1e-8), nfail);
  }
}

int cu_wave_prop_test(gkyl_wave_prop **slvr, const int ndim)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));

  for (int d=0; d<ndim; ++d) {
    ker_cu_wave_prop_test<<<1, 1>>>(slvr[d]->on_dev, d, nfail_dev);
    checkCuda(cudaGetLastError());
  }

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

