#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_iso_euler_mixture.h>
#include <gkyl_wv_iso_euler_mixture_priv.h>

void
test_iso_euler_mixture_twocomponent_basic()
{
  double vt1 = 1.0;
  double vt2 = 10.0;

  double *vt_s = gkyl_malloc(sizeof(double[2]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(2, vt_s, false);

  TEST_CHECK( iso_euler_mixture->num_equations == 7 );
  TEST_CHECK( iso_euler_mixture->num_waves == 2 );

  double alpha1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3;
  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2);
  double vt_total = (alpha1 * vt1) + ((1.0 - alpha1) * vt2);

  double q[7] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, rho_total * alpha1, alpha1 * rho1,
    (1.0 - alpha1) * rho2 };

  double prims[7];
  gkyl_iso_euler_mixture_prim_vars(2, vt_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) ); 
  TEST_CHECK( gkyl_compare(prims[4], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[5], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho2, 1e-16) );

  double fluxes[3][7] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + (rho_total * (vt_total * vt_total)), rho_total * (vx_total * vy_total),
     rho_total * (vx_total * vz_total), rho_total * vx_total * alpha1, alpha1 * (vx_total * rho1), (1.0 - alpha1) * (vx_total * rho2) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + (rho_total * (vt_total * vt_total)),
     rho_total * (vy_total * vz_total), rho_total * vy_total * alpha1, alpha1 * (vy_total * rho1), (1.0 - alpha1) * (vy_total * rho2) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total),
     (rho_total * (vz_total * vz_total)) + (rho_total * (vt_total * vt_total)), rho_total * vz_total * alpha1, alpha1 * (vz_total * rho1),
     (1.0 - alpha1) * (vz_total * rho2) },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[7], flux_local[7], flux[7];
  for (int d = 0; d < 3; d++) {
    iso_euler_mixture->rotate_to_local_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_mixture_flux(2, vt_s, q_local, flux_local);
    iso_euler_mixture->rotate_to_global_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[7], q_g[7];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[7], q1[7];
    iso_euler_mixture->cons_to_riem(iso_euler_mixture, q_local, q_local, w1);
    iso_euler_mixture->riem_to_cons(iso_euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_threecomponent_basic()
{
  double vt1 = 1.0;
  double vt2 = 10.0;
  double vt3 = 100.0;

  double *vt_s = gkyl_malloc(sizeof(double[3]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;
  vt_s[2] = vt3;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(3, vt_s, false);

  TEST_CHECK( iso_euler_mixture->num_equations == 9 );
  TEST_CHECK( iso_euler_mixture->num_waves == 2 );

  double alpha1 = 0.5, alpha2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3;
  double rho_total = (alpha1 * rho1) + (alpha2 * rho2) + (1.0 - (alpha1 + alpha2)) * rho3;
  double vt_total = (alpha1 * vt1) + (alpha2 * vt2) + (1.0 - (alpha1 + alpha2)) * vt3;

  double q[9] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, rho_total * alpha1, rho_total * alpha2,
    alpha1 * rho1, alpha2 * rho2, (1.0 - (alpha1 + alpha2)) * rho3 };

  double prims[9];
  gkyl_iso_euler_mixture_prim_vars(3, vt_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) ); 
  TEST_CHECK( gkyl_compare(prims[4], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[5], alpha2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], rho3, 1e-16) );

  double fluxes[3][9] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + (rho_total * (vt_total * vt_total)), rho_total * (vx_total * vy_total),
     rho_total * (vx_total * vz_total), rho_total * vx_total * alpha1, rho_total * vx_total * alpha2, alpha1 * (vx_total * rho1),
     alpha2 * (vx_total * rho2), (1.0 - (alpha1 + alpha2)) * (vx_total * rho3) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + (rho_total * (vt_total * vt_total)),
     rho_total * (vy_total * vz_total), rho_total * vy_total * alpha1, rho_total * vy_total * alpha2, alpha1 * (vy_total * rho1),
     alpha2 * (vy_total * rho2), (1.0 - (alpha1 + alpha2)) * (vy_total * rho3) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total),
     (rho_total * (vz_total * vz_total)) + (rho_total * (vt_total * vt_total)), rho_total * vz_total * alpha1, rho_total * vz_total * alpha2,
     alpha1 * (vz_total * rho1), alpha2 * (vz_total * rho2), (1.0 - (alpha1 + alpha2)) * (vz_total * rho3) },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[9], flux_local[9], flux[9];
  for (int d = 0; d < 3; d++) {
    iso_euler_mixture->rotate_to_local_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_mixture_flux(3, vt_s, q_local, flux_local);
    iso_euler_mixture->rotate_to_global_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[9], q_g[9];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[9], q1[9];
    iso_euler_mixture->cons_to_riem(iso_euler_mixture, q_local, q_local, w1);
    iso_euler_mixture->riem_to_cons(iso_euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

TEST_LIST = {
  { "iso_euler_mixture_twocomponent_basic", test_iso_euler_mixture_twocomponent_basic },
  { "iso_euler_mixture_threecomponent_basic", test_iso_euler_mixture_threecomponent_basic },
  { NULL, NULL },
};