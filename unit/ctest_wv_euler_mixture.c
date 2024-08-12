#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_euler_mixture_priv.h>

void
test_euler_mixture_twocomponent_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(2, gas_gamma_s, false);

  TEST_CHECK( euler_mixture->num_equations == 8 );
  TEST_CHECK( euler_mixture->num_waves == 2 );

  double alpha1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (alpha1 * E1) + ((1.0 - alpha1) * E2);

  double q[8] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * alpha1,
    alpha1 * rho1, (1.0 - alpha1) * rho2 };

  double prims[8];
  gkyl_euler_mixture_prim_vars(2, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed mixture pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical mixture rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho2, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed mixture pressure.

  double fluxes[3][8] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * alpha1, alpha1 * (vx_total * rho1), (1.0 - alpha1) * (vx_total * rho2) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * alpha1, alpha1 * (vy_total * rho1), (1.0 - alpha1) * (vy_total * rho2) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * alpha1, alpha1 * (vz_total * rho1), (1.0 - alpha1) * (vz_total * rho2) },
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

  double q_local[8], flux_local[8], flux[8];
  for (int d = 0; d < 3; d++) {
    euler_mixture->rotate_to_local_func(euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_mixture_flux(2, gas_gamma_s, q_local, flux_local);
    euler_mixture->rotate_to_global_func(euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[8], q_g[8];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[8], q1[8];
    euler_mixture->cons_to_riem(euler_mixture, q_local, q_local, w1);
    euler_mixture->riem_to_cons(euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_threecomponent_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(3, gas_gamma_s, false);

  TEST_CHECK( euler_mixture->num_equations == 10 );
  TEST_CHECK( euler_mixture->num_waves == 2 );

  double alpha1 = 0.5, alpha2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (alpha1 * rho1) + (alpha2 * rho2) + ((1.0 - (alpha1 + alpha2)) * rho3);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E3 = (p_total / (gas_gamma3 - 1.0)) + (0.5 * rho3 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (alpha1 * E1) + (alpha2 * E2) + ((1.0 - (alpha1 + alpha2)) * rho3);

  double q[10] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * alpha1,
    rho_total * alpha2, alpha1 * rho1, alpha2 * rho2, (1.0 - (alpha1 + alpha2)) * rho3 };

  double prims[10];
  gkyl_euler_mixture_prim_vars(3, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed mixture pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical mixture rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], alpha2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[9], rho3, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed mixture pressure.

  double fluxes[3][10] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * alpha1, rho_total * vx_total * alpha2, alpha1 * (vx_total * rho1),
     alpha2 * (vx_total * rho2), (1.0 - (alpha1 + alpha2)) * (vx_total * rho3) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * alpha1, rho_total * vy_total * alpha2, alpha1 * (vy_total * rho1),
     alpha2 * (vy_total * rho2), (1.0 - (alpha1 + alpha2)) * (vy_total * rho3) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * alpha1, rho_total * vz_total * alpha2, alpha1 * (vz_total * rho1),
     alpha2 * (vz_total * rho2), (1.0 - (alpha1 + alpha2)) * (vz_total * rho3) },
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

  double q_local[10], flux_local[10], flux[10];
  for (int d = 0; d < 3; d++) {
    euler_mixture->rotate_to_local_func(euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_mixture_flux(3, gas_gamma_s, q_local, flux_local);
    euler_mixture->rotate_to_global_func(euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[10], q_g[10];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[10], q1[10];
    euler_mixture->cons_to_riem(euler_mixture, q_local, q_local, w1);
    euler_mixture->riem_to_cons(euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

TEST_LIST = {
  { "euler_mixture_twocomponent_basic", test_euler_mixture_twocomponent_basic },
  { "euler_mixture_threecomponent_basic", test_euler_mixture_threecomponent_basic },
  { NULL, NULL },
};