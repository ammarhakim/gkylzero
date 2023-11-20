
#include <acutest.h>
#include <gkyl_wv_cold_sr_fluid.h>
#include <math.h>

// TODO:
// 1. Add u-avg output in roe-averaged-velocity routine
// 2. check these statisfy momentum conservation 
// 3. If so, setup a static test for comparison
// 4. Satisfied entropy conditions
// 5. verify Vxl > vxr (conditions to compute U-hat)


bool 
local_entropy_test(const double ql[4], const double qr[4], double u[3], const double c)
{
  bool cond = 1;
  for(int i=0; i<3; ++i){
      cond = cond && fmin(qr[i+1]/qr[0],ql[i+1]/ql[0]) <= u[i] && u[i] <= fmax(qr[i+1]/qr[0],ql[i+1]/ql[0]);
      cond = cond && (NU_cons_func( ql, qr, i, u, c) <= 1.e-10);
  }
  return cond;
}

void 
test_1x3v()
{

  // CASE 1: 

  // Give constants for the test:
  const double ql[4] = {1.0,15.0,0.01,-5.0};
  const double qr[4] = {0.1,0.5,0.01,25.0};
  const double c = 1.0;
  double u[3];

  // Call the roe-averaged velocity routine
  double v_avg = compute_sr_roe_averaged_velocity_via_ridders(ql, qr, c, u);

  // Compute comparison with conservation laws
  double dens_cons = dens_cons_func( ql, qr, u, c );
  double NU_cons[3];
  for(int i=0; i<3; ++i)
    NU_cons[i] = NU_cons_func( ql, qr, i, u, c);

  // Test the data 
  if (ql[0] > 0.0 && qr[0] > 0.0){
    printf("\nv_avg = %1.16e\n",v_avg);
    printf("Ux = %1.16e\n",u[0]);
    printf("Uy = %1.16e\n",u[1]);
    printf("Uz = %1.16e\n",u[2]);
    printf("Satisfy Vx - Ux/gamma: %1.16e\n",v_avg-u[0]/sqrt(1.0+(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(c*c)));
    printf("Satisfy Density-conservation eqn.: %1.16e\n",dens_cons);
    for(int i=0; i<3; ++i)
      printf("Satisfy NU%d-conservation eqn.: %1.16e\n",i,NU_cons[i]);
    for(int i=0; i<3; ++i)
      printf("Satisfied entropy (U%d): %1.16e < %1.16e < %1.16e\n",i,
        fmin(qr[i+1]/qr[0],ql[i+1]/ql[0]),u[i],fmax(qr[i+1]/qr[0],ql[i+1]/ql[0]));

    // Test the entropy condition is satisfied
    bool pass_entropy_test = local_entropy_test(ql, qr, u, c);
    if (pass_entropy_test)
      printf("Entropy test (PASSES)\n");
    else
      printf("Entropy test (FAILS)\n");

    // Run a test on the result:
    double v_avg_ref = 5.2671432172836086e-01;
    TEST_CHECK(gkyl_compare_double(v_avg_ref, v_avg, 1e-12));
  } else {
    printf("q <= 0 on one side!\n");
  }

}

void test_1d3v_v_hat_comp() { test_1x3v(); }

TEST_LIST = {
  {"test_1d3v_v_hat_comp", test_1d3v_v_hat_comp},
  {NULL, NULL},
};