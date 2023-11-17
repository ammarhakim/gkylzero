#include <acutest.h>

#include <gkyl_wv_cold_sr_fluid.h>
 
void 
test_durand_kerner()
{

  // Test case 
  double complex roots[4];

  // Reference solutions:
  double complex roots_ref[4];

  // input the reference roots
  roots_ref[0] = 0.0;
  roots_ref[1] = 0.0;
  roots_ref[2] = 0.0;
  roots_ref[3] = 0.0;

  //Left and right states
  double rhor = 9.9846560007158514e-01;
  double urx = 4.3872434904750417e+00;
  double ury = -1.5862931899624997e+00;
  double urz = 0.0000000000000000e+00;
  double rhol = 1.0000000000000000e+00;
  double ulx = 4.3521512798835111e+00;
  double uly = -1.5777826086635633e+00;
  double ulz = 0.0000000000000000e+00;

  // compute gammas
  double gammal = sqrt(1.0 + ulx*ulx + uly*uly + ulz*ulz);
  double gammar = sqrt(1.0 + urx*urx + ury*ury + urz*urz);

  // compute velocities
  double vLx = ulx/gammal;
  double vLy = uly/gammal;
  double vLz = ulz/gammal;
  double vRx = urx/gammar;
  double vRy = ury/gammar;
  double vRz = urz/gammar;

  // Hand off the structured data
  struct cold_sr csr = {
    .rhol = rhol,
    .rhor = rhor,
    .vl = { vLx, vLy, vLz },
    .vr = { vRx, vRy, vRz }
  };

  // Run the durand_kerner method to find the roots
  durand_kerner_method(&csr, roots);

  // Print the real and imaginary components with %1.16e format, for computed/reference
  printf("\n");
  for (int i = 0; i < 4; ++i) {
    printf("(Comp.) Root %d: Real = %1.16e, Imaginary = %1.16e\n", i+1, creal(roots[i]), cimag(roots[i]));
  }
  for (int i = 0; i < 4; ++i) {
    printf("(refr.) Root %d: Real = %1.16e, Imaginary = %1.16e\n", i+1, creal(roots_ref[i]), cimag(roots_ref[i]));
  }
}

void test_dk() { test_durand_kerner(); }

TEST_LIST = {
  {"test_dk", test_dk},
  {NULL, NULL},
};