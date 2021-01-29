#pragma once

#include <math.h>
#include <stdlib.h>

#include <gkyl_real_type.h>

// Maximum order ordinate/weight data
static const int gkyl_gauss_max = 8;

// Ordinates
static const gkyl_real gkyl_gauss_ordinates_1[] =
{ 0.0 };
static const gkyl_real gkyl_gauss_ordinates_2[] =
{ -0.5773502691896257645091, 0.5773502691896257645091 };
static const gkyl_real gkyl_gauss_ordinates_3[] =
{ -0.7745966692414833770359, 0, 0.7745966692414833770359 };
static const gkyl_real gkyl_gauss_ordinates_4[] =
{ -0.8611363115940525752239, -0.3399810435848562648027, 0.3399810435848562648027, 0.8611363115940525752239 };
static const gkyl_real gkyl_gauss_ordinates_5[] =
{ -0.9061798459386639927976, -0.5384693101056830910363, 0, 0.5384693101056830910363, 0.9061798459386639927976 };
static const gkyl_real gkyl_gauss_ordinates_6[] =
{ -0.9324695142031520278123, -0.6612093864662645136614, -0.2386191860831969086305, 0.2386191860831969086305, 0.6612093864662645136614, 0.9324695142031520278123 };
static const gkyl_real gkyl_gauss_ordinates_7[] =
{ -0.9491079123427585245262, -0.7415311855993944398639, -0.4058451513773971669066, 0, 0.4058451513773971669066, 0.7415311855993944398639, 0.9491079123427585245262 };
static const gkyl_real gkyl_gauss_ordinates_8[] =
{ -0.960289856497536231684, -0.7966664774136267395916, -0.5255324099163289858177, -0.1834346424956498049395, 0.1834346424956498049395, 0.525532409916328985818, 0.796666477413626739592, 0.9602898564975362316836 };

// Weights
static const gkyl_real gkyl_gauss_weights_1[] =
{ 2.0 };
static const gkyl_real gkyl_gauss_weights_2[] =
{ 1.0, 1.0 };
static const gkyl_real gkyl_gauss_weights_3[] =
{ 0.5555555555555555555556, 0.888888888888888888889, 0.555555555555555555556 };
static const gkyl_real gkyl_gauss_weights_4[] =
{ 0.3478548451374538573731, 0.6521451548625461426269, 0.652145154862546142627, 0.3478548451374538573731 };
static const gkyl_real gkyl_gauss_weights_5[] =
{ 0.2369268850561890875143, 0.4786286704993664680413, 0.568888888888888888889, 0.478628670499366468041, 0.236926885056189087514 };
static const gkyl_real gkyl_gauss_weights_6[] =
{ 0.1713244923791703450403, 0.36076157304813860757, 0.46791393457269104739, 0.46791393457269104739, 0.36076157304813860757, 0.1713244923791703450403 };
static const gkyl_real gkyl_gauss_weights_7[] =
{ 0.1294849661688696932706, 0.279705391489276667901, 0.38183005050511894495, 0.4179591836734693877552, 0.38183005050511894495, 0.279705391489276667901, 0.1294849661688696932706 };
static const gkyl_real gkyl_gauss_weights_8[] =
{ 0.1012285362903762591525, 0.222381034453374470544, 0.313706645877887287338, 0.3626837833783619829652, 0.3626837833783619829652, 0.31370664587788728734, 0.222381034453374470544, 0.1012285362903762591525 };

// gkyl_gauss_ordinates[N] are ordinates for N-point Guassian
// integration
static const gkyl_real* gkyl_gauss_ordinates[] = {
  NULL, // N=0 makes no sense,
  gkyl_gauss_ordinates_1,
  gkyl_gauss_ordinates_2,
  gkyl_gauss_ordinates_3,
  gkyl_gauss_ordinates_4,
  gkyl_gauss_ordinates_5,
  gkyl_gauss_ordinates_6,
  gkyl_gauss_ordinates_7,
  gkyl_gauss_ordinates_8
};

// gkyl_gauss_weights[N] are weights for N-point Guassian integration
static const gkyl_real* gkyl_gauss_weights[] = {
  NULL, // N=0 makes no sense,
  gkyl_gauss_weights_1,
  gkyl_gauss_weights_2,
  gkyl_gauss_weights_3,
  gkyl_gauss_weights_4,
  gkyl_gauss_weights_5,
  gkyl_gauss_weights_6,
  gkyl_gauss_weights_7,
  gkyl_gauss_weights_8
};

/**
 * Compute ordinates and weights for use in Gaussian quadrature.
 *
 * @param x1 Left coordinate of domain
 * @param x2 Right coordinate of domain
 * @param x On output, ordinates
 * @param w On output, weights
 */
void gkyl_gauleg(gkyl_real x1, gkyl_real x2,  gkyl_real x[], gkyl_real w[], int n);
