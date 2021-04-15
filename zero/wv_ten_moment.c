#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_prim_ten_moment.h>
#include <gkyl_wv_ten_moment.h>

static double sq(double x) { return x*x; }

static const int dir_u_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

static const int dir_p_shuffle[][6] = {
  {4, 5, 6, 7, 8, 9},
  {7, 8, 5, 9, 6, 4},
  {9, 6, 8, 4, 5, 7}
};

// Make indexing cleaner with the dir_shuffle
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

#define PXX dp[0]
#define PXY dp[1]
#define PXZ dp[2]
#define PYY dp[3]
#define PYZ dp[4]
#define PZZ dp[5]

/* Multiply by phi prime */
static void mulByPhiPrime(double p0, double u1, double u2, double u3, const double w[10], double out[10]) 
{ 
   out[0] = w[0]; 
   out[1] = w[0]*u1+w[1]*p0; 
   out[2] = w[0]*u2+w[2]*p0; 
   out[3] = w[0]*u3+w[3]*p0; 
   out[4] = w[0]*sq(u1)+2*w[1]*p0*u1+w[4]; 
   out[5] = w[0]*u1*u2+w[1]*p0*u2+w[2]*p0*u1+w[5]; 
   out[6] = w[0]*u1*u3+w[1]*p0*u3+w[3]*p0*u1+w[6]; 
   out[7] = w[0]*sq(u2)+2*w[2]*p0*u2+w[7]; 
   out[8] = w[0]*u2*u3+w[2]*p0*u3+w[3]*p0*u2+w[8]; 
   out[9] = w[0]*sq(u3)+2*w[3]*p0*u3+w[9]; 
} 

struct wv_ten_moment {
    struct gkyl_wv_eqn eqn; // base object
};

static void
ten_moment_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_ten_moment *ten_moment = container_of(base, struct wv_ten_moment, eqn);
  gkyl_free(ten_moment);
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_ten_moment *ten_moment = container_of(eqn, struct wv_ten_moment, eqn);
  const int *d = dir_u_shuffle[dir];
  const int *dp = dir_p_shuffle[dir];

  double vl[10], vr[10];
  gkyl_ten_moment_primitive(ql, vl);
  gkyl_ten_moment_primitive(qr, vr);

  // compute Roe averages
  double sqrl = sqrt(vl[0]), sqrr = sqrt(vr[0]);
  double sqr1 = 1/(sqrl+sqrr);
  
  double p0 = sqrl*sqrr;
  double p2s1 = sq(p0*sqr1);
  
  double u1 = (sqrl*vl[RHOU] + sqrr*vr[RHOU])*sqr1;
  double u2 = (sqrl*vl[RHOV] + sqrr*vr[RHOV])*sqr1;
  double u3 = (sqrl*vl[RHOW] + sqrr*vr[RHOW])*sqr1;
  double p11 = (sqrr*vl[PXX]+sqrl*vr[PXX])*sqr1 + 1.0/3.0*p2s1*(vr[RHOU]-vl[RHOU])*(vr[RHOU]-vl[RHOU]);
  double p12 = (sqrr*vl[PXY]+sqrl*vr[PXY])*sqr1 + 1.0/3.0*p2s1*(vr[RHOU]-vl[RHOU])*(vr[RHOV]-vl[RHOV]);
  double p13 = (sqrr*vl[PXZ]+sqrl*vr[PXZ])*sqr1 + 1.0/3.0*p2s1*(vr[RHOU]-vl[RHOU])*(vr[RHOW]-vl[RHOW]);
  double p22 = (sqrr*vl[PYY]+sqrl*vr[PYY])*sqr1 + 1.0/3.0*p2s1*(vr[RHOV]-vl[RHOV])*(vr[RHOV]-vl[RHOV]);
  double p23 = (sqrr*vl[PYZ]+sqrl*vr[PYZ])*sqr1 + 1.0/3.0*p2s1*(vr[RHOV]-vl[RHOV])*(vr[RHOW]-vl[RHOW]);
  double p33 = (sqrr*vl[PZZ]+sqrl*vr[PZZ])*sqr1 + 1.0/3.0*p2s1*(vr[RHOW]-vl[RHOW])*(vr[RHOW]-vl[RHOW]);

  // for multiplication by phi' we need to use unrotated values
  double v[4];
  v[RHOU] = u1; v[RHOV] = u2; v[RHOW] = u3;

  double phiDelta[10];

  // pre-multiply jump (delta) by phiPrime inverse: we do this as
  // jumps are in conserved variables, while left eigenvectors used
  // below are computed from primitive variables
  phiDelta[0] = delta[0];
  phiDelta[1] = delta[RHOU]/p0-(1.0*delta[0]*u1)/p0; 
  phiDelta[2] = delta[RHOV]/p0-(1.0*delta[0]*u2)/p0; 
  phiDelta[3] = delta[RHOW]/p0-(1.0*delta[0]*u3)/p0; 
  phiDelta[4] = delta[0]*sq(u1)-2.0*delta[RHOU]*u1+delta[PXX]; 
  phiDelta[5] = delta[0]*u1*u2-1.0*delta[RHOU]*u2-1.0*delta[RHOV]*u1+delta[PXY]; 
  phiDelta[6] = delta[0]*u1*u3-1.0*delta[RHOU]*u3-1.0*delta[RHOW]*u1+delta[PXZ]; 
  phiDelta[7] = delta[0]*sq(u2)-2.0*delta[RHOV]*u2+delta[PYY]; 
  phiDelta[8] = delta[0]*u2*u3-1.0*delta[RHOV]*u3-1.0*delta[RHOW]*u2+delta[PYZ]; 
  phiDelta[9] = delta[0]*sq(u3)-2.0*delta[RHOW]*u3+delta[PZZ];

  // predefine some constants
  double p11sq = sq(p11), p12sq = sq(p12), p13sq = sq(p13), p11th = sqrt(p11*p11*p11);
  double sqp0 = sqrt(p0), sqp11 = sqrt(p11);
  
  double leftProj[10];
  // project jumps on left eigenvectors [Gen from Maxima]
  leftProj[0] = (0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq-(0.5*phiDelta[4]*p12)/p11sq-(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[1] = (0.5*phiDelta[1]*sqp0*p13)/p11th-(0.5*phiDelta[4]*p13)/p11sq-(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[2] = (-(0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq)-(0.5*phiDelta[4]*p12)/p11sq+(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[3] = (-(0.5*phiDelta[1]*sqp0*p13)/p11th)-(0.5*phiDelta[4]*p13)/p11sq+(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[4] = (0.16666666666666666667*phiDelta[4])/p11sq-(0.2886751345948129*phiDelta[1]*sqp0)/p11th; 
  leftProj[5] = (0.2886751345948129*phiDelta[1]*sqp0)/p11th+(0.16666666666666666667*phiDelta[4])/p11sq; 
  leftProj[6] = phiDelta[0]-(0.3333333333333333333*phiDelta[4]*p0)/p11; 
  leftProj[7] = (-(0.3333333333333333333*phiDelta[4]*p11*p22)/p11sq)+(1.333333333333333333*phiDelta[4]*p12sq)/p11sq-(2.0*phiDelta[5]*p12)/p11+phiDelta[7]; 
  leftProj[8] = (-(0.3333333333333333333*phiDelta[4]*p11*p23)/p11sq)+(1.333333333333333333*phiDelta[4]*p12*p13)/p11sq-(1.0*phiDelta[5]*p13)/p11-(1.0*phiDelta[6]*p12)/p11+phiDelta[8]; 
  leftProj[9] = (-(0.3333333333333333333*phiDelta[4]*p11*p33)/p11sq)+(1.333333333333333333*phiDelta[4]*p13sq)/p11sq-(2.0*phiDelta[6]*p13)/p11+phiDelta[9]; 

  // compute waves and speeds
  double wv[10];

  // Wave 1: (ev 1 and 2 are repeated)
  s[0] = u1-sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[RHOU] = 0.0; 
  wv[RHOV] = -(1.0*leftProj[0]*sqp11)/sqp0; 
  wv[RHOW] = -(1.0*leftProj[1]*sqp11)/sqp0; 
  wv[PXX] = 0.0; 
  wv[PXY] = leftProj[0]*p11; 
  wv[PXZ] = leftProj[1]*p11; 
  wv[PYY] = 2.0*leftProj[0]*p12; 
  wv[PYZ] = leftProj[0]*p13+leftProj[1]*p12; 
  wv[PZZ] = 2.0*leftProj[1]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[0]);

  // Wave 2: (ev 3 and 4 are repeated)
  s[1] = u1+sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[RHOU] = 0.0; 
  wv[RHOV] = (leftProj[2]*sqp11)/sqp0; 
  wv[RHOW] = (leftProj[3]*sqp11)/sqp0; 
  wv[PXX] = 0.0; 
  wv[PXY] = leftProj[2]*p11; 
  wv[PXZ] = leftProj[3]*p11; 
  wv[PYY] = 2.0*leftProj[2]*p12; 
  wv[PYZ] = leftProj[2]*p13+leftProj[3]*p12; 
  wv[PZZ] = 2.0*leftProj[3]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[10]);
  
  // Wave 3 (ev 5)
  s[2] = u1-sqrt(3*p11/p0);
  wv[0] = leftProj[4]*p0*p11; 
  wv[RHOU] = -(1.732050807568877*leftProj[4]*p11th)/sqp0; 
  wv[RHOV] = -(1.732050807568877*leftProj[4]*sqp11*p12)/sqp0; 
  wv[RHOW] = -(1.732050807568877*leftProj[4]*sqp11*p13)/sqp0; 
  wv[PXX] = 3.0*leftProj[4]*p11sq; 
  wv[PXY] = 3.0*leftProj[4]*p11*p12; 
  wv[PXZ] = 3.0*leftProj[4]*p11*p13; 
  wv[PYY] = leftProj[4]*p11*p22+2.0*leftProj[4]*p12sq; 
  wv[PYZ] = leftProj[4]*p11*p23+2.0*leftProj[4]*p12*p13; 
  wv[PZZ] = leftProj[4]*p11*p33+2.0*leftProj[4]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[20]);

  // Wave 4 (ev 6)
  s[3] = u1+sqrt(3*p11/p0);
  wv[0] = leftProj[5]*p0*p11; 
  wv[RHOU] = (1.732050807568877*leftProj[5]*p11th)/sqp0; 
  wv[RHOV] = (1.732050807568877*leftProj[5]*sqp11*p12)/sqp0; 
  wv[RHOW] = (1.732050807568877*leftProj[5]*sqp11*p13)/sqp0; 
  wv[PXX] = 3.0*leftProj[5]*p11sq; 
  wv[PXY] = 3.0*leftProj[5]*p11*p12; 
  wv[PXZ] = 3.0*leftProj[5]*p11*p13; 
  wv[PYY] = leftProj[5]*p11*p22+2.0*leftProj[5]*p12sq; 
  wv[PYZ] = leftProj[5]*p11*p23+2.0*leftProj[5]*p12*p13; 
  wv[PZZ] = leftProj[5]*p11*p33+2.0*leftProj[5]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[30]);

  // Wave 5: (ev 7, 8, 9, 10 are repeated)
  s[4] = u1;
  wv[0] = leftProj[6]; 
  wv[RHOU] = 0.0; 
  wv[RHOV] = 0.0; 
  wv[RHOW] = 0.0; 
  wv[PXX] = 0.0; 
  wv[PXY] = 0.0; 
  wv[PXZ] = 0.0; 
  wv[PYY] = leftProj[7]; 
  wv[PYZ] = leftProj[8]; 
  wv[PZZ] = leftProj[9];

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[40]);
  
  return fabs(u1)+sqrt(3*p11/p0);
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[10], *w2 = &waves[20], *w3 = &waves[30], *w4 = &waves[40];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]), s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]), s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]);

  for (int i=0; i<10; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i] + s3m*w3[i] + s4m*w4[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i] + s3p*w3[i] + s4p*w4[i];
  }
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, int dir, const double *q)
{
  const struct wv_ten_moment *ten_moment = container_of(eqn, struct wv_ten_moment, eqn);
  return gkyl_ten_moment_max_abs_speed(dir, q);
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_new()
{
  struct wv_ten_moment *ten_moment = gkyl_malloc(sizeof(struct wv_ten_moment));

  ten_moment->eqn.num_equations = 10;
  ten_moment->eqn.num_waves = 5;
  ten_moment->eqn.waves_func = wave_roe;
  ten_moment->eqn.qfluct_func = qfluct_roe;

  ten_moment->eqn.ref_count = (struct gkyl_ref_count) { ten_moment_free, 1 };

  return &ten_moment->eqn;
}
