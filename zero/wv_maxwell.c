#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_maxwell.h>

static const int dir_shuffle[][6] = {
  {0, 1, 2, 3, 4, 5},
  {1, 2, 0, 4, 5, 3},
  {2, 0, 1, 5, 3, 4}
};

// Make indexing cleaner with the dir_shuffle
#define EX d[0]
#define EY d[1]
#define EZ d[2]
#define BX d[3]
#define BY d[4]
#define BZ d[5]

struct wv_maxwell {
  struct gkyl_wv_eqn eqn; // base object
  double c; // light speed
  double e_fact, b_fact; // electric and magnetic correction factors
};

static void
maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_maxwell *maxwell = container_of(base, struct wv_maxwell, eqn);
  gkyl_free(maxwell);
}

static void
rot_to_local_rect(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *qglobal, double *qlocal)
{
  const int *d = dir_shuffle[dir];  
  qlocal[0] = qglobal[EX];
  qlocal[1] = qglobal[EY];
  qlocal[2] = qglobal[EZ];
  qlocal[3] = qglobal[BX];
  qlocal[4] = qglobal[BY];
  qlocal[5] = qglobal[BZ];
  qlocal[6] = qglobal[6];
  qlocal[7] = qglobal[7];
}

static void
rot_to_global_rect(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *qlocal, double *qglobal)
{
  const int *d = dir_shuffle[dir];  
  qglobal[EX] = qlocal[0];
  qglobal[EY] = qlocal[1];
  qglobal[EZ] = qlocal[2];
  qglobal[BX] = qlocal[3];
  qglobal[BY] = qlocal[4];
  qglobal[BZ] = qlocal[5];
  qglobal[6] = qlocal[6];
  qglobal[7] = qlocal[7];
}

static inline void
rot_to_local(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Rotate E to local coordinates
  qlocal[0] = qglobal[0]*norm[0] + qglobal[1]*norm[1] + qglobal[2]*norm[2];
  qlocal[1] = qglobal[0]*tau1[0] + qglobal[1]*tau1[1] + qglobal[2]*tau1[2];
  qlocal[2] = qglobal[0]*tau2[0] + qglobal[1]*tau2[1] + qglobal[2]*tau2[2];
  // Rotate B to local coordinates
  qlocal[3] = qglobal[3]*norm[0] + qglobal[4]*norm[1] + qglobal[5]*norm[2];
  qlocal[4] = qglobal[3]*tau1[0] + qglobal[4]*tau1[1] + qglobal[5]*tau1[2];
  qlocal[5] = qglobal[3]*tau2[0] + qglobal[4]*tau2[1] + qglobal[5]*tau2[2];
  // Correction potentials are scalars and unchanged
  qlocal[6] = qglobal[6];
  qlocal[7] = qglobal[7];
}

static inline void
rot_to_global(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Rotate E back to global coordinates
  qglobal[0] = qlocal[0]*norm[0] + qlocal[1]*tau1[0] + qlocal[2]*tau2[0];
  qglobal[1] = qlocal[0]*norm[1] + qlocal[1]*tau1[1] + qlocal[2]*tau2[1];
  qglobal[2] = qlocal[0]*norm[2] + qlocal[1]*tau1[2] + qlocal[2]*tau2[2];
  // Rotate B back to global coordinates
  qglobal[3] = qlocal[3]*norm[0] + qlocal[4]*tau1[0] + qlocal[5]*tau2[0];
  qglobal[4] = qlocal[3]*norm[1] + qlocal[4]*tau1[1] + qlocal[5]*tau2[1];
  qglobal[5] = qlocal[3]*norm[2] + qlocal[4]*tau1[2] + qlocal[5]*tau2[2];
  // Correction potentials are scalars and unchanged
  qglobal[6] = qlocal[6];
  qglobal[7] = qlocal[7];
}

// Waves and speeds using Roe averaging
static double
wave(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);

  double c = maxwell->c, c1 = 1/c;
  double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;
    
  // compute projections of jump
  // Note correction potentials are scalars and no dir_shuffle required
  double a1 = 0.5*(delta[3]-delta[7]*c1);
  double a2 = 0.5*(delta[3]+delta[7]*c1);
  double a3 = 0.5*(delta[0]-delta[6]*c);
  double a4 = 0.5*(delta[0]+delta[6]*c);
  double a5 = 0.5*(delta[1]-delta[5]*c);
  double a6 = 0.5*(delta[4]*c+delta[2]);
  double a7 = 0.5*(delta[5]*c+delta[1]);
  double a8 = 0.5*(delta[2]-delta[4]*c);

  // set waves to 0.0 as most entries vanish
  for (int i=0; i<8*6; ++i) waves[i] = 0.0;

  double *w = 0;

  // wave 1:
  w = &waves[0*8];
  w[3] = a1;
  w[7] = -a1*c;
  s[0] = -c*b_fact;

  // wave 2:
  w = &waves[1*8];
  w[3] = a2;
  w[7] = a2*c;
  s[1] = c*b_fact;

  // wave 3:
  w = &waves[2*8];
  w[0] = a3;
  w[6] = -a3*c1;
  s[2] = -c*e_fact;

  // wave 4:
  w = &waves[3*8];
  w[0] = a4;
  w[6] = a4*c1;
  s[3] = c*e_fact;

  // wave 5: (two waves with EV -c, -c lumped into one)
  w = &waves[4*8];
  w[1] = a5;
  w[2] = a6;
  w[4] = a6*c1;
  w[5] = -a5*c1;
  s[4] = -c;

  // wave 6: (two waves with EV c, c lumped into one)
  w = &waves[5*8];
  w[1] = a7;
  w[2] = a8;
  w[4] = -a8*c1;
  w[5] = a7*c1;
  s[5] = c;
  
  return c;
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0*8], *w1 = &waves[1*8], *w2 = &waves[2*8];
  const double *w3 = &waves[3*8], *w4 = &waves[4*8], *w5 = &waves[5*8];
  
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  for (int i=0; i<8; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i] + s3m*w3[i] + s4m*w4[i] + s5m*w5[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i] + s3p*w3[i] + s4p*w4[i] + s5p*w5[i];
  }
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, int dir, const double *q)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  return maxwell->c;
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_new(double c, double e_fact, double b_fact)
{
  struct wv_maxwell *maxwell = gkyl_malloc(sizeof(struct wv_maxwell));

  maxwell->eqn.type = GKYL_EQN_MAXWELL;
  maxwell->eqn.num_equations = 8;
  maxwell->eqn.num_waves = 6;
  
  maxwell->c = c;
  maxwell->e_fact = e_fact;
  maxwell->b_fact = b_fact;
  
  maxwell->eqn.waves_func = wave;
  maxwell->eqn.qfluct_func = qfluct;
  maxwell->eqn.max_speed_func = max_speed;
  maxwell->eqn.rotate_to_local_func = rot_to_local_rect;
  maxwell->eqn.rotate_to_global_func = rot_to_global_rect;

  maxwell->eqn.ref_count = (struct gkyl_ref_count) { maxwell_free, 1 };

  return &maxwell->eqn;
}
