#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_em_sources.h>

// Makes indexing cleaner
static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

struct gkyl_moment_em_sources {
    struct gkyl_rect_grid grid; // grid object
    int ndim; // number of dimensions
    int nfluids; // number of fluids in multi-fluid system
    struct gkyl_moment_em_sources_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters
    double epsilon0; // permittivity of free space
    int has_pressure; // flag to update energy based on updated kinetic energy
};


gkyl_moment_em_sources*
gkyl_moment_em_sources_new(struct gkyl_moment_em_sources_inp inp)
{
  gkyl_moment_em_sources *up = gkyl_malloc(sizeof(gkyl_moment_em_sources));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->nfluids = inp.nfluids;
  for (int n=0; n<inp.nfluids; ++n) up->param[n] = inp.param[n];
  up->epsilon0 = inp.epsilon0;
  up->has_pressure = inp.has_pressure;

  return up;
}

static void
em_source_update(const gkyl_moment_em_sources *mes, double dt, double* fluids[], double* em)
{
  // based on Smithe (2007) with corrections but using Hakim (2019) notations
  // full reference in Wang, Hakim, Ng, Dong, & Germaschewski JCP 2020
  // implementation follow Wang et al. Appendix C, equation numbers referenced henceforth
  int nfluids = mes->nfluids;
  double epsilon0 = mes->epsilon0;
  double b[3];
  double Bmag = sqrt(em[BX]*em[BX] + em[BY]*em[BY] + em[BZ]*em[BZ]);
  // get magnetic field unit vector 
  if (Bmag > 0.0) {
    b[0] = em[BX]/Bmag;
    b[1] = em[BY]/Bmag;
    b[2] = em[BZ]/Bmag;
  }
  else {
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
  }
  double qbym[nfluids], Wc_dt[nfluids], wp_dt2[nfluids];
  double J[nfluids][3];
  double w02 = 0.0;
  double gam2 = 0.0;
  double delta = 0.0;
  double K[3];
  for (unsigned i=0; i < 3; ++i) K[i] = 0.0;

  for (int n=0; n < nfluids; ++n)
  {
    qbym[n] = mes->param[n].charge / mes->param[n].mass;
    double *f = fluids[n];

    J[n][0] = f[MX] * qbym[n];
    J[n][1] = f[MY] * qbym[n];
    J[n][2] = f[MZ] * qbym[n];
    // cyclotron frequency * dt
    Wc_dt[n] = qbym[n] * Bmag * dt;
    // plasma frequency^2 * dt^2
    wp_dt2[n] = f[RHO] * qbym[n] * qbym[n] * dt * dt / epsilon0;

    double tmp = 1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0);
    // w02, gam2, & delta = Eq. C.11 
    w02 += wp_dt2[n] / tmp;
    gam2 += wp_dt2[n] * Wc_dt[n] * Wc_dt[n] / tmp;
    delta += wp_dt2[n] * Wc_dt[n] / tmp;
    // K = Eq. C.9
    K[0] -= dt / tmp * (J[n][0] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[0] * (b[0]*J[n][0] + b[1]*J[n][1] + b[2]*J[n][2]) 
      - (Wc_dt[n] / 2.0) * (b[1]*J[n][2] - b[2]*J[n][1])
    );
    K[1] -= dt / tmp * (J[n][1] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[1] * (b[0]*J[n][0] + b[1]*J[n][1] + b[2]*J[n][2]) 
      - (Wc_dt[n] / 2.0) * (b[2]*J[n][0] - b[0]*J[n][2])
    );
    K[2] -= dt / tmp * (J[n][2] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[2] * (b[0]*J[n][0] + b[1]*J[n][1] + b[2]*J[n][2]) 
      - (Wc_dt[n] / 2.0) * (b[0]*J[n][1] - b[1]*J[n][0])
    );
  }
  // Delta2 (capital Delta) = Eq. C.11
  double Delta2 = delta * delta / (1.0 + w02 / 4.0);

  double F[3], F_halfK[3], Fbar[3];
  F[0] = em[EX] * epsilon0;
  F[1] = em[EY] * epsilon0;
  F[2] = em[EZ] * epsilon0;
  F_halfK[0] = F[0] + 0.5*K[0];
  F_halfK[1] = F[1] + 0.5*K[1];
  F_halfK[2] = F[2] + 0.5*K[2];
  // Fbar = Eq. C.10
  Fbar[0] = 1.0 / (1.0 + w02 / 4.0 + Delta2 / 64.0) * (F_halfK[0]
    + ((Delta2 / 64. - gam2 / 16.) / (1. + w02 / 4. + gam2 / 16.)) * b[0] * (b[0]*F_halfK[0] + b[1]*F_halfK[1] + b[2]*F_halfK[2])
    + (delta / 8. / (1. + w02 / 4.)) * (b[1]*F_halfK[2] - b[2]*F_halfK[1])
  );
  Fbar[1] = 1.0 / (1.0 + w02 / 4.0 + Delta2 / 64.0) * (F_halfK[1]
    + ((Delta2 / 64. - gam2 / 16.) / (1. + w02 / 4. + gam2 / 16.)) * b[1] * (b[0]*F_halfK[0] + b[1]*F_halfK[1] + b[2]*F_halfK[2])
    + (delta / 8. / (1. + w02 / 4.)) * (b[2]*F_halfK[0] - b[0]*F_halfK[2])
  );
  Fbar[2] = 1.0 / (1.0 + w02 / 4.0 + Delta2 / 64.0) * (F_halfK[2]
    + ((Delta2 / 64. - gam2 / 16.) / (1. + w02 / 4. + gam2 / 16.)) * b[2] * (b[0]*F_halfK[0] + b[1]*F_halfK[1] + b[2]*F_halfK[2])
    + (delta / 8. / (1. + w02 / 4.)) * (b[0]*F_halfK[1] - b[1]*F_halfK[0])
  );

  em[EX] = (2.0 * Fbar[0] - F[0])/epsilon0;
  em[EY] = (2.0 * Fbar[1] - F[1])/epsilon0;
  em[EZ] = (2.0 * Fbar[2] - F[2])/epsilon0;

  double Jstar[3], J_new[3];
  for (int n=0; n < nfluids; ++n)
  {
    double *f = fluids[n];

    // Jstar = Eq. C.7
    Jstar[0] = J[n][0] + Fbar[0] * (wp_dt2[n] / dt / 2.0);
    Jstar[1] = J[n][1] + Fbar[1] * (wp_dt2[n] / dt / 2.0);
    Jstar[2] = J[n][2] + Fbar[2] * (wp_dt2[n] / dt / 2.0);
    // J_new = 2 * Eq. C.8 - J_n
    // where Eq. C.8 is Jbar
    J_new[0] = 2.0 * (Jstar[0] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[0] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[1]*Jstar[2] - b[2]*Jstar[1])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0))
      - J[n][0];
    J_new[1] = 2.0 * (Jstar[1] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[1] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[2]*Jstar[0] - b[0]*Jstar[2])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0)) 
      - J[n][1];
    J_new[2] = 2.0 * (Jstar[2] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[2] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[0]*Jstar[1] - b[1]*Jstar[0])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0)) 
      - J[n][2];

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];
  } 
}

static void
pressure_em_source_update(const gkyl_moment_em_sources *mes, double dt, double* fluids[], double* em)
{
  // Source updater that updates the electric field and fluid momenta,
  // along with the fluid pressure Momenta update only affects kinetic
  // energy; need to separate kinetic from internal energy
  int nfluids = mes->nfluids;
  double keOld[nfluids], keNew[nfluids];
  for (int n=0; n < nfluids; ++n)
  {
    const double *f = fluids[n];
    keOld[n] = 0.5 * (f[MX]*f[MX] + f[MY]*f[MY] + f[MZ]*f[MZ]) / f[RHO];
  }

  em_source_update(mes, dt, fluids, em);

  for (int n=0; n < nfluids; ++n)
  {
    double *f = fluids[n];
    keNew[n] = 0.5 * (f[MX]*f[MX] + f[MY]*f[MY] + f[MZ]*f[MZ]) / f[RHO];
    f[ER] += keNew[n] - keOld[n];
  }
}

void
gkyl_moment_em_sources_advance(const gkyl_moment_em_sources *mes, double dt,
  const struct gkyl_range *update_range, struct gkyl_array *fluid[], struct gkyl_array *em)
{

  int ndim = mes->ndim;
  int nfluids = mes->nfluids;
  int has_pressure = mes->has_pressure;
  double *fluids[nfluids];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    
    long lidx = gkyl_range_idx(update_range, iter.idx);
    
    for (int n=0; n<nfluids; ++n)
      fluids[n] = gkyl_array_fetch(fluid[n], lidx);

    if (has_pressure)
      pressure_em_source_update(mes, dt, fluids, gkyl_array_fetch(em, lidx));
    else
      em_source_update(mes, dt, fluids, gkyl_array_fetch(em, lidx));
  }
}

void
gkyl_moment_em_sources_release(gkyl_moment_em_sources* up)
{
  free(up);
}
