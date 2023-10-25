#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_mat.h>

// Makes indexing cleaner
static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

static const unsigned P11 = 4;
static const unsigned P12 = 5;
static const unsigned P13 = 6;
static const unsigned P22 = 7;
static const unsigned P23 = 8;
static const unsigned P33 = 9;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

#define INTERNAL_ENERGY (4)

#define sq(x) ((x) * (x))

struct gkyl_moment_em_coupling {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int nfluids; // number of fluids in multi-fluid system
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters

  bool is_charged_species; 
  double epsilon0; // permittivity of free space

  bool ramp_app_E;
  double t_ramp_E; // ramp up time for external electric field
  bool ramp_app_curr;
  double t_ramp_curr; // ramp up time for applied current

  bool has_collision; // has collisions
  // normalized collision frequencies; nu_sr = nu_base[s][r] * rho_r
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool has_nT_sources;
};

// Rotate pressure tensor using magnetic field. See Wang
// et. al. 2020 for details.
static void
pressure_tensor_rotate(double qbym, double dt, const double *em, const double *ext_em,
  double prOld[6], double prRhs[6], double prOut[6])
{
  double Bx = em[BX] + ext_em[BX];
  double By = em[BY] + ext_em[BY];
  double Bz = em[BZ] + ext_em[BZ];

  double dt1 = 0.5 * dt;
  double dtsq = dt1 * dt1;
  double dt3 = dtsq * dt1;
  double dt4 = dt3 * dt1;
  double Bx2 = Bx * Bx;
  double Bx3 = Bx * Bx2;
  double Bx4 = Bx * Bx3;
  double By2 = By * By;
  double By3 = By * By2;
  double By4 = By * By3;
  double Bz2 = Bz * Bz;
  double Bz3 = Bz * Bz2;
  double Bz4 = Bz * Bz3;

  double qb2 = qbym * qbym;
  double qb3 = qbym * qb2;
  double qb4 = qbym * qb3;
  double d = 1 + 5*(Bx2 + By2 + Bz2)*dtsq*qb2 + 4*(Bx2 + By2 + Bz2)*(Bx2 + By2 + Bz2)*dt4*qb4;

  prOut[0] = 2.0*(prRhs[0] + 2*dt1*(Bz*prRhs[1] - By*prRhs[2])*qbym + dtsq*(5*Bx2*prRhs[0] + 2*Bx*(By*prRhs[1] + Bz*prRhs[2]) + Bz2*(3*prRhs[0] + 2*prRhs[3]) - 4*By*Bz*prRhs[4] + By2*(3*prRhs[0] + 2*prRhs[5]))*qb2 + 2*dt3*(4*Bx2*(Bz*prRhs[1] - By*prRhs[2]) - (By2 + Bz2)*(-(Bz*prRhs[1]) + By*prRhs[2]) - 3*Bx*(By2*prRhs[4] - Bz2*prRhs[4] + By*Bz*(-prRhs[3] + prRhs[5])))*qb3 + 2*dt4*(2*Bx4*prRhs[0] + 4*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 + Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + (By2 + Bz2)*(Bz2*(prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[5])) + Bx2*(4*By*Bz*prRhs[4] + By2*(3*prRhs[3] + prRhs[5]) + Bz2*(prRhs[3] + 3*prRhs[5])))*qb4)/d - prOld[0];
  prOut[1] =  2.0*(prRhs[1] + dt1*(Bx*prRhs[2] + Bz*(-prRhs[0] + prRhs[3]) - By*prRhs[4])*qbym + dtsq*(4*Bx2*prRhs[1] + 4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2] + Bx*(3*Bz*prRhs[4] + By*(prRhs[0] + prRhs[3] - 2*prRhs[5])))*qb2 + dt3*(4*Bx3*prRhs[2] - 2*Bx*(By2 + Bz2)*prRhs[2] + Bz3*(-prRhs[0] + prRhs[3]) - 4*By3*prRhs[4] + 2*By*Bz2*prRhs[4] - By2*Bz*(prRhs[0] - 4*prRhs[3] + 3*prRhs[5]) + Bx2*(2*By*prRhs[4] + Bz*(-4*prRhs[0] + prRhs[3] + 3*prRhs[5])))*qb3 + 2*Bx*By*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d - prOld[1];
  prOut[2] =  2.0*(prRhs[2] + dt1*(-(Bx*prRhs[1]) + Bz*prRhs[4] + By*(prRhs[0] - prRhs[5]))*qbym + dtsq*(3*By*Bz*prRhs[1] + 4*Bx2*prRhs[2] + By2*prRhs[2] + 4*Bz2*prRhs[2] + Bx*(3*By*prRhs[4] + Bz*(prRhs[0] - 2*prRhs[3] + prRhs[5])))*qb2 + dt3*(-4*Bx3*prRhs[1] + 2*Bx*(By2 + Bz2)*prRhs[1] - 2*By2*Bz*prRhs[4] + 4*Bz3*prRhs[4] + By*Bz2*(prRhs[0] + 3*prRhs[3] - 4*prRhs[5]) + By3*(prRhs[0] - prRhs[5]) - Bx2*(2*Bz*prRhs[4] + By*(-4*prRhs[0] + 3*prRhs[3] + prRhs[5])))*qb3 + 2*Bx*Bz*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d - prOld[2];
  prOut[3] =  2.0*(prRhs[3] + (-2*Bz*dt1*prRhs[1] + 2*Bx*dt1*prRhs[4])*qbym + dtsq*(2*Bx*By*prRhs[1] + 5*By2*prRhs[3] + Bz2*(2*prRhs[0] + 3*prRhs[3]) + Bz*(-4*Bx*prRhs[2] + 2*By*prRhs[4]) + Bx2*(3*prRhs[3] + 2*prRhs[5]))*qb2 + 2*dt3*(Bx2*(-(Bz*prRhs[1]) + 3*By*prRhs[2]) - Bz*(4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2]) + Bx3*prRhs[4] + Bx*(4*By2*prRhs[4] + Bz2*prRhs[4] + 3*By*Bz*(-prRhs[0] + prRhs[5])))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) + 2*Bx*(2*By2 - Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + 2*By4*prRhs[3] + Bz4*(prRhs[0] + prRhs[3]) + 4*By3*Bz*prRhs[4] - 2*By*Bz3*prRhs[4] + Bx4*(prRhs[3] + prRhs[5]) + By2*Bz2*(prRhs[0] + 3*prRhs[5]) + Bx2*(-2*By*Bz*prRhs[4] + By2*(3*prRhs[0] + prRhs[5]) + Bz2*(prRhs[0] + 2*prRhs[3] + prRhs[5])))*qb4)/d - prOld[3];
  prOut[4] =  2.0*(prRhs[4] + dt1*(By*prRhs[1] - Bz*prRhs[2] + Bx*(-prRhs[3] + prRhs[5]))*qbym + dtsq*(3*Bx*Bz*prRhs[1] + Bx2*prRhs[4] + 4*By2*prRhs[4] + 4*Bz2*prRhs[4] + By*(3*Bx*prRhs[2] + Bz*(-2*prRhs[0] + prRhs[3] + prRhs[5])))*qb2 + dt3*(4*By3*prRhs[1] - 2*By*Bz2*prRhs[1] + 2*By2*Bz*prRhs[2] - 4*Bz3*prRhs[2] + Bx2*(-2*By*prRhs[1] + 2*Bz*prRhs[2]) + Bx3*(-prRhs[3] + prRhs[5]) + Bx*(-(Bz2*(3*prRhs[0] + prRhs[3] - 4*prRhs[5])) + By2*(3*prRhs[0] - 4*prRhs[3] + prRhs[5])))*qb3 - 2*By*Bz*dt4*(-6*Bx*(By*prRhs[1] + Bz*prRhs[2]) - 6*By*Bz*prRhs[4] + Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]) + Bx2*(-2*prRhs[0] + prRhs[3] + prRhs[5]))*qb4)/d - prOld[4];
  prOut[5] =  2.0*(prRhs[5] + 2*dt1*(By*prRhs[2] - Bx*prRhs[4])*qbym + dtsq*(2*Bx*Bz*prRhs[2] + By*(-4*Bx*prRhs[1] + 2*Bz*prRhs[4]) + 5*Bz2*prRhs[5] + By2*(2*prRhs[0] + 3*prRhs[5]) + Bx2*(2*prRhs[3] + 3*prRhs[5]))*qb2 - 2*dt3*(Bx2*(3*Bz*prRhs[1] - By*prRhs[2]) - By*(3*By*Bz*prRhs[1] + By2*prRhs[2] + 4*Bz2*prRhs[2]) + Bx3*prRhs[4] + Bx*(3*By*Bz*(-prRhs[0] + prRhs[3]) + By2*prRhs[4] + 4*Bz2*prRhs[4]))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 - 2*Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + By2*Bz2*(prRhs[0] + 3*prRhs[3]) - 2*By3*Bz*prRhs[4] + 4*By*Bz3*prRhs[4] + 2*Bz4*prRhs[5] + By4*(prRhs[0] + prRhs[5]) + Bx4*(prRhs[3] + prRhs[5]) + Bx2*(Bz2*(3*prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[3] + 2*prRhs[5])))*qb4)/d - prOld[5];
}

// Update momentum and E field using time-centered scheme. See Wang
// et. al. 2020 for details.
static void
em_source_update(const gkyl_moment_em_coupling *mes, double tcurr, double dt, 
  double *fluids[GKYL_MAX_SPECIES], const double *app_accels[GKYL_MAX_SPECIES],
  double *em, const double *app_current, const double *ext_em)
{
  // based on Smithe (2007) with corrections but using Hakim (2019) notations
  // full reference in Wang, Hakim, Ng, Dong, & Germaschewski JCP 2020
  // implementation follow Wang et al. Appendix D, equation numbers referenced henceforth
  int nfluids = mes->nfluids;
  double epsilon0 = mes->epsilon0;
  double b[3] = { 0.0, 0.0, 0.0 };
  double Bx = em[BX] + ext_em[BX];
  double By = em[BY] + ext_em[BY];
  double Bz = em[BZ] + ext_em[BZ];
  double Bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
  // get magnetic field unit vector 
  if (Bmag > 0.0) {
    b[0] = Bx/Bmag;
    b[1] = By/Bmag;
    b[2] = Bz/Bmag;
  }
  double qbym[GKYL_MAX_SPECIES], Wc_dt[GKYL_MAX_SPECIES], wp_dt2[GKYL_MAX_SPECIES];
  double J[GKYL_MAX_SPECIES][3];
  double JOld[GKYL_MAX_SPECIES][3];
  double w02 = 0.0;
  double gam2 = 0.0;
  double delta = 0.0;
  double K[3] = { 0.0, 0.0, 0.0 };

  // Scale factors if we are ramping E fields or currents
  double scale_fac_E = 1.0;
  double scale_fac_curr = 1.0;

  for (int n=0; n < nfluids; ++n) {
    qbym[n] = mes->param[n].charge / mes->param[n].mass;
    const double *f = fluids[n];
    const double *app_accel = app_accels[n];

    JOld[n][0] = f[MX] * qbym[n];
    JOld[n][1] = f[MY] * qbym[n];
    JOld[n][2] = f[MZ] * qbym[n];

    if (mes->ramp_app_E)
      scale_fac_E = fmin(1.0, tcurr/mes->t_ramp_E);
    // Add contributions from static electric field and applied acceleration to current
    // Note: Conversion from applied acceleration to force density occurs here
    J[n][0] = JOld[n][0] + 0.5*dt*qbym[n]*f[RHO]*(qbym[n]*ext_em[EX]*scale_fac_E + app_accel[0]);
    J[n][1] = JOld[n][1] + 0.5*dt*qbym[n]*f[RHO]*(qbym[n]*ext_em[EY]*scale_fac_E + app_accel[1]);
    J[n][2] = JOld[n][2] + 0.5*dt*qbym[n]*f[RHO]*(qbym[n]*ext_em[EZ]*scale_fac_E + app_accel[2]);

    // cyclotron frequency * dt
    Wc_dt[n] = qbym[n] * Bmag * dt;
    // plasma frequency^2 * dt^2
    wp_dt2[n] = f[RHO] * qbym[n] * qbym[n] * dt * dt / epsilon0;

    double tmp = 1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0);
    // w02, gam2, & delta = Eq. D.11 
    w02 += wp_dt2[n] / tmp;
    gam2 += wp_dt2[n] * Wc_dt[n] * Wc_dt[n] / tmp;
    delta += wp_dt2[n] * Wc_dt[n] / tmp;
    // K = Eq. D.9
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
  // Delta2 (capital Delta) = Eq. D.11
  double Delta2 = delta * delta / (1.0 + w02 / 4.0);

  double FOld[3], F[3], F_halfK[3], Fbar[3];

  FOld[0] = em[EX] * epsilon0;
  FOld[1] = em[EY] * epsilon0;
  FOld[2] = em[EZ] * epsilon0;

  if (mes->ramp_app_curr)
    scale_fac_curr = fmin(1.0, tcurr/mes->t_ramp_curr);
  // Add contributions from applied currents to electric field
  // Note: Conversion from applied current to electric field units occurs here
  F[0] = FOld[0] - 0.5*dt*app_current[0]/epsilon0*scale_fac_curr;
  F[1] = FOld[1] - 0.5*dt*app_current[1]/epsilon0*scale_fac_curr;
  F[2] = FOld[2] - 0.5*dt*app_current[2]/epsilon0*scale_fac_curr;

  F_halfK[0] = F[0] + 0.5*K[0];
  F_halfK[1] = F[1] + 0.5*K[1];
  F_halfK[2] = F[2] + 0.5*K[2];
  // Fbar = Eq. D.10
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

  em[EX] = (2.0 * Fbar[0] - FOld[0])/epsilon0;
  em[EY] = (2.0 * Fbar[1] - FOld[1])/epsilon0;
  em[EZ] = (2.0 * Fbar[2] - FOld[2])/epsilon0;

  double Jstar[3], J_new[3];
  for (int n=0; n < nfluids; ++n)
  {
    double *f = fluids[n];

    // Jstar = Eq. D.7
    Jstar[0] = J[n][0] + Fbar[0] * (wp_dt2[n] / dt / 2.0);
    Jstar[1] = J[n][1] + Fbar[1] * (wp_dt2[n] / dt / 2.0);
    Jstar[2] = J[n][2] + Fbar[2] * (wp_dt2[n] / dt / 2.0);
    // J_new = 2 * Eq. D.8 - J_n
    // where Eq. D.8 is Jbar
    J_new[0] = 2.0 * (Jstar[0] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[0] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[1]*Jstar[2] - b[2]*Jstar[1])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0))
      - JOld[n][0];
    J_new[1] = 2.0 * (Jstar[1] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[1] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[2]*Jstar[0] - b[0]*Jstar[2])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0)) 
      - JOld[n][1];
    J_new[2] = 2.0 * (Jstar[2] 
      + (Wc_dt[n] * Wc_dt[n] / 4.0) * b[2] * (b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2]) 
      - (Wc_dt[n] / 2.0) * (b[0]*Jstar[1] - b[1]*Jstar[0])) / (1.0 + (Wc_dt[n] * Wc_dt[n] / 4.0)) 
      - JOld[n][2];

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];
  } 
}

// Update momentum of neutral fluid using time-centered scheme.
// Assumes no EM fields in the simulation. See Wang
// et. al. 2020 for details.
static void
neut_source_update(const gkyl_moment_em_coupling *mes, double tcurr, double dt,
  double *fluids[GKYL_MAX_SPECIES], const double *app_accels[GKYL_MAX_SPECIES])
{
  int nfluids = mes->nfluids;
  for (int n=0; n < nfluids; ++n)
  {
    double *f = fluids[n];
    const double *app_accel = app_accels[n];

    // Add contributions from applied acceleration to momentum
    // Note: Conversion from applied acceleration to force density occurs here
    f[MX] += dt*f[RHO]*app_accel[0];
    f[MY] += dt*f[RHO]*app_accel[1];
    f[MZ] += dt*f[RHO]*app_accel[2];
  }
}

static inline double
internal_energy(const double *f)
{
  return f[ER] - 0.5*(sq(f[MX])+sq(f[MY])+sq(f[MZ]))/f[RHO];
}


static void
calcNu(const gkyl_moment_em_coupling *mes,
  double **fluids,
  const double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES],
  double *nu)
{
  int nfluids = mes->nfluids;

  // TODO temperature dependence? user-defined calcNu function?
  for (int s=0; s<nfluids; ++s) {
    double *nu_s = nu + nfluids * s;
    for (int r=0; r<nfluids; ++r)
      nu_s[r] = nu_base[s][r] * fluids[r][RHO];
  }
}


static void
collision_source_update(const gkyl_moment_em_coupling *mes, double dt,
  double *fluids[GKYL_MAX_SPECIES])
{
  int nfluids = mes->nfluids;
  double nu[nfluids * nfluids];

  /* STEP 0: CALCULATE INTER-SPECIES COLLISION FREQUENCIES */
  calcNu(mes, fluids, mes->nu_base, nu);

  /* STEP 1: UPDATE VELOCITIES */
  // the lhs matrix is identical for all three (i.e., x,y,z) components
  struct gkyl_mat *lhs = gkyl_mat_new(nfluids, nfluids, 0.0);
  // rhs has 3 column vectors, each has nfluids components <-> nfluids equations
  struct gkyl_mat *rhs = gkyl_mat_new(nfluids, 3, 0.0);

  for (int s=0; s<nfluids; ++s) {
    double *fs = fluids[s];
    gkyl_mat_set(rhs, s, 0, fs[MX] / fs[RHO]);
    gkyl_mat_set(rhs, s, 1, fs[MY] / fs[RHO]);
    gkyl_mat_set(rhs, s, 2, fs[MZ] / fs[RHO]);

    gkyl_mat_set(lhs, s, s, 1.0);
    double *nu_s = nu + nfluids * s;
    for (int r=0; r<nfluids; ++r) {
      if (r==s) {
        gkyl_mat_inc(lhs, s, s, 0.5 * dt * nu_s[s]);
        continue;
      }

      double half_dt_nu_sr = 0.5 * dt * nu_s[r];
      gkyl_mat_inc(lhs, s, s, half_dt_nu_sr);
      gkyl_mat_inc(lhs, s, r, -half_dt_nu_sr);
    }
  }

  gkyl_mem_buff ipiv = gkyl_mem_buff_new(sizeof(long[nfluids]));
  // lhs and rhs will be rewritten by output, i.e., rhs is solution after solve
  bool status = gkyl_mat_linsolve_lu(lhs, rhs, gkyl_mem_buff_data(ipiv));

  /* STEP 2: UPDATE TEMPERATURE OR PRESSURE */
  if (true) {
    gkyl_mat_clear(lhs, 0.0);
    struct gkyl_mat *rhs_T = gkyl_mat_new(nfluids, 1, 0.0);

    double T[nfluids];
    for (int s=0; s<nfluids; ++s) {
      double ms = mes->param[s].mass;

      double *fs = fluids[s];
      T[s] = internal_energy(fs) / fs[RHO] * ms; // this is actually T/(gamma-1)
      gkyl_mat_set(rhs_T, s, 0, T[s]);
      gkyl_mat_set(lhs, s, s, 1.0);

      double *nu_s = nu + nfluids * s;
      for (int r=0; r<nfluids; ++r) {
        if (r==s)
          continue;

        double mr = mes->param[r].mass;
        // velocity difference due to velocities at t=n+1/2
        double du2 = sq(gkyl_mat_get(rhs,s,0)-gkyl_mat_get(rhs,r,0))
                   + sq(gkyl_mat_get(rhs,s,1)-gkyl_mat_get(rhs,r,1))
                   + sq(gkyl_mat_get(rhs,s,2)-gkyl_mat_get(rhs,r,2));
        double coeff_sr = dt * nu_s[r] * ms / (ms + mr);

        gkyl_mat_inc(rhs_T, s, 0, coeff_sr * 0.5 * mr * du2);
        gkyl_mat_inc(lhs, s, s, coeff_sr);
        gkyl_mat_inc(lhs, s, r, -coeff_sr);
      }
    }

    // Compute pressure at n+1/2 and then at n+1 using old velocities
    status = gkyl_mat_linsolve_lu(lhs, rhs_T, gkyl_mem_buff_data(ipiv));
    for (int s=0; s<nfluids; ++s) {
      double *f = fluids[s];
      f[INTERNAL_ENERGY] = (2.0 * gkyl_mat_get(rhs_T,s,0) - T[s]) * f[RHO] / mes->param[s].mass;
    }

    gkyl_mat_release(rhs_T);
  }
  else {
    for (int s=0; s<nfluids; ++s) {
      double *f = fluids[s];
      f[INTERNAL_ENERGY] = internal_energy(f);
    }
  }

  /* STEP 3: UPDATE MOMENTUM AND TOTAL ENERGY */
  for (int s=0; s<nfluids; ++s)
  {
    double *f = fluids[s];
    f[MX] = 2 * f[RHO] * gkyl_mat_get(rhs,s,0) - f[MX];
    f[MY] = 2 * f[RHO] * gkyl_mat_get(rhs,s,1) - f[MY];
    f[MZ] = 2 * f[RHO] * gkyl_mat_get(rhs,s,2) - f[MZ];
    f[ER] = f[INTERNAL_ENERGY] + 0.5*(sq(f[MX])+sq(f[MY])+sq(f[MZ]))/f[RHO];
  }

  gkyl_mem_buff_release(ipiv);
  gkyl_mat_release(lhs);
  gkyl_mat_release(rhs);
}

/***********************************************/
/* user-defined density and temperature source */
/***********************************************/

static void
nT_source_euler_update(const gkyl_moment_em_coupling *mes,
  const int s,
  const double dt,
  double *q,
  double *q1,
  const double *S)
{
  double m = mes->param[s].mass;

  double rho0 = q[RHO];
  double n0 = rho0 / m;
  double u = q[MX] / rho0;
  double v = q[MY] / rho0;
  double w = q[MZ] / rho0;
  double v2 = u*u + v*v + w*w;
  double TT0 = (q[ER] - 0.5 * rho0 * v2) / n0; // TT0=kB*T0/(gamma-1)

  double n1 = n0 + dt * S[0]; // source for number density n
  double TT1 = TT0 + dt * S[1]; // source for scaled temperature kB*T/(gamma-1)

  double rho1 = n1 * m;
  q[RHO] = rho1;
  q1[MX] = rho1 * u;
  q1[MY] = rho1 * v;
  q1[MZ] = rho1 * w;
  q1[ER] = n1 * TT1 + 0.5 * rho1 * v2;
}

void
nT_source_update(
  const gkyl_moment_em_coupling *mes,
  const double dt,
  double *qPtrs[],
  const double *sourcePtrs[])
{
  for (int s=0; s<mes->nfluids; ++s)
  {
    double *q = qPtrs[s];
    const double *S = sourcePtrs[s];
    nT_source_euler_update(mes, s, dt, q, q, S);
  }
}

/**********************************/
/* update all sources in one call */
/**********************************/


// Update momentum and E field along with appropriate pressure update.
// If isothermal Euler equations, only update momentum and E,
// if Euler equations, update energy based on new kinetic energy,
// if Ten moment equations, rotate pressure tensor around magnetic field.
// See Wang et. al. 2020 JCP for details
static void
fluid_source_update(const gkyl_moment_em_coupling *mes, double tcurr, double dt,
  double* fluids[GKYL_MAX_SPECIES], 
  const double *app_accels[GKYL_MAX_SPECIES], const double *rhs_s[GKYL_MAX_SPECIES], 
  double* em, const double* app_current, const double* ext_em,
  const double *nT_source_s[GKYL_MAX_SPECIES])
{
  int nfluids = mes->nfluids;
  double keOld[GKYL_MAX_SPECIES];
  double prInp[6], prRhs[6], prTen[GKYL_MAX_SPECIES][6];
  
  for (int n=0; n < nfluids; ++n) {
    const double *f = fluids[n];
    const double *rhs = rhs_s[n];
    // If Euler equations, store old value of the kinetic energy for later update.
    if (mes->param[n].type == GKYL_EQN_EULER) {
      keOld[n] = 0.5 * (f[MX]*f[MX] + f[MY]*f[MY] + f[MZ]*f[MZ]) / f[RHO];
    }
    // If Ten moment equations, rotate the pressure tensor based on local magnetic field.
    else if (mes->param[n].type == GKYL_EQN_TEN_MOMENT) {
      double qbym = mes->param[n].charge / mes->param[n].mass;
      prInp[0] = f[P11] - f[MX] * f[MX] / f[RHO];
      prInp[1] = f[P12] - f[MX] * f[MY] / f[RHO];
      prInp[2] = f[P13] - f[MX] * f[MZ] / f[RHO];
      prInp[3] = f[P22] - f[MY] * f[MY] / f[RHO];
      prInp[4] = f[P23] - f[MY] * f[MZ] / f[RHO];
      prInp[5] = f[P33] - f[MZ] * f[MZ] / f[RHO];

      prRhs[0] = prInp[0] + 0.5*dt*rhs[P11];
      prRhs[1] = prInp[1] + 0.5*dt*rhs[P12];
      prRhs[2] = prInp[2] + 0.5*dt*rhs[P13];
      prRhs[3] = prInp[3] + 0.5*dt*rhs[P22];
      prRhs[4] = prInp[4] + 0.5*dt*rhs[P23];
      prRhs[5] = prInp[5] + 0.5*dt*rhs[P33];

      double p = 1.0/3.0*(prInp[0] + prInp[3] + prInp[5]);
      double vth = sqrt(p/f[RHO]);
      double nu = vth*mes->param[n].k0;
      double edt = exp(nu*dt);


      // If permittivity of free-space is non-zero, EM fields exist
      // and we rotate pressure tensor around magnetic field 
      if (mes->is_charged_species)
        pressure_tensor_rotate(qbym, dt, em, ext_em, prInp, prRhs, prTen[n]);

      // Divide out integrating factor from collisional relaxation
      // (P - pI) = exp(-nu*dt)*(P - pI); where p is the isotropic pressure 
      prTen[n][0] = (prTen[n][0] - p)/edt + p;
      prTen[n][1] = prTen[n][1]/edt;
      prTen[n][2] = prTen[n][2]/edt;
      prTen[n][3] = (prTen[n][3] - p)/edt + p;
      prTen[n][4] = prTen[n][4]/edt;
      prTen[n][5] = (prTen[n][5] - p)/edt + p;
    }
  }

  // Update momentum and electric field using time-centered implicit solve.
  // If permittivity of free-space is non-zero, EM fields exist and electric
  // field is updated, otherwise just update momentum
  if (mes->is_charged_species)
    em_source_update(mes, tcurr, dt, fluids, app_accels, em, app_current, ext_em);
  else
    neut_source_update(mes, tcurr, dt, fluids, app_accels);

  for (int n=0; n < nfluids; ++n) {
    double *f = fluids[n];
    // If Euler equations, get updated energy based on new kinetic energy
    // from updated momentum (pressure should be unaffected by source terms).
    if (mes->param[n].type == GKYL_EQN_EULER) {
      f[ER] += 0.5 * (f[MX]*f[MX] + f[MY]*f[MY] + f[MZ]*f[MZ]) / f[RHO] - keOld[n];
    }
    // If Ten moment equations, update full stress tensor (pressure tensor + Reynolds stress)
    // based on rotated pressure tensor and updated momentum.
    else if (mes->param[n].type == GKYL_EQN_TEN_MOMENT) {
      f[P11] = f[MX] * f[MX] / f[RHO] + prTen[n][0];
      f[P12] = f[MX] * f[MY] / f[RHO] + prTen[n][1];
      f[P13] = f[MX] * f[MZ] / f[RHO] + prTen[n][2];
      f[P22] = f[MY] * f[MY] / f[RHO] + prTen[n][3];
      f[P23] = f[MY] * f[MZ] / f[RHO] + prTen[n][4];
      f[P33] = f[MZ] * f[MZ] / f[RHO] + prTen[n][5];
    }
  }

  if (mes->has_collision)
    collision_source_update(mes, dt, fluids);

  if (mes->has_nT_sources)
    nT_source_update(mes, dt, fluids, nT_source_s);
}

gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp)
{
  gkyl_moment_em_coupling *up = gkyl_malloc(sizeof(gkyl_moment_em_coupling));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->nfluids = inp.nfluids;
  for (int n=0; n<inp.nfluids; ++n) up->param[n] = inp.param[n];

  up->is_charged_species = false;
  up->epsilon0 = inp.epsilon0;
  if (up->epsilon0 != 0.0) 
    up->is_charged_species = true;

  up->ramp_app_E = false;
  if(inp.t_ramp_E != 0.0) {
    up->ramp_app_E = true;
    up->t_ramp_E = inp.t_ramp_E;
  }
  up->ramp_app_curr = false;
  if(inp.t_ramp_curr != 0.0) {
    up->ramp_app_curr = true;
    up->t_ramp_curr = inp.t_ramp_curr;
  }

  up->has_collision = inp.has_collision;
  if(inp.has_collision) {
    for (int s=0; s<inp.nfluids; ++s)
      for (int r=0; r<inp.nfluids; ++r)
        up->nu_base[s][r] = inp.nu_base[s][r];
  }

  up->has_nT_sources = inp.has_nT_sources;

  return up;
}

void
gkyl_moment_em_coupling_advance(const gkyl_moment_em_coupling *mes, double tcurr, double dt,
  const struct gkyl_range *update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], 
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *pr_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array *em, const struct gkyl_array *app_current, const struct gkyl_array *ext_em,
  const struct gkyl_array *nT_sources[GKYL_MAX_SPECIES])
{
  int nfluids = mes->nfluids;
  double *fluids[GKYL_MAX_SPECIES];
  const double *app_accels[GKYL_MAX_SPECIES];
  // pressure tensor rhs for gradient-based closure
  const double *pr_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_source_s[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    
    long lidx = gkyl_range_idx(update_range, iter.idx);
    
    for (int n=0; n<nfluids; ++n) {
      fluids[n] = gkyl_array_fetch(fluid[n], lidx);
      app_accels[n] = gkyl_array_cfetch(app_accel[n], lidx);
      pr_rhs_s[n] = gkyl_array_cfetch(pr_rhs[n], lidx);
      nT_source_s[n] = gkyl_array_cfetch(nT_sources[n], lidx);
    }

    fluid_source_update(mes, tcurr, dt, fluids, app_accels, pr_rhs_s, 
      gkyl_array_fetch(em, lidx),
      gkyl_array_cfetch(app_current, lidx),
      gkyl_array_cfetch(ext_em, lidx),
      nT_source_s
    );
  }
}

void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* up)
{
  gkyl_free(up);
}
