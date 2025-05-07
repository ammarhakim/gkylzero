#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>

#include <gkyl_wv_euler_rgfm_priv.h>
#include <gkyl_wv_gr_maxwell_priv.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler_tetrad_priv.h>
#include <gkyl_gr_blackhole.h>

struct gkyl_wave_prop {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_wave_limiter limiter; // limiter to use
  double cfl; // CFL number
  const struct gkyl_wv_eqn *equation; // equation object

  bool force_low_order_flux; // only use Lax flux
  bool check_inv_domain; // flag to indicate if invariant domains are checked

  enum gkyl_wave_split_type split_type; // type of splitting to use

  struct gkyl_wave_geom *geom; // geometry object
  struct gkyl_comm *comm; // communcator
  
  // data for 1D slice update
  struct gkyl_array *waves, *apdq, *amdq, *speeds, *flux2;
  // flags to indicate if fluctuations should be recomputed
  struct gkyl_array *redo_fluct;

  // some stats
  long n_calls; // number of calls to updater
  long n_bad_advance_calls; // number of calls in which positivity had to be fixed
  long n_bad_cells; // number  of cells fixed
  long n_max_bad_cells; // maximum number of cells fixed in a call
};

static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
      theta = 1.0;
      break;
    
    // ** Fully formally-verified implementation of the minmod flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_minmod_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_minmod_tvd.rkt **
    case GKYL_MIN_MOD:
      theta = fmax(0.0, fmin(1.0, r));
      break;

    // ** Partially formally-verified implementation of the superbee flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): NOT PROVEN **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_superbee_tvd.rkt **
    case GKYL_SUPERBEE:
      theta = fmax3(0.0, fmin((2.0 * r), 1.0), fmin(r, 2.0));
      break;

    // ** Partially formally-verified implementation of the van Leer flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_van_leer_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): NOT PROVEN **
    case GKYL_VAN_LEER:
      theta = ((r + fabs(r)) / (1.0 + fabs(r)));
      break;

    // ** Fully formally-verified implementation of the monotonized-centered flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_monotonized_centered_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_monotonized_centered_tvd.rkt **
    case GKYL_MONOTONIZED_CENTERED:
      theta = fmax(0.0, fmin3((2.0 * r), ((1.0 + r) / 2.0), 2.0));
      break;

    case GKYL_BEAM_WARMING:
      theta = r;
      break;

    case GKYL_ZERO:
      theta = 0;
      break;
  }
  return theta;
}

gkyl_wave_prop*
gkyl_wave_prop_new(const struct gkyl_wave_prop_inp *winp)
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(*up));

  up->grid = *(winp->grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = winp->num_up_dirs;
  for (int i=0; i<winp->num_up_dirs; ++i)
    up->update_dirs[i] = winp->update_dirs[i];

  up->limiter = winp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp->limiter;
  up->cfl = winp->cfl;
  up->equation = gkyl_wv_eqn_acquire(winp->equation);

  if (winp->comm)
    up->comm = gkyl_comm_acquire(winp->comm);
  else
    up->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  up->force_low_order_flux = winp->force_low_order_flux;
  up->check_inv_domain = winp->check_inv_domain;

  up->split_type = winp->split_type;

  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  int max_1d = 0;
  for (int d=0; d<ext_range.ndim; ++d) {
    int shape = gkyl_range_shape(&ext_range, d);
    max_1d = max_1d > shape ? max_1d : shape;
  }

  // allocate memory to store 1D slices of waves, speeds and
  // second-order correction flux
  int meqn = winp->equation->num_equations, mwaves = winp->equation->num_waves;
  up->waves = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, max_1d);
  up->apdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->amdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->speeds = gkyl_array_new(GKYL_DOUBLE, mwaves, max_1d);
  up->flux2 = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  up->redo_fluct = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  up->geom = gkyl_wave_geom_acquire(winp->geom);

  up->n_calls = up->n_bad_advance_calls = 0;
  up->n_bad_cells = up->n_max_bad_cells = 0;

  return up;
}

// some helper functions

static inline void
copy_wv_vec(int n, double * GKYL_RESTRICT out, const double * GKYL_RESTRICT inp)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT q, const double * GKYL_RESTRICT amdq_r, const double * GKYL_RESTRICT apdq_l)
{
  for (int i=0; i<meqn; ++i)
    q[i] = q[i] - dtdx*(apdq_l[i] + amdq_r[i]);
}

static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

static inline void
wave_rescale(int meqn, double fact, double *w)
{
  for (int i=0; i<meqn; ++i) w[i] *= fact; 
}

static inline void
calc_second_order_qflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

// this is the sign function for doubles
static inline int sign_double(double val) { return (0.0 < val) - (val < 0.0); }

static inline void
calc_second_order_fflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*sign_double(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

static void
limit_waves(const gkyl_wave_prop *wv, int mwaves,
  const struct gkyl_range *slice_range,
  int lower, int upper, struct gkyl_array *waves, const struct gkyl_array *speed)
{
  int meqn = wv->equation->num_equations;

  for (int mw=0; mw<mwaves; ++mw) {
    const double *wl = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower-1));
    const double *wr = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw*meqn], &wr[mw*meqn]);

    for (int i=lower; i<=upper; ++i) {
      double dotl = dotr;
      
      double * GKYL_RESTRICT wi = gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double * GKYL_RESTRICT wi1 = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i+1));
      
      double wnorm2 = wave_dot_prod(meqn, &wi[mw*meqn], &wi[mw*meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw*meqn], &wi1[mw*meqn]);

      if (wnorm2 > 0) {
        const double *s = gkyl_array_cfetch(speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw*meqn]);
      }
    }
  }
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout)
{
  wv->n_calls += 1;
  
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations;
  //  when forced to use Lax fluxes, we only have a single wave
  int mwaves = wv->force_low_order_flux ? 2 :  wv->equation->num_waves;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // delibrately a double
  
  double ql_local[meqn], qr_local[meqn];
  double fjump_local[meqn];
  double waves_local[meqn*mwaves];
  double amdq_local[meqn], apdq_local[meqn];
  double delta[meqn];

  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  double max_speed = 0.0;

  // state of the update
  enum update_state {
    WV_FIRST_SWEEP, WV_POSITIVITY_SWEEP, WV_FIN_SWEEP
  } state, next_state;

  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];

    double dtdx = dt/wv->grid.dx[dir];

    // upper/lower bounds in direction 'd'. These are edge indices
    int loidx = update_range->lower[dir]-1;
    int upidx = update_range->upper[dir]+2;

    // cell indices in 1D slice for interior cells
    int loidx_c = update_range->lower[dir];
    int upidx_c = update_range->upper[dir];

    struct gkyl_range slice_range;
    gkyl_range_init(&slice_range, 1, (int[]) { loidx }, (int[]) { upidx } );

    struct gkyl_range perp_range;
    gkyl_range_shorten_from_above(&perp_range, update_range, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    // outer loop is over perpendicular directions, inner loop over 1D
    // slice along that direction
    while (gkyl_range_iter_next(&iter)) {
      
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);

      gkyl_array_clear(wv->redo_fluct, 1.0);
      
      enum gkyl_wv_flux_type ftype = wv->force_low_order_flux ?
        GKYL_WV_LOW_ORDER_FLUX : GKYL_WV_HIGH_ORDER_FLUX;

      state = WV_FIRST_SWEEP;

      // perform 1D sweeps, fixing positivity if required
      while (state != WV_FIN_SWEEP) {

        if (state == WV_POSITIVITY_SWEEP)
          ftype = GKYL_WV_LOW_ORDER_FLUX;

        // copy previous time-step solution
        for (int i=loidx_c; i<=upidx_c; ++i) {
          idxl[dir] = i; // cell index
          long lidx = gkyl_range_idx(update_range, idxl);
          copy_wv_vec(meqn, gkyl_array_fetch(qout, lidx), gkyl_array_cfetch(qin, lidx));
        }

        for (int i=loidx; i<=upidx; ++i) {
          idxl[dir] = i-1; idxr[dir] = i;
          long sidx = gkyl_ridx(slice_range, i);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxr);
          double *s = gkyl_array_fetch(wv->speeds, sidx);
          const double *redo_fluct = gkyl_array_cfetch(wv->redo_fluct, sidx);

          if (redo_fluct[0] > 0.0) {

            // compute fluctuations and waves only if needed (this
            // prevents doing the full 1D sweep with low-order fluxes
            // on positivity violations)
            long lidx = gkyl_range_idx(update_range, idxl);
            long ridx = gkyl_range_idx(update_range, idxr);

            const double *qinl = gkyl_array_cfetch(qin, lidx);
            const double *qinr = gkyl_array_cfetch(qin, ridx);

            gkyl_wv_eqn_rotate_to_local(wv->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinl, ql_local);
            gkyl_wv_eqn_rotate_to_local(wv->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinr, qr_local);

            if (wv->split_type == GKYL_WAVE_QWAVE)
              calc_jump(meqn, ql_local, qr_local, delta);
            else
              gkyl_wv_eqn_flux_jump(wv->equation, ql_local, qr_local, delta);

            double my_max_speed = gkyl_wv_eqn_waves(wv->equation, ftype, delta,
              ql_local, qr_local, waves_local, s);
            max_speed = max_speed > my_max_speed ? max_speed : my_max_speed;

            double lenr = cg->lenr[dir];
            for (int mw=0; mw<mwaves; ++mw)
              s[mw] *= lenr; // rescale speeds

            // compute fluctuations in local coordinates
            if (wv->split_type == GKYL_WAVE_QWAVE)
              gkyl_wv_eqn_qfluct(wv->equation, ftype, ql_local, qr_local,
                waves_local, s, amdq_local, apdq_local);
            else
              gkyl_wv_eqn_ffluct(wv->equation, ftype, ql_local, qr_local,
                waves_local, s, amdq_local, apdq_local);
        
            double *waves = gkyl_array_fetch(wv->waves, sidx);
            for (int mw=0; mw<mwaves; ++mw)
              // rotate waves back
              gkyl_wv_eqn_rotate_to_global(wv->equation, 
                cg->tau1[dir], cg->tau2[dir], cg->norm[dir], &waves_local[mw*meqn], &waves[mw*meqn]
              );

            // rotate fluctuations
            double *amdq = gkyl_array_fetch(wv->amdq, sidx);
            gkyl_wv_eqn_rotate_to_global(wv->equation, 
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], amdq_local, amdq);
            
            double *apdq = gkyl_array_fetch(wv->apdq, sidx);
           gkyl_wv_eqn_rotate_to_global(wv->equation, 
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], apdq_local, apdq);
          }
          
          cfla = calc_cfla(mwaves, cfla, dtdx/cg->kappa, s);
        }

        if (cfla > cflm) // check time-step before any updates are performed
          is_cfl_violated = 1.0;

        if (is_cfl_violated > 0) {
          // we need to use this goto to jump out of this deep loop to
          // avoid potential problems with taking too large a
          // time-step. NOTE: This jump is local to a rank. An
          // all-reduce here can't be done as one may end up with a
          // hang due to missing allreduce from some ranks.
          goto outsideloop;
        }

        // compute first-order update in each cell
        for (int i=loidx_c; i<=upidx_c; ++i) { // loop is over cells
        
          idxl[dir] = i; // cell index and left-edge index
          long lidx = gkyl_range_idx(update_range, idxl);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

          calc_first_order_update(meqn, dtdx/cg->kappa,
            gkyl_array_fetch(qout, lidx), 
            gkyl_array_cfetch(wv->amdq, gkyl_ridx(slice_range, i+1)),
            gkyl_array_cfetch(wv->apdq, gkyl_ridx(slice_range, i))
          );
        }

        if (state == WV_FIRST_SWEEP) {
          // we only compute second-correction if we are in first sweep
          
          // apply limiters to waves for all edges in update range,
          // including edges that are on the range boundary
          limit_waves(wv, mwaves, &slice_range,
            update_range->lower[dir], update_range->upper[dir]+1, wv->waves, wv->speeds);

          // get the kappa in the first ghost cell on left (needed in
          // the second order flux calculation)
          idxl[dir] = update_range->lower[dir]-1;
          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
          double kappal = cg->kappa;

          gkyl_array_clear(wv->flux2, 0.0);
          // compute second-order correction fluxes at each interface:
          // note that there is one extra edge than cell
          for (int i=loidx_c; i<=upidx_c+1; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            const double *waves = gkyl_array_cfetch(wv->waves, sidx);
            const double *s = gkyl_array_cfetch(wv->speeds, sidx);
            double *flux2 = gkyl_array_fetch(wv->flux2, sidx);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
            double kappar = cg->kappa;

            if (wv->split_type == GKYL_WAVE_QWAVE) {
              for (int mw=0; mw<mwaves; ++mw)
                calc_second_order_qflux(meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);
            }
            else {
              for (int mw=0; mw<mwaves; ++mw)
                calc_second_order_fflux(meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);
            }

            kappal = kappar;
          }

          // add second correction flux to solution in each interior cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

            calc_second_order_update(meqn, dtdx/cg->kappa,
              gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i+1))
            );
          }
        }

        next_state = WV_FIN_SWEEP;
        // check invariant domains if needed
        if ( (state == WV_FIRST_SWEEP) && wv->check_inv_domain) {
          long n_bad_cells = 0;            

          gkyl_array_clear(wv->redo_fluct, 0.0); // by default no edge needs recomputing
          
          // check if invariant domains are violated, flagging edges
          // of each bad cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            idxl[dir] = i;
            const double *qt = gkyl_array_cfetch(qout, gkyl_range_idx(update_range, idxl));
            if (!gkyl_wv_eqn_check_inv(wv->equation, qt)) {

              double *redo_fluct_l = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i));
              double *redo_fluct_r = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i+1));
              // mark left and right edges so fluctuations are redone
              redo_fluct_l[0] = 1.0;
              redo_fluct_r[0] = 1.0;

              n_bad_cells += 1;
            }
          }

          if (n_bad_cells > 0) {
            // we need to resweep the 1D slice again
            next_state = WV_POSITIVITY_SWEEP;
            wv->n_bad_advance_calls += 1;
          }

          wv->n_bad_cells += n_bad_cells;
          wv->n_max_bad_cells = wv->n_max_bad_cells >  n_bad_cells ? wv->n_max_bad_cells : n_bad_cells;
        }

        if (wv->equation->type == GKYL_EQN_EULER_RGFM) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
          int num_species = euler_rgfm->num_species;
          int reinit_freq = euler_rgfm->reinit_freq;
      
          for (int i = loidx_c; i <= upidx_c; i++) {
            idxl[dir] = i;

            double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));

            double reinit_param = qnew[4 + (2 * num_species)];
      
            if (reinit_param > reinit_freq) {
              double rho_total = qnew[0];
       
              bool update_up = false;
              bool update_down = false;
              for (int j = 0; j < num_species - 1; j++) {
                if (qnew[5 + j] / rho_total >= 0.5) {
                  idxl[dir] = i - 1;
                  double *ql = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i - 2;
                  double *qll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i - 3;
                  double *qlll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  
                  double rho_total_l = ql[0];
                  double rho_total_ll = qll[0];
                  double rho_total_lll = qlll[0];

                  idxl[dir] = i + 1;
                  double *qr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i + 2;
                  double *qrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i + 3;
                  double *qrrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  
                  double rho_total_r = qr[0];
                  double rho_total_rr = qrr[0];
                  double rho_total_rrr = qrrr[0];
                  
                  if (ql[5 + j] / rho_total_l < 0.5 || qll[5 + j] / rho_total_ll < 0.5 || qlll[5 + j] / rho_total_lll < 0.5 || qr[5 + j] / rho_total_r < 0.5 ||
                      qrr[5 + j] / rho_total_rr < 0.5 || qrrr[5 + j] / rho_total_rrr < 0.5)  {
                    qnew[5 + j] = 0.99999 * rho_total;
                    qnew[4 + num_species + j] = 0.99999 * rho_total;
                    update_up = true;
                  }
                }
                
                if (qnew[5 + j] / rho_total < 0.5) {
                  idxl[dir] = i - 1;
                  double *ql = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i - 2;
                  double *qll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i - 3;
                  double *qlll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  
                  double rho_total_l = ql[0];
                  double rho_total_ll = qll[0];
                  double rho_total_lll = qlll[0];

                  idxl[dir] = i + 1;
                  double *qr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i + 2;
                  double *qrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  idxl[dir] = i + 3;
                  double *qrrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
                  
                  double rho_total_r = qr[0];
                  double rho_total_rr = qrr[0];
                  double rho_total_rrr = qrrr[0];

                  if (qr[5 + j] / rho_total_r >= 0.5 || qrr[5 + j] / rho_total_rr >= 0.5 || qrrr[5 + j] / rho_total_rrr >= 0.5 || ql[5 + j] / rho_total_l >= 0.5 ||
                      qll[5 + j] / rho_total_ll >= 0.5 || qlll[5 + j] / rho_total_lll >= 0.5) {
                    qnew[5 + j] = 0.00001 * rho_total;
                    qnew[4 + num_species + j] = 0.00001 * rho_total;
                    update_down = true;
                  }
                }
              }

              if (update_up) {
                qnew[3 + (2 * num_species)] = 0.00001 * rho_total;
              }
              if (update_down) {
                qnew[3 + (2 * num_species)] = 0.99999 * rho_total;
              }

              qnew[4 + (2 * num_species)] = 0.0;
            }
            else {
              qnew[4 + (2 * num_species)] += 1.0;
            }
          }
        }

        if (wv->equation->type == GKYL_EQN_GR_MAXWELL) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
          
          const enum gkyl_spacetime_gauge spacetime_gauge = gr_maxwell->spacetime_gauge;

          if (spacetime_gauge == GKYL_STATIC_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_maxwell->spacetime;
            int reinit_freq = gr_maxwell->reinit_freq;
            
            for (int i = loidx_c; i<= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[22];

              if (evol_param > reinit_freq) {
                double x = qnew[23];
                double y = qnew[24];
                double z = qnew[25];

                double lapse;
                double *shift = gkyl_malloc(sizeof(double[3]));
                bool in_excision_region;

                double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
                }

                spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
                spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
                spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

                spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);

                qnew[8] = lapse;
                qnew[9] = shift[0]; qnew[10] = shift[1]; qnew[11] = shift[2];

                qnew[12] = spatial_metric[0][0]; qnew[13] = spatial_metric[0][1]; qnew[14] = spatial_metric[0][2];
                qnew[15] = spatial_metric[1][0]; qnew[16] = spatial_metric[1][1]; qnew[17] = spatial_metric[1][2];
                qnew[18] = spatial_metric[2][0]; qnew[19] = spatial_metric[2][1]; qnew[20] = spatial_metric[2][2];

                if (in_excision_region) {
                  for (int i = 0; i < 22; i++) {
                    qnew[i] = 0.0;
                  }

                  qnew[21] = -1.0;
                }
                else {
                  qnew[21] = 1.0;
                }

                for (int i = 0; i < 3; i++) {
                  gkyl_free(spatial_metric[i]);
                }
                gkyl_free(spatial_metric);
                gkyl_free(shift);
                
                qnew[22] = 0.0;
              }
              else {
                qnew[22] += 1.0;
              }
            }
          }
        }

        if (wv->equation->type == GKYL_EQN_GR_EULER) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
          
          const enum gkyl_spacetime_gauge spacetime_gauge = gr_euler->spacetime_gauge;

          if (spacetime_gauge == GKYL_STATIC_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_euler->spacetime;
            int reinit_freq = gr_euler->reinit_freq;
            
            for (int i = loidx_c; i<= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[67];

              if (evol_param > reinit_freq) {
                double x = qnew[68];
                double y = qnew[69];
                double z = qnew[70];

                double lapse;
                double *shift = gkyl_malloc(sizeof(double[3]));
                bool in_excision_region;

                double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
                }

                double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
                }

                double *lapse_der = gkyl_malloc(sizeof(double[3]));
                double **shift_der = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  shift_der[i] = gkyl_malloc(sizeof(double[3]));
                }

                double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

                  for (int j = 0; j < 3; j++) {
                    spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
                  }
                }

                spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
                spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
                spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

                spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
                spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

                spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
                spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
                spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

                qnew[5] = lapse;
                qnew[6] = shift[0]; qnew[7] = shift[1]; qnew[8] = shift[2];

                qnew[9] = spatial_metric[0][0]; qnew[10] = spatial_metric[0][1]; qnew[11] = spatial_metric[0][2];
                qnew[12] = spatial_metric[1][0]; qnew[13] = spatial_metric[1][1]; qnew[14] = spatial_metric[1][2];
                qnew[15] = spatial_metric[2][0]; qnew[16] = spatial_metric[2][1]; qnew[17] = spatial_metric[2][2];

                qnew[18] = extrinsic_curvature[0][0]; qnew[19] = extrinsic_curvature[0][1]; qnew[20] = extrinsic_curvature[0][2];
                qnew[21] = extrinsic_curvature[1][0]; qnew[22] = extrinsic_curvature[1][1]; qnew[23] = extrinsic_curvature[1][2];
                qnew[24] = extrinsic_curvature[2][0]; qnew[25] = extrinsic_curvature[2][1]; qnew[26] = extrinsic_curvature[2][2];

                if (in_excision_region) {
                  qnew[27] = -1.0;
                }
                else {
                  qnew[27] = 1.0;
                }

                qnew[28] = lapse_der[0]; qnew[29] = lapse_der[1]; qnew[30] = lapse_der[2];

                qnew[31] = shift_der[0][0]; qnew[32] = shift_der[0][1]; qnew[33] = shift_der[0][2];
                qnew[34] = shift_der[1][0]; qnew[35] = shift_der[1][1]; qnew[36] = shift_der[1][2];
                qnew[37] = shift_der[2][0]; qnew[38] = shift_der[2][1]; qnew[39] = shift_der[2][2];

                qnew[40] = spatial_metric_der[0][0][0]; qnew[41] = spatial_metric_der[0][0][1]; qnew[42] = spatial_metric_der[0][0][2];
                qnew[43] = spatial_metric_der[0][1][0]; qnew[44] = spatial_metric_der[0][1][1]; qnew[45] = spatial_metric_der[0][1][2];
                qnew[46] = spatial_metric_der[0][2][0]; qnew[47] = spatial_metric_der[0][2][1]; qnew[48] = spatial_metric_der[0][2][2];

                qnew[49] = spatial_metric_der[1][0][0]; qnew[50] = spatial_metric_der[1][0][1]; qnew[51] = spatial_metric_der[1][0][2];
                qnew[52] = spatial_metric_der[1][1][0]; qnew[53] = spatial_metric_der[1][1][1]; qnew[54] = spatial_metric_der[1][1][2];
                qnew[55] = spatial_metric_der[1][2][0]; qnew[56] = spatial_metric_der[1][2][1]; qnew[57] = spatial_metric_der[1][2][2];

                qnew[58] = spatial_metric_der[2][0][0]; qnew[59] = spatial_metric_der[2][0][1]; qnew[60] = spatial_metric_der[2][0][2];
                qnew[61] = spatial_metric_der[2][1][0]; qnew[62] = spatial_metric_der[2][1][1]; qnew[63] = spatial_metric_der[2][1][2];
                qnew[64] = spatial_metric_der[2][2][0]; qnew[65] = spatial_metric_der[2][2][1]; qnew[66] = spatial_metric_der[2][2][2];

                if (in_excision_region) {
                  for (int i = 0; i < 67; i++) {
                    qnew[i] = 0.0;
                  }
                  qnew[27] = -1.0;
                }

                for (int i = 0; i < 3; i++) {
                  gkyl_free(spatial_metric[i]);
                  gkyl_free(extrinsic_curvature[i]);
                  gkyl_free(shift_der[i]);
              
                  for (int j = 0; j < 3; j++) {
                    gkyl_free(spatial_metric_der[i][j]);
                  }
                  gkyl_free(spatial_metric_der[i]);
                }
                gkyl_free(spatial_metric);
                gkyl_free(extrinsic_curvature);
                gkyl_free(shift);
                gkyl_free(lapse_der);
                gkyl_free(shift_der);
                gkyl_free(spatial_metric_der);
                
                qnew[67] = 0.0;
              }
              else {
                qnew[67] += 1.0;
              }
            }
          }
        }

        if (wv->equation->type == GKYL_EQN_GR_EULER_TETRAD) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(eqn, struct wv_gr_euler_tetrad, eqn);
          
          const enum gkyl_spacetime_gauge spacetime_gauge = gr_euler_tetrad->spacetime_gauge;

          if (spacetime_gauge == GKYL_STATIC_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_euler_tetrad->spacetime;
            int reinit_freq = gr_euler_tetrad->reinit_freq;
            
            for (int i = loidx_c; i<= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[67];

              if (evol_param > reinit_freq) {
                double x = qnew[68];
                double y = qnew[69];
                double z = qnew[70];

                double lapse;
                double *shift = gkyl_malloc(sizeof(double[3]));
                bool in_excision_region;

                double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
                }

                double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
                }

                double *lapse_der = gkyl_malloc(sizeof(double[3]));
                double **shift_der = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  shift_der[i] = gkyl_malloc(sizeof(double[3]));
                }

                double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

                  for (int j = 0; j < 3; j++) {
                    spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
                  }
                }

                spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
                spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
                spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

                spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
                spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

                spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
                spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
                spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

                qnew[5] = lapse;
                qnew[6] = shift[0]; qnew[7] = shift[1]; qnew[8] = shift[2];

                qnew[9] = spatial_metric[0][0]; qnew[10] = spatial_metric[0][1]; qnew[11] = spatial_metric[0][2];
                qnew[12] = spatial_metric[1][0]; qnew[13] = spatial_metric[1][1]; qnew[14] = spatial_metric[1][2];
                qnew[15] = spatial_metric[2][0]; qnew[16] = spatial_metric[2][1]; qnew[17] = spatial_metric[2][2];

                qnew[18] = extrinsic_curvature[0][0]; qnew[19] = extrinsic_curvature[0][1]; qnew[20] = extrinsic_curvature[0][2];
                qnew[21] = extrinsic_curvature[1][0]; qnew[22] = extrinsic_curvature[1][1]; qnew[23] = extrinsic_curvature[1][2];
                qnew[24] = extrinsic_curvature[2][0]; qnew[25] = extrinsic_curvature[2][1]; qnew[26] = extrinsic_curvature[2][2];

                if (in_excision_region) {
                  qnew[27] = -1.0;
                }
                else {
                  qnew[27] = 1.0;
                }

                qnew[28] = lapse_der[0]; qnew[29] = lapse_der[1]; qnew[30] = lapse_der[2];

                qnew[31] = shift_der[0][0]; qnew[32] = shift_der[0][1]; qnew[33] = shift_der[0][2];
                qnew[34] = shift_der[1][0]; qnew[35] = shift_der[1][1]; qnew[36] = shift_der[1][2];
                qnew[37] = shift_der[2][0]; qnew[38] = shift_der[2][1]; qnew[39] = shift_der[2][2];

                qnew[40] = spatial_metric_der[0][0][0]; qnew[41] = spatial_metric_der[0][0][1]; qnew[42] = spatial_metric_der[0][0][2];
                qnew[43] = spatial_metric_der[0][1][0]; qnew[44] = spatial_metric_der[0][1][1]; qnew[45] = spatial_metric_der[0][1][2];
                qnew[46] = spatial_metric_der[0][2][0]; qnew[47] = spatial_metric_der[0][2][1]; qnew[48] = spatial_metric_der[0][2][2];

                qnew[49] = spatial_metric_der[1][0][0]; qnew[50] = spatial_metric_der[1][0][1]; qnew[51] = spatial_metric_der[1][0][2];
                qnew[52] = spatial_metric_der[1][1][0]; qnew[53] = spatial_metric_der[1][1][1]; qnew[54] = spatial_metric_der[1][1][2];
                qnew[55] = spatial_metric_der[1][2][0]; qnew[56] = spatial_metric_der[1][2][1]; qnew[57] = spatial_metric_der[1][2][2];

                qnew[58] = spatial_metric_der[2][0][0]; qnew[59] = spatial_metric_der[2][0][1]; qnew[60] = spatial_metric_der[2][0][2];
                qnew[61] = spatial_metric_der[2][1][0]; qnew[62] = spatial_metric_der[2][1][1]; qnew[63] = spatial_metric_der[2][1][2];
                qnew[64] = spatial_metric_der[2][2][0]; qnew[65] = spatial_metric_der[2][2][1]; qnew[66] = spatial_metric_der[2][2][2];

                if (in_excision_region) {
                  for (int i = 0; i < 67; i++) {
                    qnew[i] = 0.0;
                  }
                  qnew[27] = -1.0;
                }

                for (int i = 0; i < 3; i++) {
                  gkyl_free(spatial_metric[i]);
                  gkyl_free(extrinsic_curvature[i]);
                  gkyl_free(shift_der[i]);
              
                  for (int j = 0; j < 3; j++) {
                    gkyl_free(spatial_metric_der[i][j]);
                  }
                  gkyl_free(spatial_metric_der[i]);
                }
                gkyl_free(spatial_metric);
                gkyl_free(extrinsic_curvature);
                gkyl_free(shift);
                gkyl_free(lapse_der);
                gkyl_free(shift_der);
                gkyl_free(spatial_metric_der);
                
                qnew[67] = 0.0;
              }
              else {
                qnew[67] += 1.0;
              }
            }
          }
        }

        if (wv->equation->type == GKYL_EQN_GR_ULTRA_REL_EULER) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);

          const enum gkyl_spacetime_gauge spacetime_gauge = gr_ultra_rel_euler->spacetime_gauge;

          if (spacetime_gauge == GKYL_STATIC_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler->spacetime;
            int reinit_freq = gr_ultra_rel_euler->reinit_freq;
            
            for (int i = loidx_c; i<= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[66];

              if (evol_param > reinit_freq) {
                double x = qnew[67];
                double y = qnew[68];
                double z = qnew[69];

                double lapse;
                double *shift = gkyl_malloc(sizeof(double[3]));
                bool in_excision_region;

                double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
                }

                double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
                }

                double *lapse_der = gkyl_malloc(sizeof(double[3]));
                double **shift_der = gkyl_malloc(sizeof(double*[3]));
                for (int i = 0; i < 3; i++) {
                  shift_der[i] = gkyl_malloc(sizeof(double[3]));
                }

                double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
                for (int i = 0; i < 3; i++) {
                  spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

                  for (int j = 0; j < 3; j++) {
                    spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
                  }
                }

                spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
                spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
                spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

                spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
                spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

                spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
                spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
                spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

                qnew[4] = lapse;
                qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

                qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
                qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
                qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

                qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
                qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
                qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

                if (in_excision_region) {
                  qnew[26] = -1.0;
                }
                else {
                  qnew[26] = 1.0;
                }

                qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

                qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
                qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
                qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

                qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
                qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
                qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

                qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
                qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
                qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

                qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
                qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
                qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

                if (in_excision_region) {
                  for (int i = 0; i < 66; i++) {
                    qnew[i] = 0.0;
                  }
                  qnew[26] = -1.0;
                }

                for (int i = 0; i < 3; i++) {
                  gkyl_free(spatial_metric[i]);
                  gkyl_free(extrinsic_curvature[i]);
                  gkyl_free(shift_der[i]);
              
                  for (int j = 0; j < 3; j++) {
                    gkyl_free(spatial_metric_der[i][j]);
                  }
                  gkyl_free(spatial_metric_der[i]);
                }
                gkyl_free(spatial_metric);
                gkyl_free(extrinsic_curvature);
                gkyl_free(shift);
                gkyl_free(lapse_der);
                gkyl_free(shift_der);
                gkyl_free(spatial_metric_der);

                qnew[66] = 0.0;
              }
              else {
                qnew[66] += 1.0;
              }
            }
          }
          else if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler->spacetime;
            const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

            double mass = blackhole->mass;
            double spin = blackhole->spin;

            double pos_x = blackhole->pos_x;
            double pos_y = blackhole->pos_y;
            double pos_z = blackhole->pos_z;

            for (int i = loidx_c; i <= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[66];

              if (evol_param < mass) {
                evol_param += 0.001 + (0.001 * evol_param);
              }
              qnew[66] = evol_param;

              double x = qnew[67];
              double y = qnew[68];
              double z = qnew[69];

              struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

              double lapse;
              double *shift = gkyl_malloc(sizeof(double[3]));
              bool in_excision_region;

              double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
              }

              double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
              }

              double *lapse_der = gkyl_malloc(sizeof(double[3]));
              double **shift_der = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                shift_der[i] = gkyl_malloc(sizeof(double[3]));
              }

              double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
              for (int i = 0; i < 3; i++) {
                spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

                for (int j = 0; j < 3; j++) {
                  spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
                }
              }

              new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
              new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
              new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

              new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
              new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

              new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
              new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
              new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

              qnew[4] = lapse;
              qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

              qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
              qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
              qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

              qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
              qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
              qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

              if (in_excision_region) {
                qnew[26] = -1.0;
              }
              else {
                qnew[26] = 1.0;
              }

              qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

              qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
              qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
              qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

              qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
              qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
              qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

              qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
              qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
              qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

              qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
              qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
              qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

              if (in_excision_region) {
                for (int i = 0; i < 66; i++) {
                  qnew[i] = 0.0;
                }
                qnew[26] = -1.0;
              }

              for (int i = 0; i < 3; i++) {
                gkyl_free(spatial_metric[i]);
                gkyl_free(extrinsic_curvature[i]);
                gkyl_free(shift_der[i]);
            
                for (int j = 0; j < 3; j++) {
                  gkyl_free(spatial_metric_der[i][j]);
                }
                gkyl_free(spatial_metric_der[i]);
              }
              gkyl_free(spatial_metric);
              gkyl_free(extrinsic_curvature);
              gkyl_free(shift);
              gkyl_free(lapse_der);
              gkyl_free(shift_der);
              gkyl_free(spatial_metric_der);
              gkyl_free(new_spacetime);
            }
          }
        }

        if (wv->equation->type == GKYL_EQN_GR_ULTRA_REL_EULER_TETRAD) {
          const struct gkyl_wv_eqn* eqn = wv->equation;
          const struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(eqn, struct wv_gr_ultra_rel_euler_tetrad, eqn);

          const enum gkyl_spacetime_gauge spacetime_gauge = gr_ultra_rel_euler_tetrad->spacetime_gauge;

          if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
            const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler_tetrad->spacetime;
            const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

            double mass = blackhole->mass;
            double spin = blackhole->spin;

            double pos_x = blackhole->pos_x;
            double pos_y = blackhole->pos_y;
            double pos_z = blackhole->pos_z;

            for (int i = loidx_c; i <= upidx_c; i++) {
              idxl[dir] = i;

              double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
              double evol_param = qnew[66];

              if (evol_param < mass) {
                evol_param += 0.001 + (0.001 * evol_param);
              }
              qnew[66] = evol_param;

              double x = qnew[67];
              double y = qnew[68];
              double z = qnew[69];

              struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

              double lapse;
              double *shift = gkyl_malloc(sizeof(double[3]));
              bool in_excision_region;

              double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
              }

              double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
              }

              double *lapse_der = gkyl_malloc(sizeof(double[3]));
              double **shift_der = gkyl_malloc(sizeof(double*[3]));
              for (int i = 0; i < 3; i++) {
                shift_der[i] = gkyl_malloc(sizeof(double[3]));
              }

              double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
              for (int i = 0; i < 3; i++) {
                spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

                for (int j = 0; j < 3; j++) {
                  spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
                }
              }

              new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
              new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
              new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

              new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
              new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

              new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
              new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
              new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

              qnew[4] = lapse;
              qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

              qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
              qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
              qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

              qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
              qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
              qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

              if (in_excision_region) {
                qnew[26] = -1.0;
              }
              else {
                qnew[26] = 1.0;
              }

              qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

              qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
              qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
              qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

              qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
              qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
              qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

              qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
              qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
              qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

              qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
              qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
              qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

              if (in_excision_region) {
                for (int i = 0; i < 66; i++) {
                  qnew[i] = 0.0;
                }
                qnew[26] = -1.0;
              }

              for (int i = 0; i < 3; i++) {
                gkyl_free(spatial_metric[i]);
                gkyl_free(extrinsic_curvature[i]);
                gkyl_free(shift_der[i]);
            
                for (int j = 0; j < 3; j++) {
                  gkyl_free(spatial_metric_der[i][j]);
                }
                gkyl_free(spatial_metric_der[i]);
              }
              gkyl_free(spatial_metric);
              gkyl_free(extrinsic_curvature);
              gkyl_free(shift);
              gkyl_free(lapse_der);
              gkyl_free(shift_der);
              gkyl_free(spatial_metric_der);
              gkyl_free(new_spacetime);
            }
          }
        }

        state = next_state; // change state for next sweep
        
      } // end loop over sweeps
    } // end loop over perpendicular directions
  } // end loop over directions

  outsideloop:
  ;

  // compute actual CFL, status & max-speed across all domains
  double red_vars[3] = { cfla, is_cfl_violated, max_speed };
  double red_vars_global[3] = { 0.0, 0.0, 0.0 };
  gkyl_comm_allreduce(wv->comm, GKYL_DOUBLE, GKYL_MAX, 3, &red_vars, &red_vars_global);

  cfla = red_vars_global[0];
  is_cfl_violated = red_vars_global[1];
  max_speed = red_vars_global[2];

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  if (is_cfl_violated > 0.0)
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_wave_prop_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
      .max_speed = max_speed,
    };
  
  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).

  return (struct gkyl_wave_prop_status) {
    .success = is_cfl_violated > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
    .max_speed = max_speed,
  };
}

double
gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];
      double dx = wv->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(wv->equation, q);
      max_dt = fmin(max_dt, wv->cfl*dx/maxs);
    }
    
  }

  return max_dt;
}

struct gkyl_wave_prop_stats
gkyl_wave_prop_stats(const gkyl_wave_prop *wv)
{
  return (struct gkyl_wave_prop_stats) {
    .n_calls = wv->n_calls,
    .n_bad_advance_calls = wv->n_bad_advance_calls,
    .n_bad_cells = wv->n_bad_cells,
    .n_max_bad_cells = wv->n_max_bad_cells
  };
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_array_release(up->waves);
  gkyl_array_release(up->apdq);
  gkyl_array_release(up->amdq);
  gkyl_array_release(up->speeds);
  gkyl_array_release(up->flux2);
  gkyl_array_release(up->redo_fluct);
  gkyl_comm_release(up->comm);
  
  gkyl_wave_geom_release(up->geom);
  
  gkyl_free(up);
}
