#include <float.h>
#include <math.h>

#include <gkyl_update_fsm.h>
#include <gkyl_alloc.h>

struct gkyl_update_fsm*
gkyl_update_fsm_new(int nsteps, struct gkyl_update_fsm_step steps[])
{
  struct gkyl_update_fsm *fsm;
  fsm = gkyl_malloc(sizeof *fsm);

  fsm->nsteps = nsteps;

  fsm->steps = gkyl_malloc(sizeof(struct gkyl_update_fsm_step[nsteps]));
  for (int i=0; i<nsteps; ++i)
    fsm->steps[i] = steps[i];

  return fsm;
}

struct gkyl_update_status
gkyl_update_fsm_run(struct gkyl_update_fsm *seq, double tcurr, double dt)
{
  int step = 0;
  bool success = true;
  double dt_suggested = DBL_MAX, dt_actual = dt;

  int state = GKYL_UPDATE_FSM_FIRST;
  while ((state != GKYL_UPDATE_FSM_FINISH) && (state != GKYL_UPDATE_FSM_ABORT)) {
    // run current step
    struct gkyl_update_status status =
      seq->steps[state].u(tcurr, dt_actual, seq->steps[state].ctx);

    if (step > 0 && status.dt_actual < dt_actual) {
      // if actual time-step taken after first one is smaller than
      // dt_suggested, we need to redo whole sequence
      step = 0;
      state = GKYL_UPDATE_FSM_REDO;
      dt_actual = status.dt_actual;
      dt_suggested = DBL_MAX;
    }
    else {
      step += 1;
      state = status.next_state;
      dt_actual = status.dt_actual;
      dt_suggested = fmin(dt_suggested, status.dt_suggested);
    }
  }
  
  return (struct gkyl_update_status) {
    .success = (state == GKYL_UPDATE_FSM_FINISH ? true : false),
    .dt_actual = dt_actual,
    .dt_suggested = dt_suggested,
  };
}

void
gkyl_update_fsm_release(struct gkyl_update_fsm *seq)
{
  gkyl_free(seq->steps);
  gkyl_free(seq);
}
