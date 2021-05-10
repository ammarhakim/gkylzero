#include <float.h>

#include <gkyl_update_fsm.h>
#include <gkyl_alloc.h>

struct gkyl_update_fsm*
gkyl_update_fsm_new(int nsteps, const struct gkyl_update_fsm_step *redo)
{
  struct gkyl_update_fsm *fsm;
  fsm = gkyl_malloc(sizeof *fsm);

  fsm->nsteps = nsteps;

  // we need to allocate one extra step to store REDO
  fsm->steps = gkyl_malloc(sizeof(struct gkyl_update_fsm_step[nsteps+1]));
  fsm->steps[GKYL_UPDATE_FSM_REDO].ctx = redo->ctx;
  fsm->steps[GKYL_UPDATE_FSM_REDO].u = redo->u;

  return fsm;
}

struct gkyl_update_status
gkyl_update_fsm_run(struct gkyl_update_fsm *seq, double tcurr, double dt)
{
  return (struct gkyl_update_status) {
    .success = 1,
    .next_state = GKYL_UPDATE_FSM_FINISH,
    .dt_actual = DBL_MAX,
    .dt_suggested = DBL_MAX
  };
}
