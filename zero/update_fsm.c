#include <float.h>

#include <gkyl_update_fsm.h>
#include <gkyl_alloc.h>

static struct gkyl_update_status
update_step_dummy(double tcurr, double dt, void *ctx)
{
  return (struct gkyl_update_status) {
    .success = 1,
    .next_state = GKYL_UPDATE_FSM_FINISH,
    .dt_actual = DBL_MAX,
    .dt_suggested = DBL_MAX
  };  
}

struct gkyl_update_fsm*
gkyl_update_fsm_new(int nsteps)
{
  struct gkyl_update_fsm *fsm;
  fsm = gkyl_malloc(sizeof *fsm);

  fsm->nsteps = nsteps;
  fsm->steps = gkyl_malloc(sizeof(struct gkyl_update_fsm_step[nsteps]));

  // set all steps to default to dummy steps
  for (int i=0; i<nsteps; ++i)
    fsm->steps[i] = (struct gkyl_update_fsm_step) { .u = update_step_dummy };

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
