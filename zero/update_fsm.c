#include <float.h>
#include <math.h>

#include <gkyl_update_fsm.h>
#include <gkyl_alloc.h>

struct gkyl_update_fsm*
gkyl_update_fsm_with_transitions_new(int nactions, struct gkyl_update_fsm_action actions[],
  int ntransitions, struct gkyl_update_fsm_transition *transitions,
  struct gkyl_update_fsm_action redo)
{
  struct gkyl_update_fsm *fsm;
  fsm = gkyl_malloc(sizeof *fsm);

  fsm->nactions = nactions;
  fsm->actions = gkyl_malloc(sizeof(struct gkyl_update_fsm_action[nactions]));
  for (int i=0; i<nactions; ++i)
    fsm->actions[i] = actions[i];

  fsm->ntransitions = ntransitions;
  fsm->transitions = gkyl_malloc(sizeof(struct gkyl_update_fsm_transition[ntransitions]));
  for (int i=0; i<ntransitions; ++i)
    fsm->transitions[i] = transitions[i];

  fsm->redo = redo;

  return fsm;
}

struct gkyl_update_fsm*
gkyl_update_fsm_new(int nactions, struct gkyl_update_fsm_action actions[],
  struct gkyl_update_fsm_action redo)
{
  struct gkyl_update_fsm *fsm;
  fsm = gkyl_malloc(sizeof *fsm);

  fsm->nactions = nactions;
  fsm->actions = gkyl_malloc(sizeof(struct gkyl_update_fsm_action[nactions]));
  for (int i=0; i<nactions; ++i)
    fsm->actions[i] = actions[i];

  // construct default transition list
  int nt = fsm->ntransitions = nactions;
  fsm->transitions = gkyl_malloc(sizeof(struct gkyl_update_fsm_transition[nt]));

  // following transition table ensures that each action is run in order specified
  for (int i=0; i<nt-1; ++i)
    fsm->transitions[i] = (struct gkyl_update_fsm_transition) {
      i+1, GKYL_UPDATE_FSM_REDO
    };
  fsm->transitions[nt-1] = (struct gkyl_update_fsm_transition) {
    GKYL_UPDATE_FSM_FINISH,
    GKYL_UPDATE_FSM_REDO
  };

  fsm->redo = redo;

  return fsm;
}

struct gkyl_update_fsm_status
gkyl_update_fsm_run(struct gkyl_update_fsm *seq, double tcurr, double dt)
{
  int step = 0;
  bool success = true;
  double dt_suggested = DBL_MAX, dt_actual = dt;

  int state = GKYL_UPDATE_FSM_FIRST;
  while (state != GKYL_UPDATE_FSM_FINISH) {
    // run current step
    struct gkyl_update_fsm_status status =
      seq->actions[state].u(tcurr, dt_actual, seq->actions[state].ctx);

    if (status.status == GKYL_UPDATE_FSM_STATUS_FAIL)
      break;

    state = seq->transitions[state].success_action;
    step += 1;
    dt_suggested = fmin(dt_suggested, status.dt_suggested);

    if ( (status.status == GKYL_UPDATE_FSM_STATUS_REDO) || (step > 0 && status.dt_actual < dt_actual) ) {
      // redo the whole sequence if asked to do so, or if a later step
      // takes a smaller time-step than it was asked to take
      seq->redo.u(tcurr, dt, seq->redo.ctx);
      
      step = 0;
      state = GKYL_UPDATE_FSM_FIRST;
      dt_actual = status.dt_actual;
      dt_suggested = DBL_MAX;      
    }
  }
  
  return (struct gkyl_update_fsm_status) {
    .status = (state == GKYL_UPDATE_FSM_FINISH ? GKYL_UPDATE_FSM_STATUS_SUCCESS : GKYL_UPDATE_FSM_STATUS_FAIL),
    .dt_actual = dt_actual,
    .dt_suggested = dt_suggested,
  };
}

void
gkyl_update_fsm_release(struct gkyl_update_fsm *seq)
{
  gkyl_free(seq->actions);
  gkyl_free(seq->transitions);
  gkyl_free(seq);
}
