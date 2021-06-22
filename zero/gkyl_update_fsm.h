#pragma once

enum gkyl_update_fsm_code {
  GKYL_UPDATE_FSM_STATUS_SUCCESS, // success  
  GKYL_UPDATE_FSM_STATUS_REDO, // redo 
  GKYL_UPDATE_FSM_STATUS_FAIL // catastrophic failure
};

/** Status of a step in the update process */
struct gkyl_update_fsm_status {
  enum gkyl_update_fsm_code status; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

/** 
 * Basic states of FSM. User code MUST start their state with
 * GKYL_UPDATE_FSM_FIRST for the FSM to function properly.
 */
enum gkyl_update_fsm_state {
  GKYL_UPDATE_FSM_FINISH = -1, // signals completion of update
  GKYL_UPDATE_FSM_FIRST = 0 // user code MUST have this as first state
};

// Represents an action in the FSM
struct gkyl_update_fsm_action {
  void *ctx; // closure context to pass to update method
  // function pointer to perform update step
  struct gkyl_update_fsm_status (*u)(double tcurr, double dt, void *ctx);
};

// Represent a transition of the FSM
struct gkyl_update_fsm_transition {
  int success_action; // action on success
};

// Finite-state machine
struct gkyl_update_fsm {
  int nactions; // number of actions
  struct gkyl_update_fsm_action *actions; // actions
  int ntransitions; // number of transitions
  struct gkyl_update_fsm_transition *transitions; // transitions
  struct gkyl_update_fsm_action redo; // redo action
};

/**
 * Allocate a new FSM update sequence object.
 * 
 * @param nactions Number of actions
 * @param actions List of actions
 * @param ntransitions Number of transitions
 * @param transitions List of transitions
 * @param redo Action on transition to redo state
 * @return New FSM update object
 */
struct gkyl_update_fsm* gkyl_update_fsm_new(int nactions, struct gkyl_update_fsm_action actions[],
  int ntransitions, struct gkyl_update_fsm_transition *transitions,
  struct gkyl_update_fsm_action redo);

/**
 * Given a update sequence, run it till sequence completes, returning
 * final status. Note that the input @a dt may not be actually what 
 * the system takes. The actual time-step take is returned in the
 * dt_actual field, and a suggested time-step is returned in
 * dt_suggested. Typically, if success field of returned status is false
 * then the simulation must abort (after cleanup).
 *
 * @param seq Sequence object
 * @param turr Current time
 * @param dt Suggested time-step to take.
 */
struct gkyl_update_fsm_status gkyl_update_fsm_run(struct gkyl_update_fsm *seq, double tcurr, double dt);

/**
 * Free update FSM object
 * 
 * @param seq Sequence object
 */
void gkyl_update_fsm_release(struct gkyl_update_fsm *seq);
