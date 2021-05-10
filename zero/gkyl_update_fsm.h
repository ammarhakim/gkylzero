#pragma once

#include <stdbool.h>

/** Status of a step in the update process */
struct gkyl_update_status {
    bool success; // true if update worked, false if a fatal error
    int next_state; // next state of update sequence
    double dt_actual; // actual time-step taken
    double dt_suggested; // suggested stable time-step
};

/** 
 * Basic states of FSM. User code MUST start their state with
 * GKYL_UPDATE_FSM_FIRST for the FSM to function properly.
 */
enum {
  GKYL_UPDATE_FSM_FINISH = -1, // signals completion of update sequence
  GKYL_UPDATE_FSM_REDO = 0, // redo state
  GKYL_UPDATE_FSM_FIRST // user code MUST have this as the first state number
};

struct gkyl_update_fsm_step {
    void *ctx; // closure context to pass to update method
    // function pointer to perform update step
    struct gkyl_update_status (*u)(double tcurr, double dt, void *ctx);
};

struct gkyl_update_fsm {
    int nsteps; // number of steps in sequence
    struct gkyl_update_fsm_step *steps; // steps in sequence
};

/**
 * Allocate a new FSM update sequence object. The returned object must
 * be manually populated by the caller.
 * 
 * @param nsteps Number of stepcs in update sequence
 * @param redo Step to perform a "redo" of sequence 
 */
struct gkyl_update_fsm* gkyl_update_fsm_new(int nsteps,
  struct gkyl_update_fsm_step redo);

/**
 * Given a update sequence, run it till sequence completes, returning
 * final status. Note that the input 'dt' may not be actually what a
 * given system takes. The actual time-step take is returned in the
 * dt_actual field, and a suggested time-step is returned in
 * dt_suggested. Typically, if success field of returned status is
 * 0 then the simulation must abort (after cleanup).
 *
 * @param seq Sequence of steps in update
 * @param init_state Initial state of FSM
 * @param turr Current time
 * @param dt Suggested time-step to take.
 */
struct gkyl_update_status gkyl_update_fsm_run(struct gkyl_update_fsm *seq,
  int init_state, double tcurr, double dt);

/**
 * Free update FSM object
 */
void gkyl_update_fsm_release(struct gkyl_update_fsm *seq);
