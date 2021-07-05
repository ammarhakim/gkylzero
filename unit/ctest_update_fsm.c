#include <float.h>

#include <acutest.h>
#include <gkyl_update_fsm.h>

enum seq_states {
  SRC_1 = GKYL_UPDATE_FSM_FIRST, // first step
  FLUID,
  SRC_2
};

struct seq_ctx {
  int nredo, nfluid, nsrc;
};

static void
seq_ctx_print(struct seq_ctx ctx)
{
  fprintf(stdout, "nredo: %d. nfluid: %d. nsrc: %d\n", ctx.nredo, ctx.nfluid, ctx.nsrc);
}

struct gkyl_update_fsm_status
seq_redo(double tcurr, double dt, void *ctx)
{
  //printf("seq_redo: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nredo += 1;
  
  return (struct gkyl_update_fsm_status) {
    .status = GKYL_UPDATE_FSM_STATUS_SUCCESS,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX,
  };
}

struct gkyl_update_fsm_status
seq_src_1(double tcurr, double dt, void *ctx)
{
  //printf("seq_src_1: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nsrc += 1;

  enum gkyl_update_fsm_code status = GKYL_UPDATE_FSM_STATUS_SUCCESS;
  if (dt > 100)
    status = GKYL_UPDATE_FSM_STATUS_FAIL;  

  return (struct gkyl_update_fsm_status) {
    .status = status,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX,
  };
}

struct gkyl_update_fsm_status
seq_fluid(double tcurr, double dt, void *ctx)
{
  //printf("seq_fluid: dt = %g\n", dt);
  const double max_dt = 0.1; // maximum possible time-step
  
  struct seq_ctx *sc = ctx;
  sc->nfluid += 1;

  double dt_actual = dt;

  if (dt > max_dt)
    dt_actual = max_dt;

  // take time-step of dt_actual
  
  return (struct gkyl_update_fsm_status) {
    .status = GKYL_UPDATE_FSM_STATUS_SUCCESS,
    .dt_actual = dt_actual,
    .dt_suggested = max_dt,
  };
}

struct gkyl_update_fsm_status
seq_src_2(double tcurr, double dt, void *ctx)
{
  //printf("seq_src_2: dt = %g\n", dt);

  struct seq_ctx *sc = ctx;
  sc->nsrc += 1;

  enum gkyl_update_fsm_code status = GKYL_UPDATE_FSM_STATUS_SUCCESS;
  if (dt > 100)
    status = GKYL_UPDATE_FSM_STATUS_FAIL;  

  return (struct gkyl_update_fsm_status) {
    .status = status,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX,
  };
}

void
test_seq_1()
{
  struct seq_ctx ctx = { };

  struct gkyl_update_fsm *seq = gkyl_update_fsm_new(
    3, // number of actions
    (struct gkyl_update_fsm_action[3]) { // actions in sequence
      [SRC_1] = { .ctx = &ctx, .u = seq_src_1 },
      [FLUID] = { .ctx = &ctx, .u = seq_fluid },
      [SRC_2] = { .ctx = &ctx, .u = seq_src_2 }
    },
    3, // number of transitions
    (struct gkyl_update_fsm_transition[3]) { // transition table
      [SRC_1] = { FLUID, GKYL_UPDATE_FSM_REDO },
      [FLUID] = { SRC_2, GKYL_UPDATE_FSM_REDO },
      [SRC_2] = { GKYL_UPDATE_FSM_FINISH, GKYL_UPDATE_FSM_REDO }
    },
    (struct gkyl_update_fsm_action) { .ctx = &ctx, .u = seq_redo } // redo action
  );

  TEST_CHECK( 3 == seq->nactions );
  
  struct gkyl_update_fsm_status status = gkyl_update_fsm_run(seq, 0.0, 1.0);

  TEST_CHECK( GKYL_UPDATE_FSM_STATUS_SUCCESS == status.status );
  TEST_CHECK( 1 == ctx.nredo );
  TEST_CHECK( 2 == ctx.nfluid );
  TEST_CHECK( 3 == ctx.nsrc );

  //seq_ctx_print(ctx);

  TEST_CHECK( 0.1 == status.dt_actual );
  TEST_CHECK( 0.1 == status.dt_suggested );

  status = gkyl_update_fsm_run(seq, 0.0, 200.0); // should abort
  TEST_CHECK( GKYL_UPDATE_FSM_STATUS_FAIL == status.status );

  gkyl_update_fsm_release(seq);
}

TEST_LIST = {
  { "seq_1", test_seq_1 },
  { NULL, NULL },  
};
