#include <float.h>
#include <gkyl_update_fsm.h>

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
