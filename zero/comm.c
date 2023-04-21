#include <gkyl_comm.h>

struct gkyl_comm*
gkyl_comm_acquire(const struct gkyl_comm *comm)
{
  gkyl_ref_count_inc(&comm->ref_count);
  return (struct gkyl_comm*) comm;
}

void
gkyl_comm_release(const struct gkyl_comm *comm)
{
  gkyl_ref_count_dec(&comm->ref_count);
}
