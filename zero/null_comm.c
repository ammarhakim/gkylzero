#include <gkyl_alloc.h>
#include <gkyl_null_comm.h>

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  gkyl_free(comm);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = 0;
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *data,
  void *out)
{

  return 0;
}

struct gkyl_comm*
gkyl_null_comm_new(void)
{
  struct gkyl_comm *comm = gkyl_malloc(sizeof *comm);

  comm->get_rank = get_rank;
  comm->all_reduce = all_reduce;

  comm->ref_count = gkyl_ref_count_init(comm_free);

  return comm;
}

void
gkyl_null_comm_release(struct gkyl_comm *comm)
{
  gkyl_free(comm);
}
