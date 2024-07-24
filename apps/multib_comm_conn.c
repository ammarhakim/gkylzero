#include <gkyl_alloc.h>
#include <gkyl_multib_comm_conn.h>

#include <string.h>

struct multib_comm_conn {
  struct gkyl_multib_comm_conn mcc;
};

static void
multib_comm_conn_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_multib_comm_conn *mcc = container_of(ref, struct gkyl_multib_comm_conn, ref_count);
  struct multib_comm_conn *cconn = container_of(mcc, struct multib_comm_conn, mcc);
  
  gkyl_free(cconn->mcc.comm_conn);
  gkyl_free(cconn);
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new(int num, const struct gkyl_comm_conn *comm_conn)
{
  struct multib_comm_conn *cconn = gkyl_malloc(sizeof *cconn);
  cconn->mcc.num_comm_conn = num;
  cconn->mcc.comm_conn = gkyl_malloc(sizeof(struct gkyl_comm_conn[num]));

  for (int i=0; i<num; ++i)
    memcpy(&cconn->mcc.comm_conn[i], &comm_conn[i], sizeof(struct gkyl_comm_conn));

  cconn->mcc.ref_count = gkyl_ref_count_init(multib_comm_conn_free);

  return &cconn->mcc;
}

void
gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn)
{
  gkyl_ref_count_dec(&cconn->ref_count);
}
