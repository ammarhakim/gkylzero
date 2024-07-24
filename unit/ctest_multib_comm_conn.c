#include <acutest.h>
#include <gkyl_multib_comm_conn.h>

static void
test_0(void)
{
  struct gkyl_comm_conn cclist[] = {
    { .rank = 1 },
    { .rank = 2 },
  };

  struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new(2, cclist);

  TEST_CHECK( mbcc->num_comm_conn == 2 );

  gkyl_multib_comm_conn_release(mbcc);
}

TEST_LIST = {
  { "test_0", test_0 },
  { NULL, NULL },
};
