#include <acutest.h>
#include <gkyl_null_comm.h>

void
test_1()
{
  struct gkyl_comm *comm = gkyl_null_comm_new();

  gkyl_comm_release(comm);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
