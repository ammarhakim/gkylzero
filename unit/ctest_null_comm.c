#include <acutest.h>
#include <gkyl_null_comm.h>

void
test_1()
{
  struct gkyl_comm *comm = gkyl_null_comm_new();

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  
  TEST_CHECK( rank == 0 );

  double out[3], inp[3] = { 2.0, 4.0, 8.0 };
  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MIN, 3, inp, out);

  for (int i=0; i<3; ++i)
    TEST_CHECK( out[i] == inp[i] );

  gkyl_comm_release(comm);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
