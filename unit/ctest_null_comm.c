#include <acutest.h>
#include <gkyl_null_comm.h>

void
test_1()
{
  struct gkyl_comm *comm = gkyl_null_comm_new();

  int rank;
  comm->get_rank(comm, &rank);
  
  TEST_CHECK( rank == 0 );

  double out[3], inp[3] = { 2.0, 4.0, 8.0 };
  comm->all_reduce(comm, GKYL_DOUBLE, GKYL_MIN, 3, inp, out);

  for (int i=0; i<3; ++i)
    TEST_CHECK( out[i] == inp[i] );

  gkyl_comm_release(comm);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
