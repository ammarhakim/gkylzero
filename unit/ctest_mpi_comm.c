#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <gkyl_mpi_comm.h>

void
test_1()
{
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .comm = MPI_COMM_WORLD,
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  TEST_CHECK( rank == 0 );

  int sz;
  gkyl_comm_get_size(comm, &sz);
  TEST_CHECK( sz == 1 );

  gkyl_comm_release(comm);
}

TEST_LIST = {
    {"test_1", test_1},
    {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
    {NULL, NULL},
};

#endif
