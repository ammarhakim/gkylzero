#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <math.h>
#include <mpi.h>

#include <gkyl_mpi_comm.h>

void
test_1()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  TEST_CHECK( rank == m_rank );

  int sz;
  gkyl_comm_get_size(comm, &sz);
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  TEST_CHECK( sz == m_sz );

  gkyl_comm_release(comm);
}

void
test_n2()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

  double vals[2];
  if (m_rank == 0) {
    vals[0] = 1.0;
    vals[1] = 2.0;
  }
  if (m_rank == 1) {
    vals[0] = 3.0;
    vals[1] = -1.0;
  }

  double v_max[2], v_min[2], v_sum[2];

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MAX, 2, vals, v_max);
  TEST_CHECK( v_max[0] == 3.0 );
  TEST_CHECK( v_max[1] == 2.0 );

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MIN, 2, vals, v_min);
  TEST_CHECK( v_min[0] == 1.0 );
  TEST_CHECK( v_min[1] == -1.0 );

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_SUM, 2, vals, v_sum);
  TEST_CHECK( v_sum[0] == 4.0 );
  TEST_CHECK( v_sum[1] == 1.0 );

  gkyl_comm_release(comm);
}

TEST_LIST = {
    {"test_1", test_1},
    {"test_n2", test_n2},
    {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
    {NULL, NULL},
};

#endif
