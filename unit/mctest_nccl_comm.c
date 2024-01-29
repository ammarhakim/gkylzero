#include <acutest.h>

#ifdef GKYL_HAVE_NCCL

#include <math.h>
#include <stc/cstr.h>
#include <gkyl_util.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_nccl_comm.h>

void
mpi_allreduce()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  double vals_ho[2];
  if (m_rank == 0) {
    vals_ho[0] = 1.0;
    vals_ho[1] = 2.0;
  }
  if (m_rank == 1) {
    vals_ho[0] = 3.0;
    vals_ho[1] = -1.0;
  }
  double *vals = gkyl_cu_malloc(2*sizeof(double));
  gkyl_cu_memcpy(vals, vals_ho, 2*sizeof(double), GKYL_CU_MEMCPY_H2D);

  double v_max_ho[2], v_min_ho[2], v_sum_ho[2];
  double *v_max = gkyl_cu_malloc(2*sizeof(double));
  double *v_min = gkyl_cu_malloc(2*sizeof(double));
  double *v_sum = gkyl_cu_malloc(2*sizeof(double));

  gkyl_comm_all_reduce(comm_dev, GKYL_DOUBLE, GKYL_MAX, 2, vals, v_max);
  gkyl_cu_memcpy(v_max_ho, v_max, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_max_ho[0] == 3.0 );
  TEST_CHECK( v_max_ho[1] == 2.0 );

  gkyl_comm_all_reduce(comm_dev, GKYL_DOUBLE, GKYL_MIN, 2, vals, v_min);
  gkyl_cu_memcpy(v_min_ho, v_min, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_min_ho[0] == 1.0 );
  TEST_CHECK( v_min_ho[1] == -1.0 );

  gkyl_comm_all_reduce(comm_dev, GKYL_DOUBLE, GKYL_SUM, 2, vals, v_sum);
  gkyl_cu_memcpy(v_sum_ho, v_sum, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_sum_ho[0] == 4.0 );
  TEST_CHECK( v_sum_ho[1] == 1.0 );

  gkyl_cu_free(vals);
  gkyl_cu_free(v_max);
  gkyl_cu_free(v_min);
  gkyl_cu_free(v_sum);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm_ho);
  gkyl_comm_release(comm_dev);
}

void
nccl_n2_array_isend_irecv_2d()
{
  // Test array_send and array_recv with a nonblocking comm.
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });

  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size != 2) return;
  
  int cuts[] = { comm_size, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int rank, sz;
  gkyl_comm_get_rank(comm_ho, &rank);
  gkyl_comm_get_size(comm_ho, &sz);

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  // Allocate send/recv buffers.
  double sendval = rank==0? 20005. : 30005.;
  double recvval = rank==0? 30005. : 20005.;

  // Assume the range is not decomposed.
  struct gkyl_array *arrA = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
  struct gkyl_array *arrB = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
  struct gkyl_array *recvbuff_ho = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_clear(arrA, sendval*(1-rank));
  gkyl_array_clear(arrB, sendval*rank);

  struct gkyl_array *recvbuff = rank==0? arrB : arrA;
  struct gkyl_array *sendbuff = rank==0? arrA : arrB;

  struct gkyl_comm_state *cstate = gkyl_comm_state_new(comm_dev);
  int tag = 13;
  // Communicate data from rank 0 to rank 1.
  if (rank == 1)
    gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
  if (rank == 0)
    gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);

  if (rank == 1) {
    gkyl_comm_state_wait(comm_dev, cstate);

    gkyl_array_copy(recvbuff_ho, recvbuff);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long idx = gkyl_range_idx(&range, iter.idx);
      const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
      TEST_CHECK( f[0] == recvval );
    }
  }

  // Communicate data from rank 1 to rank 0.
  if (rank == 0)
    gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
  if (rank == 1)
    gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);

  if (rank == 0) {
    gkyl_comm_state_wait(comm_dev, cstate);

    gkyl_array_copy(recvbuff_ho, recvbuff);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long idx = gkyl_range_idx(&range, iter.idx);
      const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
      TEST_CHECK( f[0] == recvval );
    }
  }

  // Communicate data between rank 0 and 1 at the same time.
  gkyl_array_clear(recvbuff, 0.);
  gkyl_comm_group_call_start(comm_dev);
  gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
  gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);
  gkyl_comm_state_wait(comm_dev, cstate);
  gkyl_comm_group_call_end(comm_dev);

  gkyl_array_copy(recvbuff_ho, recvbuff);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&range, iter.idx);
    const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
    TEST_CHECK( f[0] == recvval );
  }

  gkyl_comm_barrier(comm_dev);

  gkyl_comm_state_release(comm_dev, cstate);
  gkyl_array_release(recvbuff_ho);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm_dev);
  gkyl_comm_release(comm_ho);
}

  
TEST_LIST = {
  {"mpi_allreduce", mpi_allreduce},
  {"nccl_n2_array_isend_irecv_2d", nccl_n2_array_isend_irecv_2d},
  {NULL, NULL},
};

#else

// nothing to test if not building with NCCL
TEST_LIST = {
  {NULL, NULL},
};

#endif
