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
test_1()
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
  gkyl_array_clear(arrA, sendval*(1-rank));
  gkyl_array_clear(arrB, sendval*rank);

  struct gkyl_array *recvbuff = rank==0? arrB : arrA;
  struct gkyl_array *sendbuff = rank==0? arrA : arrB;

  struct gkyl_comm_state *cstate = gkyl_comm_state_new(comm_dev);
  int tag = 13;
  // Post irecv before send.
  if (rank == 1) 
    gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
  if (rank == 0)
    gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);
  if (rank == 1) 
    gkyl_comm_state_wait(comm_dev, cstate);

  gkyl_comm_barrier(comm_dev);

  gkyl_comm_state_release(comm_dev, cstate);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm_dev);
  gkyl_comm_release(comm_ho);
}

  
TEST_LIST = {
  {"test_1", test_1},
  {NULL, NULL},
};

#else

// nothing to test if not building with NCCL
TEST_LIST = {
  {NULL, NULL},
};

#endif
