#include <gkyl_gyrokinetic_comms.h>
#include <assert.h>

struct gkyl_rect_decomp* 
gkyl_gyrokinetic_comms_decomp_new(int cdim, const int *cells_conf, const int *cuts, bool use_mpi, FILE *iostream)
{
  // Create global range.
  struct gkyl_range global_range_conf;
  gkyl_create_global_range(cdim, cells_conf, &global_range_conf);

  // Create decomposition.
  int cuts_used[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++) {
    if (use_mpi)
      cuts_used[d] = cuts[d];
    else
      cuts_used[d] = 1;
  }
#else
  for (int d = 0; d < cdim; d++) cuts_used[d] = 1;
#endif

  int my_rank = 0, comm_size = 1;
#ifdef GKYL_HAVE_MPI
  if (use_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  }
#endif

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) ncuts *= cuts_used[d];

  if (ncuts != comm_size) {
    if (my_rank == 0)
      fprintf(iostream, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
#ifdef GKYL_HAVE_MPI
    if (use_mpi) MPI_Finalize();
#endif
    assert(false);
  }

  for (int d = 0; d < cdim - 1; d++) {
    if (cuts_used[d] > 1) {
      if (my_rank == 0)
        fprintf(iostream, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts_used[d], d);
#ifdef GKYL_HAVE_MPI
      if (use_mpi) MPI_Finalize();
#endif
      assert(false);
    }
  }

  return gkyl_rect_decomp_new_from_cuts(cdim, cuts_used, &global_range_conf);
}

struct gkyl_comm* 
gkyl_gyrokinetic_comms_new(bool use_mpi, bool use_gpu, struct gkyl_rect_decomp *decomp, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
#else
    fprintf(iostream, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = use_gpu
    }
  );
#endif

  return comm;
}

void
gyrokinetic_comms_release(struct gkyl_rect_decomp *decomp, struct gkyl_comm *comm)
{
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}
