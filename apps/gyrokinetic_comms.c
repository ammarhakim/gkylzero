#include <gkyl_gyrokinetic_comms.h>

struct gkyl_rect_decomp* 
gyrokinetic_comms_decomp_new(int cdim, const int *cells_conf, struct gkyl_app_args app_args, FILE *iostream)
{
  // Create global range.
  struct gkyl_range global_range_conf;
  gkyl_create_global_range(cdim, cells_conf, &global_range_conf);

  // Create decomposition.
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi)
      cuts[d] = app_args.cuts[d];
    else
      cuts[d] = 1;
  }
#else
  for (int d = 0; d < cdim; d++) cuts[d] = 1;
#endif

  int my_rank = 0, comm_size = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  }
#endif

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) ncuts *= cuts[d];

  if (ncuts != comm_size) {
    if (my_rank == 0)
      fprintf(iostream, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
#ifdef GKYL_HAVE_MPI
    if (app_args.use_mpi) MPI_Finalize();
#endif
    assert(false);
  }

  for (int d = 0; d < cdim - 1; d++) {
    if (cuts[d] > 1) {
      if (my_rank == 0)
        fprintf(iostream, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts[d], d);
#ifdef GKYL_HAVE_MPI
      if (app_args.use_mpi) MPI_Finalize();
#endif
      assert(false);
    }
  }

  return gkyl_rect_decomp_new_from_cuts(cdim, cuts, &global_range_conf);
}

struct gkyl_comm* 
gyrokinetic_comms_new(struct gkyl_app_args app_args, struct gkyl_rect_decomp *decomp, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm;

#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
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
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
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
