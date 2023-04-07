#pragma once

#ifdef GKYL_HAVE_MPI

#include <mpi.h>
#include <gkyl_comm.h>
#include <gkyl_rect_decomp.h>

// input to create new MPI communicator
struct gkyl_mpi_comm_inp {
  MPI_Comm comm;
  struct gkyl_rect_decomp *decomp;
};

/**
 * Return a new MPI communicator
 *
 * @param inp Input struct to use for initialiation
 * @return New MPI communicator
 */
struct gkyl_comm *gkyl_mpi_comm_new(const struct gkyl_mpi_comm_inp *inp);

#endif
