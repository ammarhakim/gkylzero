#pragma once

#ifdef GKYL_HAVE_MPI

#include <mpi.h>

#include <gkyl_comm.h>

/**
 * Return a new MPI communicator
 *
 * @param comm Base MPI communicator to use
 * @return New MPI communicator
 */
struct gkyl_comm *gkyl_mpi_comm_new(MPI_Comm comm);

#endif
