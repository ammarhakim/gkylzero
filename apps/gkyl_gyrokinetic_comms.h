#pragma once

#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

/**
 * Create a new decomposition of configuration space for the gyrokinetic app.
 *
 * @param cdim Number of configuration space dimensions.
 * @param cells_conf Number of cells in each conf-space dimension.
 * @param app_args Command line arguments passed to the app.
 * @param Pointer to place where to put error messages.
 * @return New gkyl_rect_decomp decomposition object.
 */
struct gkyl_rect_decomp* 
gkyl_gyrokinetic_comms_decomp_new(int cdim, const int *cells_conf, const int *cuts, bool use_mpi, FILE *iostream);

/**
 * Create a new communicator for the gyrokinetic app.
 *
 * @param app_args Command line argumnts passed to the app.
 * @param decomp Decomposition object.
 * @param Pointer to place where to put error messages.
 * @return New gkyl_comm communicator object.
 */
struct gkyl_comm* 
gkyl_gyrokinetic_comms_new(bool use_mpi, bool use_gpu, FILE *iostream);

/**
 * Free gyrokinetic app decomp and comm objects.
 *
 * @param decomp Decomposition object.
 * @param comm Communicator object.
 */
void
gkyl_gyrokinetic_comms_release(struct gkyl_comm *comm);
