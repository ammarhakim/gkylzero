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

#include <rt_arg_parse.h>

/**
 * Create a new decomposition of configuration space for the gyrokinetic app.
 *
 * @param cdim Number of configuration space dimensions.
 * @param cells_conf Number of cells in each conf-space dimension.
 * @param app_args Command line argumnts passed to the app.
 * @return New gkyl_rect_decomp decomposition object.
 */
struct gkyl_rect_decomp* 
gyrokinetic_comms_decomp_new(int cdim, const int *cells_conf, struct gkyl_app_args app_args);

/**
 * Create a new communicator for the gyrokinetic app.
 *
 * @param app_args Command line argumnts passed to the app.
 * @param decomp Decomposition object.
 * @return New gkyl_comm communicator object.
 */
struct gkyl_comm* 
gyrokinetic_comms_new(struct gkyl_app_args app_args, struct gkyl_rect_decomp *decomp);

/**
 * Free gyrokinetic app decomp and comm objects.
 *
 * @param decomp Decomposition object.
 * @param comm Communicator object.
 */
void
gyrokinetic_comms_release(struct gkyl_rect_decomp *decomp, struct gkyl_comm *comm);
