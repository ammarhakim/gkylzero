#pragma once

#include <gkyl_comm.h>

/**
 * Return a new "null" communicator, i.e. a communicator for a single
 * core calculation.
 *
 * @return New communicator
 */
struct gkyl_comm *gkyl_null_comm_new(void);

