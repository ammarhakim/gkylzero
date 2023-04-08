
#ifdef GKYL_HAVE_MPI

#include <gkyl_alloc.h>
#include <gkyl_elem_type.h>
#include <gkyl_mpi_comm.h>

// Mapping of Gkeyll type to MPI_Datatype
static MPI_Datatype g2_mpi_datatype[] = {
  [GKYL_INT] = MPI_INT,
  [GKYL_FLOAT] = MPI_FLOAT,
  [GKYL_DOUBLE] = MPI_DOUBLE
};

// Mapping of Gkeyll ops to MPI_Op
static MPI_Op g2_mpi_op[] = {
  [GKYL_MIN] = MPI_MIN,
  [GKYL_MAX] = MPI_MAX,
  [GKYL_SUM] = MPI_SUM
};

// Private struct wrapping MPI-specific code
struct mpi_comm {
  struct gkyl_comm base; // base communicator

  MPI_Comm mcomm; // MPI communicator to use
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition2
};

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  gkyl_rect_decomp_release(mpi->decomp);
  gkyl_free(mpi);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Comm_rank(mpi->mcomm, rank);
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Comm_size(mpi->mcomm, sz);
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);  
  int ret =
    MPI_Allreduce(inp, out, nelem, g2_mpi_datatype[type], g2_mpi_op[op], mpi->mcomm);
  return ret == MPI_SUCCESS ? 0 : 1;
}

/* static int */
/* array_sync(struct gkyl_comm *comm, struct gkyl_array *array) */
/* { */
/*   return 0; */
/* } */

static int
barrier(struct gkyl_comm *comm)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Barrier(mpi->mcomm);
  return 0;
}

struct gkyl_comm*
gkyl_mpi_comm_new(const struct gkyl_mpi_comm_inp *inp)
{
  struct mpi_comm *mpi = gkyl_malloc(sizeof *mpi);
  mpi->mcomm = inp->comm;
  mpi->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  
  mpi->base.get_rank = get_rank;
  mpi->base.get_size = get_size;
  mpi->base.barrier = barrier;
  mpi->base.all_reduce = all_reduce;

  mpi->base.ref_count = gkyl_ref_count_init(comm_free);

  return &mpi->base;
}

#endif
