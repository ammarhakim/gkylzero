#ifdef GKYL_HAVE_NCCL

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_nccl_comm.h>
#include <gkyl_util.h>

#define checkNCCL(cmd) do {                         \
  ncclResult_t res = cmd;                           \
  if (res != ncclSuccess) {                         \
    printf("Failed, NCCL error %s:%d '%s'\n",       \
        __FILE__,__LINE__,ncclGetErrorString(res)); \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)

// Private struct wrapping NCCL-specific code
struct nccl_comm {
  struct gkyl_comm base; // base communicator

  ncclComm_t ncomm; // NCCL communicator to use
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition

  cudaStream_t custream; // cuda stream for NCCL comms.
};

// Mapping of Gkeyll type to ncclDataType_t
static ncclDataType_t g2_nccl_datatype[] = {
  [GKYL_INT] = ncclInt,
  [GKYL_INT_64] = ncclInt64,
  [GKYL_FLOAT] = ncclFloat,
  [GKYL_DOUBLE] = ncclDouble,
};

struct gkyl_comm_state {
  ncclComm_t *ncomm;
  int tag;
  cudaStream_t *custream;
  int peer;
};

static struct gkyl_comm_state* comm_state_new()
{
  struct gkyl_comm_state *state = gkyl_malloc(sizeof *state);
  return state;
}

static void comm_state_release(struct gkyl_comm_state *state)
{
  gkyl_free(state);
}

static void comm_state_wait(struct gkyl_comm_state *state)
{
  ncclResult_t nstat;
  do {
    checkNCCL(ncclCommGetAsyncError(*(state->ncomm), &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(*(state->custream)));
}

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  // Finalize NCCL comm.
  ncclCommDestroy(nccl->ncomm);

  gkyl_free(nccl);
}

static int
barrier(struct gkyl_comm *comm)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  ncclResult_t nstat;
  do {
    checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(nccl->custream));
  return 0;
}

static int
array_send(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  ncclResult_t nstat = ncclSend(array->data, vol, g2_nccl_datatype[array->type], dest, nccl->ncomm, nccl->custream);
  do {
    checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(nccl->custream));
  return 1;
}

static int
array_recv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  ncclResult_t nstat = ncclRecv(array->data, vol, g2_nccl_datatype[array->type], src, nccl->ncomm, nccl->custream);
  do {
    checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(nccl->custream));
  return 1;
}

static int
array_isend(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  ncclResult_t nstat = ncclSend(array->data, vol, g2_nccl_datatype[array->type], dest, nccl->ncomm, nccl->custream);
  state->ncomm = &nccl->ncomm;
  state->tag = tag;
  state->custream = &nccl->custream;
  state->peer = dest;
  return 1;
}

static int
array_irecv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  ncclResult_t nstat = ncclRecv(array->data, vol, g2_nccl_datatype[array->type], src, nccl->ncomm, nccl->custream);
  state->ncomm = &nccl->ncomm;
  state->tag = tag;
  state->custream = &nccl->custream;
  state->peer = src;
  return 1;
}

static void
group_call_start()
{
  checkNCCL(ncclGroupStart());
}

static void
group_call_end()
{
  checkNCCL(ncclGroupEnd());
}

struct gkyl_comm*
gkyl_nccl_comm_new(const struct gkyl_nccl_comm_inp *inp)
{
  struct nccl_comm *nccl = gkyl_malloc(sizeof *nccl);
  nccl->decomp = gkyl_rect_decomp_acquire(inp->decomp);

  int rank, size;
  MPI_Comm_rank(inp->mpi_comm, &rank);
  MPI_Comm_size(inp->mpi_comm, &size);

  int num_devices[1];
  checkCuda(cudaGetDeviceCount(num_devices));

  int local_rank = rank % num_devices[0];
  checkCuda(cudaSetDevice(local_rank));

  ncclUniqueId nId;
  if (rank == 0) ncclGetUniqueId(&nId);
  MPI_Bcast(&nId, sizeof(nId), MPI_BYTE, 0, inp->mpi_comm);

  checkCuda(cudaStreamCreate(&nccl->custream));

  // Initialize NCCL comm
  ncclConfig_t config = NCCL_CONFIG_INITIALIZER;
  config.blocking = 0;  // Nonblocking (doesn't block at NCCL calls).
  ncclResult_t nstat = ncclCommInitRankConfig(&nccl->ncomm, size, nId, rank, &config);
  if (config.blocking == 0) {
    do {
      checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
    } while(nstat == ncclInProgress);
  } else {
    checkNCCL(nstat);
  }
  checkCuda(cudaStreamSynchronize(nccl->custream));
  
  nccl->base.barrier = barrier;
  nccl->base.gkyl_array_send = array_send;
  nccl->base.gkyl_array_isend = array_isend;
  nccl->base.gkyl_array_recv = array_recv;
  nccl->base.gkyl_array_irecv = array_irecv;
  nccl->base.comm_state_new = comm_state_new;
  nccl->base.comm_state_release = comm_state_release;
  nccl->base.comm_state_wait = comm_state_wait;
  nccl->base.comm_group_call_start = group_call_start;
  nccl->base.comm_group_call_end = group_call_end;
  nccl->base.ref_count = gkyl_ref_count_init(comm_free);

  return &nccl->base;
}

#endif
