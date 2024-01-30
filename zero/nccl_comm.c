#ifdef GKYL_HAVE_NCCL

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_nccl_comm.h>
#include <gkyl_util.h>

// Maximum number of recv neighbors: not sure hard-coding this is a
// good idea.
#define MAX_RECV_NEIGH 32

#define NCCL_BASE_TAG 4343
#define NCCL_BASE_PER_TAG 5353

// Some NCCL calls should return ncclSuccess when done properly,
// but others (e.g. ncclGroupEnd) may return ncclInProgress.
// If a function that should return ncclSuccess returns ncclInProgress
// for some reason, having only one check function may be a problem.
// We could create a separate check function which waits and times out
// after a set amount of time.
#define checkNCCL(cmd) do {                           \
  ncclResult_t res = cmd;                             \
  if (res != ncclSuccess  && res != ncclInProgress) { \
    printf("Failed, NCCL error %s:%d '%s'\n",         \
        __FILE__,__LINE__,ncclGetErrorString(res));   \
    exit(EXIT_FAILURE);                               \
  }                                                   \
} while(0)

// Mapping of Gkeyll type to ncclDataType_t
static ncclDataType_t g2_nccl_datatype[] = {
  [GKYL_INT] = ncclInt,
  [GKYL_INT_64] = ncclInt64,
  [GKYL_FLOAT] = ncclFloat,
  [GKYL_DOUBLE] = ncclDouble,
};

// Mapping of Gkeyll ops to ncclRedOp_t.
static ncclRedOp_t g2_nccl_op[] = {
  [GKYL_MIN] = ncclMin,
  [GKYL_MAX] = ncclMax,
  [GKYL_SUM] = ncclSum,
};

struct gkyl_comm_state {
  ncclComm_t *ncomm;
  int tag;
  cudaStream_t *custream;
  int peer;
};

// Receive data
struct comm_buff_stat {
  struct gkyl_range range;
//  MPI_Request status;
  gkyl_mem_buff buff;
};

// Private struct wrapping NCCL-specific code
struct nccl_comm {
  struct gkyl_comm base; // base communicator.

  int rank; // Process ID in this communicator.
  int size; // Size of this communicator.

  ncclComm_t ncomm; // NCCL communicator to use.
  cudaStream_t custream; // Cuda stream for NCCL comms.
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // Whether to sync corners.

  struct gkyl_rect_decomp_neigh *neigh; // neighbors of local region
  struct gkyl_rect_decomp_neigh *per_neigh[GKYL_MAX_DIM]; // periodic neighbors

  int nrecv; // number of elements in rinfo array
  struct comm_buff_stat recv[MAX_RECV_NEIGH]; // info for recv data

  int nsend; // number of elements in sinfo array
  struct comm_buff_stat send[MAX_RECV_NEIGH]; // info for send data

  MPI_Comm mpi_comm; // MPI comm this NCCL comm derives from.
};

static struct gkyl_comm_state* comm_state_new(struct gkyl_comm *comm)
{
  struct gkyl_comm_state *state = gkyl_malloc(sizeof *state);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  state->ncomm = &nccl->ncomm;
  state->custream = &nccl->custream;

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
get_rank(struct gkyl_comm *comm, int *rank)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  *rank = nccl->rank;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  *sz = nccl->size;
  return 0;
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
  return 0;
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
  return 0;
}

static int
array_isend(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  checkNCCL(ncclSend(array->data, vol, g2_nccl_datatype[array->type], dest, nccl->ncomm, nccl->custream));
  state->tag = tag;
  state->peer = dest;
  return 0;
}

static int
array_irecv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  checkNCCL(ncclRecv(array->data, vol, g2_nccl_datatype[array->type], src, nccl->ncomm, nccl->custream));
  state->tag = tag;
  state->peer = src;
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  checkNCCL(ncclAllReduce(inp, out, nelem, g2_nccl_datatype[type], g2_nccl_op[op], nccl->ncomm, nccl->custream));
  return 0;
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

static int
array_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, struct gkyl_array *array)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<nccl->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  int tag = NCCL_BASE_TAG;

  checkNCCL(ncclGroupStart());

  // post nonblocking recv to get data into ghost-cells  
  int nridx = 0;
  for (int n=0; n<nccl->neigh->num_neigh; ++n) {
    int nid = nccl->neigh->neigh[n];
    
    int isrecv = gkyl_sub_range_intersect(
      &nccl->recv[nridx].range, local_ext, &nccl->decomp->ranges[nid]);
    size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;

    if (isrecv) {
      if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
        gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);
      
      checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[nridx].buff),
        recv_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

      nridx += 1;
    }
  }

  // post non-blocking sends of skin-cell data to neighbors
  int nsidx = 0;
  for (int n=0; n<nccl->neigh->num_neigh; ++n) {
    int nid = nccl->neigh->neigh[n];
    
    struct gkyl_range neigh_ext;
    gkyl_range_extend(&neigh_ext, &nccl->decomp->ranges[nid], elo, eup);

    int issend = gkyl_sub_range_intersect(
      &nccl->send[nsidx].range, local, &neigh_ext);
    size_t send_vol = array->esznc*nccl->send[nsidx].range.volume;

    if (issend) {
      if (gkyl_mem_buff_size(nccl->send[nsidx].buff) < send_vol)
        gkyl_mem_buff_resize(nccl->send[nsidx].buff, send_vol);
      
      gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->send[nsidx].buff),
        array, &(nccl->send[nsidx].range));
      
      checkNCCL(ncclSend(gkyl_mem_buff_data(nccl->send[nsidx].buff),
        send_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

      nsidx += 1;
    }
  }
  checkNCCL(ncclGroupEnd());

  // Complete sends and recvs.
  ncclResult_t nstat;
  do {
    checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(nccl->custream));
  
  // Copy data into ghost-cells.
  for (int r=0; r<nridx; ++r) {
    int isrecv = nccl->recv[r].range.volume;
    if (isrecv) {
      gkyl_array_copy_from_buffer(array,
        gkyl_mem_buff_data(nccl->recv[r].buff),
        &(nccl->recv[r].range)
      );
    }
  }

  return 0;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, nccl->decomp);
  struct gkyl_comm *ext_comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = nccl->mpi_comm,
      .decomp = ext_decomp,
      .sync_corners = nccl->sync_corners,
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  return ext_comm;
}

static struct gkyl_comm*
split_comm(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *new_decomp)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  MPI_Comm new_mpi_comm;
  int ret = MPI_Comm_split(nccl->mpi_comm, color, nccl->rank, &new_mpi_comm);
  assert(ret == MPI_SUCCESS);

  struct gkyl_comm *newcomm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = new_mpi_comm,
      .decomp = new_decomp,
    }
  );
  return newcomm;
}

struct gkyl_comm*
gkyl_nccl_comm_new(const struct gkyl_nccl_comm_inp *inp)
{
  struct nccl_comm *nccl = gkyl_malloc(sizeof *nccl);

  nccl->mpi_comm = inp->mpi_comm;
  MPI_Comm_rank(inp->mpi_comm, &nccl->rank);
  MPI_Comm_size(inp->mpi_comm, &nccl->size);

  int num_devices[1];
  checkCuda(cudaGetDeviceCount(num_devices));

  int local_rank = nccl->rank % num_devices[0];
  checkCuda(cudaSetDevice(local_rank));

  ncclUniqueId nId;
  if (nccl->rank == 0) ncclGetUniqueId(&nId);
  MPI_Bcast(&nId, sizeof(nId), MPI_BYTE, 0, inp->mpi_comm);

  checkCuda(cudaStreamCreate(&nccl->custream));

  // Initialize NCCL comm
  ncclConfig_t config = NCCL_CONFIG_INITIALIZER;
  config.blocking = 0;  // Nonblocking (doesn't block at NCCL calls).
  ncclResult_t nstat = ncclCommInitRankConfig(&nccl->ncomm, nccl->size, nId, nccl->rank, &config);
  if (config.blocking == 0) {
    do {
      checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
    } while(nstat == ncclInProgress);
  } else {
    checkNCCL(nstat);
  }
  checkCuda(cudaStreamSynchronize(nccl->custream));

  nccl->sync_corners = inp->sync_corners;
  nccl->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  nccl->neigh = gkyl_rect_decomp_calc_neigh(nccl->decomp, inp->sync_corners, nccl->rank);
  for (int d=0; d<nccl->decomp->ndim; ++d)
    nccl->per_neigh[d] =
      gkyl_rect_decomp_calc_periodic_neigh(nccl->decomp, d, false, nccl->rank);

  nccl->nrecv = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    nccl->recv[i].buff = gkyl_mem_buff_cu_new(16);

  nccl->nsend = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    nccl->send[i].buff = gkyl_mem_buff_cu_new(16);
  
  nccl->base.get_rank = get_rank;
  nccl->base.get_size = get_size;
  nccl->base.barrier = barrier;
  nccl->base.all_reduce = all_reduce;
  nccl->base.gkyl_array_send = array_send;
  nccl->base.gkyl_array_isend = array_isend;
  nccl->base.gkyl_array_recv = array_recv;
  nccl->base.gkyl_array_irecv = array_irecv;
  nccl->base.gkyl_array_sync = array_sync;
  nccl->base.comm_state_new = comm_state_new;
  nccl->base.comm_state_release = comm_state_release;
  nccl->base.comm_state_wait = comm_state_wait;
  nccl->base.comm_group_call_start = group_call_start;
  nccl->base.comm_group_call_end = group_call_end;
  nccl->base.extend_comm = extend_comm;
  nccl->base.split_comm = split_comm;
  nccl->base.ref_count = gkyl_ref_count_init(comm_free);

  return &nccl->base;
}

#endif
