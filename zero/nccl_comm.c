#ifdef GKYL_HAVE_NCCL

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_nccl_comm.h>
#include <gkyl_util.h>

#include <errno.h>

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
#define checkNCCL(cmd) do {                            \
  ncclResult_t res = cmd;                              \
  if (res != ncclSuccess  && res != ncclInProgress) {  \
    fprintf(stderr, "Failed, NCCL error %s:%d '%s'\n", \
        __FILE__,__LINE__,ncclGetErrorString(res));    \
    exit(EXIT_FAILURE);                                \
  }                                                    \
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

// Mapping of Gkeyll type to MPI_Datatype
static MPI_Datatype g2_mpi_datatype[] = {
  [GKYL_INT] = MPI_INT,
  [GKYL_INT_64] = MPI_INT64_T,
  [GKYL_FLOAT] = MPI_FLOAT,
  [GKYL_DOUBLE] = MPI_DOUBLE
};

// Mapping of Gkeyll ops to MPI_Op
static MPI_Op g2_mpi_op[] = {
  [GKYL_MIN] = MPI_MIN,
  [GKYL_MAX] = MPI_MAX,
  [GKYL_SUM] = MPI_SUM
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
  MPI_Comm mcomm; // MPI communicator to use
  struct gkyl_comm *mpi_comm; // MPI comm this NCCL comm derives from.
  bool has_decomp; // Whether this comm is associated with a decomposition (e.g. of a range)
  cudaStream_t custream; // Cuda stream for NCCL comms.
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // Whether to sync corners.
  long local_range_offset; // Offset of the local region.

  struct gkyl_rect_decomp_neigh *neigh; // neighbors of local region
  struct gkyl_rect_decomp_neigh *per_neigh[GKYL_MAX_DIM]; // periodic neighbors

  int nrecv; // number of elements in rinfo array
  struct comm_buff_stat recv[MAX_RECV_NEIGH]; // info for recv data

  int nsend; // number of elements in sinfo array
  struct comm_buff_stat send[MAX_RECV_NEIGH]; // info for send data

  struct gkyl_range dir_edge; // for use in computing tags
  int is_on_edge[2][GKYL_MAX_DIM]; // flags to indicate if local range is on edge
  bool touches_any_edge; // true if this range touches any edge

  // buffers for for allgather
  struct comm_buff_stat allgather_buff_local; 
  struct comm_buff_stat allgather_buff_global; 
};

static struct gkyl_comm_state *
comm_state_new(struct gkyl_comm *comm)
{
  struct gkyl_comm_state *state = gkyl_malloc(sizeof *state);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  state->ncomm = &nccl->ncomm;
  state->custream = &nccl->custream;

  return state;
}

static void
comm_state_release(struct gkyl_comm_state *state)
{
  gkyl_free(state);
}

static void
comm_state_wait(struct gkyl_comm_state *state)
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

  if (nccl->has_decomp) {
    int ndim = nccl->decomp->ndim;
    gkyl_rect_decomp_release(nccl->decomp);

    gkyl_rect_decomp_neigh_release(nccl->neigh);
    for (int d=0; d<ndim; ++d)
      gkyl_rect_decomp_neigh_release(nccl->per_neigh[d]);

    for (int i=0; i<MAX_RECV_NEIGH; ++i)
      gkyl_mem_buff_release(nccl->recv[i].buff);

    for (int i=0; i<MAX_RECV_NEIGH; ++i)
      gkyl_mem_buff_release(nccl->send[i].buff);

    gkyl_mem_buff_release(nccl->allgather_buff_local.buff);
    gkyl_mem_buff_release(nccl->allgather_buff_global.buff);
  }

  // Finalize NCCL comm.
  checkCuda(cudaStreamSynchronize(nccl->custream));
  checkCuda(cudaDeviceSynchronize());
  ncclCommDestroy(nccl->ncomm);

  gkyl_comm_release(nccl->mpi_comm);
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
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_array_meta *meta,
  const struct gkyl_array *arr, const char *fname)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  return gkyl_comm_array_write(nccl->mpi_comm, grid, range, meta, arr, fname);
}

static int
array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  return gkyl_comm_array_read(nccl->mpi_comm, grid, range, arr, fname);
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
allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  checkNCCL(ncclAllReduce(inp, out, nelem, g2_nccl_datatype[type], g2_nccl_op[op], nccl->ncomm, nccl->custream));
  return 0;
}

static int
allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);  
  return gkyl_comm_allreduce(nccl->mpi_comm, type, op, nelem, inp, out);
}

static int
array_allgather(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global, 
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  assert(array_global->esznc == array_local->esznc);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  struct gkyl_range gather_range;

  int rank = nccl->rank;

  assert(local->volume == nccl->decomp->ranges[rank].volume);
  assert(global->volume == nccl->decomp->parent_range.volume);

  // potentially re-size local buffer volume
  size_t send_vol = array_local->esznc*nccl->decomp->ranges[rank].volume;
  if (gkyl_mem_buff_size(nccl->allgather_buff_local.buff) < send_vol)
    gkyl_mem_buff_resize(nccl->allgather_buff_local.buff, send_vol);

  // potentially re-size global buffer volume
  size_t buff_global_vol = array_local->esznc*nccl->decomp->parent_range.volume;
  if (gkyl_mem_buff_size(nccl->allgather_buff_global.buff) < buff_global_vol)
    gkyl_mem_buff_resize(nccl->allgather_buff_global.buff, buff_global_vol);

  // copy data to local buffer
  gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->allgather_buff_local.buff), 
    array_local, local);

  size_t nelem = array_local->esznc*nccl->decomp->ranges[rank].volume;
  // gather data into global buffer
  checkNCCL(ncclAllGather(gkyl_mem_buff_data(nccl->allgather_buff_local.buff),
                          gkyl_mem_buff_data(nccl->allgather_buff_global.buff),
			  nelem, ncclChar, nccl->ncomm, nccl->custream));

  // copy data to global array
  int idx = 0;
  for (int r=0; r<nccl->decomp->ndecomp; ++r) {
    int isrecv = gkyl_sub_range_intersect(
      &gather_range, global, &nccl->decomp->ranges[r]);
    gkyl_array_copy_from_buffer(array_global, 
      gkyl_mem_buff_data(nccl->allgather_buff_global.buff) + idx, &gather_range);
    idx += array_local->esznc*gather_range.volume;
  }

  return 0;
}

static int
array_bcast(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  assert(asend->esznc == arecv->esznc);
  assert(asend->size == arecv->size);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  size_t nelem = asend->ncomp*asend->size;

  checkNCCL(ncclBroadcast(asend->data, arecv->data, nelem, g2_nccl_datatype[asend->type],
      root, nccl->ncomm, nccl->custream));

  return 0;
}

static int
array_bcast_host(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  return gkyl_comm_array_bcast_host(nccl->mpi_comm, asend, arecv, root);
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

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  if (!nccl->touches_any_edge) return 0; // nothing to sync

  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<nccl->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int nsidx = 0;  

  int shift_sign[] = { -1, 1 };

  // post nonblocking recv to get data into ghost-cells
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      checkNCCL(ncclGroupStart());
      if (nccl->is_on_edge[e][dir]) {
        int nid = nccl->per_neigh[dir]->neigh[0];
        if (nid == nccl->rank) {
          int delta[GKYL_MAX_DIM] = { 0 };
          delta[dir] = shift_sign[e]*gkyl_range_shape(&nccl->decomp->parent_range, dir);

          struct gkyl_range neigh_shift, neigh_shift_ext;
          gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);

          int isrecv = gkyl_sub_range_intersect(
            &nccl->recv[nridx].range, local_ext, &neigh_shift);

          delta[dir] *= -1;
          gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);
          gkyl_range_extend(&neigh_shift_ext, &neigh_shift, elo, eup);
          int issend = gkyl_sub_range_intersect(
            &nccl->send[nsidx].range, local, &neigh_shift_ext);
          
          size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;
          if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
            gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

          gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->recv[nridx].buff), array, &(nccl->send[nsidx].range));
          gkyl_array_copy_from_buffer(array, gkyl_mem_buff_data(nccl->recv[nridx].buff), &(nccl->recv[nridx].range));

          nridx += 1;
          nsidx += 1;
        } else {
          int delta[GKYL_MAX_DIM] = { 0 };
          delta[dir] = shift_sign[e]*gkyl_range_shape(&nccl->decomp->parent_range, dir);

          if (nccl->per_neigh[dir]->num_neigh == 1) { // really should  be a loop
//            int nid = nccl->per_neigh[dir]->neigh[0];

            struct gkyl_range neigh_shift;
            gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);

            int isrecv = gkyl_sub_range_intersect(
              &nccl->recv[nridx].range, local_ext, &neigh_shift);
            size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;

            if (isrecv) {
              if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
                gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

              checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[nridx].buff),
                recv_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

              nridx += 1;
            }

            struct gkyl_range neigh_shift_ext;
            gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);
            gkyl_range_extend(&neigh_shift_ext, &neigh_shift, elo, eup);

            int issend = gkyl_sub_range_intersect(
              &nccl->send[nsidx].range, local, &neigh_shift_ext);
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
        }
      }
      checkNCCL(ncclGroupEnd());

      // Complete sends and recvs.
      ncclResult_t nstat;
      do {
        checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
      } while(nstat == ncclInProgress);
      checkCuda(cudaStreamSynchronize(nccl->custream));
    }
  }

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

  nccl->nrecv = nridx > nccl->nrecv ? nridx : nccl->nrecv;
  
  return 0;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);

  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, nccl->decomp);
  struct gkyl_comm *ext_comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = nccl->mcomm,
      .decomp = ext_decomp,
      .sync_corners = nccl->sync_corners,
      .device_set = 1,
      .custream = nccl->custream,
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  return ext_comm;
}

static struct gkyl_comm*
split_comm(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *new_decomp)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, base);
  MPI_Comm new_mcomm;
  int ret = MPI_Comm_split(nccl->mcomm, color, nccl->rank, &new_mcomm);
  assert(ret == MPI_SUCCESS);

  struct gkyl_comm *newcomm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = new_mcomm,
      .decomp = new_decomp,
      .device_set = 1,
      .custream = nccl->custream,
    }
  );
  return newcomm;
}

struct gkyl_comm*
gkyl_nccl_comm_new(const struct gkyl_nccl_comm_inp *inp)
{
  struct nccl_comm *nccl = gkyl_malloc(sizeof *nccl);
  nccl->mcomm = inp->mpi_comm;

  nccl->mpi_comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp)
      {
        .mpi_comm = nccl->mcomm,
        .decomp = inp->decomp, 
        .sync_corners = inp->sync_corners, 
      }
    );
  MPI_Comm_rank(nccl->mcomm, &nccl->rank);
  MPI_Comm_size(nccl->mcomm, &nccl->size);

  if (inp->device_set == 0) {
    int num_devices[1];
    checkCuda(cudaGetDeviceCount(num_devices));
  
    int local_rank = nccl->rank % num_devices[0];
    checkCuda(cudaSetDevice(local_rank));
  }

  ncclUniqueId nId;
  if (nccl->rank == 0) ncclGetUniqueId(&nId);
  MPI_Bcast((void *)&nId, sizeof(nId), MPI_BYTE, 0, nccl->mcomm);

  if (inp->custream == 0)
    checkCuda(cudaStreamCreate(&nccl->custream));
  else
    nccl->custream = inp->custream;

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

  nccl->has_decomp = false;
  // In case this nccl_comm purely an object holding an ncclComm,
  // not associated with any range nor decomposition of it.
  if (inp->decomp != 0) { 
    nccl->has_decomp = true;

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

    nccl->allgather_buff_local.buff = gkyl_mem_buff_cu_new(16);
    nccl->allgather_buff_global.buff = gkyl_mem_buff_cu_new(16);

    gkyl_range_init(&nccl->dir_edge, 2, (int[]) { 0, 0 }, (int[]) { GKYL_MAX_DIM, 2 });

    int num_touches = 0;
    for (int d=0; d<nccl->decomp->ndim; ++d) {
      nccl->is_on_edge[0][d] = gkyl_range_is_on_lower_edge(
        d, &nccl->decomp->ranges[nccl->rank], &nccl->decomp->parent_range);
      nccl->is_on_edge[1][d] = gkyl_range_is_on_upper_edge(
        d, &nccl->decomp->ranges[nccl->rank], &nccl->decomp->parent_range);
      num_touches += nccl->is_on_edge[0][d] + nccl->is_on_edge[1][d];
    }
    nccl->touches_any_edge = num_touches > 0 ? true : false;
  
    nccl->local_range_offset = gkyl_rect_decomp_calc_offset(nccl->decomp, nccl->rank);

    nccl->base.gkyl_array_allgather = array_allgather;
    nccl->base.gkyl_array_sync = array_sync;
    nccl->base.gkyl_array_per_sync = array_per_sync;
    nccl->base.gkyl_array_write = array_write;
    nccl->base.gkyl_array_read = array_read;
  }
  
  nccl->base.get_rank = get_rank;
  nccl->base.get_size = get_size;
  nccl->base.barrier = barrier;
  nccl->base.allreduce = allreduce;
  nccl->base.allreduce_host = allreduce_host;
  nccl->base.gkyl_array_send = array_send;
  nccl->base.gkyl_array_isend = array_isend;
  nccl->base.gkyl_array_recv = array_recv;
  nccl->base.gkyl_array_irecv = array_irecv;
  nccl->base.gkyl_array_bcast = array_bcast;
  nccl->base.gkyl_array_bcast_host = array_bcast_host;
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
