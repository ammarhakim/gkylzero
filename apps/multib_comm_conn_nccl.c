#include <gkyl_multib_comm_conn_priv.h>
#include <gkyl_comm_priv.h>
#include <gkyl_nccl_comm_priv.h>

int
gkyl_multib_comm_conn_array_transfer_nccl(struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  checkNCCL(ncclGroupStart());

  // post nonblocking recv to get data into ghost-cells  
  int nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];

    for (int n=0; n<mbcc_r->num_comm_conn; ++n) {
      int nid = mbcc_r->comm_conn[n].rank;
      int bid = mbcc_r->comm_conn[n].block_id;
      
      size_t recv_vol = arr_recv[bI]->esznc*mbcc_r->comm_conn[n].range.volume;

      if (recv_vol>0) {
        if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
          gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

	checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[nridx].buff),
          recv_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

        nridx += 1;
      }
    }
  }

  // post non-blocking sends of skin-cell data to neighbors
  int nsidx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bI];

    for (int n=0; n<mbcc_s->num_comm_conn; ++n) {
      int nid = mbcc_s->comm_conn[n].rank;
    
      size_t send_vol = arr_send[bI]->esznc*mbcc_s->comm_conn[n].range.volume;

      if (send_vol>0) {
        if (gkyl_mem_buff_size(nccl->send[nsidx].buff) < send_vol)
          gkyl_mem_buff_resize(nccl->send[nsidx].buff, send_vol);
        
        gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->send[nsidx].buff),
          arr_send[bI], &mbcc_s->comm_conn[n].range);

	checkNCCL(ncclSend(gkyl_mem_buff_data(nccl->send[nsidx].buff),
          send_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

        nsidx += 1;
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

  // Copy data into ghost-cells.
  nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];
    for (int n=0; n<mbcc_r->num_comm_conn; ++n) {
      int isrecv = mbcc_r->comm_conn[n].range.volume;
      if (isrecv) {
        gkyl_array_copy_from_buffer(arr_recv[bI],
          gkyl_mem_buff_data(nccl->recv[nridx].buff),
          &(mbcc_r->comm_conn[n].range)
        );
        nridx += 1;
      }
    }
  }

  return 0;
}
