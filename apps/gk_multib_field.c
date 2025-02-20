#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>
#include <gkyl_rrobin_decomp.h>
#include <assert.h>
#include <float.h>
#include <time.h>

// compute total number of ranges specified by cuts
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// Initialize multib field object
struct gk_multib_field* 
gk_multib_field_new(const struct gkyl_gyrokinetic_multib *mbinp, struct gkyl_gyrokinetic_multib_app *mbapp)
{

  struct gk_multib_field *mbf = gkyl_malloc(sizeof(struct gk_multib_field));

  mbf->info = mbinp->field;
  mbf->gkfield_id = mbf->info.gkfield_id ? mbf->info.gkfield_id : GKYL_GK_FIELD_ES;

  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbapp->comm, &my_rank);
  gkyl_comm_get_size(mbapp->comm, &num_ranks);
  int rank_list[num_ranks];

  int num_blocks = mbapp->block_topo->num_blocks;
  mbf->num_local_blocks = mbapp->num_local_blocks;
  int *local_blocks = mbapp->local_blocks;
  mbf->cdim = mbapp->block_topo->ndim;

  // Get branks (number of cuts per block)
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    branks[i] = calc_cuts(mbf->cdim, bgi->cuts);
  }
  
  // Get connected blocks
  int dir = mbf->cdim-1;
  int nconnected[num_blocks];
  int **block_list = gkyl_malloc(num_blocks*sizeof(int*));
  for (int bidx=0; bidx<num_blocks; ++bidx) {
    nconnected[bidx] = gkyl_multib_conn_get_num_connected(mbapp->block_topo, bidx, dir, 0, GKYL_CONN_ALL);
    block_list[bidx] = gkyl_malloc(nconnected[bidx]*sizeof(int));
    gkyl_multib_conn_get_connection(mbapp->block_topo, bidx, dir, 0, GKYL_CONN_ALL, block_list[bidx]);
  }

  // Construct the local and global ranges for the allgather
  mbf->multibz_ranges = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  mbf->multibz_ranges_ext = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));

  int nghost[] = { 1, 1 ,1};
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    mbf->multibz_ranges[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    mbf->multibz_ranges_ext[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    gkyl_multib_comm_conn_create_multib_ranges_in_dir(mbf->multibz_ranges_ext[bI],
        mbf->multibz_ranges[bI], nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);
  }

  // Allocate local and global arrays for charge density
  mbf->phi_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->phi_multibz_dg = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->phi_multibz_smooth = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));

  mbf->rho_c_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rho_c_multibz_dg = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rho_c_multibz_smooth = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));

  mbf->lhs_weight_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->lhs_weight_multibz = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rhs_weight_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rhs_weight_multibz = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));

  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    mbf->phi_local[bI] = gkyl_array_acquire(sbapp->field->phi_smooth);
    mbf->phi_multibz_dg[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);
    mbf->phi_multibz_smooth[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);

    mbf->rho_c_local[bI] = gkyl_array_acquire(sbapp->field->rho_c);
    mbf->rho_c_multibz_dg[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);
    mbf->rho_c_multibz_smooth[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);

    mbf->rhs_weight_local[bI] = mkarr(mbapp->use_gpu, sbapp->basis.num_basis, sbapp->local_ext.volume);
    gkyl_array_shiftc(mbf->rhs_weight_local[bI], sqrt(pow(2,mbf->cdim)), 0); // Sets weight=1.
    mbf->rhs_weight_multibz[bI] = mkarr(mbapp->use_gpu, 
      sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);
    if (mbf->cdim == 1) {
      mbf->lhs_weight_local[bI] = gkyl_array_acquire(sbapp->field->epsilon);
      mbf->lhs_weight_multibz[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);
    }
    else {
      mbf->lhs_weight_local[bI] = mkarr(mbapp->use_gpu, sbapp->basis.num_basis, sbapp->local_ext.volume);
      gkyl_array_shiftc(mbf->lhs_weight_local[bI], sqrt(pow(2,mbf->cdim)), 0); // Sets weight=1.
      mbf->lhs_weight_multibz[bI] = mkarr(mbapp->use_gpu, 
        sbapp->basis.num_basis, mbf->multibz_ranges_ext[bI]->volume);
    }
  }

  // Construct the comm_conns for the allgather
  mbf->mbcc_allgatherz_send = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbf->mbcc_allgatherz_recv = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbf->mbcc_allgatherz_send[bI] = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, 
      nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);
    mbf->mbcc_allgatherz_recv[bI] = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank,
      nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);

    for (int ns=0; ns<mbf->mbcc_allgatherz_send[bI]->num_comm_conn; ++ns) {
      // need to get the actual rank that owns this cut
      int rank_idx = mbf->mbcc_allgatherz_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, 
        mbf->mbcc_allgatherz_send[bI]->comm_conn[ns].block_id, rank_list);
      mbf->mbcc_allgatherz_send[bI]->comm_conn[ns].rank = rank_list[rank_idx];
      // Make range the local range (a subrange of local_ext)
      mbf->mbcc_allgatherz_send[bI]->comm_conn[ns].range = mbapp->singleb_apps[bI]->local;
    }
    for (int nr=0; nr<mbf->mbcc_allgatherz_recv[bI]->num_comm_conn; ++nr) {
      // need to get the actual rank that owns this cut
      int rank_idx = mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, 
        mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].block_id, rank_list);
      mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].rank = rank_list[rank_idx];
      // Make range a subrange
      gkyl_sub_range_init(&mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].range,
        mbf->multibz_ranges_ext[bI], mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].range.lower,
        mbf->mbcc_allgatherz_recv[bI]->comm_conn[nr].range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbf->mbcc_allgatherz_send[bI]);
    gkyl_multib_comm_conn_sort(mbf->mbcc_allgatherz_recv[bI]);
  }


  // Gather the weight along the magnetic field
  int stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, 
    mbapp->local_blocks, mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->lhs_weight_local, 
    mbf->lhs_weight_multibz);
  stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, 
    mbapp->local_blocks, mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->rhs_weight_local, 
    mbf->rhs_weight_multibz);

  mbf->fem_parproj = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_fem_parproj*));
  // Make the parrallel smoother
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = local_blocks[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    enum gkyl_fem_parproj_bc_type fem_parbc = mbf->info.duplicate_across_blocks? mbf->info.blocks[0].fem_parbc : mbf->info.blocks[bid].fem_parbc;
    mbf->fem_parproj[bI] = gkyl_fem_parproj_new(mbf->multibz_ranges[bI],
        &sbapp->basis, fem_parbc, mbf->lhs_weight_multibz[bI], mbf->rhs_weight_multibz[bI], mbapp->use_gpu);
  }
  
  // Set intersects for copying local rho back out after smoothing
  mbf->parent_subrangesz = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = local_blocks[bI];
    int shift[GKYL_MAX_DIM] = {0};
    for (int i=0; i<nconnected[bid]; i++) {
      if (block_list[bid][i] == bid)
        break;
      else 
        shift[dir] += gkyl_range_shape(&mbapp->decomp[block_list[bid][i]]->parent_range, dir);
    }
    struct gkyl_range shifted_parent_range;
    gkyl_range_shift(&shifted_parent_range, &mbapp->singleb_apps[bI]->global, shift);
    mbf->parent_subrangesz[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    int inter = gkyl_sub_range_intersect(mbf->parent_subrangesz[bI],
        mbf->multibz_ranges[bI], &shifted_parent_range);
  }

  // Last initialization step should be to set intersects for copying local info back out after smoothing
  mbf->block_subrangesz = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = local_blocks[bI];
    int shift[GKYL_MAX_DIM] = {0};
    for (int i=0; i<nconnected[bid]; i++) {
      if (block_list[bid][i] == bid)
        break;
      else 
        shift[dir] += gkyl_range_shape(&mbapp->decomp[block_list[bid][i]]->parent_range, dir);
    }
    struct gkyl_range shifted_block_range;
    gkyl_range_shift(&shifted_block_range, &mbapp->singleb_apps[bI]->local, shift);
    mbf->block_subrangesz[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    int inter = gkyl_sub_range_intersect(mbf->block_subrangesz[bI],
        mbf->multibz_ranges[bI], &shifted_block_range);
  }

  // Free temporary memory
  gkyl_free(branks);
  for (int bidx=0; bidx<num_blocks; ++bidx) {
    gkyl_free(block_list[bidx]);
  }
  gkyl_free(block_list);

  return mbf;
}

// Compute the electrostatic potential
void
gk_multib_field_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf, const struct gkyl_array *fin[])
{
  // Every local block calculates its charge density
  for (int bI=0; bI<mbf->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Construct fin for the local block
    const struct gkyl_array *fin_local_block[mbapp->num_species];
    int lin_idx = bI * mbapp->num_species;
    for (int i=0; i<mbapp->num_species; ++i) {
      fin_local_block[i] = fin[lin_idx+i];
    }
    // accumulate rho_c in local block
    gk_field_accumulate_rho_c(sbapp, sbapp->field, fin_local_block);
  }

  // Do the allgather of rho along the magnetic field.
  int stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->rho_c_local, mbf->rho_c_multibz_dg);
  // Do the rho smoothing on the multibz range
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    gkyl_fem_parproj_set_rhs(mbf->fem_parproj[bI], mbf->rho_c_multibz_dg[bI], mbf->rho_c_multibz_dg[bI]);
    gkyl_fem_parproj_solve(mbf->fem_parproj[bI], mbf->rho_c_multibz_smooth[bI]);
  }
  if (mbf->cdim == 1) {
    for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
      gkyl_array_copy_range_to_range(mbapp->singleb_apps[bI]->field->phi_smooth, mbf->rho_c_multibz_smooth[bI],
        &mbapp->singleb_apps[bI]->local, mbf->block_subrangesz[bI]);
    }
  }
  else {
    // Copy smooth arrays  back to apps
    for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
      gkyl_array_copy_range_to_range(mbapp->singleb_apps[bI]->field->rho_c_global_smooth,
        mbf->rho_c_multibz_smooth[bI], &mbapp->singleb_apps[bI]->global, mbf->parent_subrangesz[bI]);
    }
    // Every block Solves the poisson equation
    for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
      gkyl_deflated_fem_poisson_advance(sbapp->field->deflated_fem_poisson, 
        sbapp->field->rho_c_global_smooth, 0, sbapp->field->phi_smooth);
    }
  }

}

// Release resources for multib field
void
gk_multib_field_release(struct gk_multib_field *mbf)
{

  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    gkyl_free(mbf->multibz_ranges[bI]);
    gkyl_free(mbf->multibz_ranges_ext[bI]);
    gkyl_free(mbf->parent_subrangesz[bI]);
    gkyl_free(mbf->block_subrangesz[bI]);
    gkyl_array_release(mbf->phi_local[bI]);
    gkyl_array_release(mbf->phi_multibz_dg[bI]);
    gkyl_array_release(mbf->phi_multibz_smooth[bI]);
    gkyl_array_release(mbf->rho_c_local[bI]);
    gkyl_array_release(mbf->rho_c_multibz_dg[bI]);
    gkyl_array_release(mbf->rho_c_multibz_smooth[bI]);
    gkyl_array_release(mbf->lhs_weight_local[bI]);
    gkyl_array_release(mbf->lhs_weight_multibz[bI]);
    gkyl_array_release(mbf->rhs_weight_local[bI]);
    gkyl_array_release(mbf->rhs_weight_multibz[bI]);
    gkyl_multib_comm_conn_release(mbf->mbcc_allgatherz_send[bI]);
    gkyl_multib_comm_conn_release(mbf->mbcc_allgatherz_recv[bI]);
    gkyl_fem_parproj_release(mbf->fem_parproj[bI]);
  }

  gkyl_free(mbf->multibz_ranges);
  gkyl_free(mbf->multibz_ranges_ext);
  gkyl_free(mbf->parent_subrangesz);
  gkyl_free(mbf->block_subrangesz);
  gkyl_free(mbf->phi_local);
  gkyl_free(mbf->phi_multibz_dg);
  gkyl_free(mbf->phi_multibz_smooth);
  gkyl_free(mbf->rho_c_local);
  gkyl_free(mbf->rho_c_multibz_dg);
  gkyl_free(mbf->rho_c_multibz_smooth);
  gkyl_free(mbf->lhs_weight_local);
  gkyl_free(mbf->lhs_weight_multibz);
  gkyl_free(mbf->rhs_weight_local);
  gkyl_free(mbf->rhs_weight_multibz);
  gkyl_free(mbf->mbcc_allgatherz_send);
  gkyl_free(mbf->mbcc_allgatherz_recv);
  gkyl_free(mbf->fem_parproj);

  gkyl_free(mbf);
}
