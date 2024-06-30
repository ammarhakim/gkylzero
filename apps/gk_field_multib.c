#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_multib_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

/**
 * Count number of distinct elements in an array of ints
 * .
 * @param a input array of ints
 * @param n length of a
 * return number of unique elements in a
 */
static int count_distinct(int a[], int n)
{
   int i, j, count = 1;
   for (i = 1; i < n; i++) { // Check if a[i] is a new element
     for (j = 0; j < i; j++) {
       if (a[i] == a[j])    // Check if a[i] has already been found 
          break;            // Break if it is a duplicate
     }
     if (i == j)
       count++;     //increment the number of distinct elements
   }
   return count;
}

/**
 * Populate an output array with the distinct elements in an array of ints
 * .
 * @param a input array of ints
 * @param n length of input array
 * @param n_unique number of unique elements in a
 * @param unique_array on output contains the unique elements in a
 * return number of unique elements in a
 */
static int get_unique(int *a, int n, int *unique_array) {
   unique_array[0] = a[0]; // The first element of a is the first unique element
   int i, j, count = 1;
   for (i = 1; i < n; i++) { // Check if a[i] is a new element
     for (j = 0; j < i; j++) {
       if (a[i] == a[j])    // Check if a[i] has already been found 
          break;            // Break if it is a duplicate
     }
     if (i == j) {
       count++;     //increment the number of distinct elements
       unique_array[i] = a[i];
     }
   }
   return count;
}

static void
insert_below(int* arr, int *n, int new_val)
{
  int temp_arr[GKYL_MAX_BLOCKS] = {-1};
  for (int i = 0; i<*n; i++) {
    temp_arr[i+1] = arr[i];
  }
  temp_arr[0] = new_val;
  *n = *n+1;
  for (int i = 0; i<*n; i++) {
    arr[i] = temp_arr[i];
  }

}

static void
insert_above(int* arr, int *n, int new_val)
{
  arr[*n] = new_val;
  *n = *n+1;
}

// This function should be called in a loop over num blocks local
void
set_crossz_idxs(struct gkyl_gyrokinetic_multib_app *mba, int myidx, int* crossz_blocks, int* num_blocks){
  struct gkyl_block_topo *btopo = mba->btopo;
  struct gkyl_block_connections *conn = btopo->conn;
  int dir = 1;

  int bidx = myidx;
  crossz_blocks[0] = bidx;
  *num_blocks = 1;

  while(true) {
    if (conn[bidx].connections[dir][0].edge == GKYL_BLOCK_EDGE_PHYSICAL) {
      break;
    }
    else if (conn[bidx].connections[dir][0].edge == GKYL_BLOCK_EDGE_UPPER_POSITIVE) { 
      insert_below(crossz_blocks, num_blocks, conn[bidx].connections[dir][0].bid);
      bidx = conn[bidx].connections[dir][0].bid;
    }
  }

  bidx = myidx;
  while(true) {
    if (conn[bidx].connections[dir][1].edge == GKYL_BLOCK_EDGE_PHYSICAL) {
      break;
    }
    else if (conn[bidx].connections[dir][1].edge == GKYL_BLOCK_EDGE_LOWER_POSITIVE) { 
      insert_above(crossz_blocks, num_blocks, conn[bidx].connections[dir][1].bid);
      bidx = conn[bidx].connections[dir][1].bid;
    }
  }

}

// initialize field object
struct gk_field_multib* 
gk_field_multib_new(struct gkyl_gk_multib *inp, struct gkyl_gyrokinetic_multib_app *mba)
{

  struct gk_field_multib *mbf = gkyl_malloc(sizeof(struct gk_field_multib));

  int bc = 0; // Set block counter to 0. Assume blocks are connected along z
  int bidx = mba->block_idxs[bc];

  mbf->info = inp->field;

  mbf->gkfield_id = mbf->info.gkfield_id ? mbf->info.gkfield_id : GKYL_GK_FIELD_ES;

  int num_blocks = 0;
  for (int i =0; i < GKYL_MAX_BLOCKS; i++) mbf->crossz_block_idxs[i] = -1;
  set_crossz_idxs(mba, bidx, mbf->crossz_block_idxs, &num_blocks);
  mbf->numz_blocks = num_blocks;

  // We need first to create the cross z ranges.
  // For the two-block case:
  int lower[mba->cdim], lower_off[mba->cdim];
  int upper[mba->cdim];
  for (int d=0; d<mba->cdim-1; d++) {
    lower[d] = mba->decomp_intrab[bc]->parent_range.lower[d];
    upper[d] = mba->decomp_intrab[bc]->parent_range.upper[d];
    lower_off[d] = lower[d];
  }
  lower[mba->cdim-1] = mba->decomp_intrab[bc]->parent_range.lower[1];
  // Add up all ranges in z direction
  upper[mba->cdim-1] = 0;
  lower_off[mba->cdim-1] = 1;
  for (int i=0; i<num_blocks; i++) {
    upper[mba->cdim-1] += gkyl_range_shape(&mba->decomp_intrab[mbf->crossz_block_idxs[i]]->parent_range, mba->cdim-1);
    gkyl_range_init(&mbf->parent_range_off[mbf->crossz_block_idxs[i]], mba->cdim, lower_off, upper);
    lower_off[mba->cdim-1] += gkyl_range_shape(&mba->decomp_intrab[mbf->crossz_block_idxs[i]]->parent_range, mba->cdim-1);
  }

  struct gkyl_range crossz;
  gkyl_range_init(&crossz, mba->cdim, lower, upper);
  int nghost[mba->cdim];
  for (int d=0; d<mba->cdim; d++) nghost[d] = 1;
  gkyl_create_ranges(&crossz, nghost, &mbf->crossz_ext, &mbf->crossz);

  // Create the decomp and communicator from the mb app communicator
  // Stack all the individual ranges together
  int num_ranges_z = 0;
  for (int i=0; i<num_blocks; i++)
    num_ranges_z += mba->decomp_intrab[mbf->crossz_block_idxs[i]]->ndecomp;
  mbf->num_cuts = num_ranges_z;
  mbf->cut_ranges = gkyl_malloc(num_ranges_z * sizeof(struct gkyl_range));
  int range_count = 0;
  int num_zcells_cum = 0;
  for (int ib=0; ib<num_blocks; ib++) {
    struct gkyl_rect_decomp *decomp_i = mba->decomp_intrab[ib];
    for (int id=0; id<decomp_i->ndecomp; id++) {
      for (int d=0; d<mba->cdim-1; d++) {
        lower[d] = decomp_i->ranges[id].lower[d];
        upper[d] = decomp_i->ranges[id].upper[d];
      }
      lower[mba->cdim-1] = decomp_i->ranges[id].lower[mba->cdim-1] + num_zcells_cum;
      upper[mba->cdim-1] = decomp_i->ranges[id].upper[mba->cdim-1] + num_zcells_cum;
      gkyl_range_init(&mbf->cut_ranges[range_count], mba->cdim, lower, upper);
      range_count += 1;
    }
    num_zcells_cum += gkyl_range_shape(&decomp_i->parent_range, 1);
  }

  // Populate with ranks in each block
  int zranks_mb[num_ranges_z];
  int num_ranks_cum = 0;
  for (int i=0; i<num_blocks; i++) {
    int num_ranks = gkyl_gyrokinetic_multib_ranks_per_block(mba, mbf->crossz_block_idxs[i], &zranks_mb[num_ranks_cum]);
    num_ranks_cum += num_ranks;
  }

  int num_unique_ranks = count_distinct(zranks_mb, num_ranges_z);
  int unique_zranks_mb[num_unique_ranks];
  int num_unique_ranks2 = get_unique(zranks_mb, num_ranges_z, unique_zranks_mb);
  assert(num_unique_ranks == num_unique_ranks2);

  mbf->zcomm = gkyl_comm_create_comm(mba->comm_multib, num_unique_ranks, unique_zranks_mb, 0);

  // allocate arrays for charge density
  struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
  mbf->rho_c_global_dg = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);
  mbf->rho_c_global_smooth = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);

  // allocate arrays for electrostatic potential
  mbf->phi = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);

  mbf->phi_ho = mba->use_gpu? mkarr(false, app->confBasis.num_basis, mbf->crossz_ext.volume)
                            : gkyl_array_acquire(mbf->phi);

  mbf->fem_parproj = gkyl_fem_parproj_new(&mbf->crossz, &mbf->crossz_ext, 
    &app->confBasis, mbf->info.fem_parbc, NULL, mba->use_gpu);

  mbf->zranks = gkyl_malloc(num_ranges_z*sizeof(int));
  gkyl_comm_group_translate_ranks(mba->comm_multib, num_ranges_z, zranks_mb,  mbf->zcomm, mbf->zranks);

  mbf->unique_zranks = gkyl_malloc(num_unique_ranks*sizeof(int));
  gkyl_comm_group_translate_ranks(mba->comm_multib, num_unique_ranks, unique_zranks_mb,  mbf->zcomm, mbf->unique_zranks);
  
  // Now get the sub range intersects
  // Create global subrange we'll copy the field solver solution from (into local).
  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    int bidx = mba->block_idxs[bc];
    int intersect = gkyl_sub_range_intersect(&mbf->crossz_sub_range[bc], &mbf->crossz, &mbf->parent_range_off[bidx]);
  }

  return mbf;
}

// Compute the electrostatic potential
void
gk_field_multib_rhs(gkyl_gyrokinetic_multib_app *mba, struct gk_field_multib *mbf, int fidx)
{

  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];

    // Get fin and the field for app we own
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].distfs[fidx];
    }
    struct gk_field* field = app->field;

    // accumulate rho_c in each block
    gk_field_accumulate_rho_c(app, field, fin);
  }

  // Now we need to do two things instead of the allgather.
  // 1. Do an allgather within each block to fill and then copy the result in to the right part of
  // the crossz global rho
  // 2. Do the broadcast which will loop over num_cuts
   
  // Do 1.
  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    struct gk_field *field = app->field;
    // Each block can gather its own charge density.
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);
    // Copy into correct location of crossz rho
    gkyl_array_copy_range_to_range(mbf->rho_c_global_dg, field->rho_c_global_dg, &mbf->crossz_sub_range[bc], &app->global);
  }

  // Do 2.
  for (int ic=0; ic<mbf->num_cuts; ic++) {
    struct gkyl_range range_curr = mbf->cut_ranges[ic];
    int start_idx[mba->cdim];
    for (int d=0; d<mba->cdim; d++) start_idx[d] = range_curr.lower[d];
      int xstart_idx[mba->cdim-1];
      xstart_idx[mba->cdim-1] = start_idx[mba->cdim-1];
    for(int ix = range_curr.lower[0]; ix<=range_curr.upper[0]; ix++) {
      xstart_idx[0] = ix;
      long start_loc = gkyl_range_idx(&mbf->crossz_ext, xstart_idx);
      long send_size = range_curr.volume/gkyl_range_shape(&range_curr, 0)*mbf->rho_c_global_dg->ncomp;
      double *data_start = gkyl_array_fetch(mbf->rho_c_global_dg, start_loc);
      gkyl_comm_bcast(mbf->zcomm, data_start, send_size, mbf->zranks[ic]);
    }
  }

  // Do the smoothing on the inetrblock cross-z range
  gkyl_fem_parproj_set_rhs(mbf->fem_parproj, mbf->rho_c_global_dg, mbf->rho_c_global_dg);
  gkyl_fem_parproj_solve(mbf->fem_parproj, mbf->rho_c_global_smooth);

  int bidx = mba->block_idxs[0];

  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    struct gk_field* field = app->field ;
    // Copy inter-block cross-z smoothed charge density to intrablock global charge density per process
    gkyl_array_copy_range_to_range(field->rho_c_global_smooth, mbf->rho_c_global_smooth, &app->global, &mbf->crossz_sub_range[bc]);
    // Now call the perp solver. The perp solver already accesses its own local part of the intrablock global range.
    gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth, field->phi_smooth);
  }
  
}

// release resources for field
void
gk_field_multib_release(const gkyl_gyrokinetic_multib_app* mba, struct gk_field_multib *mbf)
{
  gkyl_free(mbf->cut_ranges);
  gkyl_free(mbf->zranks);
  gkyl_free(mbf->unique_zranks);
  gkyl_array_release(mbf->rho_c_global_dg);
  gkyl_array_release(mbf->rho_c_global_smooth);
  gkyl_array_release(mbf->phi);
  gkyl_array_release(mbf->phi_ho);

  gkyl_comm_release(mbf->zcomm);

  gkyl_fem_parproj_release(mbf->fem_parproj);

  gkyl_free(mbf);
}

