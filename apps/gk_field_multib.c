#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_multib_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static int count_distinct(int a[], int n)      //Function Definition
{
   int i, j, count = 1;
   //Traverse the array
   for (i = 1; i < n; i++)      //hold an array element
   {
      for (j = 0; j < i; j++)   
      {
         if (a[i] == a[j])    //Check for duplicate elements
         {
            break;             //If duplicate elements found then break
         }
      }
      if (i == j)
      {
         count++;     //increment the number of distinct elements
      }
   }
   return count;      //Return the number of distinct elements
}

static int get_unique(int *input_array, int n, int *unique_array) {
    // Initialize an array to keep track of unique elements
    int isUnique[n];  // Using VLA (Variable Length Array) for simplicity
    
    // Initialize all elements as unique (false)
    for (int i = 0; i < n; i++) {
        isUnique[i] = 1;
    }
    
    // Mark duplicates as non-unique
    for (int i = 0; i < n; i++) {
        if (isUnique[i]) {
            for (int j = i + 1; j < n; j++) {
                if (input_array[i] == input_array[j]) {
                    isUnique[j] = 0;  // Mark as non-unique
                }
            }
        }
    }
    
    // Copy unique elements to unique_array and count them
    int uniqueCount = 0;
    for (int i = 0; i < n; i++) {
        if (isUnique[i]) {
            unique_array[uniqueCount++] = input_array[i];
        }
    }
    
    return uniqueCount;
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
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];

    mbf->info = inp->field;

    mbf->gkfield_id = mbf->info.gkfield_id ? mbf->info.gkfield_id : GKYL_GK_FIELD_ES;

    int num_blocks = 0;
    for (int i =0; i < GKYL_MAX_BLOCKS; i++) mbf->crossz_block_idxs[i] = -1;
    set_crossz_idxs(mba, bidx, mbf->crossz_block_idxs, &num_blocks);
    mbf->numz_blocks = num_blocks;

    // We need first to create the cross z ranges.
    // For the two-block case:
    int lower[2];
    int upper[2];
    lower[0] = mba->decomp_intrab[bc]->parent_range.lower[0];
    upper[0] = mba->decomp_intrab[bc]->parent_range.upper[0];
    lower[1] = mba->decomp_intrab[bc]->parent_range.lower[1];
    // Add up all ranges in z direction
    upper[1] = 0;
    for ( int i = 0; i < num_blocks; i++) {
      upper[1] += gkyl_range_shape(&mba->decomp_intrab[mbf->crossz_block_idxs[i]]->parent_range, 1);
    }
    struct gkyl_range crossz;
    gkyl_range_init(&crossz, mba->cdim, lower, upper);
    int nghost[2] = {1,1};
    gkyl_create_ranges(&crossz, nghost, &mbf->crossz_ext, &mbf->crossz);


    printf("crossz.lower = %d %d\n", crossz.lower[0], crossz.lower[1]);
    printf("crossz.upper = %d %d\n", crossz.upper[0], crossz.upper[1]);

    // Create the decomp and communicator from the mb app communicator
    // Stack all the individual ranges together
    int num_ranges = 0;
    for ( int i = 0; i < num_blocks; i++) {
      num_ranges += mba->decomp_intrab[mbf->crossz_block_idxs[i]]->ndecomp;
    }
    mbf->num_cuts = num_ranges;
    mbf->cut_ranges = gkyl_malloc(sizeof(struct gkyl_range[num_ranges]));
    int range_count = 0;
    int num_zcells_before = 0;
    for (int ib = 0; ib < num_blocks; ib++) {
      struct gkyl_rect_decomp *decomp_i = mba->decomp_intrab[ib];
      for (int id = 0; id < decomp_i->ndecomp; id++) {
        lower[0] = decomp_i->ranges[id].lower[0];
        upper[0] = decomp_i->ranges[id].upper[0];
        lower[1] = decomp_i->ranges[id].lower[1] + num_zcells_before;
        upper[1] = decomp_i->ranges[id].upper[1] + num_zcells_before;
        gkyl_range_init(&mbf->cut_ranges[range_count], mba->cdim, lower, upper);

        printf("cut_ranges[%d].lower = %d %d\n", range_count, mbf->cut_ranges[range_count].lower[0], mbf->cut_ranges[range_count].lower[1]);
        printf("cut_ranges[%d].upper = %d %d\n", range_count, mbf->cut_ranges[range_count].upper[0], mbf->cut_ranges[range_count].upper[1]);
        range_count+=1;
      }
      num_zcells_before += gkyl_range_shape(&decomp_i->parent_range, 1);
    }

    int zranks[num_ranges];
    for (int i =0; i < num_ranges; i++) zranks[i] = -1;
    // Expected function: get_ranks(bidx, int *ranks);
    //  int
    //gyrokinetic_multib_ranks_per_block(gkyl_gyrokinetic_multib_app *mba, int bidx, int *ranks)
    //{
    //...
    //return num_ranks;
    //}
    // Populates ranks with ranks in each block
    int num_populated = 0;
    for (int i = 0; i < num_blocks; i++) {
      //printf("rank before = %d\n", zranks[num_populated]);
      int num_ranks = gkyl_gyrokinetic_multib_ranks_per_block(mba, mbf->crossz_block_idxs[i], &zranks[num_populated]);
      //printf("bidx = %d | rank after = %d\n", mbf->crossz_block_idxs[i], zranks[num_populated]);
      //printf("bidx = %d | num ranks = %d\n", bidx, num_ranks);
      //num_populated += mba->decomp_intrab[mbf->crossz_block_idxs[i]]->ndecomp; Should be same as line below
      num_populated += num_ranks;
    }

    printf("num ranges = %d\n", num_ranges);
    printf("num cuts= %d\n", mbf->num_cuts);
    printf(" the cut ranks are ");
    for(int i = 0; i < num_ranges; i++) {
      printf(" %d", zranks[i]);
    }
    printf("\n");

    int num_unique_ranks = count_distinct(zranks, num_ranges);
    int unique_zranks[num_unique_ranks];
    int num_unique_ranks2 = get_unique(zranks, num_ranges, unique_zranks);
    printf("unique ranks = %d\n", num_unique_ranks);
    //printf("unique ranks2 = %d\n", num_unique_ranks2);
    printf("the uq ranks are ");
    for(int i = 0; i < num_unique_ranks; i++){
      printf(" %d", unique_zranks[i]);
    }
    printf("\n");




    //mbf->zcomm = gkyl_comm_split_comm(mba->comm_multib, 0, mbf->zdecomp);
    //mbf->zcomm = gkyl_comm_split_comm(mba->comm_multib, 0, 0);
    //mbf->zcomm = gkyl_comm_create_group_comm(comm, num_ranks_group, ranks_group, tag, decomp_group);
    int tag = 0; // I do not know what the tag should be
    mbf->zcomm = gkyl_comm_create_group_comm(mba->comm_multib, num_unique_ranks, unique_zranks, tag, 0);

    mbf->zranks = gkyl_malloc(num_ranges*sizeof(int));
    gkyl_comm_group_translate_ranks(mba->comm_multib, num_ranges, zranks,  mbf->zcomm, mbf->zranks);

    mbf->unique_zranks = gkyl_malloc(num_unique_ranks*sizeof(int));
    gkyl_comm_group_translate_ranks(mba->comm_multib, num_unique_ranks, unique_zranks,  mbf->zcomm, mbf->unique_zranks);

    // allocate arrays for charge density
    mbf->rho_c_global_dg = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);
    mbf->rho_c_global_smooth = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);

    // allocate arrays for electrostatic potential
    mbf->phi = mkarr(mba->use_gpu, app->confBasis.num_basis, mbf->crossz_ext.volume);

    mbf->phi_host = mbf->phi;
    if (mba->use_gpu) {
      mbf->phi_host = mkarr(false, app->confBasis.num_basis, mbf->crossz_ext.volume);
    }


    mbf->fem_parproj = gkyl_fem_parproj_new(&mbf->crossz, &mbf->crossz_ext, 
      &app->confBasis, mbf->info.fem_parbc, NULL, mba->use_gpu);
  
  // Now get the sub range intersects
  // Create global subrange we'll copy the field solver solution from (into local).
  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    int intersect = gkyl_sub_range_intersect(&mbf->crossz_sub_range[bc], &mbf->crossz, &app->global);
  }


  return mbf;
}

// Compute the electrostatic potential
void
gk_field_multib_rhs(gkyl_gyrokinetic_multib_app *mba, struct gk_field_multib *mbf)
{

  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];

    // Get fin and the field for app we own
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    struct gk_field* field = app->field ;

    // accumulate rho_c in each block
    gk_field_accumulate_rho_c(app, field, fin);
    // Now gather charge density into the interblock cross-z array for smoothing in z
  }

  // Now we need to do two things instead of the allgather.
  // 1. for each local block, copy field->rho_c into the correct part of mbf->rho_c_global_dg
  // 2. do the broadcast which will loop over num_cuts
   
  // Do 1.
  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    struct gk_field* field = app->field;
    gkyl_array_copy_range_to_range(mbf->rho_c_global_dg, field->rho_c, &mbf->crossz_sub_range[bc], &app->global);
  }

  int lower[2];
  int upper[2];
  int bc = 0; // assume all are connected along z
  lower[0] = mba->decomp_intrab[bc]->parent_range.lower[0];
  upper[0] = mba->decomp_intrab[bc]->parent_range.upper[0];
  
  // Do 2.
  for (int ic = 0; ic < mbf->num_cuts; ic++) {
    //gkyl_comm_bcast(struct gkyl_comm *comm, void *data, size_t data_sz, int root)
    struct gkyl_range range_curr = mbf->cut_ranges[ic];
    int idx[2] = {range_curr.lower[0], range_curr.lower[1]};
    //printf("range_curr.lower = %d %d\n", range_curr.lower[0], range_curr.lower[1]);
    //printf("range_curr.upper = %d %d\n", range_curr.upper[0], range_curr.upper[1]);
    long start_loc = gkyl_range_idx(&mbf->crossz, idx);
    double *data_start = gkyl_array_fetch(mbf->rho_c_global_dg, start_loc);
    //printf("ic = %d, start_loc = %ld\n", ic, start_loc);
    //printf(" size = %ld\n", range_curr.volume);
    //printf(" root rank - %d\n", mbf->zranks[ic]);
    gkyl_comm_bcast(mbf->zcomm, data_start, range_curr.volume, mbf->zranks[ic]);
  }
  //int rank = 0;
  //gkyl_comm_get_rank(mba->comm_multib, &rank);
  //printf("my rank is %d. Did the broadcast\n", rank);


  // Do the smoothing on the inetrblock cross-z range
  gkyl_fem_parproj_set_rhs(mbf->fem_parproj, mbf->rho_c_global_dg, mbf->rho_c_global_dg);
  gkyl_fem_parproj_solve(mbf->fem_parproj, mbf->rho_c_global_smooth);
  //printf("my rank is %d. Did the smoothing\n", rank);

  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bc];
    struct gk_field* field = app->field ;
    // Copy inter-block cross-z smoothed charge density to intrablock global charge density per process
    gkyl_array_copy_range_to_range(field->rho_c_global_smooth, mbf->rho_c_global_smooth, &app->global, &mbf->crossz_sub_range[bc]);
    // Now call the perp solver. The perp solver already accesses its own local part of the intrablock global range.
    gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth, field->phi_smooth);
  }
  //printf("my rank is %d. Did the solve\n", rank);
  
}

// release resources for field
void
gk_field_multib_release(const gkyl_gyrokinetic_multib_app* mba, struct gk_field_multib *mbf)
{
  gkyl_array_release(mbf->rho_c_global_dg);
  gkyl_array_release(mbf->rho_c_global_smooth);
  gkyl_array_release(mbf->phi);
  gkyl_free(mbf->cut_ranges);

  gkyl_comm_release(mbf->zcomm);

  gkyl_fem_parproj_release(mbf->fem_parproj);

  if (mba->use_gpu) {
    gkyl_array_release(mbf->phi_host);
  }

  gkyl_free(mbf);
}

