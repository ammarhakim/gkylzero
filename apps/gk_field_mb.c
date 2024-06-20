#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_mb_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

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
set_crossz_idxs(struct gkyl_gyrokinetic_mb_app *mb_app, struct gkyl_gyrokinetic_app *app, int myidx){
  struct gkyl_block_topo *btopo = mb_app->btopo;
  struct gkyl_block_connections *conn = btopo->conn;
  int dir = 1;

  int crossz_blocks[GKYL_MAX_BLOCKS] = {-1};

  int bidx = myidx;
  crossz_blocks[0] = bidx;
  int num_blocks = 1;
  printf("my idx = %d\n", bidx);

  while(true) {
    if (conn[bidx].connections[dir][0].edge == GKYL_BLOCK_EDGE_PHYSICAL) {
      break;
    }
    else if (conn[bidx].connections[dir][0].edge == GKYL_BLOCK_EDGE_UPPER_POSITIVE) { 
      insert_below(crossz_blocks, &num_blocks, conn[bidx].connections[dir][0].bid);
      bidx = conn[bidx].connections[dir][0].bid;
      printf(" Looking below, num_blocks = %d\n", num_blocks);
    }
  }

  bidx = myidx;
  while(true) {
    if (conn[bidx].connections[dir][1].edge == GKYL_BLOCK_EDGE_PHYSICAL) {
      break;
    }
    else if (conn[bidx].connections[dir][1].edge == GKYL_BLOCK_EDGE_LOWER_POSITIVE) { 
      insert_above(crossz_blocks, &num_blocks, conn[bidx].connections[dir][1].bid);
      bidx = conn[bidx].connections[dir][1].bid;
      printf(" Looking above, num_blocks = %d\n", num_blocks);
    }
  }

  printf("num blocks = %d\n block order = ", num_blocks);
  for (int i = 0; i < num_blocks; i++) {
    printf(" %d ", crossz_blocks[i]);
  }
  printf("\n");
}

// initialize field object
struct gk_field_mb* 
gk_field_mb_new(struct gkyl_gk_mb *gk_mb, struct gkyl_gyrokinetic_mb_app *mb_app, struct gkyl_gyrokinetic_app *app, int bidx)
{

  struct gk_field_mb *f = gkyl_malloc(sizeof(struct gk_field_mb));

  f->info = gk_mb->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  set_crossz_idxs(mb_app, app, bidx);

  // We need first to create the global z ranges.
  // For the two-block case:
  int lower[2];
  int upper[2];
  lower[0] = mb_app->decomp_intrab[0]->parent_range.lower[0];
  upper[0] = mb_app->decomp_intrab[0]->parent_range.upper[0];
  lower[1] = mb_app->decomp_intrab[0]->parent_range.lower[1];
  upper[1] = gkyl_range_shape(&mb_app->decomp_intrab[0]->parent_range, 1) + gkyl_range_shape(&mb_app->decomp_intrab[1]->parent_range, 1); // Need to add up all ranges connected in z. So would need one more [2] in for 3 blocks
  struct gkyl_range crossz;
  gkyl_range_init(&crossz, mb_app->cdim, lower, upper);
  int nghost[2] = {1,1};
  gkyl_create_ranges(&crossz, nghost, &f->crossz_ext, &f->crossz);

  // Create the decomp and communicator from the mb app communicator
  int cuts[2] = {1,2}; // For only two blocks. Would be 1,3 for 3 blocks along z
  struct gkyl_rect_decomp *zdecomp = gkyl_rect_decomp_new_from_cuts(mb_app->cdim, cuts, &f->crossz);
  f->zcomm = gkyl_comm_split_comm(mb_app->comm_mb, 0, zdecomp); // Would have different colors for other 
                                                               // block groups like 11-12, 7-8-9
                                                               // but just 0 for now


  // Now get the sub range intersects
  // Create global subrange we'll copy the field solver solution from (into local).
  int rank;
  gkyl_comm_get_rank(f->zcomm, &rank);
  int intersect = gkyl_sub_range_intersect(&f->crossz_sub_range, &f->crossz, &app->global);

  // allocate arrays for charge density
  f->rho_c_global_dg = mkarr(mb_app->use_gpu, app->confBasis.num_basis, f->crossz_ext.volume);
  f->rho_c_global_smooth = mkarr(mb_app->use_gpu, app->confBasis.num_basis, f->crossz_ext.volume);

  // allocate arrays for electrostatic potential
  f->phi = mkarr(mb_app->use_gpu, app->confBasis.num_basis, f->crossz_ext.volume);

  f->phi_host = f->phi;
  if (mb_app->use_gpu) {
    f->phi_host = mkarr(false, app->confBasis.num_basis, f->crossz_ext.volume);
  }


  f->fem_parproj = gkyl_fem_parproj_new(&f->crossz, &f->crossz_ext, 
    &app->confBasis, f->info.fem_parbc, NULL, mb_app->use_gpu);


  return f;
}

// Compute the electrostatic potential
void
gk_field_mb_rhs(gkyl_gyrokinetic_mb_app *mb_app, struct gk_field_mb *field_mb, struct gkyl_gyrokinetic_app *app)
{
    // Get rank and figure out which app is ours
    //int rank;
    //gkyl_comm_get_rank(field_mb->zcomm, &rank);
    //struct gkyl_gyrokinetic_app* app = mb_app->blocks[rank];

    // Get fin and the field for app we own
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    struct gk_field* field = app->field ;

    // accumulate rho_c in each block
    gk_field_accumulate_rho_c(app, field, fin);
    // Each block can gather its own charge density.
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);

    // Now gather charge density into global interblock array for smoothing in z
    gkyl_comm_array_allgather(field_mb->zcomm, &app->global, &field_mb->crossz, field->rho_c_global_dg, field_mb->rho_c_global_dg);
    // Do the smoothing on the inetrblock global z range
    gkyl_fem_parproj_set_rhs(field_mb->fem_parproj, field_mb->rho_c_global_dg, field_mb->rho_c_global_dg);
    gkyl_fem_parproj_solve(field_mb->fem_parproj, field_mb->rho_c_global_smooth);
    // Copy inter-block globally (in z) smoothed charge density to intrablock global charge density per process
    gkyl_array_copy_range_to_range(field->rho_c_global_smooth, field_mb->rho_c_global_smooth, &app->global, &field_mb->crossz_sub_range);
    // Now call the perp solver. The perp solver already accesses its own local part of the intra-block global range.
    gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth, field->phi_smooth);
}

// release resources for field
void
gk_field_mb_release(const gkyl_gyrokinetic_mb_app* mb_app, struct gk_field_mb *f)
{
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi);

  gkyl_fem_parproj_release(f->fem_parproj);

  if (mb_app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_free(f);
}

