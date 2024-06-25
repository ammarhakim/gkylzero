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
set_crossz_idxs(struct gkyl_gyrokinetic_mb_app *mb_app, struct gkyl_gyrokinetic_app *app, int myidx, int* crossz_blocks, int* num_blocks){
  struct gkyl_block_topo *btopo = mb_app->btopo;
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
      printf(" Looking below, num_blocks = %d\n", *num_blocks);
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
      printf(" Looking above, num_blocks = %d\n", *num_blocks);
    }
  }

}

// initialize field object
struct gk_field_mb* 
gk_field_mb_new(struct gkyl_gk_mb *gk_mb, struct gkyl_gyrokinetic_mb_app *mb_app, struct gkyl_gyrokinetic_app *app, int bidx)
{

  struct gk_field_mb *f = gkyl_malloc(sizeof(struct gk_field_mb));

  f->info = gk_mb->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  int num_blocks = 0;
  int crossz_blocks[GKYL_MAX_BLOCKS] = {-1};
  set_crossz_idxs(mb_app, app, bidx, crossz_blocks, &num_blocks);

  printf("my idx = %d\n", bidx);
  printf("num blocks = %d\n block order = ", num_blocks);
  for (int i = 0; i < num_blocks; i++) {
    printf(" %d ", crossz_blocks[i]);
  }
  printf("\n");

  // We need first to create the cross z ranges.
  // For the two-block case:
  int lower[2];
  int upper[2];
  lower[0] = mb_app->decomp_intrab[0]->parent_range.lower[0];
  upper[0] = mb_app->decomp_intrab[0]->parent_range.upper[0];
  lower[1] = mb_app->decomp_intrab[0]->parent_range.lower[1];
  // Add up all ranges in z direction
  upper[1] = 0;
  for ( int i = 0; i < num_blocks; i++) {
    upper[1] += gkyl_range_shape(&mb_app->decomp_intrab[crossz_blocks[i]]->parent_range, 1);
  }
  printf("upper[1] = %d\n", upper[1]);
  struct gkyl_range crossz;
  gkyl_range_init(&crossz, mb_app->cdim, lower, upper);
  int nghost[2] = {1,1};
  gkyl_create_ranges(&crossz, nghost, &f->crossz_ext, &f->crossz);

  // Create the decomp and communicator from the mb app communicator
  // Stack all the individual ranges together
  int num_ranges = 0;
  for ( int i = 0; i < num_blocks; i++) {
    num_ranges += mb_app->decomp_intrab[crossz_blocks[i]]->ndecomp;
  }
  printf("num ranges = %d\n", num_ranges);
  struct gkyl_range ranges[num_ranges];
  int range_count = 0;
  int num_zcells_before = 0;
  for (int ib = 0; ib < num_blocks; ib++) {
    //printf("ib = %d\n", ib);
    struct gkyl_rect_decomp *decomp_i = mb_app->decomp_intrab[ib];
    for (int id = 0; id < decomp_i->ndecomp; id++) {
      lower[0] = decomp_i->ranges[id].lower[0];
      upper[0] = decomp_i->ranges[id].upper[0];
      lower[1] = decomp_i->ranges[id].lower[1] + num_zcells_before;
      upper[1] = decomp_i->ranges[id].upper[1] + num_zcells_before;
      gkyl_range_init(&ranges[range_count], mb_app->cdim, lower, upper);
      range_count+=1;
    }
    num_zcells_before += gkyl_range_shape(&decomp_i->parent_range, 1);
  }

  f->zdecomp = gkyl_rect_decomp_new_from_ranges(mb_app->cdim, ranges, num_ranges, &f->crossz);

  for (int i = 0; i<num_ranges; i++) {
    printf("For range %d, we have z in [%d, %d]\n", i, f->zdecomp->ranges[i].lower[1], f->zdecomp->ranges[i].upper[1]);
  }


  printf("For the decomp ndecomp is = %d\n", f->zdecomp->ndecomp);
  f->zcomm = gkyl_comm_split_comm(mb_app->comm_mb, 0, f->zdecomp);

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
    int rank;
    gkyl_comm_get_rank(field_mb->zcomm, &rank);

    // Get fin and the field for app we own
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    struct gk_field* field = app->field ;

    // accumulate rho_c in each block
    gk_field_accumulate_rho_c(app, field, fin);

    // Now gather charge density into global interblock array for smoothing in z
    gkyl_comm_array_allgather(field_mb->zcomm, &field_mb->zdecomp->ranges[rank], &field_mb->crossz, field->rho_c, field_mb->rho_c_global_dg);
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

