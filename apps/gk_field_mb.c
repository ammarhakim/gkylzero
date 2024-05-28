#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_mb_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>


// initialize field object
struct gk_field_mb* 
gk_field_mb_new(struct gkyl_gk_mb *gk_mb, struct gkyl_gyrokinetic_mb_app *mb_app)
{

  struct gk_field_mb *f = gkyl_malloc(sizeof(struct gk_field_mb));

  f->info = gk_mb->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  // We need a bunch of ranges
  // Handle only the 3 block case at first
  struct gkyl_gyrokinetic_app **apps = mb_app->blocks;

  // We need first to create the global z ranges.
  int lower[2];
  int upper[2];
  lower[0] = apps[0]->grid.lower[0];
  upper[0] = apps[0]->grid.upper[0];
  lower[1] = apps[0]->grid.lower[1];
  upper[1] = apps[0]->grid.cells[1] + apps[1]->grid.cells[1] + apps[2]->grid.cells[1];
  struct gkyl_range globalz;
  gkyl_range_init(&globalz, mb_app->cdim, lower, upper);
  int nghost[2] = {1,1};
  gkyl_create_ranges(&globalz, nghost, &f->globalz_ext, &f->globalz);

  // Create the decomp and communicator from the mb app communicator
  int cuts[2] = {0,3}; // Sticking to blocks 2,3,4 for now
  struct gkyl_rect_decomp *zdecomp = gkyl_rect_decomp_new_from_cuts(mb_app->cdim, cuts, &f->globalz);
  f->zcomm = gkyl_comm_split_comm(mb_app->comm, 0, zdecomp); // Would have different colors for other block groups like 11-12, 7-8-9
                                                             // Just 0 for now

  // Now get the sub range intersects
  // Create global subrange we'll copy the field solver solution from (into local).
  int rank;
  gkyl_comm_get_rank(f->zcomm, &rank);
  int intersect = gkyl_sub_range_intersect(&f->globalz_sub_range, &f->globalz, &apps[rank]->local);

  // allocate arrays for charge density
  f->rho_c_global_dg = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, f->globalz_ext.volume);
  f->rho_c_global_smooth = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, f->globalz_ext.volume);

  // allocate arrays for electrostatic potential
  // global phi (only used in 1x simulations)
  f->phi = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, f->globalz_ext.volume);

  f->phi_host = f->phi;
  if (mb_app->use_gpu) {
    f->phi_host = mkarr(false, mb_app->confBasis.num_basis, f->globalz_ext.volume);
  }


  f->fem_parproj = gkyl_fem_parproj_new(&f->globalz, &f->globalz_ext, 
    &mb_app->confBasis, f->info.fem_parbc, NULL, mb_app->use_gpu);


  return f;
}

// Compute the electrostatic potential
void
gk_field_mb_rhs(gkyl_gyrokinetic_mb_app *mb_app, struct gk_field_mb *field_mb)
{
    // Get rank and figure out which app is ours
    int rank;
    gkyl_comm_get_rank(field_mb->zcomm, &rank);
    struct gkyl_gyrokinetic_app* app = mb_app->blocks[rank];

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

    // Now gather charge density into global cross-block array for smoothing in z
    gkyl_comm_array_allgather(field_mb->zcomm, &app->local, &field_mb->globalz, field->rho_c_global_dg, field_mb->rho_c_global_dg);
    // Do the smoothing
    gkyl_fem_parproj_set_rhs(field_mb->fem_parproj, field_mb->rho_c_global_dg, field_mb->rho_c_global_dg);
    gkyl_fem_parproj_solve(field_mb->fem_parproj, field_mb->rho_c_global_smooth);
    // Copy globally (in z) smoothed potential to local potential per process
    gkyl_array_copy_range_to_range(field->rho_c_global_smooth, field_mb->rho_c_global_smooth, &app->local, &field_mb->globalz_sub_range);
    // Now call the perp solver
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

