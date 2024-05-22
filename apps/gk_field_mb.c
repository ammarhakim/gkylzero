#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_mb_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

// For the 1x case this is what needs to happen:
// Each block's field object can gather its charge density into the field->rho_c_global_dg
// This mb_field will gather all of those into mb_field->rho_c_global_dg which is defined on the true global (cross-block) range
// mb_field needs a smoother defined on the true global range
// mb_field will smooth/multiply by kperpsqrho_s and place the result in mb_field->phi_fem defined on the global range
// mb_field will copy this back into the individual blocks field->phi_smooth arrays

// Memory needed by the mb field:
// - rho_c_global_dg on cross-block range
// - rho_c_global_smooth on cross-block range
// - phi_fem_flobal defined on the cross block range
// - fem parproj updater defined on cross block range

// initialize field object
struct gk_field_mb* 
gk_field_mb_new(struct gkyl_gk_mb *gk_mb, struct gkyl_gyrokinetic_mb_app *mb_app)
{
  struct gk_field_mb *f = gkyl_malloc(sizeof(struct gk_field_mb));

  f->info = gk_mb->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  // allocate arrays for charge density
  f->rho_c_global_dg = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, mb_app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, mb_app->global_ext.volume);

  // allocate arrays for electrostatic potential
  // global phi (only used in 1x simulations)
  f->phi_fem = mkarr(mb_app->use_gpu, mb_app->confBasis.num_basis, mb_app->global_ext.volume);
  f->phi_host = mkarr(false, mb_app->confBasis.num_basis, mb_app->global_ext.volume);


  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&f->global_sub_range, &mb_app->global, &mb_app->local);
  if (mb_app->cdim == 1) {
    // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
    f->weight = mkarr(false, mb_app->confBasis.num_basis, mb_app->global_ext.volume); // fem_parproj expects weight on host
    // Here I need to gather all the weights from each block into one global weight
    int rank;
    gkyl_comm_get_rank(mb_app->comm, &rank);
    struct gk_field* bfield = mb_app->blocks[rank]->field ;
    gkyl_comm_array_allgather(mb_app->comm, &mb_app->local, &mb_app->global, bfield->weight, f->weight);
  }
  f->fem_parproj = gkyl_fem_parproj_new(&mb_app->global, &mb_app->global_ext, 
    &mb_app->confBasis, f->info.fem_parbc, f->weight, mb_app->use_gpu);


  return f;
}

// Compute the electrostatic potential
void
gk_field_mb_rhs(gkyl_gyrokinetic_mb_app *mb_app, struct gk_field_mb *field)
{
    // gather charge density into global array for smoothing in z
    int rank;
    gkyl_comm_get_rank(mb_app->comm, &rank);
    struct gk_field* bfield = mb_app->blocks[rank]->field ;
    gkyl_comm_array_allgather(mb_app->comm, &mb_app->local, &mb_app->global, bfield->rho_c, field->rho_c_global_dg);
    if (mb_app->cdim == 1) {
      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
      gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);
      // copy globally smoothed potential to local potential per process for update
      gkyl_array_copy_range_to_range(bfield->phi_smooth, field->phi_fem, &mb_app->local, &field->global_sub_range);
    }
}

// release resources for field
void
gk_field_mb_release(const gkyl_gyrokinetic_mb_app* mb_app, struct gk_field_mb *f)
{
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);

  if (mb_app->cdim == 1) {
    gkyl_array_release(f->weight);
  }

  gkyl_fem_parproj_release(f->fem_parproj);

  if (mb_app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_free(f);
}

