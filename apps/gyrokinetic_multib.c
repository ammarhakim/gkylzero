#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>

#include <mpack.h>

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(struct gkyl_gyrokinetic_multib *mbinp)
{
  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  
  int num_blocks = mbapp->num_blocks = gkyl_block_geom_num_blocks(mbinp->block_geom);
  mbapp->app = gkyl_malloc(sizeof(struct gkyl_gyrokinetic_app*[num_blocks]));

  // Fill each block's app with the needed single block inputs
  for (int i=0; i<num_blocks; ++i) {
    struct gkyl_gk app_pb = {};
    // Fetch top-level parameters
    int cdim = app_pb.cdim = mbinp->cdim; 
    int vdim = app_pb.vdim = mbinp->vdim;
    int num_species = app_pb.num_species = mbinp->num_species;
    int num_neut_species = app_pb.num_neut_species = mbinp->num_neut_species; 

    app_pb.poly_order = mbinp->poly_order;
    app_pb.basis_type = mbinp->basis_type;
    app_pb.use_gpu = mbinp->use_gpu;
    app_pb.cfl_frac = mbinp->cfl_frac;

    const struct gkyl_block_geom_info *info = gkyl_block_geom_get_block(mbinp->block_geom, i);
    // Fetch the configuration space extents and geometry from the block info  
    for (int c=0; c<cdim; ++c) {
      app_pb.lower[c] = info->lower[c];
      app_pb.upper[c] = info->upper[c];
      app_pb.cells[c] = info->cells[c];
    }
    app_pb.geometry = info->geometry; 

    // Fetch the per-block field info
    struct gkyl_gyrokinetic_field field = {};
    app_pb.field = field; 

    // Fetch the per-block species info
    for (int s=0; s<num_species; ++s) {
      struct gkyl_gyrokinetic_species species = {};
      strcpy(species.name, mbinp->species[s].name);
      species.charge = mbinp->species[s].charge;
      species.mass = mbinp->species[s].mass;
      species.gkmodel_id = mbinp->species[s].gkmodel_id;
      species.no_by = mbinp->species[s].no_by;
      species.enforce_positivity = mbinp->species[s].enforce_positivity;

      // Velocity space information
      for (int v=0; v<vdim; ++v) {
        species.lower[v] = mbinp->species[s].lower[v];
        species.upper[v] = mbinp->species[s].upper[v];
        species.cells[v] = mbinp->species[s].cells[v];
      }
      species.mapc2p = mbinp->species[s].mapc2p;

      // Species physics modules
      species.collisions = mbinp->species[s].collisions;
      species.diffusion = mbinp->species[s].diffusion;
      species.radiation = mbinp->species[s].radiation;
      species.react = mbinp->species[s].react;
      species.react_neut = mbinp->species[s].react_neut;

      // Species diagnostics
      species.num_diag_moments = mbinp->species[s].num_diag_moments;
      for (int n=0; n<species.num_diag_moments; ++n) {
        strcpy(species.diag_moments[n], mbinp->species[s].diag_moments[n]);
      }

      // Species inputs which can differ block-to-block
      if (mbinp->species[s].duplicate_across_blocks) {
        species.projection = mbinp->species[s].blocks[0].projection;
        species.source = mbinp->species[s].blocks[0].source;
        species.polarization_density = mbinp->species[s].blocks[0].polarization_density;
      }
      else {
        species.projection = mbinp->species[s].blocks[i].projection;
        species.source = mbinp->species[s].blocks[i].source;
        species.polarization_density = mbinp->species[s].blocks[i].polarization_density;
      }

      // Assign the block's app's species object to the constructed species object
      app_pb.species[s] = species;
    }

    // Fetch the per-block neutral species info
    for (int ns=0; ns<num_neut_species; ++ns) {
      struct gkyl_gyrokinetic_neut_species neut_species = {};
      strcpy(neut_species.name, mbinp->neut_species[ns].name);
      neut_species.mass = mbinp->neut_species[ns].mass;
      neut_species.is_static = mbinp->neut_species[ns].is_static;

      // Velocity space information
      for (int v=0; v<vdim; ++v) {
        neut_species.lower[v] = mbinp->neut_species[ns].lower[v];
        neut_species.upper[v] = mbinp->neut_species[ns].upper[v];
        neut_species.cells[v] = mbinp->neut_species[ns].cells[v];
      }
      neut_species.mapc2p = mbinp->neut_species[ns].mapc2p;

      // Neutral species physics modules
      neut_species.react_neut = mbinp->neut_species[ns].react_neut;

      // Neutral species diagnostics
      neut_species.num_diag_moments = mbinp->neut_species[ns].num_diag_moments;
      for (int n=0; n<neut_species.num_diag_moments; ++n) {
        strcpy(neut_species.diag_moments[n], mbinp->neut_species[ns].diag_moments[n]);
      }

      // Neutral species inputs which can differ block-to-block
      if (mbinp->neut_species[ns].duplicate_across_blocks) {
        neut_species.projection = mbinp->neut_species[ns].blocks[0].projection;
        neut_species.source = mbinp->neut_species[ns].blocks[0].source;
      }
      else {
        neut_species.projection = mbinp->neut_species[ns].blocks[i].projection;
        neut_species.source = mbinp->neut_species[ns].blocks[i].source;
      }

      // Assign the block's app's neutral species object to the constructed neutral species object
      app_pb.neut_species[ns] = neut_species;
    }

    mbapp->app[i] = gkyl_gyrokinetic_app_new(&app_pb);     
  }

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {
  };
  
  return mbapp;
}


void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  // TO DO
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}


void gkyl_gyrokinetic_multib_app_calc_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_integrated_neut_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}


void gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_coll_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_drag(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_emissivity(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_integrated_moms(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  // TO DO
}


void gkyl_gyrokinetic_multib_app_write_iz_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_iz_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_source_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_max_corr_status(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_read_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  // TO DO
}

struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  // TO DO
  return (struct gkyl_update_status) { };
}

struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  return app->stat;
}

void gkyl_gyrokinetic_multib_app_species_ktm_rhs(gkyl_gyrokinetic_multib_app* app, int update_vol_term)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  gkyl_free(mbapp->app);
  gkyl_free(mbapp);
}
