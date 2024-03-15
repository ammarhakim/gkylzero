#pragma once

#include <gkyl_amr_block_priv.h>

// Definitions of private structs and APIs attached to these objects, for use in the coupled block AMR subsystem.

// Block-structured data (including grid, range and boundary condition information) for the coupled five-moment equations.
struct five_moment_block_data {
  struct gkyl_rect_grid grid;

  gkyl_fv_proj *fv_proj_elc;
  gkyl_fv_proj *fv_proj_ion;
  gkyl_fv_proj *fv_proj_maxwell;

  struct gkyl_range ext_range;
  struct gkyl_range range;

  struct gkyl_array *fdup_elc;
  struct gkyl_array *fdup_ion;
  struct gkyl_array *fdup_maxwell;

  struct gkyl_array *f_elc[9];
  struct gkyl_array *f_ion[9];
  struct gkyl_array *f_maxwell[9];

  struct gkyl_array *app_accel_elc;
  struct gkyl_array *app_accel_ion;

  struct gkyl_array *rhs_source_elc;
  struct gkyl_array *rhs_source_ion;

  struct gkyl_array *nT_src_elc;
  struct gkyl_array *nT_src_ion;

  struct gkyl_array *app_current;
  struct gkyl_array *ext_em;

  struct gkyl_wave_geom *geom;

  struct gkyl_wv_eqn *euler_elc;
  struct gkyl_wv_eqn *euler_ion;
  struct gkyl_wv_eqn *maxwell;

  gkyl_wave_prop *slvr_elc[2];
  gkyl_wave_prop *slvr_ion[2];
  gkyl_wave_prop *slvr_maxwell[2];

  struct skin_ghost_ranges skin_ghost;

  struct gkyl_array *bc_buffer_elc;
  struct gkyl_array *bc_buffer_ion;
  struct gkyl_array *bc_buffer_maxwell;

  struct gkyl_wv_apply_bc *lower_bc_elc[2];
  struct gkyl_wv_apply_bc *lower_bc_ion[2];
  struct gkyl_wv_apply_bc *lower_bc_maxwell[2];

  struct gkyl_wv_apply_bc *upper_bc_elc[2];
  struct gkyl_wv_apply_bc *upper_bc_ion[2];
  struct gkyl_wv_apply_bc *upper_bc_maxwell[2];

  gkyl_moment_em_coupling *src_slvr;
};

// Job pool information context for updating block-structured data for the coupled five-moment equations using threads.
struct five_moment_update_block_ctx {
  const struct five_moment_block_data *bdata;
  int bidx;
  int dir;
  double t_curr;
  double dt;

  struct gkyl_wave_prop_status stat_elc;
  struct gkyl_wave_prop_status stat_ion;
  struct gkyl_wave_prop_status stat_maxwell;
};

// Context for copying job pool information for the block-structured, coupled five-moment equations.
struct five_moment_copy_job_ctx {
  int bidx;

  const struct gkyl_array *inp_elc;
  const struct gkyl_array *inp_ion;
  const struct gkyl_array *inp_maxwell;

  struct gkyl_array *out_elc;
  struct gkyl_array *out_ion;
  struct gkyl_array *out_maxwell;
};

/**
* Boundary condition function for applying wall boundary conditions for the coupled five-moment equations.
*
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
static void five_moment_wall_bc(double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying wall boundary conditions for the Maxwell equations.
*
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
static void maxwell_wall_bc(double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Initialize updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @param conn Topology/connectivity data for the block hierarchy.
*/
void five_moment_block_bc_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn);

/**
* Release updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
*/
void five_moment_block_bc_updaters_release(struct five_moment_block_data* bdata);

/**
* Apply both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @param tm Simulation time at which the boundary conditions are applied.
* @param fld Output array.
*/
void five_moment_block_bc_updaters_apply(const struct five_moment_block_data* bdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array *fld_ion, struct gkyl_array* fld_maxwell);

/**
* Synchronize all blocks in the block hierarchy by applying all appropriate physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param btopo Topology/connectivity information for the block hierarchy.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param fld Output array.
*/
void five_moment_sync_blocks(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Write block-structured simulation data for the coupled five-moment equations onto disk.
*
* @param file_nm_elc File name schema to use for the electron simulation output.
* @param file_nm_ion File name schema to use for the ion simulation output.
* @param file_nm_maxwell File name schema to use for the Maxwell field simulation output.
* @param bdata Block-structured data for the coupled five-moment equations.
*/
void five_moment_block_data_write(const char* file_nm_elc, const char* file_nm_ion, const char* file_nm_maxwell, const struct five_moment_block_data* bdata);