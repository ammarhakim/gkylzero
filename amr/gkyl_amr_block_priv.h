#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_block_topo.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_null_pool.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_thread_pool.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_maxwell.h>
#include <rt_arg_parse.h>

#include <thpool.h>

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[2];
  struct gkyl_range lower_ghost[2];
  
  struct gkyl_range upper_skin[2];
  struct gkyl_range upper_ghost[2];
};

struct euler_block_data {
  struct gkyl_rect_grid grid;
  gkyl_fv_proj *fv_proj;
  
  struct gkyl_range ext_range, range;
  struct gkyl_array *fdup;
  struct gkyl_array *f[9];
  
  struct gkyl_wave_geom *geom;
  
  struct gkyl_wv_eqn *euler;
  gkyl_wave_prop *slvr[2];
  
  struct skin_ghost_ranges skin_ghost;
  struct gkyl_array *bc_buffer;
  
  struct gkyl_wv_apply_bc *lower_bc[2];
  struct gkyl_wv_apply_bc *upper_bc[2];
};

struct euler_update_block_ctx {
  const struct euler_block_data *bdata;
  int bidx;
  int dir;
  double t_curr;
  double dt;
  struct gkyl_wave_prop_status stat;
};

struct copy_job_ctx {
  int bidx;
  const struct gkyl_array *inp;
  struct gkyl_array *out;
};

struct sim_stats {
  int nfail;
};

void skin_ghost_ranges_init(struct skin_ghost_ranges* sgr, const struct gkyl_range* parent, const int* ghost);

void euler_block_bc_updaters_init(struct euler_block_data* bdata, const struct gkyl_block_connections* conn);

void euler_block_bc_updaters_release(struct euler_block_data* bdata);

void euler_block_bc_updaters_apply(const struct euler_block_data* bdata, double tm, struct gkyl_array* fld);

void euler_sync_blocks(const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], struct gkyl_array* fld[]);

double euler_block_data_max_dt(const struct euler_block_data* bdata);

void euler_update_block_job_func(void *ctx);

struct gkyl_update_status euler_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], double t_curr, double dt);

void euler_init_job_func(void* ctx);

void copy_job_func(void* ctx);

struct gkyl_update_status euler_update(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], double t_curr, double dt0, struct sim_stats* stats);

struct gkyl_block_topo* create_block_topo();