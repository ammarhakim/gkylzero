#include <gkyl_amr_block_priv.h>

void
skin_ghost_ranges_init(struct skin_ghost_ranges* sgr, const struct gkyl_range* parent, const int* ghost)
{
  int ndim = parent -> ndim;
  
  for (int d = 0; d < ndim; d++)
  {
    gkyl_skin_ghost_ranges(&sgr -> lower_skin[d], &sgr -> lower_ghost[d], d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr -> upper_skin[d], &sgr -> upper_ghost[d], d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

static void
euler_transmissive_bc(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++)
  {
    ghost[i] = skin[i];
  }
}

void
euler_block_bc_updaters_init(struct euler_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[] = { 2, 2, 2, 2, 2, 2, 2, 2, 2 };
  
  for (int d = 0; d < 2; d++)
  {
    bdata -> lower_bc[d] = bdata -> upper_bc[d] = 0;
    
    if (conn -> connections[d][0].edge == GKYL_PHYSICAL)
    {
      bdata -> lower_bc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler, bdata -> geom, d, GKYL_LOWER_EDGE, nghost, euler_transmissive_bc, 0);
    }
    
    if (conn -> connections[d][1].edge == GKYL_PHYSICAL)
    {
      bdata -> upper_bc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler, bdata -> geom, d, GKYL_UPPER_EDGE, nghost, euler_transmissive_bc, 0);
    }
  }
  
  skin_ghost_ranges_init(&bdata -> skin_ghost, &bdata -> ext_range, nghost);
  long buff_sz = 0;
  
  for (int d = 0; d < 2; d++)
  {
    long vol = bdata -> skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol)
    {
      buff_sz = vol;
    }
  }
  
  bdata -> bc_buffer = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
}

void
euler_block_bc_updaters_release(struct euler_block_data* bdata)
{
  for (int d = 0; d < 2; d++)
  {
    if (bdata -> lower_bc[d])
    {
      gkyl_wv_apply_bc_release(bdata -> lower_bc[d]);
    }
    if (bdata -> upper_bc[d])
    {
      gkyl_wv_apply_bc_release(bdata -> upper_bc[d]);
    }
  }
  
  gkyl_array_release(bdata -> bc_buffer);
}

void
euler_block_bc_updaters_apply(const struct euler_block_data* bdata, double tm, struct gkyl_array* fld)
{
  for (int d = 0; d < 2; d++)
  {
    if (bdata -> lower_bc[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc[d], tm, &bdata -> range, fld);
    }
    if (bdata -> upper_bc[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc[d], tm, &bdata -> range, fld);
    }
  }
}

void
euler_sync_blocks(const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], struct gkyl_array* fld[])
{
  long num_blocks = btopo -> num_blocks;
  long ndim = btopo -> ndim;

  for (int i = 0; i < num_blocks; i++)
  {
    for (int d = 0; d < ndim; d++)
    {
      const struct gkyl_target_edge *te = btopo -> conn[i].connections[d];

      if (te[0].edge != GKYL_PHYSICAL)
      {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer -> data, fld[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        if (te[0].edge == GKYL_UPPER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }

      if (te[1].edge != GKYL_PHYSICAL)
      {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer -> data, fld[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        if (te[1].edge == GKYL_UPPER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }
    }
  }
}

void
euler_block_data_write(const char* fileNm, const struct euler_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata -> grid, &bdata -> range, bdata -> f[0], fileNm);
}

double
euler_block_data_max_dt(const struct euler_block_data* bdata)
{
  double dt = DBL_MAX;

  for (int d = 0; d < 2; d++)
  {
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata -> slvr[d], &bdata -> range, bdata -> f[0]));
  }

  return dt;
}

void
euler_update_block_job_func(void* ctx)
{
  struct euler_update_block_ctx *ub_ctx = ctx;
  const struct euler_block_data *bdata = ub_ctx -> bdata;

  int d = ub_ctx -> dir;
  double t_curr = ub_ctx -> t_curr;
  double dt = ub_ctx -> dt;

  ub_ctx -> stat = gkyl_wave_prop_advance(bdata -> slvr[d], t_curr, dt, &bdata -> range, bdata -> f[d], bdata -> f[d + 1]);

  euler_block_bc_updaters_apply(bdata, t_curr, bdata -> f[d + 1]);
}

struct gkyl_update_status
euler_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], double t_curr, double dt)
{
  long num_blocks = btopo -> num_blocks;
  long ndim = btopo -> ndim;

  double dt_suggested = DBL_MAX;

  for (int d = 0; d < ndim; d++)
  {
    struct euler_update_block_ctx euler_block_ctx[num_blocks];

    for (int i = 0; i < num_blocks; i++)
    {
      euler_block_ctx[i] = (struct euler_update_block_ctx) {
        .bdata = &bdata[i],
        .t_curr = t_curr,
        .dir = d,
        .dt = dt,
        .bidx = i,
      };
    }

    for (int i = 0; i < num_blocks; i++)
    {
      gkyl_job_pool_add_work(job_pool, euler_update_block_job_func, &euler_block_ctx[i]);
    }
    gkyl_job_pool_wait(job_pool);

    struct gkyl_array *fld[num_blocks];

    for (int i = 0; i < num_blocks; i++)
    {
      if (euler_block_ctx[i].stat.success == false)
      {
        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = euler_block_ctx[i].stat.dt_suggested
        };
      }

      dt_suggested = fmin(dt_suggested, euler_block_ctx[i].stat.dt_suggested);
      fld[i] = bdata[i].f[d + 1];
    }

    euler_sync_blocks(btopo, bdata, fld);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

void
euler_init_job_func(void* ctx)
{
  struct euler_block_data *bdata = ctx;
  gkyl_fv_proj_advance(bdata -> fv_proj, 0.0, &bdata -> ext_range, bdata -> f[0]);
}

void
euler_copy_job_func(void* ctx)
{
  struct euler_copy_job_ctx *j_ctx = ctx;
  gkyl_array_copy(j_ctx -> out, j_ctx -> inp);
}

struct gkyl_update_status euler_update(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[],
    double t_curr, double dt0, struct sim_stats* stats)
{
  long num_blocks = btopo -> num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FLUID_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct euler_copy_job_ctx euler_copy_ctx[num_blocks];
  double dt = dt0;

  while (state != UPDATE_DONE)
  {
    if (state == PRE_UPDATE)
    {
      state = FLUID_UPDATE;

      for (int i = 0; i < num_blocks; i++)
      {
        euler_copy_ctx[i] = (struct euler_copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].f[0],
          .out = bdata[i].fdup
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, euler_copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
    else if (state == FLUID_UPDATE)
    {
      state = POST_UPDATE;

      struct gkyl_update_status s = euler_update_all_blocks(job_pool, btopo, bdata, t_curr, dt);

      if (!s.success)
      {
        stats -> nfail += 1;
        dt = s.dt_suggested;
        state = UPDATE_REDO;
      }
      else
      {
        dt_suggested = fmin(dt_suggested, s.dt_suggested);
      }
    }
    else if (state == POST_UPDATE)
    {
      state = UPDATE_DONE;

      for (int i = 0; i < num_blocks; i++)
      {
        euler_copy_ctx[i] = (struct euler_copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].f[2],
          .out = bdata[i].f[0]
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, euler_copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
    else if (state == UPDATE_REDO)
    {
      state = PRE_UPDATE;

      for (int i = 0; i < num_blocks; i++)
      {
        euler_copy_ctx[i] = (struct euler_copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].fdup,
          .out = bdata[i].f[0]
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, euler_copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

void
euler_write_sol(const char* fbase, int num_blocks, const struct euler_block_data bdata[])
{
  for (int i = 0; i < num_blocks; i++)
  {
    const char *fmt = "%s_b%d.gkyl";
    int sz = snprintf(0, 0, fmt, fbase, i);
    char fileNm[sz + 1];
    
    snprintf(fileNm, sizeof fileNm, fmt, fbase, i);
    euler_block_data_write(fileNm, &bdata[i]);
  }
}

double
euler_max_dt(int num_blocks, const struct euler_block_data bdata[])
{
  double dt = DBL_MAX;
  
  for (int i = 0; i < num_blocks; i++)
  {
    dt = fmin(dt, euler_block_data_max_dt(&bdata[i]));
  }
  
  return dt;
}

void
five_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[] = { 2, 2, 2, 2, 2, 2, 2, 2, 2 };
  
  for (int d = 0; d < 2; d++)
  {
    bdata -> lower_bc_elc[d] = bdata -> upper_bc_elc[d] = 0;
    bdata -> lower_bc_ion[d] = bdata -> upper_bc_ion[d] = 0;
    bdata -> lower_bc_maxwell[d] = bdata -> upper_bc_maxwell[d] = 0;
    
    if (conn -> connections[d][0].edge == GKYL_PHYSICAL)
    {
      bdata -> lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_elc, bdata -> geom, d, GKYL_LOWER_EDGE, nghost, bdata -> euler_elc -> wall_bc_func, 0);
      bdata -> lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_ion, bdata -> geom, d, GKYL_LOWER_EDGE, nghost, bdata -> euler_ion -> wall_bc_func, 0);
      bdata -> lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> maxwell, bdata -> geom, d, GKYL_LOWER_EDGE, nghost, bdata -> maxwell -> wall_bc_func, 0);
    }
    
    if (conn -> connections[d][1].edge == GKYL_PHYSICAL)
    {
      bdata -> upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_elc, bdata -> geom, d, GKYL_UPPER_EDGE, nghost, bdata -> euler_elc -> wall_bc_func, 0);
      bdata -> upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_ion, bdata -> geom, d, GKYL_UPPER_EDGE, nghost, bdata -> euler_ion -> wall_bc_func, 0);
      bdata -> upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> maxwell, bdata -> geom, d, GKYL_UPPER_EDGE, nghost, bdata -> maxwell -> wall_bc_func, 0);
    }
  }
  
  skin_ghost_ranges_init(&bdata -> skin_ghost, &bdata -> ext_range, nghost);
  long buff_sz = 0;
  
  for (int d = 0; d < 2; d++)
  {
    long vol = bdata -> skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol)
    {
      buff_sz = vol;
    }
  }
  
  bdata -> bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata -> bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata -> bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
five_moment_block_bc_updaters_release(struct five_moment_block_data* bdata)
{
  for (int d = 0; d < 2; d++)
  {
    if (bdata -> lower_bc_elc[d])
    {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_elc[d]);
    }
    if (bdata -> upper_bc_elc[d])
    {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_elc[d]);
    }

    if (bdata -> lower_bc_ion[d])
    {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_ion[d]);
    }
    if (bdata -> upper_bc_ion[d])
    {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_ion[d]);
    }

    if (bdata -> lower_bc_maxwell[d])
    {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_maxwell[d]);
    }
    if (bdata -> upper_bc_maxwell[d])
    {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_maxwell[d]);
    }
  }
  
  gkyl_array_release(bdata -> bc_buffer_elc);
  gkyl_array_release(bdata -> bc_buffer_ion);
  gkyl_array_release(bdata -> bc_buffer_maxwell);
}

void
five_moment_block_apply_periodic_bc(const struct five_moment_block_data* bdata, int dir, struct gkyl_array* fld_elc, struct gkyl_array* fld_ion, struct gkyl_array* fld_maxwell)
{
  gkyl_array_copy_to_buffer(bdata -> bc_buffer_elc -> data, fld_elc, &(bdata -> skin_ghost.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld_elc, bdata -> bc_buffer_elc -> data, &(bdata -> skin_ghost.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(bdata -> bc_buffer_elc -> data, fld_elc, &(bdata -> skin_ghost.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld_elc, bdata -> bc_buffer_elc -> data, &(bdata -> skin_ghost.lower_ghost[dir]));

  gkyl_array_copy_to_buffer(bdata -> bc_buffer_ion -> data, fld_ion, &(bdata -> skin_ghost.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld_ion, bdata -> bc_buffer_ion -> data, &(bdata -> skin_ghost.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(bdata -> bc_buffer_ion -> data, fld_ion, &(bdata -> skin_ghost.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld_ion, bdata -> bc_buffer_ion -> data, &(bdata -> skin_ghost.lower_ghost[dir]));

  gkyl_array_copy_to_buffer(bdata -> bc_buffer_maxwell -> data, fld_maxwell, &(bdata -> skin_ghost.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld_maxwell, bdata -> bc_buffer_maxwell -> data, &(bdata -> skin_ghost.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(bdata -> bc_buffer_maxwell -> data, fld_maxwell, &(bdata -> skin_ghost.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld_maxwell, bdata -> bc_buffer_maxwell -> data, &(bdata -> skin_ghost.lower_ghost[dir]));
}

void
five_moment_block_bc_updaters_apply(const struct five_moment_block_data* bdata, double tm, struct gkyl_array* fld_elc, struct gkyl_array* fld_ion, struct gkyl_array* fld_maxwell)
{
  five_moment_block_apply_periodic_bc(bdata, 0, fld_elc, fld_ion, fld_maxwell);

  for (int d = 1; d < 2; d++)
  {
    if (bdata -> lower_bc_elc[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_elc[d], tm, &bdata -> range, fld_elc);
    }
    if (bdata -> upper_bc_elc[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_elc[d], tm, &bdata -> range, fld_elc);
    }

    if (bdata -> lower_bc_ion[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_ion[d], tm, &bdata -> range, fld_ion);
    }
    if (bdata -> upper_bc_ion[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_ion[d], tm, &bdata -> range, fld_ion);
    }

    if (bdata -> lower_bc_maxwell[d])
    {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_maxwell[d], tm, &bdata -> range, fld_maxwell);
    }
    else
    {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_maxwell[d], tm, &bdata -> range, fld_maxwell);
    }
  }
}

void
five_moment_sync_blocks(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[], struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  long num_blocks = btopo -> num_blocks;
  long ndim = btopo -> ndim;

  for (int i = 0; i < num_blocks; i++)
  {
    for (int d = 0; d < ndim; d++)
    {
      const struct gkyl_target_edge *te = btopo -> conn[i].connections[d];

      if (te[0].edge != GKYL_PHYSICAL)
      {
        struct gkyl_array *bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array *bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array *bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc -> data, fld_elc[i], &(bdata[i].skin_ghost.lower_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_ion -> data, fld_ion[i], &(bdata[i].skin_ghost.lower_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_maxwell -> data, fld_maxwell[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        if (te[0].edge == GKYL_UPPER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }

      if (te[1].edge != GKYL_PHYSICAL)
      {
        struct gkyl_array *bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array *bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array *bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc -> data, fld_elc[i], &(bdata[i].skin_ghost.upper_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_ion -> data, fld_ion[i], &(bdata[i].skin_ghost.upper_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_maxwell -> data, fld_maxwell[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        if (te[1].edge == GKYL_UPPER_POSITIVE)
        {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }
    }
  }
}

void
five_moment_block_data_write_elc(const char* fileNmElc, const struct five_moment_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata -> grid, &bdata -> range, bdata -> f_elc[0], fileNmElc);
}

void
five_moment_block_data_write_ion(const char* fileNmIon, const struct five_moment_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata -> grid, &bdata -> range, bdata -> f_ion[0], fileNmIon);
}

void
five_moment_block_data_write_maxwell(const char* fileNmMaxwell, const struct five_moment_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata -> grid, &bdata -> range, bdata -> f_maxwell[0], fileNmMaxwell);
}

double
five_moment_block_data_max_dt(const struct five_moment_block_data* bdata)
{
  double dt = DBL_MAX;

  for (int d = 0; d < 2; d++)
  {
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata -> slvr_elc[d], &bdata -> range, bdata -> f_elc[0]));
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata -> slvr_ion[d], &bdata -> range, bdata -> f_ion[0]));
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata -> slvr_maxwell[d], &bdata -> range, bdata -> f_maxwell[0]));
  }

  return dt;
}

void
five_moment_update_block_job_func(void* ctx)
{
  struct five_moment_update_block_ctx *ub_ctx = ctx;
  const struct five_moment_block_data *bdata = ub_ctx -> bdata;

  int d = ub_ctx -> dir;
  double t_curr = ub_ctx -> t_curr;
  double dt = ub_ctx -> dt;

  ub_ctx -> stat_elc = gkyl_wave_prop_advance(bdata -> slvr_elc[d], t_curr, dt, &bdata -> range, bdata -> f_elc[d], bdata -> f_elc[d + 1]);
  ub_ctx -> stat_ion = gkyl_wave_prop_advance(bdata -> slvr_ion[d], t_curr, dt, &bdata -> range, bdata -> f_ion[d], bdata -> f_ion[d + 1]);
  ub_ctx -> stat_maxwell = gkyl_wave_prop_advance(bdata -> slvr_maxwell[d], t_curr, dt, &bdata -> range, bdata -> f_maxwell[d], bdata -> f_maxwell[d + 1]);

  five_moment_block_bc_updaters_apply(bdata, t_curr, bdata -> f_elc[d + 1], bdata -> f_ion[d + 1], bdata -> f_maxwell[d + 1]);
}

void
five_moment_update_block_job_func_source(void* ctx, int nstrang)
{
  struct five_moment_update_block_ctx *ub_ctx = ctx;
  const struct five_moment_block_data *bdata = ub_ctx -> bdata;

  int d = ub_ctx -> dir;
  double t_curr = ub_ctx -> t_curr;
  double dt = ub_ctx -> dt;

  struct gkyl_array *species[2];

  species[0] = bdata -> f_elc[nstrang];
  species[1] = bdata -> f_ion[nstrang];

  const struct gkyl_array *app_accel[2];
  const struct gkyl_array *rhs_source[2];
  const struct gkyl_array *nT_source[2];

  app_accel[0] = bdata -> app_accel_elc;
  app_accel[1] = bdata -> app_accel_ion;
  
  rhs_source[0] = bdata -> rhs_source_elc;
  rhs_source[1] = bdata -> rhs_source_ion;

  nT_source[0] = bdata -> nT_source_elc;
  nT_source[1] = bdata -> nT_source_ion;

  gkyl_moment_em_coupling_advance(bdata -> src_slvr, t_curr, dt, &bdata -> range, species, app_accel, rhs_source, bdata -> f_maxwell[nstrang], bdata -> app_current, bdata -> ext_em, nT_source);

  five_moment_block_bc_updaters_apply(bdata, t_curr, bdata -> f_elc[nstrang], bdata -> f_ion[nstrang], bdata -> f_maxwell[nstrang]);
}

struct gkyl_update_status
five_moment_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[], double t_curr, double dt)
{
  long num_blocks = btopo -> num_blocks;
  long ndim = btopo -> ndim;

  double dt_suggested = DBL_MAX;

  for (int d = 0; d < ndim; d++)
  {
    struct five_moment_update_block_ctx five_moment_block_ctx[num_blocks];

    for (int i = 0; i < num_blocks; i++)
    {
      five_moment_block_ctx[i] = (struct five_moment_update_block_ctx) {
        .bdata = &bdata[i],
        .t_curr = t_curr,
        .dir = d,
        .dt = dt,
        .bidx = i,
      };
    }

    for (int i = 0; i < num_blocks; i++)
    {
      gkyl_job_pool_add_work(job_pool, five_moment_update_block_job_func, &five_moment_block_ctx[i]);
    }
    gkyl_job_pool_wait(job_pool);

    struct gkyl_array *fld_elc[num_blocks];
    struct gkyl_array *fld_ion[num_blocks];
    struct gkyl_array *fld_maxwell[num_blocks];

    for (int i = 0; i < num_blocks; i++)
    {
      if (five_moment_block_ctx[i].stat_elc.success == false || five_moment_block_ctx[i].stat_ion.success == false || five_moment_block_ctx[i].stat_maxwell.success == false)
      {
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_elc.dt_suggested);
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_ion.dt_suggested);
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_maxwell.dt_suggested);

        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = dt_suggested
        };
      }

      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_elc.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_ion.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_maxwell.dt_suggested);

      fld_elc[i] = bdata[i].f_elc[d + 1];
      fld_ion[i] = bdata[i].f_ion[d + 1];
      fld_maxwell[i] = bdata[i].f_maxwell[d + 1];
    }

    five_moment_sync_blocks(btopo, bdata, fld_elc, fld_ion, fld_maxwell);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

void
five_moment_update_all_blocks_sources(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[], double t_curr, double dt, int nstrang)
{
  int num_blocks = btopo -> num_blocks;
  double dt_suggested = DBL_MAX;

  struct five_moment_update_block_ctx five_moment_block_ctx[num_blocks];

  for (int i = 0; i < num_blocks; i++)
  {
    five_moment_block_ctx[i] = (struct five_moment_update_block_ctx) {
      .bdata = &bdata[i],
      .t_curr = t_curr,
      .dir = 0,
      .dt = dt,
      .bidx = i,
    };
  }

  /*
  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_job_pool_add_work(job_pool, five_moment_update_block_job_func_source, &five_moment_block_ctx[i]);
  }
  gkyl_job_pool_wait(job_pool);
  */
  for (int i = 0; i < num_blocks; i++)
  {
    five_moment_update_block_job_func_source(&five_moment_block_ctx[i], nstrang);
  }

  struct gkyl_array *fld_elc[num_blocks];
  struct gkyl_array *fld_ion[num_blocks];
  struct gkyl_array *fld_maxwell[num_blocks];

  for (int i = 0; i < num_blocks; i++)
  {
    fld_elc[i] = bdata[i].f_elc[nstrang];
    fld_ion[i] = bdata[i].f_ion[nstrang];
    fld_maxwell[i] = bdata[i].f_maxwell[nstrang];
  }

  five_moment_sync_blocks(btopo, bdata, fld_elc, fld_ion, fld_maxwell);
}

void
five_moment_init_job_func(void* ctx)
{
  struct five_moment_block_data *bdata = ctx;

  gkyl_fv_proj_advance(bdata -> fv_proj_elc, 0.0, &bdata -> ext_range, bdata -> f_elc[0]);
  gkyl_fv_proj_advance(bdata -> fv_proj_ion, 0.0, &bdata -> ext_range, bdata -> f_ion[0]);
  gkyl_fv_proj_advance(bdata -> fv_proj_maxwell, 0.0, &bdata -> ext_range, bdata -> f_maxwell[0]);
}

void
five_moment_copy_job_func(void* ctx)
{
  struct five_moment_copy_job_ctx *j_ctx = ctx;

  gkyl_array_copy(j_ctx -> out_elc, j_ctx -> inp_elc);
  gkyl_array_copy(j_ctx -> out_ion, j_ctx -> inp_ion);
  gkyl_array_copy(j_ctx -> out_maxwell, j_ctx -> inp_maxwell);
}

struct gkyl_update_status five_moment_update(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[],
    double t_curr, double dt0, struct sim_stats* stats)
{
  long num_blocks = btopo -> num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    BLOCK_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct five_moment_copy_job_ctx five_moment_copy_ctx[num_blocks];
  double dt = dt0;

  while (state != UPDATE_DONE)
  {
    if (state == PRE_UPDATE)
    {
      state = FIRST_COUPLING_UPDATE;

      for (int i = 0; i < num_blocks; i++)
      {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].f_elc[0],
          .inp_ion = bdata[i].f_ion[0],
          .inp_maxwell = bdata[i].f_maxwell[0],
          .out_elc = bdata[i].fdup_elc,
          .out_ion = bdata[i].fdup_ion,
          .out_maxwell = bdata[i].fdup_maxwell,
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
    else if (state == FIRST_COUPLING_UPDATE)
    {
      state = BLOCK_UPDATE;

      five_moment_update_all_blocks_sources(btopo, bdata, t_curr, 0.5 * dt, 0);
    }
    else if (state == BLOCK_UPDATE)
    {
      state = SECOND_COUPLING_UPDATE;

      struct gkyl_update_status s = five_moment_update_all_blocks(job_pool, btopo, bdata, t_curr, dt);

      if (!s.success)
      {
        stats -> nfail += 1;
        dt = s.dt_suggested;
        state = UPDATE_REDO;
      }
      else
      {
        dt_suggested = fmin(dt_suggested, s.dt_suggested);
      }
    }
    else if (state == SECOND_COUPLING_UPDATE)
    {
      state = POST_UPDATE;

      five_moment_update_all_blocks_sources(btopo, bdata, t_curr, 0.5 * dt, 2);
    }
    else if (state == POST_UPDATE)
    {
      state = UPDATE_DONE;

      for (int i = 0; i < num_blocks; i++)
      {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].f_elc[2],
          .out_elc = bdata[i].f_elc[0],
          .inp_ion = bdata[i].f_ion[2],
          .out_ion = bdata[i].f_ion[0],
          .inp_maxwell = bdata[i].f_maxwell[2],
          .out_maxwell = bdata[i].f_maxwell[0],
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
    else if (state == UPDATE_REDO)
    {
      state = PRE_UPDATE;

      for (int i = 0; i < num_blocks; i++)
      {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].fdup_elc,
          .out_elc = bdata[i].f_elc[0],
          .inp_ion = bdata[i].fdup_ion,
          .out_ion = bdata[i].f_ion[0],
          .inp_maxwell = bdata[i].fdup_maxwell,
          .out_maxwell = bdata[i].f_maxwell[0],
        };
      }

      for (int i = 0; i < num_blocks; i++)
      {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
    }
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

void
five_moment_write_sol(const char* fbase, int num_blocks, const struct five_moment_block_data bdata[])
{
  for (int i = 0; i < num_blocks; i++)
  {
    const char *fmt_elc = "%s_elc_b%d.gkyl";
    int sz_elc = snprintf(0, 0, fmt_elc, fbase, i);
    char fileNm_elc[sz_elc + 1];
    
    snprintf(fileNm_elc, sizeof fileNm_elc, fmt_elc, fbase, i);
    five_moment_block_data_write_elc(fileNm_elc, &bdata[i]);

    const char *fmt_ion = "%s_ion_b%d.gkyl";
    int sz_ion = snprintf(0, 0, fmt_ion, fbase, i);
    char fileNm_ion[sz_ion + 1];
    
    snprintf(fileNm_ion, sizeof fileNm_ion, fmt_ion, fbase, i);
    five_moment_block_data_write_ion(fileNm_ion, &bdata[i]);

    const char *fmt_maxwell = "%s_maxwell_b%d.gkyl";
    int sz_maxwell = snprintf(0, 0, fmt_maxwell, fbase, i);
    char fileNm_maxwell[sz_maxwell + 1];
    
    snprintf(fileNm_maxwell, sizeof fileNm_maxwell, fmt_maxwell, fbase, i);
    five_moment_block_data_write_maxwell(fileNm_maxwell, &bdata[i]);
  }
}

double
five_moment_max_dt(int num_blocks, const struct five_moment_block_data bdata[])
{
  double dt = DBL_MAX;
  
  for (int i = 0; i < num_blocks; i++)
  {
    dt = fmin(dt, five_moment_block_data_max_dt(&bdata[i]));
  }
  
  return dt;
}

struct gkyl_block_topo*
create_block_topo()
{
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 9);
  
  btopo -> conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  btopo -> conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 4, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } }
  };

  btopo -> conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE  }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } }
  };

  btopo -> conn[3] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } }
  };

  btopo -> conn[4] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  btopo -> conn[5] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 8, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  btopo -> conn[6] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  btopo -> conn[7] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 8, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  btopo -> conn[8] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE } }
  };

  return btopo;
}