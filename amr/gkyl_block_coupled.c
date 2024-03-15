#include <gkyl_amr_block_coupled_priv.h>

static void
five_moment_wall_bc(double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
maxwell_wall_bc(double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 8; i++) {
    if (i == 1 || i == 2 || i == 3 || i == 6) {
      ghost[i] = -skin[i];
    }
    else {
      ghost[i] = skin[i];
    }
  }
}

void
five_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[9];
  for (int i = 0; i < 9; i++) {
    nghost[i] = 2;
  }

  for (int d = 0; d < 2; d++) {
    bdata -> lower_bc_elc[d] = bdata -> upper_bc_elc[d] = 0;
    bdata -> lower_bc_ion[d] = bdata -> upper_bc_ion[d] = 0;
    bdata -> lower_bc_maxwell[d] = bdata -> upper_bc_maxwell[d] = 0;

    if (conn -> connections[d][0].edge == GKYL_PHYSICAL) {
      bdata -> lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_elc, bdata -> geom, d, GKYL_LOWER_EDGE, nghost,
        five_moment_wall_bc, 0);
      bdata -> lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_ion, bdata -> geom, d, GKYL_LOWER_EDGE, nghost,
        five_moment_wall_bc, 0);
      bdata -> lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> maxwell, bdata -> geom, d, GKYL_LOWER_EDGE, nghost,
        maxwell_wall_bc, 0);
    }

    if (conn -> connections[d][1].edge == GKYL_PHYSICAL) {
      bdata -> upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_elc, bdata -> geom, d, GKYL_UPPER_EDGE, nghost,
        five_moment_wall_bc, 0);
      bdata -> upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> euler_ion, bdata -> geom, d, GKYL_UPPER_EDGE, nghost,
        five_moment_wall_bc, 0);
      bdata -> upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata -> grid, bdata -> maxwell, bdata -> geom, d, GKYL_UPPER_EDGE, nghost,
        maxwell_wall_bc, 0);
    }
  }

  skin_ghost_ranges_init(&bdata -> skin_ghost, &bdata -> ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata -> skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }
}

void
five_moment_block_bc_updaters_release(struct five_moment_block_data* bdata)
{
  for (int d = 0; d < 2; d++) {
    if (bdata -> lower_bc_elc[d]) {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_elc[d]);
    }
    if (bdata -> lower_bc_ion[d]) {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_ion[d]);
    }
    if (bdata -> lower_bc_maxwell[d]) {
      gkyl_wv_apply_bc_release(bdata -> lower_bc_maxwell[d]);
    }

    if (bdata -> upper_bc_elc[d]) {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_elc[d]);
    }
    if (bdata -> upper_bc_ion[d]) {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_ion[d]);
    }
    if (bdata -> upper_bc_maxwell[d]) {
      gkyl_wv_apply_bc_release(bdata -> upper_bc_maxwell[d]);
    }
  }

  gkyl_array_release(bdata -> bc_buffer_elc);
  gkyl_array_release(bdata -> bc_buffer_ion);
  gkyl_array_release(bdata -> bc_buffer_maxwell);
}

void
five_moment_block_bc_updaters_apply(const struct five_moment_block_data* bdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array* fld_ion, struct gkyl_array* fld_maxwell)
{
  for (int d = 0; d < 2; d++) {
    if (bdata -> lower_bc_elc[d]) {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_elc[d], tm, &bdata -> range, fld_elc);
    }
    if (bdata -> lower_bc_ion[d]) {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_ion[d], tm, &bdata -> range, fld_ion);
    }
    if (bdata -> lower_bc_maxwell[d]) {
      gkyl_wv_apply_bc_advance(bdata -> lower_bc_maxwell[d], tm, &bdata -> range, fld_maxwell);
    }

    if (bdata -> upper_bc_elc[d]) {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_elc[d], tm, &bdata -> range, fld_elc);
    }
    if (bdata -> upper_bc_ion[d]) {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_ion[d], tm, &bdata -> range, fld_ion);
    }
    if (bdata -> upper_bc_maxwell[d]) {
      gkyl_wv_apply_bc_advance(bdata -> upper_bc_maxwell[d], tm, &bdata -> range, fld_maxwell);
    }
  }
}

void
five_moment_sync_blocks(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  int num_blocks = btopo -> num_blocks;
  int ndim = btopo -> ndim;

  for (int i = 0; i < num_blocks; i++) {
    for (int d = 0; d < ndim; d++) {
      const struct gkyl_target_edge *te = btopo -> conn[i].connections[d];

#ifdef AMR_DEBUG
      if (te[0].edge != GKYL_PHYSICAL) {
        struct gkyl_array* bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array* bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array* bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc -> data, fld_elc[i], bdata[i].skin_ghost.lower_skin[d]);
        gkyl_array_copy_to_buffer(bc_buffer_ion -> data, fld_ion[i], bdata[i].skin_ghost.lower_skin[d]);
        gkyl_array_copy_to_buffer(bc_buffer_maxwell -> data, fld_maxwell[i], bdata[i].skin_ghost.lower_skin[d]);

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (te[0].edge == GKYL_UPPER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }

      if (te[1].edge != GKYL_PHYSICAL) {
        struct gkyl_array* bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array* bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array* bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc -> data, fld_elc[i], bdata[i].skin_ghost.upper_skin[d]);
        gkyl_array_copy_to_buffer(bc_buffer_ion -> data, fld_ion[i], bdata[i].skin_ghost.upper_skin[d]);
        gkyl_array_copy_to_buffer(bc_buffer_maxwell -> data, fld_maxwell[i], bdata[i].skin_ghost.upper_skin[d]);

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (te[1].edge == GKYL_UPPER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell -> data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }
#endif
    }
  }
}