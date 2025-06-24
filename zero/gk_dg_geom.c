#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>

#include <gkyl_gk_geometry.h>
#include <gkyl_gk_dg_geom.h>

static bool
gk_dg_geom_is_cu_dev(const struct gkyl_gk_dg_geom* dgg)
{
  return GKYL_IS_CU_ALLOC(dgg->flags);
}

static void
gk_dg_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_gk_dg_geom *dgg = container_of(ref, struct gkyl_gk_dg_geom, ref_count);

  for (int d=0; d<dgg->range.ndim; ++d)
    gkyl_array_release(dgg->surf_geom[d]);

  gkyl_array_release(dgg->vol_geom);

  if (gk_dg_geom_is_cu_dev(dgg)) 
    gkyl_cu_free(dgg->on_dev); 

  gkyl_free(dgg);
}


struct gkyl_gk_dg_geom *
gkyl_gk_dg_geom_new(const struct gkyl_gk_dg_geom_inp *inp)
{
  struct gkyl_gk_dg_geom *dgg = gkyl_malloc(sizeof *dgg);

  dgg->range = *inp->range;

  int ndim = dgg->range.ndim;
  int shape[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; ++d) shape[d] = inp->nquad;

  // NOTE: surfaces are ndim-1 objects
  gkyl_range_init_from_shape(&dgg->surf_quad_range, ndim-1, shape);
  gkyl_range_init_from_shape(&dgg->vol_quad_range, ndim, shape);

  for (int d=0; d<ndim; ++d)
    dgg->surf_geom[d] = gkyl_array_new(GKYL_USER,
      sizeof(struct gkyl_dg_surf_geom[dgg->surf_quad_range.volume]), dgg->range.volume);

  dgg->vol_geom = gkyl_array_new(GKYL_USER,
    sizeof(struct gkyl_gk_dg_vol_geom[dgg->vol_quad_range.volume]), dgg->range.volume);

  
  dgg->flags = 0;
  GKYL_CLEAR_CU_ALLOC(dgg->flags);
  dgg->ref_count = gkyl_ref_count_init(gk_dg_geom_free);
  dgg->on_dev = dgg; // CPU eqn obj points to itself
  
  return dgg;
}


struct gkyl_gk_dg_geom*
gkyl_gk_dg_geom_acquire(const struct gkyl_gk_dg_geom* dgg)
{
  gkyl_ref_count_inc(&dgg->ref_count);
  return (struct gkyl_gk_dg_geom*) dgg;
}

void
gkyl_gk_dg_geom_write(const struct gkyl_gk_dg_geom* dgg, const char *fname)
{
  
}

void
gkyl_gk_dg_geom_release(const struct gkyl_gk_dg_geom *dgg)
{
  gkyl_ref_count_dec(&dgg->ref_count);
}


void
gkyl_gk_dg_geom_populate_vol(struct gkyl_dg_geom *dg_geom, struct gkyl_gk_dg_geom *gk_dg_geom, struct gk_geometry* gk_geom)
{
  int ndim = gk_geom->grid.ndim;
  // Populate volume nodes
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &gk_geom->local);
  while(gkyl_range_iter_next(&iter)){

    long loc = gkyl_range_idx(&gk_geom->local, iter.idx);
    struct gkyl_dg_vol_geom *dgv = gkyl_array_fetch(dg_geom->vol_geom, loc);
    struct gkyl_gk_dg_vol_geom *gkdgv = gkyl_array_fetch(gk_dg_geom->vol_geom, loc);

    struct gkyl_range_iter qviter;
    gkyl_range_iter_init(&qviter, &dg_geom->vol_quad_range);
    int global_nodal_idx[ndim];
    while(gkyl_range_iter_next(&qviter)){
      for (int d=0; d<ndim; ++d)
        global_nodal_idx[d] = (iter.idx[d]-1)*2 + qviter.idx[d];
      long qvloc = gkyl_range_idx(&dg_geom->vol_quad_range, qviter.idx);
      long global_loc = gkyl_range_idx(&gk_geom->nrange_int, global_nodal_idx);

      // set Jc
      const double *global_val = gkyl_array_cfetch(gk_geom->geo_int.jacobgeo_nodal, global_loc);
      dgv[qvloc].Jc = global_val[0];

      // set tangents
      global_val = gkyl_array_cfetch(gk_geom->geo_int.dxdz_nodal, global_loc);
      dgv[qvloc].tang[0].x[0] = global_val[0];
      dgv[qvloc].tang[0].x[1] = global_val[1];
      dgv[qvloc].tang[0].x[2] = global_val[2];
      dgv[qvloc].tang[1].x[0] = global_val[3];
      dgv[qvloc].tang[1].x[1] = global_val[4];
      dgv[qvloc].tang[1].x[2] = global_val[5];
      dgv[qvloc].tang[2].x[0] = global_val[6];
      dgv[qvloc].tang[2].x[1] = global_val[7];
      dgv[qvloc].tang[2].x[2] = global_val[8];

      // set duals
      global_val = gkyl_array_cfetch(gk_geom->geo_int.dzdx_nodal, global_loc);
      dgv[qvloc].dual[0].x[0] = global_val[0];
      dgv[qvloc].dual[0].x[1] = global_val[1];
      dgv[qvloc].dual[0].x[2] = global_val[2];
      dgv[qvloc].dual[1].x[0] = global_val[3];
      dgv[qvloc].dual[1].x[1] = global_val[4];
      dgv[qvloc].dual[1].x[2] = global_val[5];
      dgv[qvloc].dual[2].x[0] = global_val[6];
      dgv[qvloc].dual[2].x[1] = global_val[7];
      dgv[qvloc].dual[2].x[2] = global_val[8];

      // set bmag
      global_val = gkyl_array_cfetch(gk_geom->geo_int.bmag_nodal, global_loc);
      gkdgv[qvloc].bmag = global_val[0];

      // set bmag
      global_val = gkyl_array_cfetch(gk_geom->geo_int.bmag_nodal, global_loc);
      gkdgv[qvloc].bmag = global_val[0];

      // set B3 = e^3 \dot B
      global_val = gkyl_array_cfetch(gk_geom->geo_int.B3_nodal, global_loc);
      gkdgv[qvloc].B3= global_val[0];

      // set e^i \dot curl(bhat)
      global_val = gkyl_array_cfetch(gk_geom->geo_int.dualcurlbhat_nodal, global_loc);
      gkdgv[qvloc].dualcurlbhat.x[0] = global_val[0];
      gkdgv[qvloc].dualcurlbhat.x[1] = global_val[1];
      gkdgv[qvloc].dualcurlbhat.x[2] = global_val[2];
    }
  }
}

void
gkyl_gk_dg_geom_populate_surf(struct gkyl_dg_geom *dg_geom, struct gkyl_gk_dg_geom *gk_dg_geom, struct gk_geometry* gk_geom)
{
  int ndim = gk_geom->grid.ndim;
  // Populate surface nodes
  for (int dir=0; dir<ndim; dir++) {

    struct gkyl_range local_ext_in_dir;
    int lower[3] = {gk_geom->local.lower[0], gk_geom->local.lower[1], gk_geom->local.lower[2]};
    int upper[3] = {gk_geom->local.upper[0], gk_geom->local.upper[1], gk_geom->local.upper[2]};
    upper[dir]+=1;
    gkyl_sub_range_init(&local_ext_in_dir, &gk_geom->local_ext, lower, upper);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local_ext_in_dir);
    while(gkyl_range_iter_next(&iter)){

      long loc = gkyl_range_idx(&local_ext_in_dir, iter.idx);
      struct gkyl_dg_surf_geom *dgs = gkyl_array_fetch(dg_geom->surf_geom[dir], loc);
      struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_array_fetch(gk_dg_geom->surf_geom[dir], loc);


      struct gkyl_range_iter qsiter;
      gkyl_range_iter_init(&qsiter, &dg_geom->surf_quad_range);
      int global_nodal_idx[ndim];
      while(gkyl_range_iter_next(&qsiter)){
        int count = 0;
        for (int d=0; d<ndim; ++d) {
          global_nodal_idx[d] = d == dir ? iter.idx[d]-1 : (iter.idx[d]-1)*2 + qsiter.idx[count];
          if (d != dir) count+=1;
        }
        long qsloc = gkyl_range_idx(&dg_geom->surf_quad_range, qsiter.idx);
        long global_loc = gkyl_range_idx(&gk_geom->nrange_surf[dir], global_nodal_idx);

        // set lenr
        const double *global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].lenr_nodal, global_loc);
        dgs[qsloc].area_elem = global_val[0];

        // set normals
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].normals_nodal, global_loc);
        dgs[qsloc].norm.x[0] = global_val[dir*3+0];
        dgs[qsloc].norm.x[1] = global_val[dir*3+1];
        dgs[qsloc].norm.x[2] = global_val[dir*3+2];
        
        // set e3hat \dot B = B^3/sqrt(g_33}
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].B3_nodal, global_loc);
        gkdgs[qsloc].B3  = global_val[0];
        // set n \dot curl(bhat)
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].normcurlbhat_nodal, global_loc);
        gkdgs[qsloc].normcurlbhat  = global_val[dir];

        // set |B|
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].bmag_nodal, global_loc);
        gkdgs[qsloc].bmag  = global_val[0];

        // set Jacobgeo
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].jacobgeo_nodal, global_loc);
        gkdgs[qsloc].Jc = global_val[0];

        // set bhat
        global_val = gkyl_array_cfetch(gk_geom->geo_surf[dir].b_i_nodal, global_loc);
        gkdgs[qsloc].bhat.x[0]  = global_val[0];
        gkdgs[qsloc].bhat.x[1]  = global_val[1];
        gkdgs[qsloc].bhat.x[2]  = global_val[2];

      }
    }
  }
}

void
gkyl_gk_dg_geom_write_vol(struct gkyl_dg_geom *dg_geom, struct gk_geometry* gk_geom, const char *name)
{
  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM];
  double upper[GKYL_MAX_DIM];
  for (int i=0; i<gk_geom->nrange_int.ndim; ++i) {
    lower[i] = gk_geom->grid.lower[i];
    upper[i] = gk_geom->grid.upper[i];
    cells[i] = gkyl_range_shape(&gk_geom->nrange_int, i);
  }
  struct gkyl_rect_grid ngrid;
  gkyl_rect_grid_init(&ngrid, gk_geom->nrange_int.ndim, lower, upper, cells);

  struct gkyl_array* nodal_j = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  struct gkyl_array* nodal_dxdz = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_int.volume);
  struct gkyl_array* nodal_dzdx = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_int.volume);


  int ndim = gk_geom->grid.ndim;
  // Populate volume nodes
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &gk_geom->local);
  while(gkyl_range_iter_next(&iter)){

    long loc = gkyl_range_idx(&gk_geom->local, iter.idx);
    struct gkyl_dg_vol_geom *dgv = gkyl_array_fetch(dg_geom->vol_geom, loc);

    struct gkyl_range_iter qviter;
    gkyl_range_iter_init(&qviter, &dg_geom->vol_quad_range);
    int global_nodal_idx[ndim];
    while(gkyl_range_iter_next(&qviter)){
      for (int d=0; d<ndim; ++d)
        global_nodal_idx[d] = (iter.idx[d]-1)*2 + qviter.idx[d];
      long qvloc = gkyl_range_idx(&dg_geom->vol_quad_range, qviter.idx);
      long global_loc = gkyl_range_idx(&gk_geom->nrange_int, global_nodal_idx);

      // set Jc
      double *global_val = gkyl_array_fetch(nodal_j, global_loc);
      global_val[0]= dgv[qvloc].Jc ;

      // set tangents
      global_val = gkyl_array_fetch(nodal_dxdz, global_loc);
      global_val[0] = dgv[qvloc].tang[0].x[0] ;
      global_val[1] = dgv[qvloc].tang[0].x[1] ;
      global_val[2] = dgv[qvloc].tang[0].x[2] ;
      global_val[3] = dgv[qvloc].tang[1].x[0] ;
      global_val[4] = dgv[qvloc].tang[1].x[1] ;
      global_val[5] = dgv[qvloc].tang[1].x[2] ;
      global_val[6] = dgv[qvloc].tang[2].x[0] ;
      global_val[7] = dgv[qvloc].tang[2].x[1] ;
      global_val[8] = dgv[qvloc].tang[2].x[2] ;

      // set duals
      global_val = gkyl_array_fetch(nodal_dzdx, global_loc);
      global_val[0] = dgv[qvloc].dual[0].x[0];
      global_val[1] = dgv[qvloc].dual[0].x[1];
      global_val[2] = dgv[qvloc].dual[0].x[2];
      global_val[3] = dgv[qvloc].dual[1].x[0];
      global_val[4] = dgv[qvloc].dual[1].x[1];
      global_val[5] = dgv[qvloc].dual[1].x[2];
      global_val[6] = dgv[qvloc].dual[2].x[0];
      global_val[7] = dgv[qvloc].dual[2].x[1];
      global_val[8] = dgv[qvloc].dual[2].x[2];
    }
  }

  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "tang_vol");
  char fileNm[sz+1]; // ensure no buffer overflow
  sprintf(fileNm, fmt, name, "Jc_vol");
  gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_int, 0,  nodal_j, fileNm);
  sprintf(fileNm, fmt, name, "tang_vol");
  gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_int, 0,  nodal_dxdz, fileNm);
  sprintf(fileNm, fmt, name, "dual_vol");
  gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_int, 0,  nodal_dzdx, fileNm);

  gkyl_array_release(nodal_j);
  gkyl_array_release(nodal_dxdz);
  gkyl_array_release(nodal_dzdx);

}

void
gkyl_gk_dg_geom_write_surf(struct gkyl_dg_geom *dg_geom, struct gkyl_gk_dg_geom *gk_dg_geom, struct gk_geometry* gk_geom, const char *name)
{



  int ndim = gk_geom->grid.ndim;
  // Populate surface nodes
  for (int dir=0; dir<ndim; dir++) {

    int gcells[GKYL_MAX_DIM];
    double glower[GKYL_MAX_DIM];
    double gupper[GKYL_MAX_DIM];
    for (int i=0; i<ndim; ++i) {
      glower[i] = gk_geom->grid.lower[i];
      gupper[i] = gk_geom->grid.upper[i];
      gcells[i] = gkyl_range_shape(&gk_geom->nrange_surf[dir], i);
    }
    struct gkyl_rect_grid ngrid;
    gkyl_rect_grid_init(&ngrid, gk_geom->nrange_surf[dir].ndim, glower, gupper, gcells);

    struct gkyl_range local_ext_in_dir;
    int lower[3] = {gk_geom->local.lower[0], gk_geom->local.lower[1], gk_geom->local.lower[2]};
    int upper[3] = {gk_geom->local.upper[0], gk_geom->local.upper[1], gk_geom->local.upper[2]};
    upper[dir]+=1;
    gkyl_sub_range_init(&local_ext_in_dir, &gk_geom->local_ext, lower, upper);

    struct gkyl_array* nodal_lenr = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
    struct gkyl_array* nodal_normals = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
    struct gkyl_array* nodal_B3 = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
    struct gkyl_array* nodal_normcurlbhat = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
    struct gkyl_array* nodal_bmag = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
    struct gkyl_array* nodal_b_i = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);


    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local_ext_in_dir);
    while(gkyl_range_iter_next(&iter)){

      long loc = gkyl_range_idx(&local_ext_in_dir, iter.idx);
      struct gkyl_dg_surf_geom *dgs = gkyl_array_fetch(dg_geom->surf_geom[dir], loc);
      struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_array_fetch(gk_dg_geom->surf_geom[dir], loc);


      struct gkyl_range_iter qsiter;
      gkyl_range_iter_init(&qsiter, &dg_geom->surf_quad_range);
      int global_nodal_idx[ndim];
      while(gkyl_range_iter_next(&qsiter)){
        int count = 0;
        for (int d=0; d<ndim; ++d) {
          global_nodal_idx[d] = d == dir ? iter.idx[d]-1 : (iter.idx[d]-1)*2 + qsiter.idx[count];
          if (d != dir) count+=1;
        }
        long qsloc = gkyl_range_idx(&dg_geom->surf_quad_range, qsiter.idx);
        long global_loc = gkyl_range_idx(&gk_geom->nrange_surf[dir], global_nodal_idx);

        // set lenr
        double *global_val = gkyl_array_fetch(nodal_lenr, global_loc);
        global_val[0] = dgs[qsloc].area_elem ;

        // set normals
        global_val = gkyl_array_fetch(nodal_normals, global_loc);
        global_val[0] = dgs[qsloc].norm.x[0];
        global_val[1] = dgs[qsloc].norm.x[1];
        global_val[2] = dgs[qsloc].norm.x[2];

        // set e3hat \dot B = B^3/sqrt(g_33}
        global_val = gkyl_array_fetch(nodal_B3, global_loc);
        global_val[0] = gkdgs[qsloc].B3 ;
        // set n \dot curl(bhat)
        global_val = gkyl_array_fetch(nodal_normcurlbhat, global_loc);
        global_val[0] = gkdgs[qsloc].normcurlbhat ;

        // set |B|
        global_val = gkyl_array_fetch(nodal_bmag, global_loc);
        global_val[0] = gkdgs[qsloc].bmag ;

        // set bhat
        global_val = gkyl_array_fetch(nodal_b_i, global_loc);
        global_val[0] = gkdgs[qsloc].bhat.x[0];
        global_val[1] = gkdgs[qsloc].bhat.x[1];
        global_val[2] = gkdgs[qsloc].bhat.x[2];
      }
    }



    const char *fmt = "%s-%s_dir%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, name, "normals_surf", dir);
    char fileNm[sz+1]; // ensure no buffer overflow
    sprintf(fileNm, fmt, name, "lenr_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_lenr, fileNm);
    sprintf(fileNm, fmt, name, "normals_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_normals, fileNm);
    sprintf(fileNm, fmt, name, "B3_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_B3, fileNm);
    sprintf(fileNm, fmt, name, "normcurlbhat_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_normcurlbhat, fileNm);
    sprintf(fileNm, fmt, name, "bmag_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_bmag, fileNm);
    sprintf(fileNm, fmt, name, "bhat_surf", dir);
    gkyl_grid_sub_array_write(&ngrid, &gk_geom->nrange_surf[dir], 0,  nodal_b_i, fileNm);

    gkyl_array_release(nodal_lenr);
    gkyl_array_release(nodal_normals);
  }

}
