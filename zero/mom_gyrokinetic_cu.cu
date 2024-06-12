/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>
}

enum { M0, M1, M2, M2par, M2perp, M3par, M3perp, ThreeMoments, BAD };

static int
get_gk_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "M0") == 0) { // density
    mom_idx = M0;
  }
  else if (strcmp(mom, "M1") == 0) { // parallel momentum
    mom_idx = M1;
  }
  else if (strcmp(mom, "M2") == 0) { // total energy
    mom_idx = M2;
  }
  else if (strcmp(mom, "M2par") == 0) { // parallel energy
    mom_idx = M2par;
  }
  else if (strcmp(mom, "M2perp") == 0) { // perpendicular energy
    mom_idx = M2perp;
  }
  else if (strcmp(mom, "M3par") == 0) { // parallel heat flux
    mom_idx = M3par;
  }
  else if (strcmp(mom, "M3perp") == 0) { // perpendicular heat flux
    mom_idx = M3perp;
  }
  else if (strcmp(mom, "ThreeMoments") == 0) { // Zeroth (density), First (parallel momentum), 
    mom_idx = ThreeMoments;                    // and Second (total energy) computed together
  }
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

static int
gk_num_mom(int vdim, int mom_id)
{
  int num_mom = 0;
  
  switch (mom_id) {
    case M0:
      num_mom = 1;
      break;

    case M1:
      num_mom = 1;
      break;

    case M2:
      num_mom = 1;
      break;

    case M2par:
      num_mom = 1;
      break;

    case M2perp:
      num_mom = 1;
      break;

    case M3par:
      num_mom = 1;
      break;

    case M3perp:
      num_mom = 1;
      break;

    case ThreeMoments:
      num_mom = 3;
      break;      
      
    default: // can't happen
      break;
  }

  return num_mom;
}

__global__
static void
set_cu_ptrs(struct mom_type_gyrokinetic *mom_gk,
  int mom_id, enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  // choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, 
    *m2_par_kernels, *m2_perp_kernels, *m3_par_kernels, *m3_perp_kernels, *three_moments_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1_kernels = ser_m1_kernels;
      m2_kernels = ser_m2_kernels;
      m2_par_kernels = ser_m2_par_kernels;
      m2_perp_kernels = ser_m2_perp_kernels;
      m3_par_kernels = ser_m3_par_kernels;
      m3_perp_kernels = ser_m3_perp_kernels;
      three_moments_kernels = ser_three_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_id) {
    case M0:
      mom_gk->momt.kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M1:
      mom_gk->momt.kernel = m1_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M2:
      mom_gk->momt.kernel = m2_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M2par:
      mom_gk->momt.kernel = m2_par_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M2perp:
      mom_gk->momt.kernel = m2_perp_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M3par:
      mom_gk->momt.kernel = m3_par_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case M3perp:
      mom_gk->momt.kernel = m3_perp_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case ThreeMoments:
      mom_gk->momt.kernel = three_moments_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 3;
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_gyrokinetic *mom_gk = (struct mom_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gk->momt.cdim = cdim;
  mom_gk->momt.pdim = pdim;
  mom_gk->momt.poly_order = poly_order;
  mom_gk->momt.num_config = cbasis->num_basis;
  mom_gk->momt.num_phase = pbasis->num_basis;

  int mom_id = get_gk_mom_id(mom);
  if(mom_id == BAD) {
     printf("Error: requested GK moment %s not valid\n", mom);
     assert(mom_id != BAD);
  }
  mom_gk->momt.num_mom = gk_num_mom(vdim, mom_id); // number of moments

  mom_gk->mass = mass;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  mom_gk->gk_geom = geom_ho->on_dev;
  mom_gk->vel_map = vel_map_ho->on_dev;

  mom_gk->conf_range = *conf_range;

  mom_gk->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_gk->momt.flags);
  mom_gk->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // copy struct to device
  struct mom_type_gyrokinetic *mom_gk_cu = (struct mom_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_gyrokinetic));
  gkyl_cu_memcpy(mom_gk_cu, mom_gk, sizeof(struct mom_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(mom_gk_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_gk->momt.on_dev = &mom_gk_cu->momt;

  // Updater should store host pointers.
  mom_gk->gk_geom = geom_ho; 
  mom_gk->vel_map = vel_map_ho; 
  
  return &mom_gk->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_gyrokinetic* momt, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  // set kernel pointer
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      momt->momt.kernel = ser_int_mom_kernels[tblidx].kernels[poly_order];
      break;

    default:
      assert(false);
      break;    
  }
}

struct gkyl_mom_type*
gkyl_int_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_gyrokinetic *momt = (struct mom_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->momt.num_mom = vdim+2;

  momt->mass = mass;
  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  momt->gk_geom = geom_ho->on_dev;
  momt->vel_map = vel_map_ho->on_dev;

  momt->conf_range = *conf_range;

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // copy struct to device
  struct mom_type_gyrokinetic *momt_cu = (struct mom_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_gyrokinetic));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;

  // Updater should store host pointers.
  momt->gk_geom = geom_ho; 
  momt->vel_map = vel_map_ho; 
  
  return &momt->momt;
}
