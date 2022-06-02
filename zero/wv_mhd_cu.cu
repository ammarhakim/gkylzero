/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_mhd.h>    
#include <gkyl_wv_mhd_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to isothermal mhd kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_mhd_set_cu_dev_ptrs(struct wv_mhd *mhd)
{
  mhd->eqn.waves_func = wave_roe;
  mhd->eqn.qfluct_func = qfluct_roe;
  mhd->eqn.max_speed_func = max_speed;
  mhd->eqn.rotate_to_local_func = rot_to_local_rect;
  mhd->eqn.rotate_to_global_func = rot_to_global_rect;

  // mhd->eqn.wall_bc_func = mhd_wall; // TODO
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_cu_dev_new(double gas_gamma, const char *divergence_constraint)
{
  struct wv_mhd *mhd = (struct wv_mhd*) gkyl_malloc(sizeof(struct wv_mhd));

  /* STEP 1. SET PRIMITIVE DATA */
  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->gas_gamma = gas_gamma;

  if (strcmp(divergence_constraint, "none") == 0)
  {
    mhd->divergence_constraint = DIVB_NONE;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;
  }
  else if (strcmp(divergence_constraint, "eight_waves") == 0)
  {
    mhd->divergence_constraint = DIVB_EIGHT_WAVES;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;  // will merge the entropy wave and divB wave
  }
  else if (strcmp(divergence_constraint, "glm") == 0)
  {
    mhd->divergence_constraint = DIVB_GLM;
    mhd->eqn.num_equations = 9;
    mhd->eqn.num_waves = 9;
    mhd->eqn.rotate_to_local_func = rot_to_local_rect_glm;
    mhd->eqn.rotate_to_global_func = rot_to_global_rect_glm;
  }
  else {
    // Do not constrain divergence by default. TODO: Warn or throw an error
    mhd->divergence_constraint = DIVB_NONE;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;
  }

  mhd->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(mhd->eqn.flags);
  mhd->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_mhd_free);

  /* STEP 2. COPY HOST OBJECT TO DEVICE OBJECT */
  struct wv_mhd *mhd_cu = (struct wv_mhd*) gkyl_cu_malloc(sizeof(struct wv_mhd));
  gkyl_cu_memcpy(mhd_cu, mhd, sizeof(struct wv_mhd), GKYL_CU_MEMCPY_H2D);

  /* STEP 3. SET DEVICE FUNCTION POINTERS IN DEVICE OBJECT */
  wv_mhd_set_cu_dev_ptrs<<<1,1>>>(mhd_cu);

  mhd->eqn.on_dev = &mhd_cu->eqn; // CPU eqn obj points to itself
  return &mhd->eqn;
}
