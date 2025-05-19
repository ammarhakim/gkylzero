#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Maxwell equation object.
 * 
 * @param c_speed Speed of light
 * @param e_fact Factor of light-speed for electric field correction
 * @param b_fact Factor of light-speed for magnetic field correction
 * @return Pointer to Maxwell equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_maxwell_new(double c, double e_fact, double b_fact,
  struct gkyl_wv_embed_geo* embed_geo, bool use_gpu);

/**
 * Create a new Maxwell equation object that lives on NV-GPU.
 * see new() method above for documentation.
 */
struct gkyl_wv_eqn* gkyl_wv_maxwell_cu_dev_new(double c, double e_fact, double b_fact);
