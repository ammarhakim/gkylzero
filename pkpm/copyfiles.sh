#!/bin/sh

CP_CMD=git mv
RM_CMD=:
G0=..

# kernels
mkdir -p kernels/pkpm
$CP_CMD $G0/kernels/pkpm/*.h kernels/pkpm/
$CP_CMD $G0/kernels/pkpm/*.c kernels/pkpm/

$RM_CMD $G0/kernels/pkpm/*.h
$RM_CMD $G0/kernels/pkpm/*.c

# zero
mkdir -p zero
$CP_CMD $G0/zero/dg_calc_pkpm_dist_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_pkpm_dist_vars.c zero/
$CP_CMD $G0/zero/dg_calc_pkpm_em_coupling_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_pkpm_em_coupling.c zero/
$CP_CMD $G0/zero/dg_calc_pkpm_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_pkpm_vars.c zero/
$CP_CMD $G0/zero/dg_euler_pkpm_cu.cu zero/
$CP_CMD $G0/zero/dg_euler_pkpm.c zero/
$CP_CMD $G0/zero/dg_lbo_pkpm_diff_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_pkpm_diff.c zero/
$CP_CMD $G0/zero/dg_lbo_pkpm_drag_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_pkpm_drag.c zero/
$CP_CMD $G0/zero/dg_updater_lbo_pkpm.c zero/
$CP_CMD $G0/zero/dg_updater_moment_pkpm.c zero/
$CP_CMD $G0/zero/dg_updater_pkpm.c zero/
$CP_CMD $G0/zero/dg_vlasov_pkpm_cu.cu zero/
$CP_CMD $G0/zero/dg_vlasov_pkpm.c zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_dist_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_dist_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_em_coupling_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_em_coupling.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_pkpm_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_euler_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_euler_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_pkpm_diff_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_pkpm_diff.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_pkpm_drag_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_pkpm_drag.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_lbo_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_moment_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_mom_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_pkpm.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_pkpm_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_pkpm.h zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_pkpm_cu.cu zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_pkpm.c zero/
$CP_CMD $G0/zero/mom_pkpm_cu.cu zero/
$CP_CMD $G0/zero/mom_pkpm.c zero/
$CP_CMD $G0/zero/prim_lbo_pkpm_cu.cu zero/
$CP_CMD $G0/zero/prim_lbo_pkpm.c zero/
$CP_CMD $G0/zero/mom_calc_bcorr_pkpm.c zero/
$CP_CMD $G0/zero/prim_lbo_calc_pkpm.c zero/

$RM_CMD $G0/zero/dg_calc_pkpm_dist_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_pkpm_dist_vars.c
$RM_CMD $G0/zero/dg_calc_pkpm_em_coupling_cu.cu
$RM_CMD $G0/zero/dg_calc_pkpm_em_coupling.c
$RM_CMD $G0/zero/dg_calc_pkpm_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_pkpm_vars.c
$RM_CMD $G0/zero/dg_euler_pkpm_cu.cu
$RM_CMD $G0/zero/dg_euler_pkpm.c
$RM_CMD $G0/zero/dg_lbo_pkpm_diff_cu.cu
$RM_CMD $G0/zero/dg_lbo_pkpm_diff.c
$RM_CMD $G0/zero/dg_lbo_pkpm_drag_cu.cu
$RM_CMD $G0/zero/dg_lbo_pkpm_drag.c
$RM_CMD $G0/zero/dg_updater_lbo_pkpm.c
$RM_CMD $G0/zero/dg_updater_moment_pkpm.c
$RM_CMD $G0/zero/dg_updater_pkpm.c
$RM_CMD $G0/zero/dg_vlasov_pkpm_cu.cu
$RM_CMD $G0/zero/dg_vlasov_pkpm.c
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_dist_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_dist_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_em_coupling_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_em_coupling.h
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_pkpm_vars.h
$RM_CMD $G0/zero/gkyl_dg_euler_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_dg_euler_pkpm.h
$RM_CMD $G0/zero/gkyl_dg_lbo_pkpm_diff_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_pkpm_diff.h
$RM_CMD $G0/zero/gkyl_dg_lbo_pkpm_drag_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_pkpm_drag.h
$RM_CMD $G0/zero/gkyl_dg_updater_lbo_pkpm.h
$RM_CMD $G0/zero/gkyl_dg_updater_moment_pkpm.h
$RM_CMD $G0/zero/gkyl_dg_updater_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_pkpm.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_pkpm.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_pkpm.h
$RM_CMD $G0/zero/gkyl_mom_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_mom_pkpm.h
$RM_CMD $G0/zero/gkyl_prim_lbo_pkpm_priv.h
$RM_CMD $G0/zero/gkyl_prim_lbo_pkpm.h
$RM_CMD $G0/zero/mom_bcorr_lbo_pkpm_cu.cu
$RM_CMD $G0/zero/mom_bcorr_lbo_pkpm.c
$RM_CMD $G0/zero/mom_pkpm_cu.cu
$RM_CMD $G0/zero/mom_pkpm.c
$RM_CMD $G0/zero/prim_lbo_pkpm_cu.cu
$RM_CMD $G0/zero/prim_lbo_pkpm.c
$RM_CMD $G0/zero/mom_calc_bcorr_pkpm.c
$RM_CMD $G0/zero/prim_lbo_calc_pkpm.c

# app
mkdir -p apps
$CP_CMD $G0/apps/gkyl_pkpm_lw.h apps/
$CP_CMD $G0/apps/gkyl_pkpm_priv.h apps/
$CP_CMD $G0/apps/gkyl_pkpm.h apps/
$CP_CMD $G0/apps/pkpm_field.c apps/
$CP_CMD $G0/apps/pkpm_fluid_em_coupling.c apps/
$CP_CMD $G0/apps/pkpm_forward_euler.c apps/
$CP_CMD $G0/apps/pkpm_lw.c apps/
$CP_CMD $G0/apps/pkpm_species_lbo.c apps/
$CP_CMD $G0/apps/pkpm_species_moment.c apps/
$CP_CMD $G0/apps/pkpm_species.c apps/
$CP_CMD $G0/apps/pkpm_update_explicit_ssp_rk3.c apps/
$CP_CMD $G0/apps/pkpm_update_op_split.c apps/
$CP_CMD $G0/apps/pkpm.c apps/

$RM_CMD $G0/apps/gkyl_pkpm_lw.h
$RM_CMD $G0/apps/gkyl_pkpm_priv.h
$RM_CMD $G0/apps/gkyl_pkpm.h
$RM_CMD $G0/apps/pkpm_field.c
$RM_CMD $G0/apps/pkpm_fluid_em_coupling.c
$RM_CMD $G0/apps/pkpm_forward_euler.c
$RM_CMD $G0/apps/pkpm_lw.c
$RM_CMD $G0/apps/pkpm_species_lbo.c
$RM_CMD $G0/apps/pkpm_species_moment.c
$RM_CMD $G0/apps/pkpm_species.c
$RM_CMD $G0/apps/pkpm_update_explicit_ssp_rk3.c
$RM_CMD $G0/apps/pkpm_update_op_split.c
$RM_CMD $G0/apps/pkpm.c

# unit
mkdir -p unit

# C regression tests
mkdir -p creg
$CP_CMD $G0/regression/rt_pkpm_2d_travel_pulse_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_soliton_1x_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_wave_1x_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_wave_1x_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_wave_explicit_1x_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_wave_explicit_1x_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_alf_wave_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_em_advect_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_em_advect_resonant_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_es_pot_well_1x_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_es_pot_well_1x_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_es_shock_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_es_shock_reflect_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_landau_damping_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_landau_damping_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_mom_beach_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_neut_sodshock_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_neut_sodshock_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_periodic_es_shock_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_periodic_es_shock_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_periodic_neut_sodshock_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_periodic_neut_sodshock_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_sheath_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_square_relax_1x_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_square_relax_1x_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_travel_pulse_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_travel_pulse_p2.c creg/
$CP_CMD $G0/regression/rt_pkpm_wall_p1.c creg/
$CP_CMD $G0/regression/rt_pkpm_wall_p2.c creg/
$CP_CMD $G0/regression/rt_arg_parse.h creg/

$RM_CMD $G0/regression/rt_pkpm_2d_travel_pulse_p1.c
$RM_CMD $G0/regression/rt_pkpm_alf_soliton_1x_p2.c
$RM_CMD $G0/regression/rt_pkpm_alf_wave_1x_p1.c
$RM_CMD $G0/regression/rt_pkpm_alf_wave_1x_p2.c
$RM_CMD $G0/regression/rt_pkpm_alf_wave_explicit_1x_p1.c
$RM_CMD $G0/regression/rt_pkpm_alf_wave_explicit_1x_p2.c
$RM_CMD $G0/regression/rt_pkpm_alf_wave_p1.c
$RM_CMD $G0/regression/rt_pkpm_em_advect_p1.c
$RM_CMD $G0/regression/rt_pkpm_em_advect_resonant_p1.c
$RM_CMD $G0/regression/rt_pkpm_es_pot_well_1x_p1.c
$RM_CMD $G0/regression/rt_pkpm_es_pot_well_1x_p2.c
$RM_CMD $G0/regression/rt_pkpm_es_shock_p2.c
$RM_CMD $G0/regression/rt_pkpm_es_shock_reflect_p2.c
$RM_CMD $G0/regression/rt_pkpm_landau_damping_p1.c
$RM_CMD $G0/regression/rt_pkpm_landau_damping_p2.c
$RM_CMD $G0/regression/rt_pkpm_mom_beach_p2.c
$RM_CMD $G0/regression/rt_pkpm_neut_sodshock_p1.c
$RM_CMD $G0/regression/rt_pkpm_neut_sodshock_p2.c
$RM_CMD $G0/regression/rt_pkpm_periodic_es_shock_p1.c
$RM_CMD $G0/regression/rt_pkpm_periodic_es_shock_p2.c
$RM_CMD $G0/regression/rt_pkpm_periodic_neut_sodshock_p1.c
$RM_CMD $G0/regression/rt_pkpm_periodic_neut_sodshock_p2.c
$RM_CMD $G0/regression/rt_pkpm_sheath_p1.c
$RM_CMD $G0/regression/rt_pkpm_square_relax_1x_p1.c
$RM_CMD $G0/regression/rt_pkpm_square_relax_1x_p2.c
$RM_CMD $G0/regression/rt_pkpm_travel_pulse_p1.c
$RM_CMD $G0/regression/rt_pkpm_travel_pulse_p2.c
$RM_CMD $G0/regression/rt_pkpm_wall_p1.c
$RM_CMD $G0/regression/rt_pkpm_wall_p2.c

# Lua regression tests
mkdir -p luareg
$CP_CMD $G0/regression/lua/rt_pkpm_2d_travel_pulse_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_alf_wave_1x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_alf_wave_1x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_alf_wave_explicit_1x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_alf_wave_explicit_1x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_em_advect_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_em_advect_resonant_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_es_pot_well_1x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_es_pot_well_1x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_es_shock_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_es_shock_reflect_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_landau_damping_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_landau_damping_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_neut_sodshock_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_neut_sodshock_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_periodic_es_shock_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_periodic_es_shock_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_periodic_neut_sodshock_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_periodic_neut_sodshock_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_sheath_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_square_relax_1x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_square_relax_1x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_travel_pulse_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_travel_pulse_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_wall_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_pkpm_wall_p2.lua luareg/

$RM_CMD $G0/regression/lua/rt_pkpm_2d_travel_pulse_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_alf_wave_1x_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_alf_wave_1x_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_alf_wave_explicit_1x_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_alf_wave_explicit_1x_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_em_advect_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_em_advect_resonant_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_es_pot_well_1x_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_es_pot_well_1x_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_es_shock_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_es_shock_reflect_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_landau_damping_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_landau_damping_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_neut_sodshock_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_neut_sodshock_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_periodic_es_shock_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_periodic_es_shock_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_periodic_neut_sodshock_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_periodic_neut_sodshock_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_sheath_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_square_relax_1x_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_square_relax_1x_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_travel_pulse_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_travel_pulse_p2.lua
$RM_CMD $G0/regression/lua/rt_pkpm_wall_p1.lua
$RM_CMD $G0/regression/lua/rt_pkpm_wall_p2.lua