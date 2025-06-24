#!/bin/sh

CP_CMD=cp
RM_CMD=rm
G0=..

# kernels
mkdir -p kernels/advection
$CP_CMD $G0/kernels/advection/*.h kernels/advection/
$CP_CMD $G0/kernels/advection/*.c kernels/advection/
mkdir -p kernels/canonical_pb
$CP_CMD $G0/kernels/canonical_pb/*.h kernels/canonical_pb/
$CP_CMD $G0/kernels/canonical_pb/*.c kernels/canonical_pb/
mkdir -p kernels/dg_diffusion_fluid
$CP_CMD $G0/kernels/dg_diffusion_fluid/*.h kernels/dg_diffusion_fluid/
$CP_CMD $G0/kernels/dg_diffusion_fluid/*.c kernels/dg_diffusion_fluid/
mkdir -p kernels/dg_diffusion_gen
$CP_CMD $G0/kernels/dg_diffusion_gen/*.h kernels/dg_diffusion_gen/
$CP_CMD $G0/kernels/dg_diffusion_gen/*.c kernels/dg_diffusion_gen/
mkdir -p kernels/dg_diffusion_vlasov
$CP_CMD $G0/kernels/dg_diffusion_vlasov/*.h kernels/dg_diffusion_vlasov/
$CP_CMD $G0/kernels/dg_diffusion_vlasov/*.c kernels/dg_diffusion_vlasov/
mkdir -p kernels/euler
$CP_CMD $G0/kernels/euler/*.h kernels/euler/
$CP_CMD $G0/kernels/euler/*.c kernels/euler/
mkdir -p kernels/fpo
$CP_CMD $G0/kernels/fpo/*.h kernels/fpo/
$CP_CMD $G0/kernels/fpo/*.c kernels/fpo/
mkdir -p kernels/lbo_vlasov
$CP_CMD $G0/kernels/lbo_vlasov/*.h kernels/lbo_vlasov/
$CP_CMD $G0/kernels/lbo_vlasov/*.c kernels/lbo_vlasov/
mkdir -p kernels/maxwell
$CP_CMD $G0/kernels/maxwell/*.h kernels/maxwell/
$CP_CMD $G0/kernels/maxwell/*.c kernels/maxwell/
mkdir -p kernels/sr_vlasov
$CP_CMD $G0/kernels/sr_vlasov/*.h kernels/sr_vlasov/
$CP_CMD $G0/kernels/sr_vlasov/*.c kernels/sr_vlasov/
mkdir -p kernels/vlasov
$CP_CMD $G0/kernels/vlasov/*.h kernels/vlasov/
$CP_CMD $G0/kernels/vlasov/*.c kernels/vlasov/
mkdir -p kernels/vlasov_poisson
$CP_CMD $G0/kernels/vlasov_poisson/*.h kernels/vlasov_poisson/
$CP_CMD $G0/kernels/vlasov_poisson/*.c kernels/vlasov_poisson/
mkdir -p kernels/positivity_shift_vlasov
$CP_CMD $G0/kernels/positivity_shift_vlasov/*.h kernels/positivity_shift_vlasov/
$CP_CMD $G0/kernels/positivity_shift_vlasov/*.c kernels/positivity_shift_vlasov/

$RM_CMD $G0/kernels/advection/*.h
$RM_CMD $G0/kernels/advection/*.c
$RM_CMD $G0/kernels/canonical_pb/*.h
$RM_CMD $G0/kernels/canonical_pb/*.c
$RM_CMD $G0/kernels/dg_diffusion_fluid/*.h
$RM_CMD $G0/kernels/dg_diffusion_fluid/*.c
$RM_CMD $G0/kernels/dg_diffusion_gen/*.h
$RM_CMD $G0/kernels/dg_diffusion_gen/*.c
$RM_CMD $G0/kernels/dg_diffusion_vlasov/*.h
$RM_CMD $G0/kernels/dg_diffusion_vlasov/*.c
$RM_CMD $G0/kernels/euler/*.h
$RM_CMD $G0/kernels/euler/*.c
$RM_CMD $G0/kernels/fpo/*.h
$RM_CMD $G0/kernels/fpo/*.c
$RM_CMD $G0/kernels/lbo_vlasov/*.h
$RM_CMD $G0/kernels/lbo_vlasov/*.c
$RM_CMD $G0/kernels/maxwell/*.h
$RM_CMD $G0/kernels/maxwell/*.c
$RM_CMD $G0/kernels/sr_vlasov/*.h
$RM_CMD $G0/kernels/sr_vlasov/*.c
$RM_CMD $G0/kernels/vlasov/*.h
$RM_CMD $G0/kernels/vlasov/*.c
$RM_CMD $G0/kernels/vlasov_poisson/*.h
$RM_CMD $G0/kernels/vlasov_poisson/*.c
$RM_CMD $G0/kernels/positivity_shift_vlasov/*.h
$RM_CMD $G0/kernels/positivity_shift_vlasov/*.c

# zero
mkdir -p zero
$CP_CMD $G0/zero/dg_diffusion_vlasov_cu.cu zero/
$CP_CMD $G0/zero/dg_diffusion_vlasov.c zero/
$CP_CMD $G0/zero/dg_fpo_vlasov_diff_cu.cu zero/
$CP_CMD $G0/zero/dg_fpo_vlasov_diff.c zero/
$CP_CMD $G0/zero/dg_fpo_vlasov_drag_cu.cu zero/
$CP_CMD $G0/zero/dg_fpo_vlasov_drag.c zero/
$CP_CMD $G0/zero/dg_lbo_vlasov_diff_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_vlasov_diff.c zero/
$CP_CMD $G0/zero/dg_lbo_vlasov_drag_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_vlasov_drag.c zero/
$CP_CMD $G0/zero/dg_updater_bflux_vlasov.c zero/
$CP_CMD $G0/zero/dg_updater_diffusion_vlasov.c zero/
$CP_CMD $G0/zero/dg_updater_fpo_vlasov.c zero/
$CP_CMD $G0/zero/dg_updater_lbo_vlasov.c zero/
$CP_CMD $G0/zero/dg_updater_moment.c zero/
$CP_CMD $G0/zero/dg_updater_rad_vlasov.c zero/
$CP_CMD $G0/zero/dg_updater_vlasov_poisson.c zero/
$CP_CMD $G0/zero/dg_updater_vlasov.c zero/
$CP_CMD $G0/zero/dg_vlasov_cu.cu zero/
$CP_CMD $G0/zero/dg_vlasov_poisson_cu.cu zero/
$CP_CMD $G0/zero/dg_vlasov_poisson.c zero/
$CP_CMD $G0/zero/dg_vlasov_sr_cu.cu zero/
$CP_CMD $G0/zero/dg_vlasov_sr.c zero/
$CP_CMD $G0/zero/dg_vlasov.c zero/
$CP_CMD $G0/zero/gkyl_dg_canonical_pb_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_fpo_vlasov_diff_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_fpo_vlasov_diff.h zero/
$CP_CMD $G0/zero/gkyl_dg_fpo_vlasov_drag_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_fpo_vlasov_drag.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_vlasov_diff_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_vlasov_diff.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_vlasov_drag_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_vlasov_drag.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_bflux_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_bflux_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_fpo_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_lbo_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_moment.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_rad_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_rad_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_vlasov_poisson.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_vlasov_timers.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_poisson_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_poisson.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_sr_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov_sr.h zero/
$CP_CMD $G0/zero/gkyl_dg_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_mom_calc_bcorr.h zero/
$CP_CMD $G0/zero/gkyl_mom_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_vlasov_sr_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_vlasov_sr.h zero/
$CP_CMD $G0/zero/gkyl_mom_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_calc.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_vlasov_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_correct_priv.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_correct.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_moments_priv.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_moments.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_proj_on_basis_priv.h zero/
$CP_CMD $G0/zero/gkyl_vlasov_lte_proj_on_basis.h zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_vlasov_cu.cu zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_vlasov.c zero/
$CP_CMD $G0/zero/mom_calc_bcorr.c zero/
$CP_CMD $G0/zero/mom_calc_bcorr_vlasov.c zero/
$CP_CMD $G0/zero/mom_vlasov_cu.cu zero/
$CP_CMD $G0/zero/mom_vlasov_sr_cu.cu zero/
$CP_CMD $G0/zero/mom_vlasov_sr.c zero/
$CP_CMD $G0/zero/mom_vlasov.c zero/
$CP_CMD $G0/zero/prim_lbo_calc_cu.cu zero/
$CP_CMD $G0/zero/prim_lbo_calc.c zero/
$CP_CMD $G0/zero/prim_lbo_calc_vlasov.c zero/
$CP_CMD $G0/zero/prim_lbo_cross_calc_cu.cu zero/
$CP_CMD $G0/zero/prim_lbo_cross_calc.c zero/
$CP_CMD $G0/zero/prim_lbo_cross_calc_vlasov.c zero/
$CP_CMD $G0/zero/prim_lbo_vlasov_cu.cu zero/
$CP_CMD $G0/zero/prim_lbo_vlasov.c zero/
$CP_CMD $G0/zero/vlasov_lte_correct_cu.cu zero/
$CP_CMD $G0/zero/vlasov_lte_correct.c zero/
$CP_CMD $G0/zero/vlasov_lte_moments.c zero/
$CP_CMD $G0/zero/vlasov_lte_proj_on_basis_cu.cu zero/
$CP_CMD $G0/zero/vlasov_lte_proj_on_basis.c zero/
$CP_CMD $G0/zero/gkyl_mom_type.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_sr_vars.h zero/
$CP_CMD $G0/zero/gkyl_velocity_map_priv.h zero/
$CP_CMD $G0/zero/gkyl_velocity_map.h zero/
$CP_CMD $G0/zero/gkyl_dg_eqn.h zero/
$CP_CMD $G0/zero/gkyl_hyper_dg.h zero/
$CP_CMD $G0/zero/bc_basic_cu.cu zero/
$CP_CMD $G0/zero/bc_basic.c zero/
$CP_CMD $G0/zero/bc_emission_elastic_cu.cu zero/
$CP_CMD $G0/zero/bc_emission_elastic.c zero/
$CP_CMD $G0/zero/bc_emission_spectrum_cu.cu zero/
$CP_CMD $G0/zero/bc_emission_spectrum.c zero/
$CP_CMD $G0/zero/bc_emission.c zero/
$CP_CMD $G0/zero/bgk_collisions_cu.cu zero/
$CP_CMD $G0/zero/bgk_collisions.c zero/
$CP_CMD $G0/zero/dg_advection_cu.cu zero/
$CP_CMD $G0/zero/dg_advection.c zero/
$CP_CMD $G0/zero/dg_calc_canonical_pb_fluid_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_canonical_pb_fluid_vars.c zero/
$CP_CMD $G0/zero/dg_calc_canonical_pb_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_canonical_pb_vars.c zero/
$CP_CMD $G0/zero/dg_calc_em_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_em_vars.c zero/
$CP_CMD $G0/zero/dg_calc_fluid_em_coupling_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_fluid_em_coupling.c zero/
$CP_CMD $G0/zero/dg_calc_fluid_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_fluid_vars.c zero/
$CP_CMD $G0/zero/dg_calc_prim_vars.c zero/
$CP_CMD $G0/zero/dg_calc_sr_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_sr_vars.c zero/
$CP_CMD $G0/zero/dg_canonical_pb_cu.cu zero/
$CP_CMD $G0/zero/dg_canonical_pb_fluid_cu.cu zero/
$CP_CMD $G0/zero/dg_canonical_pb_fluid.c zero/
$CP_CMD $G0/zero/dg_canonical_pb.c zero/
$CP_CMD $G0/zero/dg_diffusion_fluid_cu.cu zero/
$CP_CMD $G0/zero/dg_diffusion_fluid.c zero/
$CP_CMD $G0/zero/dg_diffusion_gen_cu.cu zero/
$CP_CMD $G0/zero/dg_diffusion_gen.c zero/
$CP_CMD $G0/zero/dg_eqn.c zero/
$CP_CMD $G0/zero/dg_euler_cu.cu zero/
$CP_CMD $G0/zero/dg_euler.c zero/
$CP_CMD $G0/zero/dg_maxwell_cu.cu zero/
$CP_CMD $G0/zero/dg_maxwell.c zero/
$CP_CMD $G0/zero/dg_updater_diffusion_fluid.c zero/
$CP_CMD $G0/zero/dg_updater_diffusion_gen.c zero/
$CP_CMD $G0/zero/dg_updater_fluid.c zero/
$CP_CMD $G0/zero/emission_elastic_model.c zero/
$CP_CMD $G0/zero/emission_spectrum_model_cu.cu zero/
$CP_CMD $G0/zero/emission_spectrum_model.c zero/
$CP_CMD $G0/zero/emission_yield_model_cu.cu zero/
$CP_CMD $G0/zero/emission_yield_model.c zero/
$CP_CMD $G0/zero/ghost_surf_calc_cu.cu zero/
$CP_CMD $G0/zero/ghost_surf_calc.c zero/
$CP_CMD $G0/zero/gkyl_bc_basic_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_basic.h zero/
$CP_CMD $G0/zero/gkyl_bc_emission_elastic_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_emission_elastic.h zero/
$CP_CMD $G0/zero/gkyl_bc_emission_spectrum_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_emission_spectrum.h zero/
$CP_CMD $G0/zero/gkyl_bc_emission.h zero/
$CP_CMD $G0/zero/gkyl_bgk_collisions_priv.h zero/
$CP_CMD $G0/zero/gkyl_bgk_collisions.h zero/
$CP_CMD $G0/zero/gkyl_dg_advection_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_advection.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_canonical_pb_fluid_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_canonical_pb_fluid_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_canonical_pb_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_canonical_pb_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_em_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_em_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_fluid_em_coupling_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_fluid_em_coupling.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_fluid_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_fluid_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_prim_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_sr_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_canonical_pb_fluid_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_canonical_pb_fluid.h zero/
$CP_CMD $G0/zero/gkyl_dg_canonical_pb.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_fluid_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_fluid.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_gen_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_gen.h zero/
$CP_CMD $G0/zero/gkyl_dg_euler_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_euler.h zero/
$CP_CMD $G0/zero/gkyl_dg_maxwell_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_maxwell.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_collisions_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_fluid_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_fluid.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_gen_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_gen.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_fluid_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_fluid.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_moment_priv.h zero/
$CP_CMD $G0/zero/gkyl_emission_elastic_model.h zero/
$CP_CMD $G0/zero/gkyl_emission_spectrum_model.h zero/
$CP_CMD $G0/zero/gkyl_emission_yield_model.h zero/
$CP_CMD $G0/zero/gkyl_ghost_surf_calc_priv.h zero/
$CP_CMD $G0/zero/gkyl_ghost_surf_calc.h zero/
$CP_CMD $G0/zero/gkyl_hyper_dg_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_calc_bcorr_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_calc_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_calc.h zero/
$CP_CMD $G0/zero/gkyl_mom_canonical_pb_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_canonical_pb.h zero/
$CP_CMD $G0/zero/gkyl_prim_bgk_cross_calc.h zero/
$CP_CMD $G0/zero/gkyl_prim_cross_m0deltas_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_cross_m0deltas.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_calc_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_cross_calc_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_cross_calc.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_type.h zero/
$CP_CMD $G0/zero/gkyl_spitzer_coll_freq_priv.h zero/
$CP_CMD $G0/zero/gkyl_spitzer_coll_freq.h zero/
$CP_CMD $G0/zero/hyper_dg_cu.cu zero/
$CP_CMD $G0/zero/hyper_dg.c zero/
$CP_CMD $G0/zero/mom_calc_bcorr_cu.cu zero/
$CP_CMD $G0/zero/mom_calc_cu.cu zero/
$CP_CMD $G0/zero/mom_calc.c zero/
$CP_CMD $G0/zero/mom_canonical_pb_cu.cu zero/
$CP_CMD $G0/zero/mom_canonical_pb.c zero/
$CP_CMD $G0/zero/mom_type.c zero/
$CP_CMD $G0/zero/prim_bgk_cross_calc_cu.cu zero/
$CP_CMD $G0/zero/prim_bgk_cross_calc.c zero/
$CP_CMD $G0/zero/prim_cross_m0deltas_cu.cu zero/
$CP_CMD $G0/zero/prim_cross_m0deltas.c zero/
$CP_CMD $G0/zero/prim_lbo_type.c zero/
$CP_CMD $G0/zero/spitzer_coll_freq_cu.cu zero/
$CP_CMD $G0/zero/spitzer_coll_freq.c zero/
$CP_CMD $G0/zero/velocity_map_cu.cu zero/
$CP_CMD $G0/zero/velocity_map.c zero/
$CP_CMD $G0/zero/gkyl_positivity_shift_vlasov.h zero/
$CP_CMD $G0/zero/gkyl_positivity_shift_vlasov_priv.h zero/
$CP_CMD $G0/zero/positivity_shift_vlasov.c zero/
$CP_CMD $G0/zero/positivity_shift_vlasov_cu.cu zero/

$RM_CMD $G0/zero/dg_diffusion_vlasov_cu.cu
$RM_CMD $G0/zero/dg_diffusion_vlasov.c
$RM_CMD $G0/zero/dg_fpo_vlasov_diff_cu.cu
$RM_CMD $G0/zero/dg_fpo_vlasov_diff.c
$RM_CMD $G0/zero/dg_fpo_vlasov_drag_cu.cu
$RM_CMD $G0/zero/dg_fpo_vlasov_drag.c
$RM_CMD $G0/zero/dg_lbo_vlasov_diff_cu.cu
$RM_CMD $G0/zero/dg_lbo_vlasov_diff.c
$RM_CMD $G0/zero/dg_lbo_vlasov_drag_cu.cu
$RM_CMD $G0/zero/dg_lbo_vlasov_drag.c
$RM_CMD $G0/zero/dg_updater_bflux_vlasov.c
$RM_CMD $G0/zero/dg_updater_diffusion_vlasov.c
$RM_CMD $G0/zero/dg_updater_fpo_vlasov.c
$RM_CMD $G0/zero/dg_updater_lbo_vlasov.c
$RM_CMD $G0/zero/dg_updater_moment.c
$RM_CMD $G0/zero/dg_updater_rad_vlasov.c
$RM_CMD $G0/zero/dg_updater_vlasov_poisson.c
$RM_CMD $G0/zero/dg_updater_vlasov.c
$RM_CMD $G0/zero/dg_vlasov_cu.cu
$RM_CMD $G0/zero/dg_vlasov_poisson_cu.cu
$RM_CMD $G0/zero/dg_vlasov_poisson.c
$RM_CMD $G0/zero/dg_vlasov_sr_cu.cu
$RM_CMD $G0/zero/dg_vlasov_sr.c
$RM_CMD $G0/zero/dg_vlasov.c
$RM_CMD $G0/zero/gkyl_dg_canonical_pb_priv.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_fpo_vlasov_diff_priv.h
$RM_CMD $G0/zero/gkyl_dg_fpo_vlasov_diff.h
$RM_CMD $G0/zero/gkyl_dg_fpo_vlasov_drag_priv.h
$RM_CMD $G0/zero/gkyl_dg_fpo_vlasov_drag.h
$RM_CMD $G0/zero/gkyl_dg_lbo_vlasov_diff_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_vlasov_diff.h
$RM_CMD $G0/zero/gkyl_dg_lbo_vlasov_drag_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_vlasov_drag.h
$RM_CMD $G0/zero/gkyl_dg_updater_bflux_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_bflux_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_updater_fpo_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_updater_lbo_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_updater_moment.h
$RM_CMD $G0/zero/gkyl_dg_updater_rad_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_rad_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_updater_vlasov_poisson.h
$RM_CMD $G0/zero/gkyl_dg_updater_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_vlasov_timers.h
$RM_CMD $G0/zero/gkyl_dg_updater_vlasov.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_poisson_priv.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_poisson.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_sr_priv.h
$RM_CMD $G0/zero/gkyl_dg_vlasov_sr.h
$RM_CMD $G0/zero/gkyl_dg_vlasov.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_vlasov.h
$RM_CMD $G0/zero/gkyl_mom_calc_bcorr.h
$RM_CMD $G0/zero/gkyl_mom_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_mom_vlasov_sr_priv.h
$RM_CMD $G0/zero/gkyl_mom_vlasov_sr.h
$RM_CMD $G0/zero/gkyl_mom_vlasov.h
$RM_CMD $G0/zero/gkyl_prim_lbo_calc.h
$RM_CMD $G0/zero/gkyl_prim_lbo_vlasov_priv.h
$RM_CMD $G0/zero/gkyl_prim_lbo_vlasov.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_correct_priv.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_correct.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_moments_priv.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_moments.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_proj_on_basis_priv.h
$RM_CMD $G0/zero/gkyl_vlasov_lte_proj_on_basis.h
$RM_CMD $G0/zero/mom_bcorr_lbo_vlasov_cu.cu
$RM_CMD $G0/zero/mom_bcorr_lbo_vlasov.c
$RM_CMD $G0/zero/mom_calc_bcorr.c
$RM_CMD $G0/zero/mom_calc_bcorr_vlasov.c
$RM_CMD $G0/zero/mom_vlasov_cu.cu
$RM_CMD $G0/zero/mom_vlasov_sr_cu.cu
$RM_CMD $G0/zero/mom_vlasov_sr.c
$RM_CMD $G0/zero/mom_vlasov.c
$RM_CMD $G0/zero/prim_lbo_calc_cu.cu
$RM_CMD $G0/zero/prim_lbo_calc.c
$RM_CMD $G0/zero/prim_lbo_calc_vlasov.c
$RM_CMD $G0/zero/prim_lbo_cross_calc_cu.cu
$RM_CMD $G0/zero/prim_lbo_cross_calc.c
$RM_CMD $G0/zero/prim_lbo_cross_calc_vlasov.c
$RM_CMD $G0/zero/prim_lbo_vlasov_cu.cu
$RM_CMD $G0/zero/prim_lbo_vlasov.c
$RM_CMD $G0/zero/vlasov_lte_correct_cu.cu
$RM_CMD $G0/zero/vlasov_lte_correct.c
$RM_CMD $G0/zero/vlasov_lte_moments.c
$RM_CMD $G0/zero/vlasov_lte_proj_on_basis_cu.cu
$RM_CMD $G0/zero/vlasov_lte_proj_on_basis.c
$RM_CMD $G0/zero/gkyl_mom_type.h
$RM_CMD $G0/zero/gkyl_dg_calc_sr_vars.h
$RM_CMD $G0/zero/gkyl_velocity_map_priv.h
$RM_CMD $G0/zero/gkyl_velocity_map.h
$RM_CMD $G0/zero/gkyl_dg_eqn.h
$RM_CMD $G0/zero/gkyl_hyper_dg.h
$RM_CMD $G0/zero/bc_basic_cu.cu
$RM_CMD $G0/zero/bc_basic.c
$RM_CMD $G0/zero/bc_emission_elastic_cu.cu
$RM_CMD $G0/zero/bc_emission_elastic.c
$RM_CMD $G0/zero/bc_emission_spectrum_cu.cu
$RM_CMD $G0/zero/bc_emission_spectrum.c
$RM_CMD $G0/zero/bc_emission.c
$RM_CMD $G0/zero/bgk_collisions_cu.cu
$RM_CMD $G0/zero/bgk_collisions.c
$RM_CMD $G0/zero/dg_advection_cu.cu
$RM_CMD $G0/zero/dg_advection.c
$RM_CMD $G0/zero/dg_calc_canonical_pb_fluid_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_canonical_pb_fluid_vars.c
$RM_CMD $G0/zero/dg_calc_canonical_pb_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_canonical_pb_vars.c
$RM_CMD $G0/zero/dg_calc_em_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_em_vars.c
$RM_CMD $G0/zero/dg_calc_fluid_em_coupling_cu.cu
$RM_CMD $G0/zero/dg_calc_fluid_em_coupling.c
$RM_CMD $G0/zero/dg_calc_fluid_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_fluid_vars.c
$RM_CMD $G0/zero/dg_calc_prim_vars.c
$RM_CMD $G0/zero/dg_calc_sr_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_sr_vars.c
$RM_CMD $G0/zero/dg_canonical_pb_cu.cu
$RM_CMD $G0/zero/dg_canonical_pb_fluid_cu.cu
$RM_CMD $G0/zero/dg_canonical_pb_fluid.c
$RM_CMD $G0/zero/dg_canonical_pb.c
$RM_CMD $G0/zero/dg_diffusion_fluid_cu.cu
$RM_CMD $G0/zero/dg_diffusion_fluid.c
$RM_CMD $G0/zero/dg_diffusion_gen_cu.cu
$RM_CMD $G0/zero/dg_diffusion_gen.c
$RM_CMD $G0/zero/dg_eqn.c
$RM_CMD $G0/zero/dg_euler_cu.cu
$RM_CMD $G0/zero/dg_euler.c
$RM_CMD $G0/zero/dg_maxwell_cu.cu
$RM_CMD $G0/zero/dg_maxwell.c
$RM_CMD $G0/zero/dg_updater_diffusion_fluid.c
$RM_CMD $G0/zero/dg_updater_diffusion_gen.c
$RM_CMD $G0/zero/dg_updater_fluid.c
$RM_CMD $G0/zero/emission_elastic_model.c
$RM_CMD $G0/zero/emission_spectrum_model_cu.cu
$RM_CMD $G0/zero/emission_spectrum_model.c
$RM_CMD $G0/zero/emission_yield_model_cu.cu
$RM_CMD $G0/zero/emission_yield_model.c
$RM_CMD $G0/zero/ghost_surf_calc_cu.cu
$RM_CMD $G0/zero/ghost_surf_calc.c
$RM_CMD $G0/zero/gkyl_bc_basic_priv.h
$RM_CMD $G0/zero/gkyl_bc_basic.h
$RM_CMD $G0/zero/gkyl_bc_emission_elastic_priv.h
$RM_CMD $G0/zero/gkyl_bc_emission_elastic.h
$RM_CMD $G0/zero/gkyl_bc_emission_spectrum_priv.h
$RM_CMD $G0/zero/gkyl_bc_emission_spectrum.h
$RM_CMD $G0/zero/gkyl_bc_emission.h
$RM_CMD $G0/zero/gkyl_bgk_collisions_priv.h
$RM_CMD $G0/zero/gkyl_bgk_collisions.h
$RM_CMD $G0/zero/gkyl_dg_advection_priv.h
$RM_CMD $G0/zero/gkyl_dg_advection.h
$RM_CMD $G0/zero/gkyl_dg_calc_canonical_pb_fluid_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_canonical_pb_fluid_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_canonical_pb_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_canonical_pb_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_em_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_em_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_fluid_em_coupling_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_fluid_em_coupling.h
$RM_CMD $G0/zero/gkyl_dg_calc_fluid_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_fluid_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_prim_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_sr_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_canonical_pb_fluid_priv.h
$RM_CMD $G0/zero/gkyl_dg_canonical_pb_fluid.h
$RM_CMD $G0/zero/gkyl_dg_canonical_pb.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_fluid_priv.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_fluid.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_gen_priv.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_gen.h
$RM_CMD $G0/zero/gkyl_dg_euler_priv.h
$RM_CMD $G0/zero/gkyl_dg_euler.h
$RM_CMD $G0/zero/gkyl_dg_maxwell_priv.h
$RM_CMD $G0/zero/gkyl_dg_maxwell.h
$RM_CMD $G0/zero/gkyl_dg_updater_collisions_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_fluid_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_fluid.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_gen_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_gen.h
$RM_CMD $G0/zero/gkyl_dg_updater_fluid_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_fluid.h
$RM_CMD $G0/zero/gkyl_dg_updater_moment_priv.h
$RM_CMD $G0/zero/gkyl_emission_elastic_model.h
$RM_CMD $G0/zero/gkyl_emission_spectrum_model.h
$RM_CMD $G0/zero/gkyl_emission_yield_model.h
$RM_CMD $G0/zero/gkyl_ghost_surf_calc_priv.h
$RM_CMD $G0/zero/gkyl_ghost_surf_calc.h
$RM_CMD $G0/zero/gkyl_hyper_dg_priv.h
$RM_CMD $G0/zero/gkyl_mom_calc_bcorr_priv.h
$RM_CMD $G0/zero/gkyl_mom_calc_priv.h
$RM_CMD $G0/zero/gkyl_mom_calc.h
$RM_CMD $G0/zero/gkyl_mom_canonical_pb_priv.h
$RM_CMD $G0/zero/gkyl_mom_canonical_pb.h
$RM_CMD $G0/zero/gkyl_prim_bgk_cross_calc.h
$RM_CMD $G0/zero/gkyl_prim_cross_m0deltas_priv.h
$RM_CMD $G0/zero/gkyl_prim_cross_m0deltas.h
$RM_CMD $G0/zero/gkyl_prim_lbo_calc_priv.h
$RM_CMD $G0/zero/gkyl_prim_lbo_cross_calc_priv.h
$RM_CMD $G0/zero/gkyl_prim_lbo_cross_calc.h
$RM_CMD $G0/zero/gkyl_prim_lbo_type.h
$RM_CMD $G0/zero/gkyl_spitzer_coll_freq_priv.h
$RM_CMD $G0/zero/gkyl_spitzer_coll_freq.h
$RM_CMD $G0/zero/hyper_dg_cu.cu
$RM_CMD $G0/zero/hyper_dg.c
$RM_CMD $G0/zero/mom_calc_bcorr_cu.cu
$RM_CMD $G0/zero/mom_calc_cu.cu
$RM_CMD $G0/zero/mom_calc.c
$RM_CMD $G0/zero/mom_canonical_pb_cu.cu
$RM_CMD $G0/zero/mom_canonical_pb.c
$RM_CMD $G0/zero/mom_type.c
$RM_CMD $G0/zero/prim_bgk_cross_calc_cu.cu
$RM_CMD $G0/zero/prim_bgk_cross_calc.c
$RM_CMD $G0/zero/prim_cross_m0deltas_cu.cu
$RM_CMD $G0/zero/prim_cross_m0deltas.c
$RM_CMD $G0/zero/prim_lbo_type.c
$RM_CMD $G0/zero/spitzer_coll_freq_cu.cu
$RM_CMD $G0/zero/spitzer_coll_freq.c
$RM_CMD $G0/zero/velocity_map_cu.cu
$RM_CMD $G0/zero/velocity_map.c
$RM_CMD $G0/zero/gkyl_positivity_shift_vlasov.h
$RM_CMD $G0/zero/gkyl_positivity_shift_vlasov_priv.h
$RM_CMD $G0/zero/positivity_shift_vlasov.c
$RM_CMD $G0/zero/positivity_shift_vlasov_cu.cu

# app
mkdir -p apps
$CP_CMD $G0/apps/gkyl_vlasov_comms.h apps/
$CP_CMD $G0/apps/gkyl_vlasov_lw.h apps/
$CP_CMD $G0/apps/gkyl_vlasov_priv.h apps/
$CP_CMD $G0/apps/gkyl_vlasov.h apps/
$CP_CMD $G0/apps/vlasov_comms.c apps/
$CP_CMD $G0/apps/vlasov_forward_euler.c apps/
$CP_CMD $G0/apps/vlasov_lw.c apps/
$CP_CMD $G0/apps/vlasov_poisson_update_ssp_rk3.c apps/
$CP_CMD $G0/apps/vlasov_update_implicit_coll.c apps/
$CP_CMD $G0/apps/vlasov_update_op_split.c apps/
$CP_CMD $G0/apps/vlasov_update_ssp_rk3.c apps/
$CP_CMD $G0/apps/vlasov.c apps/
$CP_CMD $G0/apps/vm_field.c apps/
$CP_CMD $G0/apps/vm_fluid_em_coupling.c apps/
$CP_CMD $G0/apps/vm_fluid_species_source.c apps/
$CP_CMD $G0/apps/vm_fluid_species.c apps/
$CP_CMD $G0/apps/vm_species_bflux.c apps/
$CP_CMD $G0/apps/vm_species_bgk.c apps/
$CP_CMD $G0/apps/vm_species_emission.c apps/
$CP_CMD $G0/apps/vm_species_lbo.c apps/
$CP_CMD $G0/apps/vm_species_lte.c apps/
$CP_CMD $G0/apps/vm_species_moment.c apps/
$CP_CMD $G0/apps/vm_species_projection.c apps/
$CP_CMD $G0/apps/vm_species_radiation.c apps/
$CP_CMD $G0/apps/vm_species_source.c apps/
$CP_CMD $G0/apps/vm_species.c apps/
$CP_CMD $G0/apps/vp_field.c apps/

$RM_CMD $G0/apps/gkyl_vlasov_comms.h
$RM_CMD $G0/apps/gkyl_vlasov_lw.h
$RM_CMD $G0/apps/gkyl_vlasov_priv.h
$RM_CMD $G0/apps/gkyl_vlasov.h
$RM_CMD $G0/apps/vlasov_comms.c
$RM_CMD $G0/apps/vlasov_forward_euler.c
$RM_CMD $G0/apps/vlasov_lw.c
$RM_CMD $G0/apps/vlasov_poisson_update_ssp_rk3.c
$RM_CMD $G0/apps/vlasov_update_implicit_coll.c
$RM_CMD $G0/apps/vlasov_update_op_split.c
$RM_CMD $G0/apps/vlasov_update_ssp_rk3.c
$RM_CMD $G0/apps/vlasov.c
$RM_CMD $G0/apps/vm_field.c
$RM_CMD $G0/apps/vm_fluid_em_coupling.c
$RM_CMD $G0/apps/vm_fluid_species_source.c
$RM_CMD $G0/apps/vm_fluid_species.c
$RM_CMD $G0/apps/vm_species_bflux.c
$RM_CMD $G0/apps/vm_species_bgk.c
$RM_CMD $G0/apps/vm_species_emission.c
$RM_CMD $G0/apps/vm_species_lbo.c
$RM_CMD $G0/apps/vm_species_lte.c
$RM_CMD $G0/apps/vm_species_moment.c
$RM_CMD $G0/apps/vm_species_projection.c
$RM_CMD $G0/apps/vm_species_radiation.c
$RM_CMD $G0/apps/vm_species_source.c
$RM_CMD $G0/apps/vm_species.c
$RM_CMD $G0/apps/vp_field.c

# unit
mkdir -p unit
$CP_CMD $G0/unit/ctest_bc_basic.c unit/
$CP_CMD $G0/unit/ctest_canonical_pb_continuity.c unit/
$CP_CMD $G0/unit/ctest_canonical_pb_equilibrium.c unit/
$CP_CMD $G0/unit/ctest_correct_maxwellian.c unit/
$CP_CMD $G0/unit/ctest_correct_mj_integrated.c unit/
$CP_CMD $G0/unit/ctest_dg_lbo_vlasov.c unit/
$CP_CMD $G0/unit/ctest_dg_maxwell.c unit/
$CP_CMD $G0/unit/ctest_dg_vlasov.c unit/
$CP_CMD $G0/unit/ctest_hyper_dg.c unit/
$CP_CMD $G0/unit/ctest_hyper3x_dg.c unit/
$CP_CMD $G0/unit/ctest_mom_vlasov.c unit/
$CP_CMD $G0/unit/ctest_prim_vlasov.c unit/
$CP_CMD $G0/unit/ctest_proj_mj_on_basis.c unit/
$CP_CMD $G0/unit/ctest_spitzer_coll_freq.c unit/
$CP_CMD $G0/unit/ctest_dg_em_vars.c unit/
$CP_CMD $G0/unit/ctest_velocity_map.c unit/
$CP_CMD $G0/unit/ctest_positivity_shift_vlasov.c unit/

$RM_CMD $G0/unit/ctest_bc_basic.c
$RM_CMD $G0/unit/ctest_canonical_pb_continuity.c
$RM_CMD $G0/unit/ctest_canonical_pb_equilibrium.c
$RM_CMD $G0/unit/ctest_correct_maxwellian.c
$RM_CMD $G0/unit/ctest_correct_mj_integrated.c
$RM_CMD $G0/unit/ctest_dg_lbo_vlasov.c
$RM_CMD $G0/unit/ctest_dg_maxwell.c
$RM_CMD $G0/unit/ctest_dg_vlasov.c
$RM_CMD $G0/unit/ctest_hyper_dg.c
$RM_CMD $G0/unit/ctest_hyper3x_dg.c
$RM_CMD $G0/unit/ctest_mom_vlasov.c
$RM_CMD $G0/unit/ctest_prim_vlasov.c
$RM_CMD $G0/unit/ctest_proj_mj_on_basis.c
$RM_CMD $G0/unit/ctest_spitzer_coll_freq.c
$RM_CMD $G0/unit/ctest_dg_em_vars.c
$RM_CMD $G0/unit/ctest_velocity_map.c
$RM_CMD $G0/unit/ctest_positivity_shift_vlasov.c

# C regression tests
mkdir -p creg
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_ex_bgk_surf_flat_sq_ic.c creg/
$CP_CMD $G0/regression/rt_can_pb_ex_bgk_surf_flat.c creg/
$CP_CMD $G0/regression/rt_can_pb_free_streaming_surf_sphere.c creg/
$CP_CMD $G0/regression/rt_can_pb_im_bgk_surf_flat_sq_ic.c creg/
$CP_CMD $G0/regression/rt_can_pb_im_bgk_surf_flat.c creg/
$CP_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x1v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_dg_5m_mom_beach_p2.c creg/
$CP_CMD $G0/regression/rt_dg_5m_mom_beach_p3.c creg/
$CP_CMD $G0/regression/rt_dg_accel_1x1v.c creg/
$CP_CMD $G0/regression/rt_dg_euler_kh_2d.c creg/
$CP_CMD $G0/regression/rt_dg_euler_p_perturbation_p1.c creg/
$CP_CMD $G0/regression/rt_dg_euler_p_perturbation_p2.c creg/
$CP_CMD $G0/regression/rt_dg_euler_sodshock_p1.c creg/
$CP_CMD $G0/regression/rt_dg_euler_sodshock_p2.c creg/
$CP_CMD $G0/regression/rt_dg_maxwell_plane_wave_2d.c creg/
$CP_CMD $G0/regression/rt_dg_maxwell_wg_2d.c creg/
$CP_CMD $G0/regression/rt_diffusion_1x.c creg/
$CP_CMD $G0/regression/rt_diffusion_2x.c creg/
$CP_CMD $G0/regression/rt_diffusion_const_1x.c creg/
$CP_CMD $G0/regression/rt_diffusion_const_2x.c creg/
$CP_CMD $G0/regression/rt_diffusion_const_3x.c creg/
$CP_CMD $G0/regression/rt_diffusion_gen_2x.c creg/
$CP_CMD $G0/regression/rt_diffusion_gen_3x.c creg/
$CP_CMD $G0/regression/rt_diffusion4_const_1x.c creg/
$CP_CMD $G0/regression/rt_diffusion4_const_2x.c creg/
$CP_CMD $G0/regression/rt_diffusion4_const_3x.c creg/
$CP_CMD $G0/regression/rt_diffusion6_const_1x.c creg/
$CP_CMD $G0/regression/rt_diffusion6_const_2x.c creg/
$CP_CMD $G0/regression/rt_diffusion6_const_3x.c creg/
$CP_CMD $G0/regression/rt_escreen_sr.c creg/
$CP_CMD $G0/regression/rt_vlasov_advect_1x.c creg/
$CP_CMD $G0/regression/rt_vlasov_advect_2x.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x1v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_bgk_relax_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_em_advect_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_em_advect_resonant_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_emission_spectrum_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_es_pot_well.c creg/
$CP_CMD $G0/regression/rt_vlasov_es_shock_lbo_1x1v.c creg/
$CP_CMD $G0/regression/rt_vlasov_es_shock_lbo_1x3v.c creg/
$CP_CMD $G0/regression/rt_vlasov_es_shock.c creg/
$CP_CMD $G0/regression/rt_vlasov_freestream_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_freestream_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_landau_damping_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_landau_damping_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_cross_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_cross_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x1v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_lbo_relax_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x1v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x3v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_neut_lbo_wall.c creg/
$CP_CMD $G0/regression/rt_vlasov_poisson_emission_spectrum_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_sheath_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_sheath_bgk_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_bgk_relax_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_freestream.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_neut_bgk_sodshock_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_twostream_1x1v.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_twostream_1x3v.c creg/
$CP_CMD $G0/regression/rt_vlasov_sr_weibel_1x3v.c creg/
$CP_CMD $G0/regression/rt_vlasov_twostream_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_twostream_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_weibel_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_weibel_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_vlasov_weibel_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vlasov_weibel_lbo_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_landau_damping_1x1v_p1.c creg/
$CP_CMD $G0/regression/rt_vp_landau_damping_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_sheath_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_sheath_Aext_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_sheath_Bext_1x2v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_sheath_feedback_sources_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_vp_sheath_feedback_sources_bgk_1x1v_p2.c creg/
$CP_CMD $G0/regression/rt_hyper_vlasov_tm.c creg/
$CP_CMD $G0/regression/rt_vlasov_kerntm.c creg/
$CP_CMD $G0/regression/rt_vlasov_moments.c creg/
$CP_CMD $G0/regression/rt_arg_parse.h creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_hyperbolic_khi_im_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_khi_im_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_can_pb_newtonian_orbits_2x2v_p2.c creg/
$CP_CMD $G0/regression/rt_gr_can_pb_schwarzschild_bh_geodesics.c creg/

$RM_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p1.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1.c
$RM_CMD $G0/regression/rt_can_pb_ex_bgk_surf_flat_sq_ic.c
$RM_CMD $G0/regression/rt_can_pb_ex_bgk_surf_flat.c
$RM_CMD $G0/regression/rt_can_pb_free_streaming_surf_sphere.c
$RM_CMD $G0/regression/rt_can_pb_im_bgk_surf_flat_sq_ic.c
$RM_CMD $G0/regression/rt_can_pb_im_bgk_surf_flat.c
$RM_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_1x1v_p2.c
$RM_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x1v_p1.c
$RM_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x1v_p2.c
$RM_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x2v_p1.c
$RM_CMD $G0/regression/rt_can_pb_neut_bgk_sodshock_im_1x2v_p2.c
$RM_CMD $G0/regression/rt_dg_5m_mom_beach_p2.c
$RM_CMD $G0/regression/rt_dg_5m_mom_beach_p3.c
$RM_CMD $G0/regression/rt_dg_accel_1x1v.c
$RM_CMD $G0/regression/rt_dg_euler_kh_2d.c
$RM_CMD $G0/regression/rt_dg_euler_p_perturbation_p1.c
$RM_CMD $G0/regression/rt_dg_euler_p_perturbation_p2.c
$RM_CMD $G0/regression/rt_dg_euler_sodshock_p1.c
$RM_CMD $G0/regression/rt_dg_euler_sodshock_p2.c
$RM_CMD $G0/regression/rt_dg_maxwell_plane_wave_2d.c
$RM_CMD $G0/regression/rt_dg_maxwell_wg_2d.c
$RM_CMD $G0/regression/rt_diffusion_1x.c
$RM_CMD $G0/regression/rt_diffusion_2x.c
$RM_CMD $G0/regression/rt_diffusion_const_1x.c
$RM_CMD $G0/regression/rt_diffusion_const_2x.c
$RM_CMD $G0/regression/rt_diffusion_const_3x.c
$RM_CMD $G0/regression/rt_diffusion_gen_2x.c
$RM_CMD $G0/regression/rt_diffusion_gen_3x.c
$RM_CMD $G0/regression/rt_diffusion4_const_1x.c
$RM_CMD $G0/regression/rt_diffusion4_const_2x.c
$RM_CMD $G0/regression/rt_diffusion4_const_3x.c
$RM_CMD $G0/regression/rt_diffusion6_const_1x.c
$RM_CMD $G0/regression/rt_diffusion6_const_2x.c
$RM_CMD $G0/regression/rt_diffusion6_const_3x.c
$RM_CMD $G0/regression/rt_escreen_sr.c
$RM_CMD $G0/regression/rt_vlasov_advect_1x.c
$RM_CMD $G0/regression/rt_vlasov_advect_2x.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x1v_p1.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x3v_p1.c
$RM_CMD $G0/regression/rt_vlasov_bgk_relax_1x3v_p2.c
$RM_CMD $G0/regression/rt_vlasov_em_advect_1x3v_p1.c
$RM_CMD $G0/regression/rt_vlasov_em_advect_resonant_1x3v_p1.c
$RM_CMD $G0/regression/rt_vlasov_emission_spectrum_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_es_pot_well.c
$RM_CMD $G0/regression/rt_vlasov_es_shock_lbo_1x1v.c
$RM_CMD $G0/regression/rt_vlasov_es_shock_lbo_1x3v.c
$RM_CMD $G0/regression/rt_vlasov_es_shock.c
$RM_CMD $G0/regression/rt_vlasov_freestream_p1.c
$RM_CMD $G0/regression/rt_vlasov_freestream_p2.c
$RM_CMD $G0/regression/rt_vlasov_landau_damping_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_landau_damping_1x3v_p2.c
$RM_CMD $G0/regression/rt_vlasov_lbo_cross_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_lbo_cross_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x1v_p1.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x3v_p1.c
$RM_CMD $G0/regression/rt_vlasov_lbo_relax_1x3v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x1v_p1.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x2v_p1.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x3v_p1.c
$RM_CMD $G0/regression/rt_vlasov_neut_bgk_sodshock_1x3v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_lbo_sodshock_1x3v_p2.c
$RM_CMD $G0/regression/rt_vlasov_neut_lbo_wall.c
$RM_CMD $G0/regression/rt_vlasov_poisson_emission_spectrum_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_sheath_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_sheath_bgk_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_sr_bgk_relax_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_sr_freestream.c
$RM_CMD $G0/regression/rt_vlasov_sr_neut_bgk_sodshock_1x1v_p2.c
$RM_CMD $G0/regression/rt_vlasov_sr_twostream_1x1v.c
$RM_CMD $G0/regression/rt_vlasov_sr_twostream_1x3v.c
$RM_CMD $G0/regression/rt_vlasov_sr_weibel_1x3v.c
$RM_CMD $G0/regression/rt_vlasov_twostream_p1.c
$RM_CMD $G0/regression/rt_vlasov_twostream_p2.c
$RM_CMD $G0/regression/rt_vlasov_weibel_1x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_weibel_2x2v_p1.c
$RM_CMD $G0/regression/rt_vlasov_weibel_2x2v_p2.c
$RM_CMD $G0/regression/rt_vlasov_weibel_lbo_2x2v_p2.c
$RM_CMD $G0/regression/rt_vp_landau_damping_1x1v_p1.c
$RM_CMD $G0/regression/rt_vp_landau_damping_1x1v_p2.c
$RM_CMD $G0/regression/rt_vp_sheath_1x1v_p2.c
$RM_CMD $G0/regression/rt_vp_sheath_Aext_1x2v_p2.c
$RM_CMD $G0/regression/rt_vp_sheath_Bext_1x2v_p2.c
$RM_CMD $G0/regression/rt_vp_sheath_feedback_sources_1x1v_p2.c
$RM_CMD $G0/regression/rt_vp_sheath_feedback_sources_bgk_1x1v_p2.c
$RM_CMD $G0/regression/rt_hyper_vlasov_tm.c
$RM_CMD $G0/regression/rt_vlasov_kerntm.c
$RM_CMD $G0/regression/rt_vlasov_moments.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_hyperbolic_khi_im_2x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_bgk_surf_sphere_khi_im_2x2v_p2.c
$RM_CMD $G0/regression/rt_can_pb_newtonian_orbits_2x2v_p2.c
$RM_CMD $G0/regression/rt_gr_can_pb_schwarzschild_bh_geodesics.c

# Lua regression tests
mkdir -p luareg
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_ex_bgk_surf_flat_sq_ic.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_ex_bgk_surf_flat.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_free_streaming_surf_sphere.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_im_bgk_surf_flat_sq_ic.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_im_bgk_surf_flat.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x1v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_5m_mom_beach_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_5m_mom_beach_p3.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_advect_1x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_advect_1x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_advect_2x_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_advect_2x_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_euler_p_perturbation_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_euler_p_perturbation_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_euler_sodshock_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_euler_sodshock_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_hasegawa_mima_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_hasegawa_mima_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_hasegawa_wakatani_low_adiabat_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_hasegawa_wakatani_moderate_adiabat_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_incompress_euler_double_shear_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_incompress_euler_double_shear_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_incompress_euler_vortex_waltz_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_incompress_euler_vortex_waltz_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_maxwell_plane_wave_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_maxwell_wg_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_modified_hasegawa_wakatani_low_adiabat_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_dg_modified_hasegawa_wakatani_moderate_adiabat_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x1v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_buneman_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_em_advect_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_em_advect_resonant_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_es_pot_well.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_es_shock_lbo_1x1v.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_es_shock_lbo_1x3v.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_es_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_freestream_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_freestream_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_landau_damping_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_landau_damping_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_cross_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_cross_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x1v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x1v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x3v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x3v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_neut_lbo_wall.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sheath_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sheath_bgk_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_bgk_relax_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_freestream.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_neut_bgk_sodshock_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_twostream_1x1v.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_twostream_1x3v.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_sr_weibel_1x3v.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_twostream_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_twostream_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_weibel_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_weibel_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_weibel_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vlasov_weibel_lbo_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_buneman_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_landau_damping_1x1v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_landau_damping_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_sheath_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_sheath_Aext_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_sheath_Bext_1x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_sheath_feedback_sources_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_vp_sheath_feedback_sources_bgk_1x1v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_hyperbolic_khi_im_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_khi_im_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_distorted_mesh_ic_2x2v_p2.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_newtonian_orbits.lua luareg/
$CP_CMD $G0/regression/lua/rt_can_pb_trapped_well_diffusion_1x1v_p2.lua luareg/

$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_2x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_ex_bgk_surf_flat_sq_ic.lua
$RM_CMD $G0/regression/lua/rt_can_pb_ex_bgk_surf_flat.lua
$RM_CMD $G0/regression/lua/rt_can_pb_free_streaming_surf_sphere.lua
$RM_CMD $G0/regression/lua/rt_can_pb_im_bgk_surf_flat_sq_ic.lua
$RM_CMD $G0/regression/lua/rt_can_pb_im_bgk_surf_flat.lua
$RM_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x1v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_can_pb_neut_bgk_sodshock_im_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_5m_mom_beach_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_5m_mom_beach_p3.lua
$RM_CMD $G0/regression/lua/rt_dg_advect_1x_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_advect_1x_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_advect_2x_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_advect_2x_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_euler_p_perturbation_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_euler_p_perturbation_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_euler_sodshock_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_euler_sodshock_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_hasegawa_mima_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_hasegawa_mima_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_hasegawa_wakatani_low_adiabat_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_hasegawa_wakatani_moderate_adiabat_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_incompress_euler_double_shear_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_incompress_euler_double_shear_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_incompress_euler_vortex_waltz_p1.lua
$RM_CMD $G0/regression/lua/rt_dg_incompress_euler_vortex_waltz_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_maxwell_plane_wave_2d.lua
$RM_CMD $G0/regression/lua/rt_dg_maxwell_wg_2d.lua
$RM_CMD $G0/regression/lua/rt_dg_modified_hasegawa_wakatani_low_adiabat_p2.lua
$RM_CMD $G0/regression/lua/rt_dg_modified_hasegawa_wakatani_moderate_adiabat_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x1v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_bgk_relax_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_buneman_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_em_advect_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_em_advect_resonant_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_es_pot_well.lua
$RM_CMD $G0/regression/lua/rt_vlasov_es_shock_lbo_1x1v.lua
$RM_CMD $G0/regression/lua/rt_vlasov_es_shock_lbo_1x3v.lua
$RM_CMD $G0/regression/lua/rt_vlasov_es_shock.lua
$RM_CMD $G0/regression/lua/rt_vlasov_freestream_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_freestream_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_landau_damping_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_landau_damping_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_cross_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_cross_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x1v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_lbo_relax_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x1v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x3v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_bgk_sodshock_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_lbo_sodshock_1x3v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_neut_lbo_wall.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sheath_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sheath_bgk_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_bgk_relax_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_freestream.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_neut_bgk_sodshock_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_twostream_1x1v.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_twostream_1x3v.lua
$RM_CMD $G0/regression/lua/rt_vlasov_sr_weibel_1x3v.lua
$RM_CMD $G0/regression/lua/rt_vlasov_twostream_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_twostream_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_weibel_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_weibel_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_vlasov_weibel_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vlasov_weibel_lbo_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_buneman_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_landau_damping_1x1v_p1.lua
$RM_CMD $G0/regression/lua/rt_vp_landau_damping_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_sheath_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_sheath_Aext_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_sheath_Bext_1x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_sheath_feedback_sources_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_vp_sheath_feedback_sources_bgk_1x1v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_hyperbolic_khi_im_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_bgk_surf_sphere_khi_im_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_distorted_mesh_ic_2x2v_p2.lua
$RM_CMD $G0/regression/lua/rt_can_pb_newtonian_orbits.lua
$RM_CMD $G0/regression/lua/rt_can_pb_trapped_well_diffusion_1x1v_p2.lua