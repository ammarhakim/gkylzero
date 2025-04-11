#!/bin/sh

CP_CMD=git mv
RM_CMD=:
G0=..

# data
mkdir -p data/regression
$CP_CMD $G0/data/regression/euler_riem_2d_hllc-euler_0.gkyl data/regression/

$RM_CMD $G0/data/regression/euler_riem_2d_hllc-euler_0.gkyl

# proofs
mkdir -p proofs/finite_volume
$CP_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_minmod_symmetry.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_minmod_tvd.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_monotonized_centered_symmetry.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_monotonized_centered_tvd.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_superbee_tvd.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_limiter_van_leer_symmetry.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_linear_advection_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_linear_advection_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_cfl_stability.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_local_lipschitz.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_flux_conservation.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_strict_hyperbolicity.rkt proofs/finite_volume/
$CP_CMD $G0/proofs/prover_core.rkt proofs/

$RM_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_inviscid_burgers_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_x_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_minmod_symmetry.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_minmod_tvd.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_monotonized_centered_symmetry.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_monotonized_centered_tvd.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_superbee_tvd.rkt
$RM_CMD $G0/proofs/finite_volume/proof_limiter_van_leer_symmetry.rkt
$RM_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_linear_advection_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_linear_advection_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_linear_advection_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_cfl_stability.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_local_lipschitz.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_flux_conservation.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_hyperbolicity.rkt
$RM_CMD $G0/proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_strict_hyperbolicity.rkt
$RM_CMD $G0/proofs/prover_core.rkt

# kernels
mkdir -p kernels/fem_poisson
$CP_CMD $G0/kernels/fem_poisson/*.h kernels/fem_poisson/
$CP_CMD $G0/kernels/fem_poisson/*.c kernels/fem_poisson/

$RM_CMD $G0/kernels/fem_poisson/*.h
$RM_CMD $G0/kernels/fem_poisson/*.c

# zero
mkdir -p zero
$CP_CMD $G0/zero/gkyl_gr_blackhole.h  zero/
$CP_CMD $G0/zero/gkyl_gr_minkowski.h  zero/
$CP_CMD $G0/zero/gkyl_gr_neutronstar.h  zero/
$CP_CMD $G0/zero/gkyl_gr_spacetime.h  zero/
$CP_CMD $G0/zero/gkyl_gr_spacetime_diff.h  zero/
$CP_CMD $G0/zero/gkyl_kep_scheme.h zero/
$CP_CMD $G0/zero/gkyl_mhd_src.h zero/
$CP_CMD $G0/zero/gkyl_moment_braginskii.h zero/
$CP_CMD $G0/zero/gkyl_moment_em_coupling.h zero/
$CP_CMD $G0/zero/gkyl_moment_em_coupling_priv.h zero/
$CP_CMD $G0/zero/gkyl_moment_non_ideal_priv.h zero/
$CP_CMD $G0/zero/gkyl_moment_prim_mhd.h zero/
$CP_CMD $G0/zero/gkyl_moment_prim_sr_euler.h zero/
$CP_CMD $G0/zero/gkyl_mp_scheme.h zero/
$CP_CMD $G0/zero/gkyl_sources_explicit_priv.h zero/
$CP_CMD $G0/zero/gkyl_sources_implicit_priv.h zero/
$CP_CMD $G0/zero/gkyl_ten_moment_grad_closure.h zero/
$CP_CMD $G0/zero/gkyl_wave_geom.h zero/
$CP_CMD $G0/zero/gkyl_wave_geom_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wave_prop.h  zero/
$CP_CMD $G0/zero/gkyl_wv_advect.h  zero/
$CP_CMD $G0/zero/gkyl_wv_apply_bc.h  zero/
$CP_CMD $G0/zero/gkyl_wv_burgers.h  zero/
$CP_CMD $G0/zero/gkyl_wv_canonical_pb_fluid.h  zero/
$CP_CMD $G0/zero/gkyl_wv_canonical_pb_fluid_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_coldfluid.h  zero/
$CP_CMD $G0/zero/gkyl_wv_eqn.h  zero/
$CP_CMD $G0/zero/gkyl_wv_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_euler_mixture.h  zero/
$CP_CMD $G0/zero/gkyl_wv_euler_mixture_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_euler_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_euler_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_euler_tetrad.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_euler_tetrad_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_maxwell.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_maxwell_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_maxwell_tetrad.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_maxwell_tetrad_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_medium.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_medium_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_tetrad.h  zero/
$CP_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_tetrad_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_iso_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_iso_euler_mixture.h  zero/
$CP_CMD $G0/zero/gkyl_wv_iso_euler_mixture_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_maxwell.h  zero/
$CP_CMD $G0/zero/gkyl_wv_maxwell_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_mhd.h  zero/
$CP_CMD $G0/zero/gkyl_wv_reactive_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_reactive_euler_priv.h  zero/
$CP_CMD $G0/zero/gkyl_wv_sr_euler.h  zero/
$CP_CMD $G0/zero/gkyl_wv_ten_moment.h  zero/
$CP_CMD $G0/zero/gkyl_wv_ten_moment_priv.h  zero/
$CP_CMD $G0/zero/gr_blackhole.c  zero/
$CP_CMD $G0/zero/gr_minkowski.c  zero/
$CP_CMD $G0/zero/gr_neutronstar.c  zero/
$CP_CMD $G0/zero/gr_spacetime.c  zero/
$CP_CMD $G0/zero/gr_spacetime_diff.c  zero/
$CP_CMD $G0/zero/kep_scheme.c zero/
$CP_CMD $G0/zero/mhd_src.c zero/
$CP_CMD $G0/zero/moment_braginskii.c zero/
$CP_CMD $G0/zero/moment_em_coupling.c zero/
$CP_CMD $G0/zero/moment_prim_mhd.c zero/
$CP_CMD $G0/zero/moment_prim_sr_euler.c zero/
$CP_CMD $G0/zero/mp_scheme.c zero/
$CP_CMD $G0/zero/sources_explicit.c zero/
$CP_CMD $G0/zero/sources_implicit.c zero/
$CP_CMD $G0/zero/ten_moment_grad_closure.c zero/
$CP_CMD $G0/zero/wave_geom.c  zero/
$CP_CMD $G0/zero/wave_geom_cu.cu  zero/
$CP_CMD $G0/zero/wave_prop.c  zero/
$CP_CMD $G0/zero/wv_advect.c  zero/
$CP_CMD $G0/zero/wv_apply_bc.c  zero/
$CP_CMD $G0/zero/wv_burgers.c  zero/
$CP_CMD $G0/zero/wv_canonical_pb_fluid.c  zero/
$CP_CMD $G0/zero/wv_coldfluid.c  zero/
$CP_CMD $G0/zero/wv_eqn.c  zero/
$CP_CMD $G0/zero/wv_euler.c  zero/
$CP_CMD $G0/zero/wv_euler_cu.cu  zero/
$CP_CMD $G0/zero/wv_euler_mixture.c  zero/
$CP_CMD $G0/zero/wv_gr_euler.c  zero/
$CP_CMD $G0/zero/wv_gr_euler_tetrad.c  zero/
$CP_CMD $G0/zero/wv_gr_maxwell.c  zero/
$CP_CMD $G0/zero/wv_gr_maxwell_tetrad.c  zero/
$CP_CMD $G0/zero/wv_gr_medium.c  zero/
$CP_CMD $G0/zero/wv_gr_ultra_rel_euler.c  zero/
$CP_CMD $G0/zero/wv_gr_ultra_rel_euler_tetrad.c  zero/
$CP_CMD $G0/zero/wv_iso_euler.c  zero/
$CP_CMD $G0/zero/wv_iso_euler_mixture.c  zero/
$CP_CMD $G0/zero/wv_maxwell.c  zero/
$CP_CMD $G0/zero/wv_maxwell_cu.cu  zero/
$CP_CMD $G0/zero/wv_mhd.c  zero/
$CP_CMD $G0/zero/wv_reactive_euler.c  zero/
$CP_CMD $G0/zero/wv_sr_euler.c  zero/
$CP_CMD $G0/zero/wv_ten_moment.c  zero/
$CP_CMD $G0/zero/wv_ten_moment_cu.cu zero/
$CP_CMD $G0/zero/gkyl_fem_poisson_bctype.h zero/
$CP_CMD $G0/zero/gkyl_fem_poisson_priv.h zero/
$CP_CMD $G0/zero/gkyl_fem_poisson.h zero/
$CP_CMD $G0/zero/fem_poisson_cu.cu zero/
$CP_CMD $G0/zero/fem_poisson.c zero/
$CP_CMD $G0/zero/gkyl_wv_advect_priv.h zero/
$CP_CMD $G0/zero/gkyl_wv_burgers_priv.h zero/
$CP_CMD $G0/zero/gkyl_wv_iso_euler_priv.h zero/

$RM_CMD $G0/zero/gkyl_gr_blackhole.h
$RM_CMD $G0/zero/gkyl_gr_minkowski.h
$RM_CMD $G0/zero/gkyl_gr_neutronstar.h
$RM_CMD $G0/zero/gkyl_gr_spacetime.h
$RM_CMD $G0/zero/gkyl_gr_spacetime_diff.h
$RM_CMD $G0/zero/gkyl_kep_scheme.h
$RM_CMD $G0/zero/gkyl_mhd_src.h
$RM_CMD $G0/zero/gkyl_moment_braginskii.h
$RM_CMD $G0/zero/gkyl_moment_em_coupling.h
$RM_CMD $G0/zero/gkyl_moment_em_coupling_priv.h
$RM_CMD $G0/zero/gkyl_moment_non_ideal_priv.h
$RM_CMD $G0/zero/gkyl_moment_prim_mhd.h
$RM_CMD $G0/zero/gkyl_moment_prim_sr_euler.h
$RM_CMD $G0/zero/gkyl_mp_scheme.h
$RM_CMD $G0/zero/gkyl_sources_explicit_priv.h
$RM_CMD $G0/zero/gkyl_sources_implicit_priv.h
$RM_CMD $G0/zero/gkyl_ten_moment_grad_closure.h
$RM_CMD $G0/zero/gkyl_wave_geom.h
$RM_CMD $G0/zero/gkyl_wave_geom_priv.h
$RM_CMD $G0/zero/gkyl_wave_prop.h
$RM_CMD $G0/zero/gkyl_wv_advect.h
$RM_CMD $G0/zero/gkyl_wv_apply_bc.h
$RM_CMD $G0/zero/gkyl_wv_burgers.h
$RM_CMD $G0/zero/gkyl_wv_canonical_pb_fluid.h
$RM_CMD $G0/zero/gkyl_wv_canonical_pb_fluid_priv.h
$RM_CMD $G0/zero/gkyl_wv_coldfluid.h
$RM_CMD $G0/zero/gkyl_wv_eqn.h
$RM_CMD $G0/zero/gkyl_wv_euler.h
$RM_CMD $G0/zero/gkyl_wv_euler_mixture.h
$RM_CMD $G0/zero/gkyl_wv_euler_mixture_priv.h
$RM_CMD $G0/zero/gkyl_wv_euler_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_euler.h
$RM_CMD $G0/zero/gkyl_wv_gr_euler_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_euler_tetrad.h
$RM_CMD $G0/zero/gkyl_wv_gr_euler_tetrad_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_maxwell.h
$RM_CMD $G0/zero/gkyl_wv_gr_maxwell_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_maxwell_tetrad.h
$RM_CMD $G0/zero/gkyl_wv_gr_maxwell_tetrad_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_medium.h
$RM_CMD $G0/zero/gkyl_wv_gr_medium_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler.h
$RM_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_priv.h
$RM_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_tetrad.h
$RM_CMD $G0/zero/gkyl_wv_gr_ultra_rel_euler_tetrad_priv.h
$RM_CMD $G0/zero/gkyl_wv_iso_euler.h
$RM_CMD $G0/zero/gkyl_wv_iso_euler_mixture.h
$RM_CMD $G0/zero/gkyl_wv_iso_euler_mixture_priv.h
$RM_CMD $G0/zero/gkyl_wv_maxwell.h
$RM_CMD $G0/zero/gkyl_wv_maxwell_priv.h
$RM_CMD $G0/zero/gkyl_wv_mhd.h
$RM_CMD $G0/zero/gkyl_wv_reactive_euler.h
$RM_CMD $G0/zero/gkyl_wv_reactive_euler_priv.h
$RM_CMD $G0/zero/gkyl_wv_sr_euler.h
$RM_CMD $G0/zero/gkyl_wv_ten_moment.h
$RM_CMD $G0/zero/gkyl_wv_ten_moment_priv.h
$RM_CMD $G0/zero/gr_blackhole.c
$RM_CMD $G0/zero/gr_minkowski.c
$RM_CMD $G0/zero/gr_neutronstar.c
$RM_CMD $G0/zero/gr_spacetime.c
$RM_CMD $G0/zero/gr_spacetime_diff.c
$RM_CMD $G0/zero/kep_scheme.c
$RM_CMD $G0/zero/mhd_src.c
$RM_CMD $G0/zero/moment_braginskii.c
$RM_CMD $G0/zero/moment_em_coupling.c
$RM_CMD $G0/zero/moment_prim_mhd.c
$RM_CMD $G0/zero/moment_prim_sr_euler.c
$RM_CMD $G0/zero/mp_scheme.c
$RM_CMD $G0/zero/sources_explicit.c
$RM_CMD $G0/zero/sources_implicit.c
$RM_CMD $G0/zero/ten_moment_grad_closure.c
$RM_CMD $G0/zero/wave_geom.c
$RM_CMD $G0/zero/wave_geom_cu.cu
$RM_CMD $G0/zero/wave_prop.c
$RM_CMD $G0/zero/wv_advect.c
$RM_CMD $G0/zero/wv_apply_bc.c
$RM_CMD $G0/zero/wv_burgers.c
$RM_CMD $G0/zero/wv_canonical_pb_fluid.c
$RM_CMD $G0/zero/wv_coldfluid.c
$RM_CMD $G0/zero/wv_eqn.c
$RM_CMD $G0/zero/wv_euler.c
$RM_CMD $G0/zero/wv_euler_cu.cu
$RM_CMD $G0/zero/wv_euler_mixture.c
$RM_CMD $G0/zero/wv_gr_euler.c
$RM_CMD $G0/zero/wv_gr_euler_tetrad.c
$RM_CMD $G0/zero/wv_gr_maxwell.c
$RM_CMD $G0/zero/wv_gr_maxwell_tetrad.c
$RM_CMD $G0/zero/wv_gr_medium.c
$RM_CMD $G0/zero/wv_gr_ultra_rel_euler.c
$RM_CMD $G0/zero/wv_gr_ultra_rel_euler_tetrad.c
$RM_CMD $G0/zero/wv_iso_euler.c
$RM_CMD $G0/zero/wv_iso_euler_mixture.c
$RM_CMD $G0/zero/wv_maxwell.c
$RM_CMD $G0/zero/wv_maxwell_cu.cu
$RM_CMD $G0/zero/wv_mhd.c
$RM_CMD $G0/zero/wv_reactive_euler.c
$RM_CMD $G0/zero/wv_sr_euler.c
$RM_CMD $G0/zero/wv_ten_moment.c
$RM_CMD $G0/zero/wv_ten_moment_cu.cu
$RM_CMD $G0/zero/gkyl_fem_poisson_bctype.h
$RM_CMD $G0/zero/gkyl_fem_poisson_priv.h
$RM_CMD $G0/zero/gkyl_fem_poisson.h
$RM_CMD $G0/zero/fem_poisson_cu.cu
$RM_CMD $G0/zero/fem_poisson.c
$RM_CMD $G0/zero/gkyl_wv_advect_priv.h
$RM_CMD $G0/zero/gkyl_wv_burgers_priv.h
$RM_CMD $G0/zero/gkyl_wv_iso_euler_priv.h

# app
mkdir -p apps
$CP_CMD $G0/apps/gkyl_moment.h apps/
$CP_CMD $G0/apps/gkyl_moment_lw.h apps/
$CP_CMD $G0/apps/gkyl_moment_multib.h apps/
$CP_CMD $G0/apps/gkyl_moment_multib_priv.h apps/
$CP_CMD $G0/apps/gkyl_moment_priv.h apps/
$CP_CMD $G0/apps/mom_coupling.c apps/
$CP_CMD $G0/apps/mom_field.c apps/
$CP_CMD $G0/apps/mom_priv.c apps/
$CP_CMD $G0/apps/mom_species.c apps/
$CP_CMD $G0/apps/mom_update_one_step.c apps/
$CP_CMD $G0/apps/mom_update_ssp_rk.c apps/
$CP_CMD $G0/apps/moment.c apps/
$CP_CMD $G0/apps/moment_lw.c apps/
$CP_CMD $G0/apps/moment_multib.c apps/

$RM_CMD $G0/apps/gkyl_moment.h
$RM_CMD $G0/apps/gkyl_moment_lw.h
$RM_CMD $G0/apps/gkyl_moment_multib.h
$RM_CMD $G0/apps/gkyl_moment_multib_priv.h
$RM_CMD $G0/apps/gkyl_moment_priv.h
$RM_CMD $G0/apps/mom_coupling.c
$RM_CMD $G0/apps/mom_field.c
$RM_CMD $G0/apps/mom_priv.c
$RM_CMD $G0/apps/mom_species.c
$RM_CMD $G0/apps/mom_update_one_step.c
$RM_CMD $G0/apps/mom_update_ssp_rk.c
$RM_CMD $G0/apps/moment.c
$RM_CMD $G0/apps/moment_lw.c
$RM_CMD $G0/apps/moment_multib.c

# AMR
mkdir -p amr
$CP_CMD $G0/amr/amr_block_coupled.c amr/
$CP_CMD $G0/amr/amr_block.c amr/
$CP_CMD $G0/amr/amr_core_euler_mixture.c amr/
$CP_CMD $G0/amr/amr_core_euler.c amr/
$CP_CMD $G0/amr/amr_core_five_moment.c amr/
$CP_CMD $G0/amr/amr_core_gr_euler.c amr/
$CP_CMD $G0/amr/amr_core_ten_moment.c amr/
$CP_CMD $G0/amr/amr_patch_coupled.c amr/
$CP_CMD $G0/amr/amr_patch.c amr/
$CP_CMD $G0/amr/gkyl_amr_block_coupled_priv.h amr/
$CP_CMD $G0/amr/gkyl_amr_block_priv.h amr/
$CP_CMD $G0/amr/gkyl_amr_core.h amr/
$CP_CMD $G0/amr/gkyl_amr_patch_coupled_priv.h amr/
$CP_CMD $G0/amr/gkyl_amr_patch_priv.h amr/

$RM_CMD $G0/amr/amr_block_coupled.c
$RM_CMD $G0/amr/amr_block.c
$RM_CMD $G0/amr/amr_core_euler_mixture.c
$RM_CMD $G0/amr/amr_core_euler.c
$RM_CMD $G0/amr/amr_core_five_moment.c
$RM_CMD $G0/amr/amr_core_gr_euler.c
$RM_CMD $G0/amr/amr_core_ten_moment.c
$RM_CMD $G0/amr/amr_patch_coupled.c
$RM_CMD $G0/amr/amr_patch.c
$RM_CMD $G0/amr/gkyl_amr_block_coupled_priv.h
$RM_CMD $G0/amr/gkyl_amr_block_priv.h
$RM_CMD $G0/amr/gkyl_amr_core.h
$RM_CMD $G0/amr/gkyl_amr_patch_coupled_priv.h
$RM_CMD $G0/amr/gkyl_amr_patch_priv.h

# unit
mkdir -p unit
$CP_CMD $G0/unit/ctest_wave_geom.c unit/
$CP_CMD $G0/unit/ctest_wave_geom_helpers.c unit/
$CP_CMD $G0/unit/ctest_wv_apply_bc.c unit/
$CP_CMD $G0/unit/ctest_wv_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_euler_mixture.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_euler_tetrad.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_maxwell.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_maxwell_tetrad.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_medium.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_ultra_rel_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_gr_ultra_rel_euler_tetrad.c unit/
$CP_CMD $G0/unit/ctest_wv_iso_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_iso_euler_mixture.c unit/
$CP_CMD $G0/unit/ctest_wv_maxwell.c unit/
$CP_CMD $G0/unit/ctest_wv_mhd.c unit/
$CP_CMD $G0/unit/ctest_wv_reactive_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_sr_euler.c unit/
$CP_CMD $G0/unit/ctest_wv_ten_moment.c unit/
$CP_CMD $G0/unit/ctest_gr_spacetime.c unit/
$CP_CMD $G0/unit/ctest_wv_euler_cu.cu unit/
$CP_CMD $G0/unit/ctest_wv_maxwell_cu.cu unit/
$CP_CMD $G0/unit/ctest_wv_ten_moment_cu.cu unit/
$CP_CMD $G0/unit/ctest_wave_geom_cu.cu unit/
$CP_CMD $G0/unit/ctest_fem_poisson.c unit/
$CP_CMD $G0/unit/ctest_fem_poisson_vareps.c unit/
$CP_CMD $G0/unit/ctest_fem_helmholtz.c unit/

$RM_CMD $G0/unit/ctest_wave_geom.c
$RM_CMD $G0/unit/ctest_wave_geom_helpers.c
$RM_CMD $G0/unit/ctest_wv_apply_bc.c
$RM_CMD $G0/unit/ctest_wv_euler.c
$RM_CMD $G0/unit/ctest_wv_euler_mixture.c
$RM_CMD $G0/unit/ctest_wv_gr_euler.c
$RM_CMD $G0/unit/ctest_wv_gr_euler_tetrad.c
$RM_CMD $G0/unit/ctest_wv_gr_maxwell.c
$RM_CMD $G0/unit/ctest_wv_gr_maxwell_tetrad.c
$RM_CMD $G0/unit/ctest_wv_gr_medium.c
$RM_CMD $G0/unit/ctest_wv_gr_ultra_rel_euler.c
$RM_CMD $G0/unit/ctest_wv_gr_ultra_rel_euler_tetrad.c
$RM_CMD $G0/unit/ctest_wv_iso_euler.c
$RM_CMD $G0/unit/ctest_wv_iso_euler_mixture.c
$RM_CMD $G0/unit/ctest_wv_maxwell.c
$RM_CMD $G0/unit/ctest_wv_mhd.c
$RM_CMD $G0/unit/ctest_wv_reactive_euler.c
$RM_CMD $G0/unit/ctest_wv_sr_euler.c
$RM_CMD $G0/unit/ctest_wv_ten_moment.c
$RM_CMD $G0/unit/ctest_gr_spacetime.c
$RM_CMD $G0/unit/ctest_wv_euler_cu.cu
$RM_CMD $G0/unit/ctest_wv_maxwell_cu.cu
$RM_CMD $G0/unit/ctest_wv_ten_moment_cu.cu
$RM_CMD $G0/unit/ctest_wave_geom_cu.cu
$RM_CMD $G0/unit/ctest_fem_poisson.c
$RM_CMD $G0/unit/ctest_fem_poisson_vareps.c
$RM_CMD $G0/unit/ctest_fem_helmholtz.c

# C regression tests
mkdir -p creg
$CP_CMD $G0/regression/rt_10m_burch.c creg/
$CP_CMD $G0/regression/rt_10m_burch_grad_closure.c creg/
$CP_CMD $G0/regression/rt_10m_expanding.c creg/
$CP_CMD $G0/regression/rt_10m_expanding_axi_sodshock.c creg/
$CP_CMD $G0/regression/rt_10m_expanding_sodshock.c creg/
$CP_CMD $G0/regression/rt_10m_gem.c creg/
$CP_CMD $G0/regression/rt_10m_gem_grad_closure.c creg/
$CP_CMD $G0/regression/rt_10m_lhdi.c creg/
$CP_CMD $G0/regression/rt_10m_lhdi_grad_closure.c creg/
$CP_CMD $G0/regression/rt_10m_par_firehose.c creg/
$CP_CMD $G0/regression/rt_10m_par_firehose_grad_closure.c creg/
$CP_CMD $G0/regression/rt_10m_riem.c creg/
$CP_CMD $G0/regression/rt_10m_riem_grad_closure.c creg/
$CP_CMD $G0/regression/rt_10m_sodshock.c creg/
$CP_CMD $G0/regression/rt_10m_sodshock_lax.c creg/
$CP_CMD $G0/regression/rt_5m_burch.c creg/
$CP_CMD $G0/regression/rt_5m_elc_heat_flux.c creg/
$CP_CMD $G0/regression/rt_5m_em_advect.c creg/
$CP_CMD $G0/regression/rt_5m_em_advect_resonant.c creg/
$CP_CMD $G0/regression/rt_5m_expanding.c creg/
$CP_CMD $G0/regression/rt_5m_expanding_axi_sodshock.c creg/
$CP_CMD $G0/regression/rt_5m_expanding_sodshock.c creg/
$CP_CMD $G0/regression/rt_5m_friction.c creg/
$CP_CMD $G0/regression/rt_5m_gem.c creg/
$CP_CMD $G0/regression/rt_5m_hartmann.c creg/
$CP_CMD $G0/regression/rt_5m_ion_heat_flux.c creg/
$CP_CMD $G0/regression/rt_5m_mom_beach.c creg/
$CP_CMD $G0/regression/rt_5m_riem.c creg/
$CP_CMD $G0/regression/rt_5m_rt.c creg/
$CP_CMD $G0/regression/rt_advect_wv.c creg/
$CP_CMD $G0/regression/rt_advect_wv_mp.c creg/
$CP_CMD $G0/regression/rt_burgers_shock.c creg/
$CP_CMD $G0/regression/rt_burgers_shock_mp.c creg/
$CP_CMD $G0/regression/rt_coldfluid_beach.c creg/
$CP_CMD $G0/regression/rt_coldfluid_clouda.c creg/
$CP_CMD $G0/regression/rt_coldfluid_em_coupling.c creg/
$CP_CMD $G0/regression/rt_euler_axi_sodshock.c creg/
$CP_CMD $G0/regression/rt_euler_axi_vac_riem.c creg/
$CP_CMD $G0/regression/rt_euler_bump_in_channel.c creg/
$CP_CMD $G0/regression/rt_euler_c2p_sodshock.c creg/
$CP_CMD $G0/regression/rt_euler_cart_axi_sodshock.c creg/
$CP_CMD $G0/regression/rt_euler_kh_2d.c creg/
$CP_CMD $G0/regression/rt_euler_mixture_fedkiw_shock.c creg/
$CP_CMD $G0/regression/rt_euler_mixture_shock_bubble.c creg/
$CP_CMD $G0/regression/rt_euler_noh_1d.c creg/
$CP_CMD $G0/regression/rt_euler_p_perturbation.c creg/
$CP_CMD $G0/regression/rt_euler_reflect_2d.c creg/
$CP_CMD $G0/regression/rt_euler_riem_2d_hll.c creg/
$CP_CMD $G0/regression/rt_euler_riem_2d_hllc.c creg/
$CP_CMD $G0/regression/rt_euler_riem_2d_lax.c creg/
$CP_CMD $G0/regression/rt_euler_riem_2d_roe.c creg/
$CP_CMD $G0/regression/rt_euler_riem_3d.c creg/
$CP_CMD $G0/regression/rt_euler_rt.c creg/
$CP_CMD $G0/regression/rt_euler_sodshock.c creg/
$CP_CMD $G0/regression/rt_euler_sodshock_lax.c creg/
$CP_CMD $G0/regression/rt_euler_sodshock_mp.c creg/
$CP_CMD $G0/regression/rt_euler_superwedge.c creg/
$CP_CMD $G0/regression/rt_euler_vac.c creg/
$CP_CMD $G0/regression/rt_euler_vac_riem_1d.c creg/
$CP_CMD $G0/regression/rt_euler_wave_2d_kep.c creg/
$CP_CMD $G0/regression/rt_euler_wave_2d_mp.c creg/
$CP_CMD $G0/regression/rt_euler_wave_2d_wv.c creg/
$CP_CMD $G0/regression/rt_euler_wedge_sodshock.c creg/
$CP_CMD $G0/regression/rt_gr_bhl_spinning.c creg/
$CP_CMD $G0/regression/rt_gr_bhl_spinning_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_bhl_static.c creg/
$CP_CMD $G0/regression/rt_gr_bhl_static_neutronstar.c creg/
$CP_CMD $G0/regression/rt_gr_bhl_static_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_blackhole_spinning.c creg/
$CP_CMD $G0/regression/rt_gr_blackhole_static.c creg/
$CP_CMD $G0/regression/rt_gr_bz_monopole_fast.c creg/
$CP_CMD $G0/regression/rt_gr_bz_monopole_fast_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_bz_monopole_slow.c creg/
$CP_CMD $G0/regression/rt_gr_bz_monopole_slow_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_current_sheet.c creg/
$CP_CMD $G0/regression/rt_gr_current_sheet_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_einstein_plane_shock.c creg/
$CP_CMD $G0/regression/rt_gr_kh_2d.c creg/
$CP_CMD $G0/regression/rt_gr_mild_shock.c creg/
$CP_CMD $G0/regression/rt_gr_mild_shock_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_perturbed_density.c creg/
$CP_CMD $G0/regression/rt_gr_quadrants_2d.c creg/
$CP_CMD $G0/regression/rt_gr_strong_blast.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_bhl_spinning.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_bhl_spinning_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_bhl_static.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_bhl_static_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_shock.c creg/
$CP_CMD $G0/regression/rt_gr_ultra_rel_shock_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_wald_magnetosphere_spinning.c creg/
$CP_CMD $G0/regression/rt_gr_wald_magnetosphere_spinning_tetrad.c creg/
$CP_CMD $G0/regression/rt_gr_wald_magnetosphere_static.c creg/
$CP_CMD $G0/regression/rt_gr_wald_magnetosphere_static_neutronstar.c creg/
$CP_CMD $G0/regression/rt_gr_wald_magnetosphere_static_tetrad.c creg/
$CP_CMD $G0/regression/rt_iso_euler_friction.c creg/
$CP_CMD $G0/regression/rt_iso_euler_hartmann.c creg/
$CP_CMD $G0/regression/rt_iso_euler_mixture_fedkiw_shock.c creg/
$CP_CMD $G0/regression/rt_iso_euler_mixture_shock_bubble.c creg/
$CP_CMD $G0/regression/rt_iso_euler_sodshock.c creg/
$CP_CMD $G0/regression/rt_iso_euler_sodshock_lax.c creg/
$CP_CMD $G0/regression/rt_iso_gem.c creg/
$CP_CMD $G0/regression/rt_maxwell_annular_wg.c creg/
$CP_CMD $G0/regression/rt_maxwell_annulus.c creg/
$CP_CMD $G0/regression/rt_maxwell_axi_wg.c creg/
$CP_CMD $G0/regression/rt_maxwell_expanding.c creg/
$CP_CMD $G0/regression/rt_maxwell_expanding_2d.c creg/
$CP_CMD $G0/regression/rt_maxwell_plane_wave_1d.c creg/
$CP_CMD $G0/regression/rt_maxwell_plane_wave_1d_mp.c creg/
$CP_CMD $G0/regression/rt_maxwell_plane_wave_2d.c creg/
$CP_CMD $G0/regression/rt_maxwell_plane_wave_2d_mp.c creg/
$CP_CMD $G0/regression/rt_maxwell_reflect_2d.c creg/
$CP_CMD $G0/regression/rt_maxwell_wg_2d.c creg/
$CP_CMD $G0/regression/rt_mhd_brio_wu.c creg/
$CP_CMD $G0/regression/rt_mhd_ot.c creg/
$CP_CMD $G0/regression/rt_mhd_rj2.c creg/
$CP_CMD $G0/regression/rt_reactive_euler_detonation.c creg/
$CP_CMD $G0/regression/rt_sr_euler_KH_2d.c creg/
$CP_CMD $G0/regression/rt_sr_euler_riem_2d.c creg/
$CP_CMD $G0/regression/rt_sr_euler_sodshock.c creg/
$CP_CMD $G0/regression/rt_euler_multiblock.c creg/
$CP_CMD $G0/regression/rt_multib_euler_2d.c creg/
$CP_CMD $G0/regression/rt_arg_parse.h creg/

$RM_CMD $G0/regression/rt_10m_burch.c
$RM_CMD $G0/regression/rt_10m_burch_grad_closure.c
$RM_CMD $G0/regression/rt_10m_expanding.c
$RM_CMD $G0/regression/rt_10m_expanding_axi_sodshock.c
$RM_CMD $G0/regression/rt_10m_expanding_sodshock.c
$RM_CMD $G0/regression/rt_10m_gem.c
$RM_CMD $G0/regression/rt_10m_gem_grad_closure.c
$RM_CMD $G0/regression/rt_10m_lhdi.c
$RM_CMD $G0/regression/rt_10m_lhdi_grad_closure.c
$RM_CMD $G0/regression/rt_10m_par_firehose.c
$RM_CMD $G0/regression/rt_10m_par_firehose_grad_closure.c
$RM_CMD $G0/regression/rt_10m_riem.c
$RM_CMD $G0/regression/rt_10m_riem_grad_closure.c
$RM_CMD $G0/regression/rt_10m_sodshock.c
$RM_CMD $G0/regression/rt_10m_sodshock_lax.c
$RM_CMD $G0/regression/rt_5m_burch.c
$RM_CMD $G0/regression/rt_5m_elc_heat_flux.c
$RM_CMD $G0/regression/rt_5m_em_advect.c
$RM_CMD $G0/regression/rt_5m_em_advect_resonant.c
$RM_CMD $G0/regression/rt_5m_expanding.c
$RM_CMD $G0/regression/rt_5m_expanding_axi_sodshock.c
$RM_CMD $G0/regression/rt_5m_expanding_sodshock.c
$RM_CMD $G0/regression/rt_5m_friction.c
$RM_CMD $G0/regression/rt_5m_gem.c
$RM_CMD $G0/regression/rt_5m_hartmann.c
$RM_CMD $G0/regression/rt_5m_ion_heat_flux.c
$RM_CMD $G0/regression/rt_5m_mom_beach.c
$RM_CMD $G0/regression/rt_5m_riem.c
$RM_CMD $G0/regression/rt_5m_rt.c
$RM_CMD $G0/regression/rt_advect_wv.c
$RM_CMD $G0/regression/rt_advect_wv_mp.c
$RM_CMD $G0/regression/rt_burgers_shock.c
$RM_CMD $G0/regression/rt_burgers_shock_mp.c
$RM_CMD $G0/regression/rt_coldfluid_beach.c
$RM_CMD $G0/regression/rt_coldfluid_clouda.c
$RM_CMD $G0/regression/rt_coldfluid_em_coupling.c
$RM_CMD $G0/regression/rt_euler_axi_sodshock.c
$RM_CMD $G0/regression/rt_euler_axi_vac_riem.c
$RM_CMD $G0/regression/rt_euler_bump_in_channel.c
$RM_CMD $G0/regression/rt_euler_c2p_sodshock.c
$RM_CMD $G0/regression/rt_euler_cart_axi_sodshock.c
$RM_CMD $G0/regression/rt_euler_kh_2d.c
$RM_CMD $G0/regression/rt_euler_mixture_fedkiw_shock.c
$RM_CMD $G0/regression/rt_euler_mixture_shock_bubble.c
$RM_CMD $G0/regression/rt_euler_noh_1d.c
$RM_CMD $G0/regression/rt_euler_p_perturbation.c
$RM_CMD $G0/regression/rt_euler_reflect_2d.c
$RM_CMD $G0/regression/rt_euler_riem_2d_hll.c
$RM_CMD $G0/regression/rt_euler_riem_2d_hllc.c
$RM_CMD $G0/regression/rt_euler_riem_2d_lax.c
$RM_CMD $G0/regression/rt_euler_riem_2d_roe.c
$RM_CMD $G0/regression/rt_euler_riem_3d.c
$RM_CMD $G0/regression/rt_euler_rt.c
$RM_CMD $G0/regression/rt_euler_sodshock.c
$RM_CMD $G0/regression/rt_euler_sodshock_lax.c
$RM_CMD $G0/regression/rt_euler_sodshock_mp.c
$RM_CMD $G0/regression/rt_euler_superwedge.c
$RM_CMD $G0/regression/rt_euler_vac.c
$RM_CMD $G0/regression/rt_euler_vac_riem_1d.c
$RM_CMD $G0/regression/rt_euler_wave_2d_kep.c
$RM_CMD $G0/regression/rt_euler_wave_2d_mp.c
$RM_CMD $G0/regression/rt_euler_wave_2d_wv.c
$RM_CMD $G0/regression/rt_euler_wedge_sodshock.c
$RM_CMD $G0/regression/rt_gr_bhl_spinning.c
$RM_CMD $G0/regression/rt_gr_bhl_spinning_tetrad.c
$RM_CMD $G0/regression/rt_gr_bhl_static.c
$RM_CMD $G0/regression/rt_gr_bhl_static_neutronstar.c
$RM_CMD $G0/regression/rt_gr_bhl_static_tetrad.c
$RM_CMD $G0/regression/rt_gr_blackhole_spinning.c
$RM_CMD $G0/regression/rt_gr_blackhole_static.c
$RM_CMD $G0/regression/rt_gr_bz_monopole_fast.c
$RM_CMD $G0/regression/rt_gr_bz_monopole_fast_tetrad.c
$RM_CMD $G0/regression/rt_gr_bz_monopole_slow.c
$RM_CMD $G0/regression/rt_gr_bz_monopole_slow_tetrad.c
$RM_CMD $G0/regression/rt_gr_current_sheet.c
$RM_CMD $G0/regression/rt_gr_current_sheet_tetrad.c
$RM_CMD $G0/regression/rt_gr_einstein_plane_shock.c
$RM_CMD $G0/regression/rt_gr_kh_2d.c
$RM_CMD $G0/regression/rt_gr_mild_shock.c
$RM_CMD $G0/regression/rt_gr_mild_shock_tetrad.c
$RM_CMD $G0/regression/rt_gr_perturbed_density.c
$RM_CMD $G0/regression/rt_gr_quadrants_2d.c
$RM_CMD $G0/regression/rt_gr_strong_blast.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_bhl_spinning.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_bhl_spinning_tetrad.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_bhl_static.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_bhl_static_tetrad.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_shock.c
$RM_CMD $G0/regression/rt_gr_ultra_rel_shock_tetrad.c
$RM_CMD $G0/regression/rt_gr_wald_magnetosphere_spinning.c
$RM_CMD $G0/regression/rt_gr_wald_magnetosphere_spinning_tetrad.c
$RM_CMD $G0/regression/rt_gr_wald_magnetosphere_static.c
$RM_CMD $G0/regression/rt_gr_wald_magnetosphere_static_neutronstar.c
$RM_CMD $G0/regression/rt_gr_wald_magnetosphere_static_tetrad.c
$RM_CMD $G0/regression/rt_iso_euler_friction.c
$RM_CMD $G0/regression/rt_iso_euler_hartmann.c
$RM_CMD $G0/regression/rt_iso_euler_mixture_fedkiw_shock.c
$RM_CMD $G0/regression/rt_iso_euler_mixture_shock_bubble.c
$RM_CMD $G0/regression/rt_iso_euler_sodshock.c
$RM_CMD $G0/regression/rt_iso_euler_sodshock_lax.c
$RM_CMD $G0/regression/rt_iso_gem.c
$RM_CMD $G0/regression/rt_maxwell_annular_wg.c
$RM_CMD $G0/regression/rt_maxwell_annulus.c
$RM_CMD $G0/regression/rt_maxwell_axi_wg.c
$RM_CMD $G0/regression/rt_maxwell_expanding.c
$RM_CMD $G0/regression/rt_maxwell_expanding_2d.c
$RM_CMD $G0/regression/rt_maxwell_plane_wave_1d.c
$RM_CMD $G0/regression/rt_maxwell_plane_wave_1d_mp.c
$RM_CMD $G0/regression/rt_maxwell_plane_wave_2d.c
$RM_CMD $G0/regression/rt_maxwell_plane_wave_2d_mp.c
$RM_CMD $G0/regression/rt_maxwell_reflect_2d.c
$RM_CMD $G0/regression/rt_maxwell_wg_2d.c
$RM_CMD $G0/regression/rt_mhd_brio_wu.c
$RM_CMD $G0/regression/rt_mhd_ot.c
$RM_CMD $G0/regression/rt_mhd_rj2.c
$RM_CMD $G0/regression/rt_reactive_euler_detonation.c
$RM_CMD $G0/regression/rt_sr_euler_KH_2d.c
$RM_CMD $G0/regression/rt_sr_euler_riem_2d.c
$RM_CMD $G0/regression/rt_sr_euler_sodshock.c
$RM_CMD $G0/regression/rt_euler_multiblock.c
$RM_CMD $G0/regression/rt_multib_euler_2d.c

# Lua regression tests
mkdir -p luareg
$CP_CMD $G0/regression/lua/rt_10m_burch.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_burch_grad_closure.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_expanding.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_expanding_axi_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_expanding_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_gem.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_gem_grad_closure.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_lhdi.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_lhdi_grad_closure.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_par_firehose.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_par_firehose_grad_closure.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_riem.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_riem_grad_closure.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_10m_sodshock_lax.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_burch.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_em_advect.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_em_advect_resonant.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_expanding.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_expanding_axi_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_expanding_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_gem.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_hartmann.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_mom_beach.lua luareg/
$CP_CMD $G0/regression/lua/rt_5m_riem.lua luareg/
$CP_CMD $G0/regression/lua/rt_advect_wv.lua luareg/
$CP_CMD $G0/regression/lua/rt_advect_wv_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_burgers_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_burgers_shock_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_axi_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_axi_vac_riem.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_bump_in_channel.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_c2p_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_cart_axi_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_mixture_fedkiw_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_mixture_shock_bubble.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_noh_1d.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_p_perturbation.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_reflect_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_riem_2d_hll.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_riem_2d_hllc.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_riem_2d_lax.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_riem_2d_roe.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_riem_3d.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_rt.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_sodshock_lax.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_sodshock_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_superwedge.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_vac.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_vac_riem_1d.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_wave_2d_kep.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_wave_2d_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_wave_2d_wv.lua luareg/
$CP_CMD $G0/regression/lua/rt_euler_wedge_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bhl_spinning.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bhl_spinning_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bhl_static.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bhl_static_neutronstar.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bhl_static_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_blackhole_spinning.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_blackhole_static.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bz_monopole_fast.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bz_monopole_fast_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bz_monopole_slow.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_bz_monopole_slow_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_current_sheet.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_current_sheet_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_einstein_plane_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_kh_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_mild_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_mild_shock_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_perturbed_density.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_quadrants_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_strong_blast.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_spinning.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_spinning_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_static.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_static_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_ultra_rel_shock_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_spinning.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_spinning_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static_neutronstar.lua luareg/
$CP_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static_tetrad.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_euler_hartmann.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_euler_mixture_fedkiw_shock.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_euler_mixture_shock_bubble.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_euler_sodshock.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_euler_sodshock_lax.lua luareg/
$CP_CMD $G0/regression/lua/rt_iso_gem.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_annulus.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_expanding.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_expanding_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_plane_wave_1d.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_plane_wave_1d_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_plane_wave_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_plane_wave_2d_mp.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_reflect_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_maxwell_wg_2d.lua luareg/
$CP_CMD $G0/regression/lua/rt_mhd_brio_wu.lua luareg/
$CP_CMD $G0/regression/lua/rt_mhd_ot.lua luareg/
$CP_CMD $G0/regression/lua/rt_mhd_rj2.lua luareg/
$CP_CMD $G0/regression/lua/rt_reactive_euler_detonation.lua luareg/

$RM_CMD $G0/regression/lua/rt_10m_burch.lua
$RM_CMD $G0/regression/lua/rt_10m_burch_grad_closure.lua
$RM_CMD $G0/regression/lua/rt_10m_expanding.lua
$RM_CMD $G0/regression/lua/rt_10m_expanding_axi_sodshock.lua
$RM_CMD $G0/regression/lua/rt_10m_expanding_sodshock.lua
$RM_CMD $G0/regression/lua/rt_10m_gem.lua
$RM_CMD $G0/regression/lua/rt_10m_gem_grad_closure.lua
$RM_CMD $G0/regression/lua/rt_10m_lhdi.lua
$RM_CMD $G0/regression/lua/rt_10m_lhdi_grad_closure.lua
$RM_CMD $G0/regression/lua/rt_10m_par_firehose.lua
$RM_CMD $G0/regression/lua/rt_10m_par_firehose_grad_closure.lua
$RM_CMD $G0/regression/lua/rt_10m_riem.lua
$RM_CMD $G0/regression/lua/rt_10m_riem_grad_closure.lua
$RM_CMD $G0/regression/lua/rt_10m_sodshock.lua
$RM_CMD $G0/regression/lua/rt_10m_sodshock_lax.lua
$RM_CMD $G0/regression/lua/rt_5m_burch.lua
$RM_CMD $G0/regression/lua/rt_5m_em_advect.lua
$RM_CMD $G0/regression/lua/rt_5m_em_advect_resonant.lua
$RM_CMD $G0/regression/lua/rt_5m_expanding.lua
$RM_CMD $G0/regression/lua/rt_5m_expanding_axi_sodshock.lua
$RM_CMD $G0/regression/lua/rt_5m_expanding_sodshock.lua
$RM_CMD $G0/regression/lua/rt_5m_gem.lua
$RM_CMD $G0/regression/lua/rt_5m_hartmann.lua
$RM_CMD $G0/regression/lua/rt_5m_mom_beach.lua
$RM_CMD $G0/regression/lua/rt_5m_riem.lua
$RM_CMD $G0/regression/lua/rt_advect_wv.lua
$RM_CMD $G0/regression/lua/rt_advect_wv_mp.lua
$RM_CMD $G0/regression/lua/rt_burgers_shock.lua
$RM_CMD $G0/regression/lua/rt_burgers_shock_mp.lua
$RM_CMD $G0/regression/lua/rt_euler_axi_sodshock.lua
$RM_CMD $G0/regression/lua/rt_euler_axi_vac_riem.lua
$RM_CMD $G0/regression/lua/rt_euler_bump_in_channel.lua
$RM_CMD $G0/regression/lua/rt_euler_c2p_sodshock.lua
$RM_CMD $G0/regression/lua/rt_euler_cart_axi_sodshock.lua
$RM_CMD $G0/regression/lua/rt_euler_mixture_fedkiw_shock.lua
$RM_CMD $G0/regression/lua/rt_euler_mixture_shock_bubble.lua
$RM_CMD $G0/regression/lua/rt_euler_noh_1d.lua
$RM_CMD $G0/regression/lua/rt_euler_p_perturbation.lua
$RM_CMD $G0/regression/lua/rt_euler_reflect_2d.lua
$RM_CMD $G0/regression/lua/rt_euler_riem_2d_hll.lua
$RM_CMD $G0/regression/lua/rt_euler_riem_2d_hllc.lua
$RM_CMD $G0/regression/lua/rt_euler_riem_2d_lax.lua
$RM_CMD $G0/regression/lua/rt_euler_riem_2d_roe.lua
$RM_CMD $G0/regression/lua/rt_euler_riem_3d.lua
$RM_CMD $G0/regression/lua/rt_euler_rt.lua
$RM_CMD $G0/regression/lua/rt_euler_sodshock.lua
$RM_CMD $G0/regression/lua/rt_euler_sodshock_lax.lua
$RM_CMD $G0/regression/lua/rt_euler_sodshock_mp.lua
$RM_CMD $G0/regression/lua/rt_euler_superwedge.lua
$RM_CMD $G0/regression/lua/rt_euler_vac.lua
$RM_CMD $G0/regression/lua/rt_euler_vac_riem_1d.lua
$RM_CMD $G0/regression/lua/rt_euler_wave_2d_kep.lua
$RM_CMD $G0/regression/lua/rt_euler_wave_2d_mp.lua
$RM_CMD $G0/regression/lua/rt_euler_wave_2d_wv.lua
$RM_CMD $G0/regression/lua/rt_euler_wedge_sodshock.lua
$RM_CMD $G0/regression/lua/rt_gr_bhl_spinning.lua
$RM_CMD $G0/regression/lua/rt_gr_bhl_spinning_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_bhl_static.lua
$RM_CMD $G0/regression/lua/rt_gr_bhl_static_neutronstar.lua
$RM_CMD $G0/regression/lua/rt_gr_bhl_static_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_blackhole_spinning.lua
$RM_CMD $G0/regression/lua/rt_gr_blackhole_static.lua
$RM_CMD $G0/regression/lua/rt_gr_bz_monopole_fast.lua
$RM_CMD $G0/regression/lua/rt_gr_bz_monopole_fast_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_bz_monopole_slow.lua
$RM_CMD $G0/regression/lua/rt_gr_bz_monopole_slow_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_current_sheet.lua
$RM_CMD $G0/regression/lua/rt_gr_current_sheet_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_einstein_plane_shock.lua
$RM_CMD $G0/regression/lua/rt_gr_kh_2d.lua
$RM_CMD $G0/regression/lua/rt_gr_mild_shock.lua
$RM_CMD $G0/regression/lua/rt_gr_mild_shock_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_perturbed_density.lua
$RM_CMD $G0/regression/lua/rt_gr_quadrants_2d.lua
$RM_CMD $G0/regression/lua/rt_gr_strong_blast.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_spinning.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_spinning_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_static.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_bhl_static_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_shock.lua
$RM_CMD $G0/regression/lua/rt_gr_ultra_rel_shock_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_spinning.lua
$RM_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_spinning_tetrad.lua
$RM_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static.lua
$RM_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static_neutronstar.lua
$RM_CMD $G0/regression/lua/rt_gr_wald_magnetosphere_static_tetrad.lua
$RM_CMD $G0/regression/lua/rt_iso_euler_hartmann.lua
$RM_CMD $G0/regression/lua/rt_iso_euler_mixture_fedkiw_shock.lua
$RM_CMD $G0/regression/lua/rt_iso_euler_mixture_shock_bubble.lua
$RM_CMD $G0/regression/lua/rt_iso_euler_sodshock.lua
$RM_CMD $G0/regression/lua/rt_iso_euler_sodshock_lax.lua
$RM_CMD $G0/regression/lua/rt_iso_gem.lua
$RM_CMD $G0/regression/lua/rt_maxwell_annulus.lua
$RM_CMD $G0/regression/lua/rt_maxwell_expanding.lua
$RM_CMD $G0/regression/lua/rt_maxwell_expanding_2d.lua
$RM_CMD $G0/regression/lua/rt_maxwell_plane_wave_1d.lua
$RM_CMD $G0/regression/lua/rt_maxwell_plane_wave_1d_mp.lua
$RM_CMD $G0/regression/lua/rt_maxwell_plane_wave_2d.lua
$RM_CMD $G0/regression/lua/rt_maxwell_plane_wave_2d_mp.lua
$RM_CMD $G0/regression/lua/rt_maxwell_reflect_2d.lua
$RM_CMD $G0/regression/lua/rt_maxwell_wg_2d.lua
$RM_CMD $G0/regression/lua/rt_mhd_brio_wu.lua
$RM_CMD $G0/regression/lua/rt_mhd_ot.lua
$RM_CMD $G0/regression/lua/rt_mhd_rj2.lua
$RM_CMD $G0/regression/lua/rt_reactive_euler_detonation.lua

# AMR C regression tests
mkdir -p amr_creg
$CP_CMD $G0/amr_regression/rt_amr_5m_gem_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_5m_gem_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_5m_riem_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_5m_riem_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_10m_gem_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_10m_gem_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_10m_riem_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_cart_axi_sodshock_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_cart_axi_sodshock_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_mixture_fedkiw_shock_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_mixture_fedkiw_shock_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_mixture_shock_bubble_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_mixture_shock_bubble_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_riem_2d_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_riem_2d_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_shock_bubble_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_shock_bubble_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_sodshock_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_euler_sodshock_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_bhl_spinning_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_bhl_spinning_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_bhl_static_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_bhl_static_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_blackhole_spinning_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_blackhole_spinning_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_blackhole_static_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_blackhole_static_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_mild_shock_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_mild_shock_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_perturbed_density_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_perturbed_density_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_quadrants_2d_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_quadrants_2d_l2.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_strong_blast_l1.c amr_creg/
$CP_CMD $G0/amr_regression/rt_amr_gr_strong_blast_l2.c amr_creg/

$RM_CMD $G0/amr_regression/rt_amr_5m_gem_l1.c
$RM_CMD $G0/amr_regression/rt_amr_5m_gem_l2.c
$RM_CMD $G0/amr_regression/rt_amr_5m_riem_l1.c
$RM_CMD $G0/amr_regression/rt_amr_5m_riem_l2.c
$RM_CMD $G0/amr_regression/rt_amr_10m_gem_l1.c
$RM_CMD $G0/amr_regression/rt_amr_10m_gem_l2.c
$RM_CMD $G0/amr_regression/rt_amr_10m_riem_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_cart_axi_sodshock_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_cart_axi_sodshock_l2.c
$RM_CMD $G0/amr_regression/rt_amr_euler_mixture_fedkiw_shock_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_mixture_fedkiw_shock_l2.c
$RM_CMD $G0/amr_regression/rt_amr_euler_mixture_shock_bubble_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_mixture_shock_bubble_l2.c
$RM_CMD $G0/amr_regression/rt_amr_euler_riem_2d_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_riem_2d_l2.c
$RM_CMD $G0/amr_regression/rt_amr_euler_shock_bubble_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_shock_bubble_l2.c
$RM_CMD $G0/amr_regression/rt_amr_euler_sodshock_l1.c
$RM_CMD $G0/amr_regression/rt_amr_euler_sodshock_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_bhl_spinning_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_bhl_spinning_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_bhl_static_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_bhl_static_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_blackhole_spinning_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_blackhole_spinning_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_blackhole_static_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_blackhole_static_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_mild_shock_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_mild_shock_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_perturbed_density_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_perturbed_density_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_quadrants_2d_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_quadrants_2d_l2.c
$RM_CMD $G0/amr_regression/rt_amr_gr_strong_blast_l1.c
$RM_CMD $G0/amr_regression/rt_amr_gr_strong_blast_l2.c