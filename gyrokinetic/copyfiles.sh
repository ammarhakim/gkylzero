#!/bin/sh

CP_CMD=cp
RM_CMD=rm
G0=..

# data
mkdir -p data/adas
$CP_CMD $G0/data/adas/adas_to_numpy.py data/adas/
$CP_CMD $G0/data/adas/adf11.py data/adas/
$CP_CMD $G0/data/adas/download_adas.py data/adas/
$CP_CMD $G0/data/adas/radiation_fit_parameters.txt data/adas/
$CP_CMD $G0/data/adas/read_adas.c data/adas/
$CP_CMD $G0/data/adas/read_adas.h data/adas/
$CP_CMD $G0/data/adas/read_radiation.py data/adas/
$CP_CMD $G0/data/adas/README.md data/adas/

$RM_CMD $G0/data/adas/adas_to_numpy.py
$RM_CMD $G0/data/adas/adf11.py
$RM_CMD $G0/data/adas/download_adas.py
$RM_CMD $G0/data/adas/radiation_fit_parameters.txt
$RM_CMD $G0/data/adas/read_adas.c
$RM_CMD $G0/data/adas/read_adas.h
$RM_CMD $G0/data/adas/read_radiation.py
$RM_CMD $G0/data/adas/README.md

mkdir -p data/eqdsk
$CP_CMD $G0/data/eqdsk/asdex.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/cerfon.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/elliptical.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/mast.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/single_coil.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/solovev.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/step.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/straight_cylinder.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/tcv_upper_SN.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/tcv.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/wham_hires.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/wham.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_cerfon.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_elliptical.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_single_coil.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_solovev.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_straight_cylinder.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/ltx_miller.geqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/write_efit_ltx_miller.py data/eqdsk/
$CP_CMD $G0/data/eqdsk/LTX_103955_03.eqdsk data/eqdsk/
$CP_CMD $G0/data/eqdsk/README data/eqdsk/

$RM_CMD $G0/data/eqdsk/asdex.geqdsk
$RM_CMD $G0/data/eqdsk/cerfon.geqdsk
$RM_CMD $G0/data/eqdsk/elliptical.geqdsk
$RM_CMD $G0/data/eqdsk/mast.geqdsk
$RM_CMD $G0/data/eqdsk/single_coil.geqdsk
$RM_CMD $G0/data/eqdsk/solovev.geqdsk
$RM_CMD $G0/data/eqdsk/step.geqdsk
$RM_CMD $G0/data/eqdsk/straight_cylinder.geqdsk
$RM_CMD $G0/data/eqdsk/tcv_upper_SN.geqdsk
$RM_CMD $G0/data/eqdsk/tcv.geqdsk
$RM_CMD $G0/data/eqdsk/wham_hires.geqdsk
$RM_CMD $G0/data/eqdsk/wham.geqdsk
$RM_CMD $G0/data/eqdsk/write_efit_cerfon.py
$RM_CMD $G0/data/eqdsk/write_efit_elliptical.py
$RM_CMD $G0/data/eqdsk/write_efit_single_coil.py
$RM_CMD $G0/data/eqdsk/write_efit_solovev.py
$RM_CMD $G0/data/eqdsk/write_efit_straight_cylinder.py
$RM_CMD $G0/data/eqdsk/ltx_miller.geqdsk
$RM_CMD $G0/data/eqdsk/write_efit_ltx_miller.py
$RM_CMD $G0/data/eqdsk/LTX_103955_03.eqdsk
$RM_CMD $G0/data/eqdsk/README

# kernels
mkdir -p kernels/ambi_bolt_potential
$CP_CMD $G0/kernels/ambi_bolt_potential/*.h kernels/ambi_bolt_potential/
$CP_CMD $G0/kernels/ambi_bolt_potential/*.c kernels/ambi_bolt_potential/
mkdir -p kernels/bgk_gyrokinetic
$CP_CMD $G0/kernels/bgk_gyrokinetic/*.h kernels/bgk_gyrokinetic/
$CP_CMD $G0/kernels/bgk_gyrokinetic/*.c kernels/bgk_gyrokinetic/
mkdir -p kernels/deflate_geo
$CP_CMD $G0/kernels/deflate_geo/*.h kernels/deflate_geo/
$CP_CMD $G0/kernels/deflate_geo/*.c kernels/deflate_geo/
mkdir -p kernels/deflate_surf
$CP_CMD $G0/kernels/deflate_surf/*.h kernels/deflate_surf/
$CP_CMD $G0/kernels/deflate_surf/*.c kernels/deflate_surf/
mkdir -p kernels/derived_geo
$CP_CMD $G0/kernels/derived_geo/*.h kernels/derived_geo/
$CP_CMD $G0/kernels/derived_geo/*.c kernels/derived_geo/
mkdir -p kernels/dg_diffusion_gyrokinetic
$CP_CMD $G0/kernels/dg_diffusion_gyrokinetic/*.h kernels/dg_diffusion_gyrokinetic/
$CP_CMD $G0/kernels/dg_diffusion_gyrokinetic/*.c kernels/dg_diffusion_gyrokinetic/
mkdir -p kernels/fem_parproj
$CP_CMD $G0/kernels/fem_parproj/*.h kernels/fem_parproj/
$CP_CMD $G0/kernels/fem_parproj/*.c kernels/fem_parproj/
mkdir -p kernels/fem_poisson_perp
$CP_CMD $G0/kernels/fem_poisson_perp/*.h kernels/fem_poisson_perp/
$CP_CMD $G0/kernels/fem_poisson_perp/*.c kernels/fem_poisson_perp/
mkdir -p kernels/gyrokinetic
$CP_CMD $G0/kernels/gyrokinetic/*.h kernels/gyrokinetic/
$CP_CMD $G0/kernels/gyrokinetic/*.c kernels/gyrokinetic/
mkdir -p kernels/gyrokinetic_pol_density
$CP_CMD $G0/kernels/gyrokinetic_pol_density/*.h kernels/gyrokinetic_pol_density/
$CP_CMD $G0/kernels/gyrokinetic_pol_density/*.c kernels/gyrokinetic_pol_density/
mkdir -p kernels/inflate_surf
$CP_CMD $G0/kernels/inflate_surf/*.h kernels/inflate_surf/
$CP_CMD $G0/kernels/inflate_surf/*.c kernels/inflate_surf/
mkdir -p kernels/lbo_gyrokinetic
$CP_CMD $G0/kernels/lbo_gyrokinetic/*.h kernels/lbo_gyrokinetic/
$CP_CMD $G0/kernels/lbo_gyrokinetic/*.c kernels/lbo_gyrokinetic/
mkdir -p kernels/neutral
$CP_CMD $G0/kernels/neutral/*.h kernels/neutral/
$CP_CMD $G0/kernels/neutral/*.c kernels/neutral/
mkdir -p kernels/positivity_shift_gyrokinetic
$CP_CMD $G0/kernels/positivity_shift_gyrokinetic/*.h kernels/positivity_shift_gyrokinetic/
$CP_CMD $G0/kernels/positivity_shift_gyrokinetic/*.c kernels/positivity_shift_gyrokinetic/
mkdir -p kernels/rad
$CP_CMD $G0/kernels/rad/*.h kernels/rad/
$CP_CMD $G0/kernels/rad/*.c kernels/rad/
mkdir -p kernels/translate_dim
$CP_CMD $G0/kernels/translate_dim/*.h kernels/translate_dim/
$CP_CMD $G0/kernels/translate_dim/*.c kernels/translate_dim/
mkdir -p kernels/twistshift
$CP_CMD $G0/kernels/twistshift/*.h kernels/twistshift/
$CP_CMD $G0/kernels/twistshift/*.c kernels/twistshift/

$RM_CMD $G0/kernels/ambi_bolt_potential/*.h
$RM_CMD $G0/kernels/ambi_bolt_potential/*.c
$RM_CMD $G0/kernels/bgk_gyrokinetic/*.h
$RM_CMD $G0/kernels/bgk_gyrokinetic/*.c
$RM_CMD $G0/kernels/deflate_geo/*.h
$RM_CMD $G0/kernels/deflate_geo/*.c
$RM_CMD $G0/kernels/deflate_surf/*.h
$RM_CMD $G0/kernels/deflate_surf/*.c
$RM_CMD $G0/kernels/derived_geo/*.h
$RM_CMD $G0/kernels/derived_geo/*.c
$RM_CMD $G0/kernels/dg_diffusion_gyrokinetic/*.h
$RM_CMD $G0/kernels/dg_diffusion_gyrokinetic/*.c
$RM_CMD $G0/kernels/fem_parproj/*.h
$RM_CMD $G0/kernels/fem_parproj/*.c
$RM_CMD $G0/kernels/fem_poisson_perp/*.h
$RM_CMD $G0/kernels/fem_poisson_perp/*.c
$RM_CMD $G0/kernels/gyrokinetic/*.h
$RM_CMD $G0/kernels/gyrokinetic/*.c
$RM_CMD $G0/kernels/gyrokinetic_pol_density/*.h
$RM_CMD $G0/kernels/gyrokinetic_pol_density/*.c
$RM_CMD $G0/kernels/inflate_surf/*.h
$RM_CMD $G0/kernels/inflate_surf/*.c
$RM_CMD $G0/kernels/lbo_gyrokinetic/*.h
$RM_CMD $G0/kernels/lbo_gyrokinetic/*.c
$RM_CMD $G0/kernels/neutral/*.h
$RM_CMD $G0/kernels/neutral/*.c
$RM_CMD $G0/kernels/positivity_shift_gyrokinetic/*.h
$RM_CMD $G0/kernels/positivity_shift_gyrokinetic/*.c
$RM_CMD $G0/kernels/rad/*.h
$RM_CMD $G0/kernels/rad/*.c
$RM_CMD $G0/kernels/translate_dim/*.h
$RM_CMD $G0/kernels/translate_dim/*.c
$RM_CMD $G0/kernels/twistshift/*.h
$RM_CMD $G0/kernels/twistshift/*.c

# zero
mkdir -p zero
$CP_CMD $G0/zero/ambi_bolt_potential_cu.cu zero/
$CP_CMD $G0/zero/ambi_bolt_potential.c zero/
$CP_CMD $G0/zero/bc_block_tensor.c zero/
$CP_CMD $G0/zero/bc_sheath_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/bc_sheath_gyrokinetic.c zero/
$CP_CMD $G0/zero/bc_twistshift_cu.cu zero/
$CP_CMD $G0/zero/bc_twistshift.c zero/
$CP_CMD $G0/zero/boundary_flux_cu.cu zero/
$CP_CMD $G0/zero/boundary_flux.c zero/
$CP_CMD $G0/zero/calc_bmag.c zero/
$CP_CMD $G0/zero/calc_derived_geo.c zero/
$CP_CMD $G0/zero/calc_metric.c zero/
$CP_CMD $G0/zero/deflate_geo.c zero/
$CP_CMD $G0/zero/deflate_zsurf_cu.cu zero/
$CP_CMD $G0/zero/deflate_zsurf.c zero/
$CP_CMD $G0/zero/deflated_dg_bin_ops.c zero/
$CP_CMD $G0/zero/deflated_fem_poisson.c zero/
$CP_CMD $G0/zero/dg_calc_gk_rad_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_gk_rad_vars.c zero/
$CP_CMD $G0/zero/dg_calc_gyrokinetic_vars_cu.cu zero/
$CP_CMD $G0/zero/dg_calc_gyrokinetic_vars.c zero/
$CP_CMD $G0/zero/dg_cx_cu.cu zero/
$CP_CMD $G0/zero/dg_cx.c zero/
$CP_CMD $G0/zero/dg_diffusion_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/dg_diffusion_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/dg_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_iz_cu.cu zero/
$CP_CMD $G0/zero/dg_iz.c zero/
$CP_CMD $G0/zero/dg_lbo_gyrokinetic_diff_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_gyrokinetic_diff.c zero/
$CP_CMD $G0/zero/dg_lbo_gyrokinetic_drag_cu.cu zero/
$CP_CMD $G0/zero/dg_lbo_gyrokinetic_drag.c zero/
$CP_CMD $G0/zero/dg_rad_gyrokinetic_drag_cu.cu zero/
$CP_CMD $G0/zero/dg_rad_gyrokinetic_drag.c zero/
$CP_CMD $G0/zero/dg_recomb_cu.cu zero/
$CP_CMD $G0/zero/dg_recomb.c zero/
$CP_CMD $G0/zero/dg_updater_bflux_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_updater_diffusion_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_updater_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_updater_lbo_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_updater_moment_gyrokinetic.c zero/
$CP_CMD $G0/zero/dg_updater_rad_gyrokinetic.c zero/
$CP_CMD $G0/zero/efit_utils.c zero/
$CP_CMD $G0/zero/efit.c zero/
$CP_CMD $G0/zero/fem_parproj_cu.cu zero/
$CP_CMD $G0/zero/fem_parproj.c zero/
$CP_CMD $G0/zero/fem_poisson_perp_cu.cu zero/
$CP_CMD $G0/zero/fem_poisson_perp.c zero/
$CP_CMD $G0/zero/gk_geometry_cu.cu zero/
$CP_CMD $G0/zero/gk_geometry_mapc2p.c zero/
$CP_CMD $G0/zero/gk_geometry_mirror.c zero/
$CP_CMD $G0/zero/gk_geometry_tok.c zero/
$CP_CMD $G0/zero/gk_geometry.c zero/
$CP_CMD $G0/zero/gk_maxwellian_correct_cu.cu zero/
$CP_CMD $G0/zero/gk_maxwellian_correct.c zero/
$CP_CMD $G0/zero/gk_maxwellian_moments.c zero/
$CP_CMD $G0/zero/gk_maxwellian_proj_on_basis_cu.cu zero/
$CP_CMD $G0/zero/gk_maxwellian_proj_on_basis.c zero/
$CP_CMD $G0/zero/gkgeom.c zero/
$CP_CMD $G0/zero/gkyl_ambi_bolt_potential_priv.h zero/
$CP_CMD $G0/zero/gkyl_ambi_bolt_potential.h zero/
$CP_CMD $G0/zero/gkyl_bc_block_tensor_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_block_tensor.h zero/
$CP_CMD $G0/zero/gkyl_bc_sheath_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_sheath_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_bc_twistshift_priv.h zero/
$CP_CMD $G0/zero/gkyl_bc_twistshift.h zero/
$CP_CMD $G0/zero/gkyl_boundary_flux_priv.h zero/
$CP_CMD $G0/zero/gkyl_boundary_flux.h zero/
$CP_CMD $G0/zero/gkyl_calc_bmag_priv.h zero/
$CP_CMD $G0/zero/gkyl_calc_bmag.h zero/
$CP_CMD $G0/zero/gkyl_calc_derived_geo_priv.h zero/
$CP_CMD $G0/zero/gkyl_calc_derived_geo.h zero/
$CP_CMD $G0/zero/gkyl_calc_metric_priv.h zero/
$CP_CMD $G0/zero/gkyl_calc_metric.h zero/
$CP_CMD $G0/zero/gkyl_deflate_geo_priv.h zero/
$CP_CMD $G0/zero/gkyl_deflate_geo.h zero/
$CP_CMD $G0/zero/gkyl_deflate_zsurf_priv.h zero/
$CP_CMD $G0/zero/gkyl_deflate_zsurf.h zero/
$CP_CMD $G0/zero/gkyl_deflated_dg_bin_ops_priv.h zero/
$CP_CMD $G0/zero/gkyl_deflated_dg_bin_ops.h zero/
$CP_CMD $G0/zero/gkyl_deflated_fem_poisson_priv.h zero/
$CP_CMD $G0/zero/gkyl_deflated_fem_poisson.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gk_rad_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gk_rad_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gyrokinetic_vars_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gyrokinetic_vars.h zero/
$CP_CMD $G0/zero/gkyl_dg_cx_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_cx.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_diffusion_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_iz_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_iz.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_diff_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_diff.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_drag_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_drag.h zero/
$CP_CMD $G0/zero/gkyl_dg_rad_gyrokinetic_drag_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_rad_gyrokinetic_drag.h zero/
$CP_CMD $G0/zero/gkyl_dg_recomb_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_recomb.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_bflux_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_bflux_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_diffusion_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_lbo_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_moment_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_dg_updater_rad_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_efit_priv.h zero/
$CP_CMD $G0/zero/gkyl_efit.h zero/
$CP_CMD $G0/zero/gkyl_fem_parproj_priv.h zero/
$CP_CMD $G0/zero/gkyl_fem_parproj.h zero/
$CP_CMD $G0/zero/gkyl_fem_poisson_perp_priv.h zero/
$CP_CMD $G0/zero/gkyl_fem_poisson_perp.h zero/
$CP_CMD $G0/zero/gkyl_gk_geometry_mapc2p.h zero/
$CP_CMD $G0/zero/gkyl_gk_geometry_mirror.h zero/
$CP_CMD $G0/zero/gkyl_gk_geometry_tok.h zero/
$CP_CMD $G0/zero/gkyl_gk_geometry.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_correct_priv.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_correct.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_moments_priv.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_moments.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_proj_on_basis_priv.h zero/
$CP_CMD $G0/zero/gkyl_gk_maxwellian_proj_on_basis.h zero/
$CP_CMD $G0/zero/gkyl_gkgeom.h zero/
$CP_CMD $G0/zero/gkyl_gyrokinetic_cross_prim_moms_bgk_priv.h zero/
$CP_CMD $G0/zero/gkyl_gyrokinetic_cross_prim_moms_bgk.h zero/
$CP_CMD $G0/zero/gkyl_gyrokinetic_pol_density_priv.h zero/
$CP_CMD $G0/zero/gkyl_gyrokinetic_pol_density.h zero/
$CP_CMD $G0/zero/gkyl_mirror_geo_priv.h zero/
$CP_CMD $G0/zero/gkyl_mirror_geo.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_bcorr_lbo_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_mom_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_mom_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_position_map_priv.h zero/
$CP_CMD $G0/zero/gkyl_position_map.h zero/
$CP_CMD $G0/zero/gkyl_positivity_shift_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_positivity_shift_gyrokinetic.h zero/
$CP_CMD $G0/zero/gkyl_proj_powsqrt_on_basis_priv.h zero/
$CP_CMD $G0/zero/gkyl_proj_powsqrt_on_basis.h zero/
$CP_CMD $G0/zero/gkyl_radiation_read.h zero/
$CP_CMD $G0/zero/gkyl_rescale_ghost_jacf_priv.h zero/
$CP_CMD $G0/zero/gkyl_rescale_ghost_jacf.h zero/
$CP_CMD $G0/zero/gkyl_tok_calc_derived_geo_priv.h zero/
$CP_CMD $G0/zero/gkyl_tok_calc_derived_geo.h zero/
$CP_CMD $G0/zero/gkyl_tok_geo_priv.h zero/
$CP_CMD $G0/zero/gkyl_tok_geo.h zero/
$CP_CMD $G0/zero/gyrokinetic_cross_prim_moms_bgk_cu.cu zero/
$CP_CMD $G0/zero/gyrokinetic_cross_prim_moms_bgk.c zero/
$CP_CMD $G0/zero/gyrokinetic_pol_density_cu.cu zero/
$CP_CMD $G0/zero/gyrokinetic_pol_density.c zero/
$CP_CMD $G0/zero/mirror_geo_utils.c zero/
$CP_CMD $G0/zero/mirror_geo.c zero/
$CP_CMD $G0/zero/mom_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/mom_gyrokinetic.c zero/
$CP_CMD $G0/zero/position_map.c zero/
$CP_CMD $G0/zero/positivity_shift_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/positivity_shift_gyrokinetic.c zero/
$CP_CMD $G0/zero/proj_powsqrt_on_basis_cu.cu zero/
$CP_CMD $G0/zero/proj_powsqrt_on_basis.c zero/
$CP_CMD $G0/zero/radiation_read.c zero/
$CP_CMD $G0/zero/rescale_ghost_jacf_cu.cu zero/
$CP_CMD $G0/zero/rescale_ghost_jacf.c zero/
$CP_CMD $G0/zero/tok_calc_derived_geo.c zero/
$CP_CMD $G0/zero/tok_geo_utils.c zero/
$CP_CMD $G0/zero/tok_geo.c zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_gyrokinetic_priv.h zero/
$CP_CMD $G0/zero/gkyl_prim_lbo_gyrokinetic.h zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/mom_bcorr_lbo_gyrokinetic.c zero/
$CP_CMD $G0/zero/prim_lbo_gyrokinetic_cu.cu zero/
$CP_CMD $G0/zero/prim_lbo_gyrokinetic.c zero/
$CP_CMD $G0/zero/mom_calc_bcorr_gyrokinetic.c zero/
$CP_CMD $G0/zero/prim_lbo_calc_gyrokinetic.c zero/
$CP_CMD $G0/zero/prim_lbo_cross_calc_gyrokinetic.c zero/
$CP_CMD $G0/zero/gkyl_translate_dim_priv.h zero/
$CP_CMD $G0/zero/gkyl_translate_dim.h zero/
$CP_CMD $G0/zero/translate_dim_cu.cu zero/
$CP_CMD $G0/zero/translate_dim.c zero/
$CP_CMD $G0/zero/dg_calc_gk_neut_hamil.c zero/
$CP_CMD $G0/zero/dg_calc_gk_neut_hamil_cu.cu zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gk_neut_hamil.h zero/
$CP_CMD $G0/zero/gkyl_dg_calc_gk_neut_hamil_priv.h zero/

$RM_CMD $G0/zero/ambi_bolt_potential_cu.cu
$RM_CMD $G0/zero/ambi_bolt_potential.c
$RM_CMD $G0/zero/bc_block_tensor.c
$RM_CMD $G0/zero/bc_sheath_gyrokinetic_cu.cu
$RM_CMD $G0/zero/bc_sheath_gyrokinetic.c
$RM_CMD $G0/zero/bc_twistshift_cu.cu
$RM_CMD $G0/zero/bc_twistshift.c
$RM_CMD $G0/zero/boundary_flux_cu.cu
$RM_CMD $G0/zero/boundary_flux.c
$RM_CMD $G0/zero/calc_bmag.c
$RM_CMD $G0/zero/calc_derived_geo.c
$RM_CMD $G0/zero/calc_metric.c
$RM_CMD $G0/zero/deflate_geo.c
$RM_CMD $G0/zero/deflate_zsurf_cu.cu
$RM_CMD $G0/zero/deflate_zsurf.c
$RM_CMD $G0/zero/deflated_dg_bin_ops.c
$RM_CMD $G0/zero/deflated_fem_poisson.c
$RM_CMD $G0/zero/dg_calc_gk_rad_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_gk_rad_vars.c
$RM_CMD $G0/zero/dg_calc_gyrokinetic_vars_cu.cu
$RM_CMD $G0/zero/dg_calc_gyrokinetic_vars.c
$RM_CMD $G0/zero/dg_cx_cu.cu
$RM_CMD $G0/zero/dg_cx.c
$RM_CMD $G0/zero/dg_diffusion_gyrokinetic_cu.cu
$RM_CMD $G0/zero/dg_diffusion_gyrokinetic.c
$RM_CMD $G0/zero/dg_gyrokinetic_cu.cu
$RM_CMD $G0/zero/dg_gyrokinetic.c
$RM_CMD $G0/zero/dg_iz_cu.cu
$RM_CMD $G0/zero/dg_iz.c
$RM_CMD $G0/zero/dg_lbo_gyrokinetic_diff_cu.cu
$RM_CMD $G0/zero/dg_lbo_gyrokinetic_diff.c
$RM_CMD $G0/zero/dg_lbo_gyrokinetic_drag_cu.cu
$RM_CMD $G0/zero/dg_lbo_gyrokinetic_drag.c
$RM_CMD $G0/zero/dg_rad_gyrokinetic_drag_cu.cu
$RM_CMD $G0/zero/dg_rad_gyrokinetic_drag.c
$RM_CMD $G0/zero/dg_recomb_cu.cu
$RM_CMD $G0/zero/dg_recomb.c
$RM_CMD $G0/zero/dg_updater_bflux_gyrokinetic.c
$RM_CMD $G0/zero/dg_updater_diffusion_gyrokinetic.c
$RM_CMD $G0/zero/dg_updater_gyrokinetic.c
$RM_CMD $G0/zero/dg_updater_lbo_gyrokinetic.c
$RM_CMD $G0/zero/dg_updater_moment_gyrokinetic.c
$RM_CMD $G0/zero/dg_updater_rad_gyrokinetic.c
$RM_CMD $G0/zero/efit_utils.c
$RM_CMD $G0/zero/efit.c
$RM_CMD $G0/zero/fem_parproj_cu.cu
$RM_CMD $G0/zero/fem_parproj.c
$RM_CMD $G0/zero/fem_poisson_perp_cu.cu
$RM_CMD $G0/zero/fem_poisson_perp.c
$RM_CMD $G0/zero/gk_geometry_cu.cu
$RM_CMD $G0/zero/gk_geometry_mapc2p.c
$RM_CMD $G0/zero/gk_geometry_mirror.c
$RM_CMD $G0/zero/gk_geometry_tok.c
$RM_CMD $G0/zero/gk_geometry.c
$RM_CMD $G0/zero/gk_maxwellian_correct_cu.cu
$RM_CMD $G0/zero/gk_maxwellian_correct.c
$RM_CMD $G0/zero/gk_maxwellian_moments.c
$RM_CMD $G0/zero/gk_maxwellian_proj_on_basis_cu.cu
$RM_CMD $G0/zero/gk_maxwellian_proj_on_basis.c
$RM_CMD $G0/zero/gkgeom.c
$RM_CMD $G0/zero/gkyl_ambi_bolt_potential_priv.h
$RM_CMD $G0/zero/gkyl_ambi_bolt_potential.h
$RM_CMD $G0/zero/gkyl_bc_block_tensor_priv.h
$RM_CMD $G0/zero/gkyl_bc_block_tensor.h
$RM_CMD $G0/zero/gkyl_bc_sheath_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_bc_sheath_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_bc_twistshift_priv.h
$RM_CMD $G0/zero/gkyl_bc_twistshift.h
$RM_CMD $G0/zero/gkyl_boundary_flux_priv.h
$RM_CMD $G0/zero/gkyl_boundary_flux.h
$RM_CMD $G0/zero/gkyl_calc_bmag_priv.h
$RM_CMD $G0/zero/gkyl_calc_bmag.h
$RM_CMD $G0/zero/gkyl_calc_derived_geo_priv.h
$RM_CMD $G0/zero/gkyl_calc_derived_geo.h
$RM_CMD $G0/zero/gkyl_calc_metric_priv.h
$RM_CMD $G0/zero/gkyl_calc_metric.h
$RM_CMD $G0/zero/gkyl_deflate_geo_priv.h
$RM_CMD $G0/zero/gkyl_deflate_geo.h
$RM_CMD $G0/zero/gkyl_deflate_zsurf_priv.h
$RM_CMD $G0/zero/gkyl_deflate_zsurf.h
$RM_CMD $G0/zero/gkyl_deflated_dg_bin_ops_priv.h
$RM_CMD $G0/zero/gkyl_deflated_dg_bin_ops.h
$RM_CMD $G0/zero/gkyl_deflated_fem_poisson_priv.h
$RM_CMD $G0/zero/gkyl_deflated_fem_poisson.h
$RM_CMD $G0/zero/gkyl_dg_calc_gk_rad_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_gk_rad_vars.h
$RM_CMD $G0/zero/gkyl_dg_calc_gyrokinetic_vars_priv.h
$RM_CMD $G0/zero/gkyl_dg_calc_gyrokinetic_vars.h
$RM_CMD $G0/zero/gkyl_dg_cx_priv.h
$RM_CMD $G0/zero/gkyl_dg_cx.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_dg_diffusion_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_dg_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_iz_priv.h
$RM_CMD $G0/zero/gkyl_dg_iz.h
$RM_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_diff_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_diff.h
$RM_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_drag_priv.h
$RM_CMD $G0/zero/gkyl_dg_lbo_gyrokinetic_drag.h
$RM_CMD $G0/zero/gkyl_dg_rad_gyrokinetic_drag_priv.h
$RM_CMD $G0/zero/gkyl_dg_rad_gyrokinetic_drag.h
$RM_CMD $G0/zero/gkyl_dg_recomb_priv.h
$RM_CMD $G0/zero/gkyl_dg_recomb.h
$RM_CMD $G0/zero/gkyl_dg_updater_bflux_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_bflux_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_diffusion_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_updater_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_dg_updater_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_updater_lbo_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_updater_moment_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_dg_updater_rad_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_efit_priv.h
$RM_CMD $G0/zero/gkyl_efit.h
$RM_CMD $G0/zero/gkyl_fem_parproj_priv.h
$RM_CMD $G0/zero/gkyl_fem_parproj.h
$RM_CMD $G0/zero/gkyl_fem_poisson_perp_priv.h
$RM_CMD $G0/zero/gkyl_fem_poisson_perp.h
$RM_CMD $G0/zero/gkyl_gk_geometry_mapc2p.h
$RM_CMD $G0/zero/gkyl_gk_geometry_mirror.h
$RM_CMD $G0/zero/gkyl_gk_geometry_tok.h
$RM_CMD $G0/zero/gkyl_gk_geometry.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_correct_priv.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_correct.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_moments_priv.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_moments.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_proj_on_basis_priv.h
$RM_CMD $G0/zero/gkyl_gk_maxwellian_proj_on_basis.h
$RM_CMD $G0/zero/gkyl_gkgeom.h
$RM_CMD $G0/zero/gkyl_gyrokinetic_cross_prim_moms_bgk_priv.h
$RM_CMD $G0/zero/gkyl_gyrokinetic_cross_prim_moms_bgk.h
$RM_CMD $G0/zero/gkyl_gyrokinetic_pol_density_priv.h
$RM_CMD $G0/zero/gkyl_gyrokinetic_pol_density.h
$RM_CMD $G0/zero/gkyl_mirror_geo_priv.h
$RM_CMD $G0/zero/gkyl_mirror_geo.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_mom_bcorr_lbo_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_mom_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_mom_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_position_map_priv.h
$RM_CMD $G0/zero/gkyl_position_map.h
$RM_CMD $G0/zero/gkyl_positivity_shift_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_positivity_shift_gyrokinetic.h
$RM_CMD $G0/zero/gkyl_proj_powsqrt_on_basis_priv.h
$RM_CMD $G0/zero/gkyl_proj_powsqrt_on_basis.h
$RM_CMD $G0/zero/gkyl_radiation_read.h
$RM_CMD $G0/zero/gkyl_rescale_ghost_jacf_priv.h
$RM_CMD $G0/zero/gkyl_rescale_ghost_jacf.h
$RM_CMD $G0/zero/gkyl_tok_calc_derived_geo_priv.h
$RM_CMD $G0/zero/gkyl_tok_calc_derived_geo.h
$RM_CMD $G0/zero/gkyl_tok_geo_priv.h
$RM_CMD $G0/zero/gkyl_tok_geo.h
$RM_CMD $G0/zero/gyrokinetic_cross_prim_moms_bgk_cu.cu
$RM_CMD $G0/zero/gyrokinetic_cross_prim_moms_bgk.c
$RM_CMD $G0/zero/gyrokinetic_pol_density_cu.cu
$RM_CMD $G0/zero/gyrokinetic_pol_density.c
$RM_CMD $G0/zero/mirror_geo_utils.c
$RM_CMD $G0/zero/mirror_geo.c
$RM_CMD $G0/zero/mom_gyrokinetic_cu.cu
$RM_CMD $G0/zero/mom_gyrokinetic.c
$RM_CMD $G0/zero/position_map.c
$RM_CMD $G0/zero/positivity_shift_gyrokinetic_cu.cu
$RM_CMD $G0/zero/positivity_shift_gyrokinetic.c
$RM_CMD $G0/zero/proj_powsqrt_on_basis_cu.cu
$RM_CMD $G0/zero/proj_powsqrt_on_basis.c
$RM_CMD $G0/zero/radiation_read.c
$RM_CMD $G0/zero/rescale_ghost_jacf_cu.cu
$RM_CMD $G0/zero/rescale_ghost_jacf.c
$RM_CMD $G0/zero/tok_calc_derived_geo.c
$RM_CMD $G0/zero/tok_geo_utils.c
$RM_CMD $G0/zero/tok_geo.c
$RM_CMD $G0/zero/gkyl_prim_lbo_gyrokinetic_priv.h
$RM_CMD $G0/zero/gkyl_prim_lbo_gyrokinetic.h
$RM_CMD $G0/zero/mom_bcorr_lbo_gyrokinetic_cu.cu
$RM_CMD $G0/zero/mom_bcorr_lbo_gyrokinetic.c
$RM_CMD $G0/zero/prim_lbo_gyrokinetic_cu.cu
$RM_CMD $G0/zero/prim_lbo_gyrokinetic.c
$RM_CMD $G0/zero/mom_calc_bcorr_gyrokinetic.c
$RM_CMD $G0/zero/prim_lbo_calc_gyrokinetic.c
$RM_CMD $G0/zero/prim_lbo_cross_calc_gyrokinetic.c
$RM_CMD $G0/zero/gkyl_translate_dim_priv.h
$RM_CMD $G0/zero/gkyl_translate_dim.h
$RM_CMD $G0/zero/translate_dim_cu.cu
$RM_CMD $G0/zero/translate_dim.c
$RM_CMD $G0/zero/dg_calc_gk_neut_hamil.c
$RM_CMD $G0/zero/dg_calc_gk_neut_hamil_cu.cu
$RM_CMD $G0/zero/gkyl_dg_calc_gk_neut_hamil.h
$RM_CMD $G0/zero/gkyl_dg_calc_gk_neut_hamil_priv.h

# app
mkdir -p apps
$CP_CMD $G0/apps/block_gk_geom.c apps/
$CP_CMD $G0/apps/gk_field.c apps/
$CP_CMD $G0/apps/gk_multib_field.c apps/
$CP_CMD $G0/apps/gk_neut_species_bgk.c apps/
$CP_CMD $G0/apps/gk_neut_species_lte.c apps/
$CP_CMD $G0/apps/gk_neut_species_moment.c apps/
$CP_CMD $G0/apps/gk_neut_species_projection.c apps/
$CP_CMD $G0/apps/gk_neut_species_react.c apps/
$CP_CMD $G0/apps/gk_neut_species_source.c apps/
$CP_CMD $G0/apps/gk_neut_species.c apps/
$CP_CMD $G0/apps/gk_neut_species_bflux.c apps/
$CP_CMD $G0/apps/gk_neut_species_recycle.c apps/
$CP_CMD $G0/apps/gk_species_bflux.c apps/
$CP_CMD $G0/apps/gk_species_bgk.c apps/
$CP_CMD $G0/apps/gk_species_lbo.c apps/
$CP_CMD $G0/apps/gk_species_lte.c apps/
$CP_CMD $G0/apps/gk_species_moment.c apps/
$CP_CMD $G0/apps/gk_species_projection.c apps/
$CP_CMD $G0/apps/gk_species_radiation.c apps/
$CP_CMD $G0/apps/gk_species_react.c apps/
$CP_CMD $G0/apps/gk_species_source.c apps/
$CP_CMD $G0/apps/gk_species.c apps/
$CP_CMD $G0/apps/gkyl_gk_block_geom.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic_comms.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic_lw.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic_multib_priv.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic_multib.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic_priv.h apps/
$CP_CMD $G0/apps/gkyl_gyrokinetic.h apps/
$CP_CMD $G0/apps/gkyl_multib_conn.h apps/
$CP_CMD $G0/apps/gyrokinetic_comms.c apps/
$CP_CMD $G0/apps/gyrokinetic_lw.c apps/
$CP_CMD $G0/apps/gyrokinetic_multib_update_ssp_rk3.c apps/
$CP_CMD $G0/apps/gyrokinetic_multib.c apps/
$CP_CMD $G0/apps/gyrokinetic_update_implicit_coll.c apps/
$CP_CMD $G0/apps/gyrokinetic_update_op_split.c apps/
$CP_CMD $G0/apps/gyrokinetic_update_ssp_rk3.c apps/
$CP_CMD $G0/apps/gyrokinetic.c apps/
$CP_CMD $G0/apps/multib_conn.c apps/

$RM_CMD $G0/apps/block_gk_geom.c
$RM_CMD $G0/apps/gk_field.c
$RM_CMD $G0/apps/gk_multib_field.c
$RM_CMD $G0/apps/gk_neut_species_bgk.c
$RM_CMD $G0/apps/gk_neut_species_lte.c
$RM_CMD $G0/apps/gk_neut_species_moment.c
$RM_CMD $G0/apps/gk_neut_species_projection.c
$RM_CMD $G0/apps/gk_neut_species_react.c
$RM_CMD $G0/apps/gk_neut_species_source.c
$RM_CMD $G0/apps/gk_neut_species.c
$RM_CMD $G0/apps/gk_neut_species_bflux.c
$RM_CMD $G0/apps/gk_neut_species_recycle.c
$RM_CMD $G0/apps/gk_species_bflux.c
$RM_CMD $G0/apps/gk_species_bgk.c
$RM_CMD $G0/apps/gk_species_lbo.c
$RM_CMD $G0/apps/gk_species_lte.c
$RM_CMD $G0/apps/gk_species_moment.c
$RM_CMD $G0/apps/gk_species_projection.c
$RM_CMD $G0/apps/gk_species_radiation.c
$RM_CMD $G0/apps/gk_species_react.c
$RM_CMD $G0/apps/gk_species_source.c
$RM_CMD $G0/apps/gk_species.c
$RM_CMD $G0/apps/gkyl_gk_block_geom.h
$RM_CMD $G0/apps/gkyl_gyrokinetic_comms.h
$RM_CMD $G0/apps/gkyl_gyrokinetic_lw.h
$RM_CMD $G0/apps/gkyl_gyrokinetic_multib_priv.h
$RM_CMD $G0/apps/gkyl_gyrokinetic_multib.h
$RM_CMD $G0/apps/gkyl_gyrokinetic_priv.h
$RM_CMD $G0/apps/gkyl_gyrokinetic.h
$RM_CMD $G0/apps/gkyl_multib_conn.h
$RM_CMD $G0/apps/gyrokinetic_comms.c
$RM_CMD $G0/apps/gyrokinetic_lw.c
$RM_CMD $G0/apps/gyrokinetic_multib_update_ssp_rk3.c
$RM_CMD $G0/apps/gyrokinetic_multib.c
$RM_CMD $G0/apps/gyrokinetic_update_implicit_coll.c
$RM_CMD $G0/apps/gyrokinetic_update_op_split.c
$RM_CMD $G0/apps/gyrokinetic_update_ssp_rk3.c
$RM_CMD $G0/apps/gyrokinetic.c
$RM_CMD $G0/apps/multib_conn.c

# unit
mkdir -p unit
$CP_CMD $G0/unit/ctest_ambi_bolt_potential.c unit/
$CP_CMD $G0/unit/ctest_asdex.c unit/
$CP_CMD $G0/unit/ctest_bc_sheath_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_bc_twistshift.c unit/
$CP_CMD $G0/unit/ctest_block_tensor.c unit/
$CP_CMD $G0/unit/ctest_cerfon.c unit/
$CP_CMD $G0/unit/ctest_coll_cx.c unit/
$CP_CMD $G0/unit/ctest_coll_iz.c unit/
$CP_CMD $G0/unit/ctest_coll_recomb.c unit/
$CP_CMD $G0/unit/ctest_correct_maxwellian_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_deflate_zsurf.c unit/
$CP_CMD $G0/unit/ctest_deflated_dg_bin_ops.c unit/
$CP_CMD $G0/unit/ctest_deflated_fem_poisson.c unit/
$CP_CMD $G0/unit/ctest_dg_gyrokinetic_kern_tm.c unit/
$CP_CMD $G0/unit/ctest_dg_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_dg_interpolate.c unit/
$CP_CMD $G0/unit/ctest_dg_rad_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_efit.c unit/
$CP_CMD $G0/unit/ctest_gkneut_hamil.c unit/
$CP_CMD $G0/unit/ctest_fem_parproj.c unit/
$CP_CMD $G0/unit/ctest_fem_poisson_perp.c unit/
$CP_CMD $G0/unit/ctest_gk_geometry_mapc2p.c unit/
$CP_CMD $G0/unit/ctest_gk_geometry_mirror.c unit/
$CP_CMD $G0/unit/ctest_gk_geometry_tok.c unit/
$CP_CMD $G0/unit/ctest_gkgeom.c unit/
$CP_CMD $G0/unit/ctest_gyrokinetic_cross_prim_moms_bgk.c unit/
$CP_CMD $G0/unit/ctest_gyrokinetic_pol_density.c unit/
$CP_CMD $G0/unit/ctest_integrated_moms.c unit/
$CP_CMD $G0/unit/ctest_mom_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_nodal_ops.c unit/
$CP_CMD $G0/unit/ctest_position_map.c unit/
$CP_CMD $G0/unit/ctest_positivity_shift_gyrokinetic.c unit/
$CP_CMD $G0/unit/ctest_proj_gk_bimaxwellian_on_basis.c unit/
$CP_CMD $G0/unit/ctest_proj_gk_maxwellian_on_basis.c unit/
$CP_CMD $G0/unit/ctest_proj_powsqrt_on_basis.c unit/
$CP_CMD $G0/unit/ctest_rescale_ghost_jacf.c unit/
$CP_CMD $G0/unit/ctest_step_compare.c unit/
$CP_CMD $G0/unit/ctest_step_outboard.c unit/
$CP_CMD $G0/unit/ctest_time_roots.c unit/
$CP_CMD $G0/unit/mctest_multib_sync.c unit/
$CP_CMD $G0/unit/mctest_multib_allgather.c unit/
$CP_CMD $G0/unit/ctest_translate_dim.c unit/
$CP_CMD $G0/unit/ctest_ltx_miller.c unit/

$RM_CMD $G0/unit/ctest_ambi_bolt_potential.c
$RM_CMD $G0/unit/ctest_asdex.c
$RM_CMD $G0/unit/ctest_bc_sheath_gyrokinetic.c
$RM_CMD $G0/unit/ctest_bc_twistshift.c
$RM_CMD $G0/unit/ctest_block_tensor.c
$RM_CMD $G0/unit/ctest_cerfon.c
$RM_CMD $G0/unit/ctest_coll_cx.c
$RM_CMD $G0/unit/ctest_coll_iz.c
$RM_CMD $G0/unit/ctest_coll_recomb.c
$RM_CMD $G0/unit/ctest_correct_maxwellian_gyrokinetic.c
$RM_CMD $G0/unit/ctest_deflate_zsurf.c
$RM_CMD $G0/unit/ctest_deflated_dg_bin_ops.c
$RM_CMD $G0/unit/ctest_deflated_fem_poisson.c
$RM_CMD $G0/unit/ctest_dg_gyrokinetic_kern_tm.c
$RM_CMD $G0/unit/ctest_dg_gyrokinetic.c
$RM_CMD $G0/unit/ctest_dg_interpolate.c
$RM_CMD $G0/unit/ctest_dg_rad_gyrokinetic.c
$RM_CMD $G0/unit/ctest_efit.c
$RM_CMD $G0/unit/ctest_gkneut_hamil.c
$RM_CMD $G0/unit/ctest_fem_parproj.c
$RM_CMD $G0/unit/ctest_fem_poisson_perp.c
$RM_CMD $G0/unit/ctest_gk_geometry_mapc2p.c
$RM_CMD $G0/unit/ctest_gk_geometry_mirror.c
$RM_CMD $G0/unit/ctest_gk_geometry_tok.c
$RM_CMD $G0/unit/ctest_gkgeom.c
$RM_CMD $G0/unit/ctest_gyrokinetic_cross_prim_moms_bgk.c
$RM_CMD $G0/unit/ctest_gyrokinetic_pol_density.c
$RM_CMD $G0/unit/ctest_integrated_moms.c
$RM_CMD $G0/unit/ctest_mom_gyrokinetic.c
$RM_CMD $G0/unit/ctest_nodal_ops.c
$RM_CMD $G0/unit/ctest_position_map.c
$RM_CMD $G0/unit/ctest_positivity_shift_gyrokinetic.c
$RM_CMD $G0/unit/ctest_proj_gk_bimaxwellian_on_basis.c
$RM_CMD $G0/unit/ctest_proj_gk_maxwellian_on_basis.c
$RM_CMD $G0/unit/ctest_proj_powsqrt_on_basis.c
$RM_CMD $G0/unit/ctest_rescale_ghost_jacf.c
$RM_CMD $G0/unit/ctest_step_compare.c
$RM_CMD $G0/unit/ctest_step_outboard.c
$RM_CMD $G0/unit/ctest_time_roots.c
$RM_CMD $G0/unit/mctest_multib_sync.c
$RM_CMD $G0/unit/mctest_multib_allgather.c
$RM_CMD $G0/unit/ctest_translate_dim.c
$RM_CMD $G0/unit/ctest_ltx_miller.c

# C regression tests
mkdir -p creg
$CP_CMD $G0/regression/rt_gk_ar_react_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_asdex_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_cross_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_im_asdex_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_im_cross_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_im_periodic_sod_shock_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_periodic_sodshock_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_relax_bimaxwellian_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_bgk_relax_bimaxwellian_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_d3d_iwl_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_d3d_iwl_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_ion_sound_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_ion_sound_adiabatic_elc_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_ion_sound_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lapd_cart_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lapd_cyl_3x2v_p1_nonuniformr.c creg/
$CP_CMD $G0/regression/rt_gk_lapd_cyl_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_cross_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_relax_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_relax_bimaxwellian_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_relax_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_lbo_relax_varnu_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_li_react_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_li_react_nonuniformv_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_ltx_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_ltx_boltz_elc_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_mdpx_cart_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_mirror_boltz_elc_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_mirror_kinetic_elc_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_multib_sheath_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_multib_slab_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_multib_step_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_multib_step_sol_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_multib_step_sol_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_nozzle_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_nozzle_half_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_rad_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_rad_low_Te_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_rad_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_1x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_3x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_bgk_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_cx_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_flr_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_neut_sheath_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformv_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformv_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformv_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformx_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformx_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_sheath_nonuniformx_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_slab_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_solovev_out_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_static_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_step_1x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_step_3x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_step_nonuniformx_1x2v_p1_cons_numeric.c creg/
$CP_CMD $G0/regression/rt_gk_step_nonuniformx_1x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_step_nonuniformx_2x2v_p1_out.c creg/
$CP_CMD $G0/regression/rt_gk_step_nonuniformx_3x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_gk_step_out_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_1x2v_p1_static_field.c creg/
$CP_CMD $G0/regression/rt_gk_wham_1x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_2xIC_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_3x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_boltz_elc_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_1x2v_p1_numeric.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_1x2v_p1_polynomial.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_2x2v_p1_numeric.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_2x2v_p1_polynomial.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_3x2v_p1_numeric.c creg/
$CP_CMD $G0/regression/rt_gk_wham_nonuniformx_3x2v_p1_polynomial.c creg/
$CP_CMD $G0/regression/rt_gkgeom.c creg/
$CP_CMD $G0/regression/rt_gk_ltx_iwl_2x2v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_neut_recycle_1x3v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_neut_step_2x3v_p1.c creg/
$CP_CMD $G0/regression/rt_gk_step_2x2v_p1_cons.c creg/
$CP_CMD $G0/regression/rt_arg_parse.h creg/

$RM_CMD $G0/regression/rt_gk_ar_react_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_asdex_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_cross_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_im_asdex_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_im_cross_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_im_periodic_sod_shock_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_periodic_sodshock_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_relax_bimaxwellian_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_bgk_relax_bimaxwellian_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_d3d_iwl_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_d3d_iwl_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_ion_sound_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_ion_sound_adiabatic_elc_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_ion_sound_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lapd_cart_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lapd_cyl_3x2v_p1_nonuniformr.c
$RM_CMD $G0/regression/rt_gk_lapd_cyl_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_cross_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_relax_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_relax_bimaxwellian_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_relax_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_lbo_relax_varnu_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_li_react_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_li_react_nonuniformv_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_ltx_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_ltx_boltz_elc_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_mdpx_cart_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_mirror_boltz_elc_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_mirror_kinetic_elc_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_multib_sheath_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_multib_slab_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_multib_step_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_multib_step_sol_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_multib_step_sol_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_nozzle_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_nozzle_half_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_rad_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_rad_low_Te_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_rad_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_1x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_sheath_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_3x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_sheath_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_bgk_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_cx_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_flr_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_neut_sheath_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformv_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformv_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformv_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformx_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformx_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_sheath_nonuniformx_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_slab_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_solovev_out_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_static_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_step_1x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_step_3x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_step_nonuniformx_1x2v_p1_cons_numeric.c
$RM_CMD $G0/regression/rt_gk_step_nonuniformx_1x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_step_nonuniformx_2x2v_p1_out.c
$RM_CMD $G0/regression/rt_gk_step_nonuniformx_3x2v_p1_cons.c
$RM_CMD $G0/regression/rt_gk_step_out_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_1x2v_p1_static_field.c
$RM_CMD $G0/regression/rt_gk_wham_1x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_2xIC_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_3x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_boltz_elc_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_1x2v_p1_numeric.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_1x2v_p1_polynomial.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_2x2v_p1_numeric.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_2x2v_p1_polynomial.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_3x2v_p1_numeric.c
$RM_CMD $G0/regression/rt_gk_wham_nonuniformx_3x2v_p1_polynomial.c
$RM_CMD $G0/regression/rt_gkgeom.c
$RM_CMD $G0/regression/rt_gk_ltx_iwl_2x2v_p1.c
$RM_CMD $G0/regression/rt_gk_neut_recycle_1x3v_p1.c
$RM_CMD $G0/regression/rt_gk_neut_step_2x3v_p1.c
$RM_CMD $G0/regression/rt_gk_step_2x2v_p1_cons.c

# Lua regression tests
mkdir -p luareg
$CP_CMD $G0/regression/lua/rt_gk_ar_react_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_bgk_cross_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_bgk_periodic_sodshock_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_bgk_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_bgk_relax_bimaxwellian_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_bgk_relax_bimaxwellian_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_ion_sound_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_ion_sound_adiabatic_elc_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_ion_sound_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_cross_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_relax_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_relax_bimaxwellian_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_relax_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_lbo_relax_varnu_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_li_react_3x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_li_react_nonuniformv_3x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_rad_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_rad_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_1x2v_p1_cons.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_3x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_bgk_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_3x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_1x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_2x2v_p1.lua luareg/
$CP_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_3x2v_p1.lua luareg/

$RM_CMD $G0/regression/lua/rt_gk_ar_react_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_bgk_cross_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_bgk_periodic_sodshock_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_bgk_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_bgk_relax_bimaxwellian_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_bgk_relax_bimaxwellian_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_ion_sound_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_ion_sound_adiabatic_elc_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_ion_sound_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_cross_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_relax_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_relax_bimaxwellian_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_relax_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_lbo_relax_varnu_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_li_react_3x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_li_react_nonuniformv_3x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_rad_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_rad_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_1x2v_p1_cons.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_3x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_bgk_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformv_3x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_1x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_2x2v_p1.lua
$RM_CMD $G0/regression/lua/rt_gk_sheath_nonuniformx_3x2v_p1.lua