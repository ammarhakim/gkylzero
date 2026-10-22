"""Parameters and exact physical endpoint surfaces of rt_gk_d3d_iwl_3x2v_p1.c."""
geometry = dict(R_axis=1.6, Z_axis=0.013055028, B_axis=2.0*1.6486461/1.6,
                R_lcfs=2.17, kappa=1.35, delta=0.4, a_shift=0.0)
R_omp_bounds = (2.17-5*0.15/8, 2.17+3*0.15/8)


def q(R_omp):
  # Safety factor as a function of outboard major radius [m].
  return ((49.46395467479657*R_omp-260.79513158768754)*R_omp
          + 458.42618139184754)*R_omp-267.63441353752336
