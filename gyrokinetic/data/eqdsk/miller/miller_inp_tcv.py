"""Exact surfaces of rt_gk_tcv_iwl_adapt_source_3x2v_p1.c, including shift."""
import math

geometry = dict(R_axis=0.8867856264, Z_axis=0.1414361745, B_axis=1.4,
                R_lcfs=1.0870056099999, kappa=1.4, delta=-0.38, a_shift=0.5)
offset = geometry["R_lcfs"]-geometry["R_axis"]
a = 2*offset/(1+math.sqrt(1-2*geometry["a_shift"]*offset/geometry["R_axis"]))


def R_omp(r):
  # Map minor radius r [m] to outboard midplane major radius [m],
  # including the quadratic Shafranov shift of the flux-surface center.
  return geometry["R_axis"]-geometry["a_shift"]*r*r/(2*geometry["R_axis"])+r


R_omp_bounds = (R_omp(a-0.04), R_omp(a+0.08))


def q(R_omp):
  # Safety factor as a function of outboard major radius [m].
  return ((497.3420166252413*R_omp-1408.736172826569)*R_omp
          + 1331.4134861681464)*R_omp-419.00692601227627
