#pragma once

struct gkyl_mom_updater_lbo_vlasov {
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator

  struct gkyl_prim_lbo_type *coll_prim; // LBO primitive moments type
  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
};
