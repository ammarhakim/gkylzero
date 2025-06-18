#include <gkyl_array.h>
#include <gkyl_mirror_geo_dg.h>
#include <gkyl_gk_geometry.h>

void
gk_geometry_mirror_array_acquire(struct gk_geometry *up, struct gkyl_mirror_geo_dg *geo_dg)
{
    up->geo_corn.mc2p = gkyl_array_acquire(geo_dg->mapc2p);
    up->geo_corn.mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_corn.mc2nu_pos = gkyl_array_acquire(geo_dg->mc2nu_pos);
    up->geo_corn.mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.dxdz = gkyl_array_acquire(geo_dg->tang);
    up->geo_int.dzdx = gkyl_array_acquire(geo_dg->dual);
    up->geo_int.dualmag = gkyl_array_acquire(geo_dg->dualmag);
    up->geo_int.normals = gkyl_array_acquire(geo_dg->normals);
    up->geo_int.g_ij = gkyl_array_acquire(geo_dg->metric_covar);
    up->geo_int.g_ij_neut = gkyl_array_acquire(geo_dg->metric_covar_neut);
    up->geo_int.gij = gkyl_array_acquire(geo_dg->metric_contr);
    up->geo_int.gij_neut = gkyl_array_acquire(geo_dg->metric_contr_neut);
    up->geo_int.gxxj = gkyl_array_acquire(geo_dg->gxxj);
    up->geo_int.gxyj = gkyl_array_acquire(geo_dg->gxyj);
    up->geo_int.gyyj = gkyl_array_acquire(geo_dg->gyyj);
    up->geo_int.gxzj = gkyl_array_acquire(geo_dg->gxzj);
    up->geo_int.jacobgeo = gkyl_array_acquire(geo_dg->Jc);
    up->geo_int.jacobgeo_inv = gkyl_array_acquire(geo_dg->Jc_inv);
    up->geo_int.jacobtot = gkyl_array_acquire(geo_dg->JB);
    up->geo_int.jacobtot_inv = gkyl_array_acquire(geo_dg->JB_inv);
    up->geo_int.b_i = gkyl_array_acquire(geo_dg->b_covar);
    up->geo_int.bcart = gkyl_array_acquire(geo_dg->b_cart);
    up->geo_corn.bmag = gkyl_array_acquire(geo_dg->Bmag);
    up->geo_int.bmag_inv = gkyl_array_acquire(geo_dg->Bmag_inv);
    up->geo_int.bmag_inv_sq = gkyl_array_acquire(geo_dg->Bmag_inv_sq);
    up->geo_int.cmag = gkyl_array_acquire(geo_dg->C);
    up->geo_int.eps2 = gkyl_array_acquire(geo_dg->eps2);
}

void 
gk_geometry_mirror_array_new(struct gk_geometry *up)
{
    up->geo_corn.mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_corn.mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_corn.mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_corn.mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->geo_corn.bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->geo_int.cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->geo_int.eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
}

void
gk_geometry_mirror_array_copy(struct gk_geometry *up, struct gkyl_mirror_geo_dg *mirror_geo_dg)
{
    struct gkyl_range global_sub_range;
    int intersect = gkyl_sub_range_intersect(&global_sub_range, &up->global, &up->local);
    gkyl_array_copy_range_to_range(up->geo_corn.mc2p, mirror_geo_dg->mapc2p, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_corn.mc2nu_pos, mirror_geo_dg->mc2nu_pos, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_corn.bmag, mirror_geo_dg->Bmag, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.g_ij, mirror_geo_dg->metric_covar, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.g_ij_neut, mirror_geo_dg->metric_covar_neut, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.dxdz, mirror_geo_dg->tang, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.dzdx, mirror_geo_dg->dual, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.dualmag, mirror_geo_dg->dualmag, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.normals, mirror_geo_dg->normals, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.jacobgeo, mirror_geo_dg->Jc, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.jacobgeo_inv, mirror_geo_dg->Jc_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gij, mirror_geo_dg->metric_contr, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gij_neut, mirror_geo_dg->metric_contr_neut, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.b_i, mirror_geo_dg->b_covar, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.bcart, mirror_geo_dg->b_cart, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.cmag, mirror_geo_dg->C, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.jacobtot, mirror_geo_dg->JB, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.jacobtot_inv, mirror_geo_dg->JB_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.bmag_inv, mirror_geo_dg->Bmag_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.bmag_inv_sq, mirror_geo_dg->Bmag_inv_sq, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gxxj, mirror_geo_dg->gxxj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gxyj, mirror_geo_dg->gxyj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gyyj, mirror_geo_dg->gyyj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.gxzj, mirror_geo_dg->gxzj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->geo_int.eps2, mirror_geo_dg->eps2, &up->local, &global_sub_range);
}
