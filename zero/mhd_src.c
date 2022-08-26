#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mhd_src.h>


// Makes indexing cleaner
#define DN (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)
#define BX (5)
#define BY (6)
#define BZ (7)
#define PSI_GLM (8)

gkyl_mhd_src *gkyl_mhd_src_new(struct gkyl_mhd_src_inp inp,
                               const struct gkyl_range *local_ext) {
  gkyl_mhd_src *up = gkyl_malloc(sizeof(gkyl_mhd_src));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->divergence_constraint = inp.divergence_constraint;
  up->glm_ch = inp.glm_ch;
  up->glm_alpha = inp.glm_alpha;
  up->dxyz_min = inp.dxyz_min;

  if (up->divergence_constraint == GKYL_MHD_DIVB_EIGHT_WAVES) {
    up->divB_array = gkyl_array_new(GKYL_DOUBLE, 1, local_ext->volume);
  }

  return up;
}

void gkyl_mhd_src_advance(const gkyl_mhd_src *up, double dt,
                          const struct gkyl_range *update_range,
                          struct gkyl_array *q_array,
                          const struct gkyl_array *acc_array) {
  int div_type = up->divergence_constraint;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  // compute divB, B_dot_gradPsi if needed
  if (div_type == GKYL_MHD_DIVB_EIGHT_WAVES) {
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(update_range, iter.idx);

      double *q = gkyl_array_fetch(q_array, lidx);
      double *divB = gkyl_array_fetch(up->divB_array, lidx);

      int idxl[3], idxr[3];
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxl);
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxr);
      
      double my_divB = 0.0;
      for(int d=0; d<up->grid.ndim; ++d) {
        idxl[d]--;
        idxr[d]++;

        lidx = gkyl_range_idx(update_range, idxl);
        double *ql = gkyl_array_fetch(q_array, lidx);
        lidx = gkyl_range_idx(update_range, idxr);
        double *qr = gkyl_array_fetch(q_array, lidx);

        double delta = 2.0 * up->grid.dx[d];
        int Bn = BX + d;
        my_divB += (qr[Bn] - ql[Bn]) / delta;

        idxl[d]++;
        idxr[d]--;
      }
      divB[0] = my_divB;
    }
  }

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);

    double *q = gkyl_array_fetch(q_array, lidx);
    const double *acc = gkyl_array_cfetch(acc_array, lidx);

    if (div_type == GKYL_MHD_DIVB_EIGHT_WAVES) {
      // Powell et al JCP 1999, doi:10.1006/jcph.1999.6299 Equation (25)
      double divB = ((double *)gkyl_array_fetch(up->divB_array, lidx))[0];

      double ux = q[MX] / q[DN], uy = q[MY] / q[DN], uz = q[MZ] / q[DN];
      double Bx = q[BX], By = q[BY], Bz = q[BZ];
      q[MX] -= dt * divB * Bx;
      q[MY] -= dt * divB * By;
      q[MZ] -= dt * divB * Bz;
      q[BX] -= dt * divB * ux;
      q[BY] -= dt * divB * uy;
      q[BZ] -= dt * divB * uz;
      q[ER] -= dt * divB * (ux*Bx + uy*By + uz*Bz);
    } else if (div_type == GKYL_MHD_DIVB_GLM) {
      double ch = up->glm_ch;
      double alpha = up->glm_alpha;

      // Mignone & Tzeferacos, JCP (2010) 229, 2117, Equation (27).
      double rate = alpha * ch / up->dxyz_min;
      q[PSI_GLM] *= exp(- rate * dt);
    }
  }
}

void gkyl_mhd_src_release(gkyl_mhd_src *up) {
  if (up->divergence_constraint == GKYL_MHD_DIVB_EIGHT_WAVES) {
    gkyl_array_release(up->divB_array);
  }
  gkyl_free(up); }
