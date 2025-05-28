#include <gkyl_alloc.h>
#include <gkyl_math.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_geo_gen *
gkyl_mirror_geo_gen_inew(const struct gkyl_mirror_geo_gen_inp *inp)
{
  struct gkyl_mirror_geo_gen *geo = gkyl_malloc(sizeof *geo);

  enum { NPSI, NZ };

  // Outline the grid to loop over
  int nc[2]; // Node count
  nc[NPSI] = inp->comp_grid->cells[0]+1;
  nc[NZ] = inp->comp_grid->cells[2]+1;
  long nctot = nc[NPSI]*nc[NZ];

  struct gkyl_range node_rng;
  gkyl_range_init_from_shape(&node_rng, 2, nc);

  // Compute geometry quantities at node locations
  bool status = true;
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_geo_gen_geom), nctot);

  // Loop over all node locations
  for (int iz=0; iz<nc[NZ]; ++iz) {
    for (int ipsi=0; ipsi<nc[NPSI]; ++ipsi) {
      // Find our location in the array
      int idx[2] = { ipsi, iz };
      long loc = gkyl_range_idx(&node_rng, idx);

      // Access relevant array component
      struct gkyl_mirror_geo_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc);
      const struct gkyl_mirror_grid_gen_geom *grid = gkyl_array_fetch(inp->mirror_grid->nodes_geom, loc);
      const double *node_rz = gkyl_array_fetch(inp->mirror_grid->nodes_rz, loc);
      const double *node_psi = gkyl_array_fetch(inp->mirror_grid->nodes_psi, loc);

      // Copy node locations
      for (int i=0; i<2; ++i)
        g->rz_coord[i] = node_rz[i];
      g->psi = node_psi[0];
      
      // Dual and tangent vectors in cartesian coordinates
      for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
          g->tang[i] = gkyl_vec3_polar_con_to_cart(g->rz_coord[0], 0.0, grid->tang[i]);
          g->dual[i] = gkyl_vec3_polar_con_to_cart(g->rz_coord[0], 0.0, grid->dual[i]);
        }
      }

      // Magnatude of dual vectors
      for (int i=0; i<3; ++i)
        g->dualmag[i] = gkyl_vec3_len(g->dual[i]);

      // Normal vectors
      for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
          g->normal[i].x[j] = g->dual[i].x[j] / gkyl_vec3_len(g->dual[i]);
        }
      }

      // Determine metric tensor
      int count = 0;
      for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
          if (i > j)
            continue;
          g->metric_contr[count] = gkyl_vec3_dot(g->dual[i], g->dual[j]);
          g->metric_covar[count] = gkyl_vec3_dot(g->tang[i], g->tang[j]);
          // For neutrals, we assume the metric is the same for now
          g->metric_contr_neut[count] = g->metric_contr[count];
          g->metric_covar_neut[count] = g->metric_covar[count];
          count++;
        }
      }

      // Determine covariant magnetic field quantities
      g->B_covar.x[0] = grid->B.x[0];
      g->B_covar.x[1] = grid->B.x[1] * g->rz_coord[0] * g->rz_coord[0]; // B_1 = B^1 * r^2
      g->B_covar.x[2] = grid->B.x[2];
        
      // Determine Cartesian components of magnetic field vector
      g->B_cart = gkyl_vec3_polar_con_to_cart(g->rz_coord[0], 0.0, grid->B);

      g->Bmag = gkyl_vec3_len(g->B_cart);
      g->Bmag_inv = 1.0 / g->Bmag;
      g->Bmag_inv_sq = g->Bmag_inv * g->Bmag_inv;

      // Determine scalar quantities combining J and B
      g->Jc = grid->Jc;
      g->Jc_inv = 1.0 / g->Jc;
      g->JB = g->Jc * g->Bmag;
      g->JB_inv = 1.0 / g->JB;

      // Compute cmag
      g->C = g->JB / sqrt(g->metric_covar[5]); // g_33 is the last element in the covariant metric tensor
      g->eps2 = g->Jc * g->metric_contr[5] - g->JB / g->metric_covar[5];
    }
  }
  return geo;
}

void 
gkyl_mirror_geo_gen_release(struct gkyl_mirror_geo_gen *geo) 
{
  gkyl_array_release(geo->nodes_geom);
  gkyl_free(geo);

}
