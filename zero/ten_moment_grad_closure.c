#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_ten_moment_grad_closure_priv.h>

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

gkyl_ten_moment_grad_closure*
gkyl_ten_moment_grad_closure_new(struct gkyl_ten_moment_grad_closure_inp inp)
{
  gkyl_ten_moment_grad_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_grad_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;
  up->omega = inp.omega;
  up->cfl = inp.cfl;
  up->mag = inp.mag;

  if (inp.comm)
    up->comm = gkyl_comm_acquire(inp.comm);
  else
    up->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  grad_closure_calc_q_choose(up);
  
  return up;
}

struct gkyl_ten_moment_grad_closure_status
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces,
  const struct gkyl_range *heat_flux_range, const struct gkyl_range *update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, double dt, struct gkyl_array *heat_flux,
  struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };

  double cfla = 0.0, cfl = gces->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // deliberately a double

  gkyl_array_clear(rhs, 0.0);

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(update_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(heat_flux_range, offsets_centers);

  const double* fluid_d[sz[ndim-1]];
  const double* em_tot_d[sz[ndim-1]];
  double *rhs_d[sz[ndim-1]];

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_range);
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_range, iter_vertex.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
      rhs_d[i] = gkyl_array_fetch(rhs, linc_center + offsets_vertices[i]);
    }

    cfla = gces->calc_q(gces, fluid_d, em_tot_d, gkyl_array_fetch(cflrate, linc_center),
      cfla, dt, rhs_d);
  }

  if (cfla > cflm)
    is_cfl_violated = 1.0;

  // compute actual CFL, status & max-speed across all domains
  double red_vars[2] = { cfla, is_cfl_violated };
  double red_vars_global[2] = { 0.0, 0.0 };
  gkyl_comm_allreduce(gces->comm, GKYL_DOUBLE, GKYL_MAX, 2, &red_vars, &red_vars_global);

  cfla = red_vars_global[0];
  is_cfl_violated = red_vars_global[1];

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  if (is_cfl_violated > 0.0)
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_ten_moment_grad_closure_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
    };

  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).
  return (struct gkyl_ten_moment_grad_closure_status) {
    .success = is_cfl_violated > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
  };

}

void
gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure* up)
{
  gkyl_comm_release(up->comm);
  free(up);
}
