#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>


// Object type
typedef struct gkyl_deflated_fem_poisson gkyl_deflated_fem_poisson;

struct gkyl_deflated_fem_poisson* gkyl_deflated_fem_poisson_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, struct gkyl_range local, struct gkyl_range local_ext, 
  struct gkyl_array *epsilon, struct gkyl_array *kSq, struct gkyl_poisson_bc poisson_bc, bool use_gpu);


void gkyl_deflated_fem_poisson_advance(struct gkyl_deflated_fem_poisson* up, struct gkyl_array *field, struct gkyl_array* phi);


void gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson* up);
