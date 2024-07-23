#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

void
gkyl_wv_maxwell_free(const struct gkyl_ref_count *ref)
{ 
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_maxwell *maxwell = container_of(base->on_dev, struct wv_maxwell, eqn);
    gkyl_cu_free(maxwell);
  }
  
  struct wv_maxwell *maxwell = container_of(base, struct wv_maxwell, eqn);
  gkyl_free(maxwell);  
}

static inline void
maxwell_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 8; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_new(double c, double e_fact, double b_fact, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_wv_maxwell_cu_dev_new(c, e_fact, b_fact);
  } 
#endif    
  struct wv_maxwell *maxwell = gkyl_malloc(sizeof(struct wv_maxwell));

  maxwell->eqn.type = GKYL_EQN_MAXWELL;
  maxwell->eqn.num_equations = 8;  
  maxwell->eqn.num_waves = 6;
  maxwell->eqn.num_diag = 6; // Ex^2, Ey^2, Ez^2, Bx^2, By^2, Bz^2
  
  maxwell->c = c;
  maxwell->e_fact = e_fact;
  maxwell->b_fact = b_fact;
  
  maxwell->eqn.waves_func = wave;
  maxwell->eqn.qfluct_func = qfluct;

  maxwell->eqn.flux_jump = flux_jump;
  maxwell->eqn.check_inv_func = check_inv;
  maxwell->eqn.max_speed_func = max_speed;
  maxwell->eqn.rotate_to_local_func = rot_to_local;
  maxwell->eqn.rotate_to_global_func = rot_to_global;

  maxwell->eqn.cons_to_riem = cons_to_riem;
  maxwell->eqn.riem_to_cons = riem_to_cons;

  maxwell->eqn.wall_bc_func = maxwell_wall;

  maxwell->eqn.cons_to_diag = maxwell_cons_to_diag;

  maxwell->eqn.source_func = maxwell_source;

  maxwell->eqn.ref_count = gkyl_ref_count_init(maxwell_free);
  maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_maxwell_free);
  maxwell->eqn.on_dev = &maxwell->eqn; // CPU eqn obj points to itself

  return &maxwell->eqn;
}
