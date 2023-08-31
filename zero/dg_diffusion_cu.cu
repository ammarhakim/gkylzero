/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion.h>    
#include <gkyl_dg_diffusion_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__ void static
dg_diffusion_set_cu_dev_ptrs(struct dg_diffusion* diffusion, 
  enum gkyl_basis_type b_type, int cdim, int poly_order, 
  int order, enum gkyl_diffusion_id diffusion_id)
{
  const gkyl_dg_diffusion_vol_kern_list* vol_kernels;
  const gkyl_dg_diffusion_surf_kern_list* surf_x_kernels;
  const gkyl_dg_diffusion_surf_kern_list* surf_y_kernels;
  const gkyl_dg_diffusion_surf_kern_list* surf_z_kernels; 

  // Because our scheme uses 1-cell recovery, the volume kernel
  // *only* returns the CFL rate and is thus independent of basis
  // type (e.g., serendipity vs. tensor) or equation system (e.g., Euler vs. PKPM)
  if (order == 4) 
    vol_kernels = ser_vol4_kernels;
  else if (order == 6)
    vol_kernels = ser_vol6_kernels;
  else 
    vol_kernels = ser_vol_kernels;

  // Isotropic diffusion
  if (diffusion_id == GKYL_ISO_DIFFUSION) {
    diffusion->eqn.num_equations = 1;
    if (order == 4) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf4_x_kernels;
          surf_y_kernels = ser_surf4_y_kernels;
          surf_z_kernels = ser_surf4_z_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf4_x_kernels;
          surf_y_kernels = ten_surf4_y_kernels;
          surf_z_kernels = ten_surf4_z_kernels;
          break;

        default:
          assert(false);
          break;    
      } 
    }
    else if (order == 6) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf6_x_kernels;
          surf_y_kernels = ser_surf6_y_kernels;
          surf_z_kernels = ser_surf6_z_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf6_x_kernels;
          surf_y_kernels = ten_surf6_y_kernels;
          surf_z_kernels = ten_surf6_z_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
    else {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf_x_kernels;
          surf_y_kernels = ser_surf_y_kernels;
          surf_z_kernels = ser_surf_z_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf_x_kernels;
          surf_y_kernels = ten_surf_y_kernels;
          surf_z_kernels = ten_surf_z_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
  } // PKPM isotropic diffusion (3 components: rhoux, rhouy, rhouz)
  else if (diffusion_id == GKYL_PKPM_DIFFUSION) {
    diffusion->eqn.num_equations = 3;
    if (order == 4) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf4_x_pkpm_kernels;
          surf_y_kernels = ser_surf4_y_pkpm_kernels;
          surf_z_kernels = ser_surf4_z_pkpm_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf4_x_pkpm_kernels;
          surf_y_kernels = ten_surf4_y_pkpm_kernels;
          surf_z_kernels = ten_surf4_z_pkpm_kernels;
          break;

        default:
          assert(false);
          break;    
      } 
    }
    else if (order == 6) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf6_x_pkpm_kernels;
          surf_y_kernels = ser_surf6_y_pkpm_kernels;
          surf_z_kernels = ser_surf6_z_pkpm_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf6_x_pkpm_kernels;
          surf_y_kernels = ten_surf6_y_pkpm_kernels;
          surf_z_kernels = ten_surf6_z_pkpm_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
    else {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf_x_pkpm_kernels;
          surf_y_kernels = ser_surf_y_pkpm_kernels;
          surf_z_kernels = ser_surf_z_pkpm_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf_x_pkpm_kernels;
          surf_y_kernels = ten_surf_y_pkpm_kernels;
          surf_z_kernels = ten_surf_z_pkpm_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
  } // Isothermal Euler isotropic diffusion (4 components: rho, rhoux, rhouy, rhouz)
  else if (diffusion_id == GKYL_ISO_EULER_DIFFUSION) {
    diffusion->eqn.num_equations = 4;
    if (order == 4) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf4_x_iso_euler_kernels;
          surf_y_kernels = ser_surf4_y_iso_euler_kernels;
          surf_z_kernels = ser_surf4_z_iso_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf4_x_iso_euler_kernels;
          surf_y_kernels = ten_surf4_y_iso_euler_kernels;
          surf_z_kernels = ten_surf4_z_iso_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      } 
    }
    else if (order == 6) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf6_x_iso_euler_kernels;
          surf_y_kernels = ser_surf6_y_iso_euler_kernels;
          surf_z_kernels = ser_surf6_z_iso_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf6_x_iso_euler_kernels;
          surf_y_kernels = ten_surf6_y_iso_euler_kernels;
          surf_z_kernels = ten_surf6_z_iso_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
    else {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf_x_iso_euler_kernels;
          surf_y_kernels = ser_surf_y_iso_euler_kernels;
          surf_z_kernels = ser_surf_z_iso_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf_x_iso_euler_kernels;
          surf_y_kernels = ten_surf_y_iso_euler_kernels;
          surf_z_kernels = ten_surf_z_iso_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
  } // grad^2 Euler isotropic diffusion (5 components: rho, rhoux, rhouy, rhouz, Energy)
  else if (diffusion_id == GKYL_EULER_DIFFUSION) {
    diffusion->eqn.num_equations = 5;
    if (order == 4) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf4_x_euler_kernels;
          surf_y_kernels = ser_surf4_y_euler_kernels;
          surf_z_kernels = ser_surf4_z_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf4_x_euler_kernels;
          surf_y_kernels = ten_surf4_y_euler_kernels;
          surf_z_kernels = ten_surf4_z_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      } 
    }
    else if (order == 6) {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf6_x_euler_kernels;
          surf_y_kernels = ser_surf6_y_euler_kernels;
          surf_z_kernels = ser_surf6_z_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf6_x_euler_kernels;
          surf_y_kernels = ten_surf6_y_euler_kernels;
          surf_z_kernels = ten_surf6_z_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
    else {
      switch (b_type) {
        case GKYL_BASIS_MODAL_SERENDIPITY:
          surf_x_kernels = ser_surf_x_euler_kernels;
          surf_y_kernels = ser_surf_y_euler_kernels;
          surf_z_kernels = ser_surf_z_euler_kernels;
          break;

        case GKYL_BASIS_MODAL_TENSOR:
          surf_x_kernels = ten_surf_x_euler_kernels;
          surf_y_kernels = ten_surf_y_euler_kernels;
          surf_z_kernels = ten_surf_z_euler_kernels;
          break;

        default:
          assert(false);
          break;    
      }       
    }
  }
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;
  diffusion->eqn.vol_term = CK(vol_kernels, cdim, poly_order);
  diffusion->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    diffusion->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    diffusion->surf[2] = CK(surf_z_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_cu_dev_new(const struct gkyl_basis* cbasis, double D, int order, enum gkyl_diffusion_id diffusion_id)
{
  struct dg_diffusion* diffusion = (struct dg_diffusion*) gkyl_malloc(sizeof(struct dg_diffusion));

  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_free);

  diffusion->D = D;

  // copy the host struct to device struct
  struct dg_diffusion* diffusion_cu = (struct dg_diffusion*) gkyl_cu_malloc(sizeof(struct dg_diffusion));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(struct dg_diffusion), GKYL_CU_MEMCPY_H2D);
  dg_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order, order, diffusion_id);

  // set parent on_dev pointer
  diffusion->eqn.on_dev = &diffusion_cu->eqn;

  return &diffusion->eqn;
}
