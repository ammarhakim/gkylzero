#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_util.h>
#include <gkyl_maxwell_kernels.h>

// Types for various kernels
typedef double (*maxwell_vol_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double *q, double* restrict out);

typedef double (*maxwell_surf_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double tau,
  const double *ql, const double *qc, const double *qr, double* restrict out);

// Volume kernel list
static struct { maxwell_vol_t kernels[3]; } vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_ser_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_ser_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
static struct { maxwell_surf_t kernels[3]; } surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
static struct { maxwell_surf_t kernels[3]; } surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
static struct { maxwell_surf_t kernels[3]; } surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL }, // 2
};

struct dg_maxwell {
    struct gkyl_dg_eqn eqn; // Base object    
    gkyl_maxwell_inp maxwell_data; // Parameters needed by kernels
    maxwell_vol_t vol; // pointer to volume kernel
    maxwell_surf_t surf[3]; // pointers to surface kernels
};

static void
maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_maxwell *maxwell = container_of(base, struct dg_maxwell, eqn);
  free(maxwell);
}

static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* restrict qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell->vol(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  double maxsOld, const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* restrict qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  double maxs = maxwell->maxwell_data.c;
  return maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
    maxs, qInL, qInC, qInR, qRhsOut);
}

static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  double maxsOld, const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* restrict qRhsOut)
{
  return 0;
}

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,polyOrder) lst[cdim-1].kernels[polyOrder]

struct gkyl_dg_eqn*
gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  struct dg_maxwell *maxwell = gkyl_malloc(sizeof(struct dg_maxwell));

  int cdim = cbasis->ndim;
  int polyOrder = cbasis->polyOrder;

  maxwell->eqn.num_equations = 8;
  maxwell->eqn.vol_term = vol;
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  maxwell->maxwell_data.c = lightSpeed;
  maxwell->maxwell_data.chi = lightSpeed*elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = lightSpeed*mgnErrorSpeedFactor;

  maxwell->vol =  CK(vol_kernels, cdim, polyOrder);
  assert(maxwell->vol);

  maxwell->surf[0] = CK(surf_x_kernels, cdim, polyOrder);
  if (cdim>1)
    maxwell->surf[1] = CK(surf_y_kernels, cdim, polyOrder);
  if (cdim>2)
    maxwell->surf[2] = CK(surf_z_kernels, cdim, polyOrder);

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(maxwell->surf[i]);

  // set reference counter
  maxwell->eqn.ref_count = (struct gkyl_ref_count) { maxwell_free, 1 };
  
  return &maxwell->eqn;
}
