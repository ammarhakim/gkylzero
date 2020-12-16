#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_util.h>

#include <kernels/vlasov/gkyl_vlasov_kernels.h>

// Types for various kernels
typedef double (*vlasov_vol_t)(const double *w, const double *dxv,
  const double *EM, const double *f, double* out);

typedef void (*vlasov_stream_surf_t)(const double *wl, const double *wr,
  const double *dxvl, const double *dxvr,
  const double *fl, const double *fr, double* outl, double* outr);

typedef double (*vlasov_accel_surf_t)(const double *wl, const double *wr,
  const double *dxvl, const double *dxvr,
  const double amax, const double *EM,
  const double *fl, const double *fr, double* outl, double* outr);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// Volume kernel list
static struct { vlasov_vol_t kernels[3]; } vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_vol_1x1v_p1, vlasov_vol_1x1v_p2 }, // 0
  { NULL, vlasov_vol_1x2v_p1, vlasov_vol_1x2v_p2 }, // 1
  { NULL, vlasov_vol_1x3v_p1, vlasov_vol_1x3v_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_vol_2x2v_p1, vlasov_vol_2x2v_p2 }, // 3
  { NULL, vlasov_vol_2x3v_p1, vlasov_vol_2x3v_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_vol_3x3v_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
static struct { vlasov_stream_surf_t kernels[3]; } stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surf_1x1v_x_p1, vlasov_surf_1x1v_x_p2 }, // 0
  { NULL, vlasov_surf_1x2v_x_p1, vlasov_surf_1x2v_x_p2 }, // 1
  { NULL, vlasov_surf_1x3v_x_p1, vlasov_surf_1x3v_x_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surf_2x2v_x_p1, vlasov_surf_2x2v_x_p2 }, // 3
  { NULL, vlasov_surf_2x3v_x_p1, vlasov_surf_2x3v_x_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_x_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
static struct { vlasov_stream_surf_t kernels[3]; } stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surf_2x2v_y_p1, vlasov_surf_2x2v_y_p2 }, // 3
  { NULL, vlasov_surf_2x3v_y_p1, vlasov_surf_2x3v_y_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_y_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
static struct { vlasov_stream_surf_t kernels[3]; } stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_z_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
static struct { vlasov_accel_surf_t kernels[3]; } accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surf_1x1v_vx_p1, vlasov_surf_1x1v_vx_p2 }, // 0
  { NULL, vlasov_surf_1x2v_vx_p1, vlasov_surf_1x2v_vx_p2 }, // 1
  { NULL, vlasov_surf_1x3v_vx_p1, vlasov_surf_1x3v_vx_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surf_2x2v_vx_p1, vlasov_surf_2x2v_vx_p2 }, // 3
  { NULL, vlasov_surf_2x3v_vx_p1, vlasov_surf_2x3v_vx_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_vx_p1, NULL                   }, // 5
};

// Accel_Surf_Vx_Kernels surface kernel list: vy-direction
static struct { vlasov_accel_surf_t kernels[3]; } accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_surf_1x2v_vy_p1, vlasov_surf_1x2v_vy_p2 }, // 1
  { NULL, vlasov_surf_1x3v_vy_p1, vlasov_surf_1x3v_vy_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surf_2x2v_vy_p1, vlasov_surf_2x2v_vy_p2 }, // 3
  { NULL, vlasov_surf_2x3v_vy_p1, vlasov_surf_2x3v_vy_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_vy_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
static struct { vlasov_accel_surf_t kernels[3]; } accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_surf_1x3v_vz_p1, vlasov_surf_1x3v_vz_p2}, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_surf_2x3v_vz_p1, vlasov_surf_2x3v_vz_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surf_3x3v_vz_p1, NULL }, // 5
};

struct dg_vlasov {
    struct gkyl_dg_eqn eqn; // Base object
    int cdim; // Config-space dimensions
    int pdim; // Phase-space dimensions
    vlasov_vol_t vol; // Volume kernel
    vlasov_stream_surf_t stream_surf[3]; // Surface terms for streaming
    vlasov_accel_surf_t accel_surf[3]; // Surface terms for acceleration

    struct gkyl_range conf_range; // configuration space range
    struct gkyl_array *qmem; // Pointer to q/m*EM field
};

static void
vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_vlasov *vlasov = container_of(base, struct dg_vlasov, eqn);
  gkyl_free(vlasov);
}

void
gkyl_vlasov_set_qmem(const struct gkyl_dg_eqn *eqn, struct gkyl_array *qmem)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->qmem = qmem;
}

static double
vol(const struct gkyl_dg_eqn *eqn, 
  const double*  xc, const double*  dx, const int* idx, const double* qIn, double *qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  assert(vlasov->qmem);
  
  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov->vol(xc, dx, gkyl_array_fetch(vlasov->qmem, cidx),
    qIn, qRhsOut);
}

static double
surf(const struct gkyl_dg_eqn *eqn, int dir,
  const double*  xcL, const double*  xcR, const double*  dxL, const double* dxR,
  double maxsOld, const int*  idxL, const int*  idxR,
  const double* qInL, const double*  qInR, double *qRhsOutL, double *qRhsOutR)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  double amax = 0.0;
  if (dir < vlasov->cdim) {
    vlasov->stream_surf[dir]
      (xcL, xcR, dxL, dxR, qInL, qInR, qRhsOutL, qRhsOutR);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxL);
    amax = vlasov->accel_surf[dir-vlasov->cdim]
      (xcL, xcR, dxL, dxR, maxsOld, gkyl_array_fetch(vlasov->qmem, cidx),
        qInL, qInR, qRhsOutL, qRhsOutR);
  }
  
  return amax;
}

static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcL, const double*  xcR, const double*  dxL, const double*  dxR,
  double maxsOld, const int*  idxL, const int*  idxR,
  const double* qInL, const double* qInR, double *qRhsOutL, double *qRhsOutR)
{
  return 0;
}

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vdim,polyOrder) lst[cv_index[cdim].vdim[vdim]].kernels[polyOrder]

struct gkyl_dg_eqn*
gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov *vlasov = gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int polyOrder = cbasis->polyOrder;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;

  vlasov->vol = CK(vol_kernels,cdim,vdim,polyOrder);

  vlasov->stream_surf[0] = CK(stream_surf_x_kernels,cdim,vdim,polyOrder);
  if (cdim>1)
    vlasov->stream_surf[1] = CK(stream_surf_y_kernels,cdim,vdim,polyOrder);
  if (cdim>2)
    vlasov->stream_surf[2] = CK(stream_surf_z_kernels,cdim,vdim,polyOrder);

  vlasov->accel_surf[0] = CK(accel_surf_vx_kernels,cdim,vdim,polyOrder);
  if (vdim>1)
    vlasov->accel_surf[1] = CK(accel_surf_vy_kernels,cdim,vdim,polyOrder);
  if (vdim>2)
    vlasov->accel_surf[2] = CK(accel_surf_vz_kernels,cdim,vdim,polyOrder);

  // ensure non-NULL pointers
  assert(vlasov->vol);
  for (int i=0; i<cdim; ++i) assert(vlasov->stream_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov->accel_surf[i]);

  vlasov->qmem = 0; 
  vlasov->conf_range = *conf_range;

  // set reference counter
  vlasov->eqn.ref_count = (struct gkyl_ref_count) { vlasov_free, 1 };
  
  return &vlasov->eqn;
}
