#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_sr_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_vlasov_sr {
  struct gkyl_mom_type momt;
  double v_thresh; // Threshold velocity for upper/lower half-plane moments. 
  struct gkyl_range conf_range; // Configuration-space range.
  struct gkyl_range vel_range; // Velocity-space range.
  struct gkyl_mom_vlasov_sr_auxfields auxfields; // Auxiliary fields.
};

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_vlasov_sr_mom_kern_list;

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x1v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x1v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x2v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x2v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x3v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x1v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x1v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x2v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x2v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x3v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_3x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x1v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x1v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x2v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x2v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x3v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x1v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x1v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x2v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x2v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x3v_ser_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_3x3v_ser_p1(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_M0_1x1v_ser_p1, kernel_vlasov_sr_M0_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_M0_1x2v_ser_p1, kernel_vlasov_sr_M0_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_M0_1x3v_ser_p1, kernel_vlasov_sr_M0_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_M0_2x1v_ser_p1, kernel_vlasov_sr_M0_2x1v_ser_p2 }, // 3  
  { NULL, kernel_vlasov_sr_M0_2x2v_ser_p1, kernel_vlasov_sr_M0_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_M0_2x3v_ser_p1, kernel_vlasov_sr_M0_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_M0_3x3v_ser_p1, NULL                            }, // 6
};

// M1i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_m1i_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_M1i_1x1v_ser_p1, kernel_vlasov_sr_M1i_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_M1i_1x2v_ser_p1, kernel_vlasov_sr_M1i_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_M1i_1x3v_ser_p1, kernel_vlasov_sr_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_M1i_2x1v_ser_p1, kernel_vlasov_sr_M1i_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_M1i_2x2v_ser_p1, kernel_vlasov_sr_M1i_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_M1i_2x3v_ser_p1, kernel_vlasov_sr_M1i_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_M1i_3x3v_ser_p1, NULL                             }, // 6
};

// M2 kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_M2_1x1v_ser_p1, kernel_vlasov_sr_M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_M2_1x2v_ser_p1, kernel_vlasov_sr_M2_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_M2_1x3v_ser_p1, kernel_vlasov_sr_M2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_M2_2x1v_ser_p1, kernel_vlasov_sr_M2_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_M2_2x2v_ser_p1, kernel_vlasov_sr_M2_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_M2_2x3v_ser_p1, kernel_vlasov_sr_M2_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_M2_3x3v_ser_p1, NULL                            }, // 6
};

// M3i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_m3i_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_M3i_1x1v_ser_p1, kernel_vlasov_sr_M3i_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_M3i_1x2v_ser_p1, kernel_vlasov_sr_M3i_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_M3i_1x3v_ser_p1, kernel_vlasov_sr_M3i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_M3i_2x1v_ser_p1, kernel_vlasov_sr_M3i_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_M3i_2x2v_ser_p1, kernel_vlasov_sr_M3i_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_M3i_2x3v_ser_p1, kernel_vlasov_sr_M3i_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_M3i_3x3v_ser_p1, NULL                             }, // 6
};

// Ni = (M0, M1i) kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_Ni_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_Ni_1x1v_ser_p1, kernel_vlasov_sr_Ni_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_Ni_1x2v_ser_p1, kernel_vlasov_sr_Ni_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_Ni_1x3v_ser_p1, kernel_vlasov_sr_Ni_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_Ni_2x1v_ser_p1, kernel_vlasov_sr_Ni_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_Ni_2x2v_ser_p1, kernel_vlasov_sr_Ni_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_Ni_2x3v_ser_p1, kernel_vlasov_sr_Ni_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_Ni_3x3v_ser_p1, NULL                            }, // 6
};

// Tij = (M2, M3i (vdim components), Stress tensor (vdim*(vdim+1))/2 components)) kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_Tij_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_Tij_1x1v_ser_p1, kernel_vlasov_sr_Tij_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_Tij_1x2v_ser_p1, kernel_vlasov_sr_Tij_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_Tij_1x3v_ser_p1, kernel_vlasov_sr_Tij_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_Tij_2x1v_ser_p1, kernel_vlasov_sr_Tij_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_Tij_2x2v_ser_p1, kernel_vlasov_sr_Tij_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_Tij_2x3v_ser_p1, kernel_vlasov_sr_Tij_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_Tij_3x3v_ser_p1, NULL                             }, // 6
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_int_mom_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_int_mom_1x1v_ser_p1, kernel_vlasov_sr_int_mom_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_int_mom_1x2v_ser_p1, kernel_vlasov_sr_int_mom_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_int_mom_1x3v_ser_p1, kernel_vlasov_sr_int_mom_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_int_mom_2x1v_ser_p1, kernel_vlasov_sr_int_mom_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_int_mom_2x2v_ser_p1, kernel_vlasov_sr_int_mom_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_int_mom_2x3v_ser_p1, kernel_vlasov_sr_int_mom_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_int_mom_3x3v_ser_p1, NULL                                 }, // 6
};

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_upper_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_upper_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_upper_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_upper_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_upper_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_upper_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_upper_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_upper_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_lower_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_lower_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_lower_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_lower_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_lower_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_lower_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_lower_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_lower_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x1v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x1v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x2v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x2v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x3v_ser_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_3x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_3x3v_ser_p1(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

//
// Serendipity basis kernels with a mapped momentum (four-velocity)-space grid
//

// M0 upper half-plane kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_vmap_m0_upper_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vmap_M0_upper_1x1v_ser_p1, kernel_vlasov_sr_vmap_M0_upper_1x1v_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vmap_M0_upper_2x1v_ser_p1, kernel_vlasov_sr_vmap_M0_upper_2x1v_ser_p2 }, // 3
  { NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M0 lower half-plane kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_vmap_m0_lower_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vmap_M0_lower_1x1v_ser_p1, kernel_vlasov_sr_vmap_M0_lower_1x1v_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vmap_M0_lower_2x1v_ser_p1, kernel_vlasov_sr_vmap_M0_lower_2x1v_ser_p2 }, // 3
  { NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M1i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_vmap_m1i_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vmap_M1i_1x1v_ser_p1, kernel_vlasov_sr_vmap_M1i_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_vmap_M1i_1x2v_ser_p1, kernel_vlasov_sr_vmap_M1i_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_vmap_M1i_1x3v_ser_p1, kernel_vlasov_sr_vmap_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vmap_M1i_2x1v_ser_p1, kernel_vlasov_sr_vmap_M1i_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_vmap_M1i_2x2v_ser_p1, kernel_vlasov_sr_vmap_M1i_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_vmap_M1i_2x3v_ser_p1, kernel_vlasov_sr_vmap_M1i_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_vmap_M1i_3x3v_ser_p1, NULL                                  }, // 6
};

// M3i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_vmap_m3i_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vmap_M3i_1x1v_ser_p1, kernel_vlasov_sr_vmap_M3i_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_vmap_M3i_1x2v_ser_p1, kernel_vlasov_sr_vmap_M3i_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_vmap_M3i_1x3v_ser_p1, kernel_vlasov_sr_vmap_M3i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vmap_M3i_2x1v_ser_p1, kernel_vlasov_sr_vmap_M3i_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_vmap_M3i_2x2v_ser_p1, kernel_vlasov_sr_vmap_M3i_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_vmap_M3i_2x3v_ser_p1, kernel_vlasov_sr_vmap_M3i_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_vmap_M3i_3x3v_ser_p1, NULL                                  }, // 6
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list ser_vmap_int_mom_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vmap_int_mom_1x1v_ser_p1, kernel_vlasov_sr_vmap_int_mom_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_vmap_int_mom_1x2v_ser_p1, kernel_vlasov_sr_vmap_int_mom_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_vmap_int_mom_1x3v_ser_p1, kernel_vlasov_sr_vmap_int_mom_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vmap_int_mom_2x1v_ser_p1, kernel_vlasov_sr_vmap_int_mom_2x1v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_vmap_int_mom_2x2v_ser_p1, kernel_vlasov_sr_vmap_int_mom_2x2v_ser_p2 }, // 4
  { NULL, kernel_vlasov_sr_vmap_int_mom_2x3v_ser_p1, kernel_vlasov_sr_vmap_int_mom_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, kernel_vlasov_sr_vmap_int_mom_3x3v_ser_p1, NULL                                      }, // 6
};

//
// Tensor kernels for moment computation
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x1v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x2v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_1x3v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x2v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M0_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M0_2x3v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M1i_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M1i_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M2_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_M2_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x1v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x2v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_1x3v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x2v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_M3i_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  return vlasov_sr_M3i_2x3v_tensor_p2(xc, dx, idx, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Ni_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Ni_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_Tij_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_Tij_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_int_mom_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_int_mom_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

//
// Tensor basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_m0_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_M0_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_M0_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_M0_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_M0_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_M0_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M1i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_m1i_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_M1i_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_M1i_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_M1i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_M1i_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_M1i_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M2 kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_m2_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_M2_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_M2_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_M2_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_M2_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_M2_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M3i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_m3i_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_M3i_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_M3i_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_M3i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_M3i_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_M3i_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// Ni = (M0, M1i) kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_Ni_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_Ni_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_Ni_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_Ni_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_Ni_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_Ni_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// Tij = (M2, M3i (vdim components), Stress tensor (vdim*(vdim+1))/2 components)) kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_Tij_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_Tij_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_Tij_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_Tij_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_Tij_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_Tij_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_int_mom_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_int_mom_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_int_mom_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_int_mom_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_int_mom_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_int_mom_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

//
// Tensor kernels for moment computation on mapped momentum (four-velocity)-space grid
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_upper_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_upper_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M0_lower_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M0_lower_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    mom_vm_sr->v_thresh, f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M1i_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M1i_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_M3i_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_M3i_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    f, out);  
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x1v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_1x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x2v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

GKYL_CU_DH
static void
kernel_vlasov_sr_vmap_int_mom_2x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return vlasov_sr_vmap_int_mom_2x3v_tensor_p2(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.vmap, vidx),
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.gamma, vidx),
    f, out);
}

//
// Tensor basis kernels with a mapped momentum (four-velocity)-space grid
//

// M0 upper half-plane kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_vmap_m0_upper_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_vmap_M0_upper_1x1v_tensor_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M0 lower half-plane kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_vmap_m0_lower_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_vmap_M0_lower_1x1v_tensor_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M1i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_vmap_m1i_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_vmap_M1i_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_vmap_M1i_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_vmap_M1i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_vmap_M1i_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_vmap_M1i_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// M3i kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_vmap_m3i_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_vmap_M3i_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_vmap_M3i_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_vmap_M3i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_vmap_M3i_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_vmap_M3i_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_vlasov_sr_mom_kern_list tensor_vmap_int_mom_kernels[] = {
  // 1x kernels
  { NULL, NULL, kernel_vlasov_sr_vmap_int_mom_1x1v_tensor_p2 }, // 0
  { NULL, NULL, kernel_vlasov_sr_vmap_int_mom_1x2v_tensor_p2 }, // 1
  { NULL, NULL, kernel_vlasov_sr_vmap_int_mom_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, kernel_vlasov_sr_vmap_int_mom_2x2v_tensor_p2 }, // 4
  { NULL, NULL, kernel_vlasov_sr_vmap_int_mom_2x3v_tensor_p2 }, // 5
  // 3x kernels
  { NULL, NULL, NULL }, // 6
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_vm_sr_free(const struct gkyl_ref_count *ref);
