#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_vmap_M0_upper_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *vmap, double v_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *p0 = &vmap[0]; 
  double f_nodes[6] = {0.0};
  if (0.6324555320336759*p0[2]-0.9486832980505137*p0[1]+0.7071067811865475*p0[0] > v_thresh) { 
    f_nodes[0] = (-0.447213595499958*f[5])+0.4472135954999579*f[4]+0.6708203932499369*f[3]-0.6708203932499369*f[2]-0.5*f[1]+0.5*f[0]; 
    f_nodes[3] = 0.447213595499958*f[5]+0.4472135954999579*f[4]-0.6708203932499369*(f[3]+f[2])+0.5*(f[1]+f[0]); 
  } 
  if (0.7071067811865475*p0[0]-0.7905694150420947*p0[2] > v_thresh) { 
    f_nodes[1] = 0.5590169943749476*f[5]-0.5590169943749475*f[4]-0.5*f[1]+0.5*f[0]; 
    f_nodes[4] = (-0.5590169943749476*f[5])-0.5590169943749475*f[4]+0.5*(f[1]+f[0]); 
  } 
  if (0.6324555320336759*p0[2]+0.9486832980505137*p0[1]+0.7071067811865475*p0[0] > v_thresh) { 
    f_nodes[2] = (-0.447213595499958*f[5])+0.4472135954999579*f[4]-0.6708203932499369*f[3]+0.6708203932499369*f[2]-0.5*f[1]+0.5*f[0]; 
    f_nodes[5] = 0.447213595499958*f[5]+0.4472135954999579*f[4]+0.6708203932499369*(f[3]+f[2])+0.5*(f[1]+f[0]); 
  } 
  out[0] += (0.392837100659193*f_nodes[5]+0.6285393610547091*f_nodes[4]+0.392837100659193*f_nodes[3]+0.392837100659193*f_nodes[2]+0.6285393610547091*f_nodes[1]+0.392837100659193*f_nodes[0])*volFact; 
  out[1] += (0.392837100659193*f_nodes[5]+0.6285393610547091*f_nodes[4]+0.392837100659193*f_nodes[3]-0.392837100659193*f_nodes[2]-0.6285393610547091*f_nodes[1]-0.392837100659193*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M0_lower_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *vmap, double v_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *p0 = &vmap[0]; 
  double f_nodes[6] = {0.0};
  if (0.6324555320336759*p0[2]-0.9486832980505137*p0[1]+0.7071067811865475*p0[0] < -v_thresh) { 
    f_nodes[0] = (-0.447213595499958*f[5])+0.4472135954999579*f[4]+0.6708203932499369*f[3]-0.6708203932499369*f[2]-0.5*f[1]+0.5*f[0]; 
    f_nodes[3] = 0.447213595499958*f[5]+0.4472135954999579*f[4]-0.6708203932499369*(f[3]+f[2])+0.5*(f[1]+f[0]); 
  } 
  if (0.7071067811865475*p0[0]-0.7905694150420947*p0[2] < -v_thresh) { 
    f_nodes[1] = 0.5590169943749476*f[5]-0.5590169943749475*f[4]-0.5*f[1]+0.5*f[0]; 
    f_nodes[4] = (-0.5590169943749476*f[5])-0.5590169943749475*f[4]+0.5*(f[1]+f[0]); 
  } 
  if (0.6324555320336759*p0[2]+0.9486832980505137*p0[1]+0.7071067811865475*p0[0] < -v_thresh) { 
    f_nodes[2] = (-0.447213595499958*f[5])+0.4472135954999579*f[4]-0.6708203932499369*f[3]+0.6708203932499369*f[2]-0.5*f[1]+0.5*f[0]; 
    f_nodes[5] = 0.447213595499958*f[5]+0.4472135954999579*f[4]+0.6708203932499369*(f[3]+f[2])+0.5*(f[1]+f[0]); 
  } 
  out[0] += (0.392837100659193*f_nodes[5]+0.6285393610547091*f_nodes[4]+0.392837100659193*f_nodes[3]+0.392837100659193*f_nodes[2]+0.6285393610547091*f_nodes[1]+0.392837100659193*f_nodes[0])*volFact; 
  out[1] += (0.392837100659193*f_nodes[5]+0.6285393610547091*f_nodes[4]+0.392837100659193*f_nodes[3]-0.392837100659193*f_nodes[2]-0.6285393610547091*f_nodes[1]-0.392837100659193*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *jacob_vel_inv, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.449489742783178*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 

  out[0] += (p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M3i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *p0 = &vmap[0]; 
  out[0] += (p0[2]*f[4]+p0[1]*f[2]+f[0]*p0[0])*volFact; 
  out[1] += (1.0*p0[2]*f[5]+p0[1]*f[3]+p0[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *vmap, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double *p0 = &vmap[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.414213562373095*gamma[2]*f[4]+1.414213562373095*gamma[1]*f[2]+1.414213562373095*f[0]*gamma[0])*volFact; 
  out[2] += (1.414213562373095*p0[2]*f[4]+1.414213562373095*p0[1]*f[2]+1.414213562373095*f[0]*p0[0])*volFact; 
} 
