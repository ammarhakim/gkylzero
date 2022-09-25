#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_diff_boundary_surfvpar_1x1v_ser_p3(const double *w, const double *dxv, const double *nuVtSq, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuVtSq[4]:    Thermal speeds squared times collisionality. 
  // fSkin/Edge:   Distribution function in cells 
  // out:          Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[12] = {0.0}; 
  vol_incr[5] = 4.743416490252569*nuVtSq[3]*fSkin[8]*rdvSq4+4.743416490252569*nuVtSq[2]*fSkin[4]*rdvSq4+4.743416490252569*fSkin[1]*nuVtSq[1]*rdvSq4+4.743416490252569*fSkin[0]*nuVtSq[0]*rdvSq4; 
  vol_incr[7] = 4.16619044897648*nuVtSq[2]*fSkin[8]*rdvSq4+4.16619044897648*nuVtSq[3]*fSkin[4]*rdvSq4+4.242640687119286*nuVtSq[1]*fSkin[4]*rdvSq4+4.242640687119286*fSkin[1]*nuVtSq[2]*rdvSq4+4.743416490252569*fSkin[0]*nuVtSq[1]*rdvSq4+4.743416490252569*nuVtSq[0]*fSkin[1]*rdvSq4; 
  vol_incr[9] = 16.20185174601965*nuVtSq[3]*fSkin[10]*rdvSq4+16.20185174601965*nuVtSq[2]*fSkin[6]*rdvSq4+16.20185174601965*nuVtSq[1]*fSkin[3]*rdvSq4+16.20185174601965*nuVtSq[0]*fSkin[2]*rdvSq4; 
  vol_incr[11] = 14.23024947075771*nuVtSq[2]*fSkin[10]*rdvSq4+14.2302494707577*nuVtSq[3]*fSkin[6]*rdvSq4+14.49137674618944*nuVtSq[1]*fSkin[6]*rdvSq4+14.49137674618944*nuVtSq[2]*fSkin[3]*rdvSq4+16.20185174601965*nuVtSq[0]*fSkin[3]*rdvSq4+16.20185174601965*nuVtSq[1]*fSkin[2]*rdvSq4; 

  double temp_diff[12] = {0.0}; 
  double temp_edge[12] = {0.0}; 
  double diff_incr[12] = {0.0}; 
  double edge_incr[12] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6821077598838398*fSkin[9])-0.6821077598838398*fEdge[9]-1.554766015605323*fSkin[5]+1.554766015605323*fEdge[5]-1.935025511580855*fSkin[2]-1.935025511580855*fEdge[2]-1.3671875*fSkin[0]+1.3671875*fEdge[0]; 
  temp_diff[1] = (-0.6821077598838398*fSkin[11])-0.6821077598838398*fEdge[11]-1.554766015605323*fSkin[7]+1.554766015605323*fEdge[7]-1.935025511580855*fSkin[3]-1.935025511580855*fEdge[3]-1.3671875*fSkin[1]+1.3671875*fEdge[1]; 
  temp_diff[2] = (-1.508772131709791*fSkin[9])-0.8541184610018139*fEdge[9]-3.298087380754754*fSkin[5]+2.087780085064936*fEdge[5]-4.0078125*fSkin[2]-2.6953125*fEdge[2]-2.801050915365293*fSkin[0]+1.935025511580855*fEdge[0]; 
  temp_diff[3] = (-1.508772131709791*fSkin[11])-0.8541184610018139*fEdge[11]-3.298087380754754*fSkin[7]+2.087780085064936*fEdge[7]-4.0078125*fSkin[3]-2.6953125*fEdge[3]-2.801050915365293*fSkin[1]+1.935025511580855*fEdge[1]; 
  temp_diff[4] = (-1.935025511580855*fSkin[6])-1.935025511580855*fEdge[6]-1.3671875*fSkin[4]+1.3671875*fEdge[4]; 
  temp_diff[5] = (-2.792970701173145*fSkin[9])-0.2575079369875949*fEdge[9]-5.8203125*fSkin[5]+1.1328125*fEdge[5]-6.868493903039716*fSkin[2]-1.785203261142481*fEdge[2]-4.734175171112836*fSkin[0]+1.380073204863151*fEdge[0]; 
  temp_diff[6] = (-4.0078125*fSkin[6])-2.6953125*fEdge[6]-2.801050915365293*fSkin[4]+1.935025511580855*fEdge[4]; 
  temp_diff[7] = (-2.792970701173145*fSkin[11])-0.2575079369875949*fEdge[11]-5.8203125*fSkin[7]+1.1328125*fEdge[7]-6.868493903039716*fSkin[3]-1.785203261142481*fEdge[3]-4.734175171112837*fSkin[1]+1.380073204863152*fEdge[1]; 
  temp_diff[8] = (-1.935025511580855*fSkin[10])-1.935025511580855*fEdge[10]-1.3671875*fSkin[8]+1.3671875*fEdge[8]; 
  temp_diff[9] = (-4.8046875*fSkin[9])+1.1953125*fEdge[9]-9.659849020842342*fSkin[5]-1.432800572469438*fEdge[5]-11.13422688383802*fSkin[2]+0.8950343154210625*fEdge[2]-7.585865087193006*fSkin[0]-0.351388846000766*fEdge[0]; 
  temp_diff[10] = (-4.0078125*fSkin[10])-2.6953125*fEdge[10]-2.801050915365293*fSkin[8]+1.935025511580855*fEdge[8]; 
  temp_diff[11] = (-4.8046875*fSkin[11])+1.1953125*fEdge[11]-9.659849020842342*fSkin[7]-1.432800572469438*fEdge[7]-11.13422688383802*fSkin[3]+0.8950343154210625*fEdge[3]-7.585865087193007*fSkin[1]-0.351388846000766*fEdge[1]; 

  temp_edge[2] = (-2.29128784747792*fSkin[9])+1.936491673103709*fSkin[5]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[3] = (-2.29128784747792*fSkin[11])+1.936491673103709*fSkin[7]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  temp_edge[5] = 8.874119674649425*fSkin[9]-7.5*fSkin[5]+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[6] = 0.8660254037844387*fSkin[4]-1.5*fSkin[6]; 
  temp_edge[7] = 8.874119674649425*fSkin[11]-7.5*fSkin[7]+5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 
  temp_edge[9] = (-21.0*fSkin[9])+17.74823934929885*fSkin[5]-13.74772708486752*fSkin[2]+7.937253933193772*fSkin[0]; 
  temp_edge[10] = 0.8660254037844386*fSkin[8]-1.5*fSkin[10]; 
  temp_edge[11] = (-21.0*fSkin[11])+17.74823934929885*fSkin[7]-13.74772708486752*fSkin[3]+7.937253933193771*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSq[3]*temp_diff[8]+0.7071067811865475*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[1]*temp_diff[1]+0.7071067811865475*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.6210590034081186*nuVtSq[2]*temp_diff[8]+0.6210590034081186*nuVtSq[3]*temp_diff[4]+0.6324555320336759*nuVtSq[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSq[2]+0.7071067811865475*nuVtSq[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.7071067811865474*nuVtSq[3]*temp_diff[10]+0.7071067811865475*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[2]; 
  diff_incr[3] = 0.6210590034081187*nuVtSq[2]*temp_diff[10]+0.6210590034081187*nuVtSq[3]*temp_diff[6]+0.632455532033676*nuVtSq[1]*temp_diff[6]+0.6324555320336759*nuVtSq[2]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[3]+0.7071067811865475*nuVtSq[1]*temp_diff[2]; 
  diff_incr[4] = 0.421637021355784*nuVtSq[3]*temp_diff[8]+0.6210590034081186*nuVtSq[1]*temp_diff[8]+0.4517539514526256*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[0]*temp_diff[4]+0.6210590034081186*temp_diff[1]*nuVtSq[3]+0.7071067811865475*temp_diff[0]*nuVtSq[2]+0.6324555320336759*nuVtSq[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSq[1]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[5]; 
  diff_incr[6] = 0.4216370213557839*nuVtSq[3]*temp_diff[10]+0.6210590034081187*nuVtSq[1]*temp_diff[10]+0.4517539514526256*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[0]*temp_diff[6]+0.6210590034081187*nuVtSq[3]*temp_diff[3]+0.632455532033676*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSq[2]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[7]+0.7071067811865475*nuVtSq[1]*temp_diff[5]; 
  diff_incr[8] = 0.421637021355784*nuVtSq[2]*temp_diff[8]+0.7071067811865475*nuVtSq[0]*temp_diff[8]+0.421637021355784*nuVtSq[3]*temp_diff[4]+0.6210590034081186*nuVtSq[1]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSq[3]+0.6210590034081186*temp_diff[1]*nuVtSq[2]; 
  diff_incr[9] = 0.7071067811865474*nuVtSq[1]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[9]; 
  diff_incr[10] = 0.421637021355784*nuVtSq[2]*temp_diff[10]+0.7071067811865475*nuVtSq[0]*temp_diff[10]+0.4216370213557839*nuVtSq[3]*temp_diff[6]+0.6210590034081187*nuVtSq[1]*temp_diff[6]+0.6210590034081187*nuVtSq[2]*temp_diff[3]+0.7071067811865474*temp_diff[2]*nuVtSq[3]; 
  diff_incr[11] = 0.6324555320336759*nuVtSq[2]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[11]+0.7071067811865474*nuVtSq[1]*temp_diff[9]; 

  edge_incr[0] = 0.7071067811865475*nuVtSq[3]*temp_edge[8]+0.7071067811865475*nuVtSq[2]*temp_edge[4]+0.7071067811865475*nuVtSq[1]*temp_edge[1]+0.7071067811865475*nuVtSq[0]*temp_edge[0]; 
  edge_incr[1] = 0.6210590034081186*nuVtSq[2]*temp_edge[8]+0.6210590034081186*nuVtSq[3]*temp_edge[4]+0.6324555320336759*nuVtSq[1]*temp_edge[4]+0.6324555320336759*temp_edge[1]*nuVtSq[2]+0.7071067811865475*nuVtSq[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSq[1]; 
  edge_incr[2] = 0.7071067811865474*nuVtSq[3]*temp_edge[10]+0.7071067811865475*nuVtSq[2]*temp_edge[6]+0.7071067811865475*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[2]; 
  edge_incr[3] = 0.6210590034081187*nuVtSq[2]*temp_edge[10]+0.6210590034081187*nuVtSq[3]*temp_edge[6]+0.632455532033676*nuVtSq[1]*temp_edge[6]+0.6324555320336759*nuVtSq[2]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[3]+0.7071067811865475*nuVtSq[1]*temp_edge[2]; 
  edge_incr[4] = 0.421637021355784*nuVtSq[3]*temp_edge[8]+0.6210590034081186*nuVtSq[1]*temp_edge[8]+0.4517539514526256*nuVtSq[2]*temp_edge[4]+0.7071067811865475*nuVtSq[0]*temp_edge[4]+0.6210590034081186*temp_edge[1]*nuVtSq[3]+0.7071067811865475*temp_edge[0]*nuVtSq[2]+0.6324555320336759*nuVtSq[1]*temp_edge[1]; 
  edge_incr[5] = 0.7071067811865475*nuVtSq[1]*temp_edge[7]+0.7071067811865475*nuVtSq[0]*temp_edge[5]; 
  edge_incr[6] = 0.4216370213557839*nuVtSq[3]*temp_edge[10]+0.6210590034081187*nuVtSq[1]*temp_edge[10]+0.4517539514526256*nuVtSq[2]*temp_edge[6]+0.7071067811865475*nuVtSq[0]*temp_edge[6]+0.6210590034081187*nuVtSq[3]*temp_edge[3]+0.632455532033676*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[2]*temp_edge[2]; 
  edge_incr[7] = 0.6324555320336759*nuVtSq[2]*temp_edge[7]+0.7071067811865475*nuVtSq[0]*temp_edge[7]+0.7071067811865475*nuVtSq[1]*temp_edge[5]; 
  edge_incr[8] = 0.421637021355784*nuVtSq[2]*temp_edge[8]+0.7071067811865475*nuVtSq[0]*temp_edge[8]+0.421637021355784*nuVtSq[3]*temp_edge[4]+0.6210590034081186*nuVtSq[1]*temp_edge[4]+0.7071067811865475*temp_edge[0]*nuVtSq[3]+0.6210590034081186*temp_edge[1]*nuVtSq[2]; 
  edge_incr[9] = 0.7071067811865474*nuVtSq[1]*temp_edge[11]+0.7071067811865475*nuVtSq[0]*temp_edge[9]; 
  edge_incr[10] = 0.421637021355784*nuVtSq[2]*temp_edge[10]+0.7071067811865475*nuVtSq[0]*temp_edge[10]+0.4216370213557839*nuVtSq[3]*temp_edge[6]+0.6210590034081187*nuVtSq[1]*temp_edge[6]+0.6210590034081187*nuVtSq[2]*temp_edge[3]+0.7071067811865474*temp_edge[2]*nuVtSq[3]; 
  edge_incr[11] = 0.6324555320336759*nuVtSq[2]*temp_edge[11]+0.7071067811865475*nuVtSq[0]*temp_edge[11]+0.7071067811865474*nuVtSq[1]*temp_edge[9]; 


  } else { 

  temp_diff[0] = 0.6821077598838398*fSkin[9]+0.6821077598838398*fEdge[9]-1.554766015605323*fSkin[5]+1.554766015605323*fEdge[5]+1.935025511580855*fSkin[2]+1.935025511580855*fEdge[2]-1.3671875*fSkin[0]+1.3671875*fEdge[0]; 
  temp_diff[1] = 0.6821077598838398*fSkin[11]+0.6821077598838398*fEdge[11]-1.554766015605323*fSkin[7]+1.554766015605323*fEdge[7]+1.935025511580855*fSkin[3]+1.935025511580855*fEdge[3]-1.3671875*fSkin[1]+1.3671875*fEdge[1]; 
  temp_diff[2] = (-1.508772131709791*fSkin[9])-0.8541184610018139*fEdge[9]+3.298087380754754*fSkin[5]-2.087780085064936*fEdge[5]-4.0078125*fSkin[2]-2.6953125*fEdge[2]+2.801050915365293*fSkin[0]-1.935025511580855*fEdge[0]; 
  temp_diff[3] = (-1.508772131709791*fSkin[11])-0.8541184610018139*fEdge[11]+3.298087380754754*fSkin[7]-2.087780085064936*fEdge[7]-4.0078125*fSkin[3]-2.6953125*fEdge[3]+2.801050915365293*fSkin[1]-1.935025511580855*fEdge[1]; 
  temp_diff[4] = 1.935025511580855*fSkin[6]+1.935025511580855*fEdge[6]-1.3671875*fSkin[4]+1.3671875*fEdge[4]; 
  temp_diff[5] = 2.792970701173145*fSkin[9]+0.2575079369875949*fEdge[9]-5.8203125*fSkin[5]+1.1328125*fEdge[5]+6.868493903039716*fSkin[2]+1.785203261142481*fEdge[2]-4.734175171112836*fSkin[0]+1.380073204863151*fEdge[0]; 
  temp_diff[6] = (-4.0078125*fSkin[6])-2.6953125*fEdge[6]+2.801050915365293*fSkin[4]-1.935025511580855*fEdge[4]; 
  temp_diff[7] = 2.792970701173145*fSkin[11]+0.2575079369875949*fEdge[11]-5.8203125*fSkin[7]+1.1328125*fEdge[7]+6.868493903039716*fSkin[3]+1.785203261142481*fEdge[3]-4.734175171112837*fSkin[1]+1.380073204863152*fEdge[1]; 
  temp_diff[8] = 1.935025511580855*fSkin[10]+1.935025511580855*fEdge[10]-1.3671875*fSkin[8]+1.3671875*fEdge[8]; 
  temp_diff[9] = (-4.8046875*fSkin[9])+1.1953125*fEdge[9]+9.659849020842342*fSkin[5]+1.432800572469438*fEdge[5]-11.13422688383802*fSkin[2]+0.8950343154210625*fEdge[2]+7.585865087193006*fSkin[0]+0.351388846000766*fEdge[0]; 
  temp_diff[10] = (-4.0078125*fSkin[10])-2.6953125*fEdge[10]+2.801050915365293*fSkin[8]-1.935025511580855*fEdge[8]; 
  temp_diff[11] = (-4.8046875*fSkin[11])+1.1953125*fEdge[11]+9.659849020842342*fSkin[7]+1.432800572469438*fEdge[7]-11.13422688383802*fSkin[3]+0.8950343154210625*fEdge[3]+7.585865087193007*fSkin[1]+0.351388846000766*fEdge[1]; 

  temp_edge[2] = (-2.29128784747792*fSkin[9])-1.936491673103709*fSkin[5]-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  temp_edge[3] = (-2.29128784747792*fSkin[11])-1.936491673103709*fSkin[7]-1.5*fSkin[3]-0.8660254037844386*fSkin[1]; 
  temp_edge[5] = (-8.874119674649425*fSkin[9])-7.5*fSkin[5]-5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[6] = (-1.5*fSkin[6])-0.8660254037844387*fSkin[4]; 
  temp_edge[7] = (-8.874119674649425*fSkin[11])-7.5*fSkin[7]-5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 
  temp_edge[9] = (-21.0*fSkin[9])-17.74823934929885*fSkin[5]-13.74772708486752*fSkin[2]-7.937253933193772*fSkin[0]; 
  temp_edge[10] = (-1.5*fSkin[10])-0.8660254037844386*fSkin[8]; 
  temp_edge[11] = (-21.0*fSkin[11])-17.74823934929885*fSkin[7]-13.74772708486752*fSkin[3]-7.937253933193771*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSq[3]*temp_diff[8]+0.7071067811865475*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[1]*temp_diff[1]+0.7071067811865475*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.6210590034081186*nuVtSq[2]*temp_diff[8]+0.6210590034081186*nuVtSq[3]*temp_diff[4]+0.6324555320336759*nuVtSq[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSq[2]+0.7071067811865475*nuVtSq[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.7071067811865474*nuVtSq[3]*temp_diff[10]+0.7071067811865475*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[2]; 
  diff_incr[3] = 0.6210590034081187*nuVtSq[2]*temp_diff[10]+0.6210590034081187*nuVtSq[3]*temp_diff[6]+0.632455532033676*nuVtSq[1]*temp_diff[6]+0.6324555320336759*nuVtSq[2]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[3]+0.7071067811865475*nuVtSq[1]*temp_diff[2]; 
  diff_incr[4] = 0.421637021355784*nuVtSq[3]*temp_diff[8]+0.6210590034081186*nuVtSq[1]*temp_diff[8]+0.4517539514526256*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[0]*temp_diff[4]+0.6210590034081186*temp_diff[1]*nuVtSq[3]+0.7071067811865475*temp_diff[0]*nuVtSq[2]+0.6324555320336759*nuVtSq[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSq[1]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[5]; 
  diff_incr[6] = 0.4216370213557839*nuVtSq[3]*temp_diff[10]+0.6210590034081187*nuVtSq[1]*temp_diff[10]+0.4517539514526256*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[0]*temp_diff[6]+0.6210590034081187*nuVtSq[3]*temp_diff[3]+0.632455532033676*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSq[2]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[7]+0.7071067811865475*nuVtSq[1]*temp_diff[5]; 
  diff_incr[8] = 0.421637021355784*nuVtSq[2]*temp_diff[8]+0.7071067811865475*nuVtSq[0]*temp_diff[8]+0.421637021355784*nuVtSq[3]*temp_diff[4]+0.6210590034081186*nuVtSq[1]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSq[3]+0.6210590034081186*temp_diff[1]*nuVtSq[2]; 
  diff_incr[9] = 0.7071067811865474*nuVtSq[1]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[9]; 
  diff_incr[10] = 0.421637021355784*nuVtSq[2]*temp_diff[10]+0.7071067811865475*nuVtSq[0]*temp_diff[10]+0.4216370213557839*nuVtSq[3]*temp_diff[6]+0.6210590034081187*nuVtSq[1]*temp_diff[6]+0.6210590034081187*nuVtSq[2]*temp_diff[3]+0.7071067811865474*temp_diff[2]*nuVtSq[3]; 
  diff_incr[11] = 0.6324555320336759*nuVtSq[2]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[11]+0.7071067811865474*nuVtSq[1]*temp_diff[9]; 

  edge_incr[0] = 0.7071067811865475*nuVtSq[3]*temp_edge[8]+0.7071067811865475*nuVtSq[2]*temp_edge[4]+0.7071067811865475*nuVtSq[1]*temp_edge[1]+0.7071067811865475*nuVtSq[0]*temp_edge[0]; 
  edge_incr[1] = 0.6210590034081186*nuVtSq[2]*temp_edge[8]+0.6210590034081186*nuVtSq[3]*temp_edge[4]+0.6324555320336759*nuVtSq[1]*temp_edge[4]+0.6324555320336759*temp_edge[1]*nuVtSq[2]+0.7071067811865475*nuVtSq[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSq[1]; 
  edge_incr[2] = 0.7071067811865474*nuVtSq[3]*temp_edge[10]+0.7071067811865475*nuVtSq[2]*temp_edge[6]+0.7071067811865475*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[2]; 
  edge_incr[3] = 0.6210590034081187*nuVtSq[2]*temp_edge[10]+0.6210590034081187*nuVtSq[3]*temp_edge[6]+0.632455532033676*nuVtSq[1]*temp_edge[6]+0.6324555320336759*nuVtSq[2]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[3]+0.7071067811865475*nuVtSq[1]*temp_edge[2]; 
  edge_incr[4] = 0.421637021355784*nuVtSq[3]*temp_edge[8]+0.6210590034081186*nuVtSq[1]*temp_edge[8]+0.4517539514526256*nuVtSq[2]*temp_edge[4]+0.7071067811865475*nuVtSq[0]*temp_edge[4]+0.6210590034081186*temp_edge[1]*nuVtSq[3]+0.7071067811865475*temp_edge[0]*nuVtSq[2]+0.6324555320336759*nuVtSq[1]*temp_edge[1]; 
  edge_incr[5] = 0.7071067811865475*nuVtSq[1]*temp_edge[7]+0.7071067811865475*nuVtSq[0]*temp_edge[5]; 
  edge_incr[6] = 0.4216370213557839*nuVtSq[3]*temp_edge[10]+0.6210590034081187*nuVtSq[1]*temp_edge[10]+0.4517539514526256*nuVtSq[2]*temp_edge[6]+0.7071067811865475*nuVtSq[0]*temp_edge[6]+0.6210590034081187*nuVtSq[3]*temp_edge[3]+0.632455532033676*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[2]*temp_edge[2]; 
  edge_incr[7] = 0.6324555320336759*nuVtSq[2]*temp_edge[7]+0.7071067811865475*nuVtSq[0]*temp_edge[7]+0.7071067811865475*nuVtSq[1]*temp_edge[5]; 
  edge_incr[8] = 0.421637021355784*nuVtSq[2]*temp_edge[8]+0.7071067811865475*nuVtSq[0]*temp_edge[8]+0.421637021355784*nuVtSq[3]*temp_edge[4]+0.6210590034081186*nuVtSq[1]*temp_edge[4]+0.7071067811865475*temp_edge[0]*nuVtSq[3]+0.6210590034081186*temp_edge[1]*nuVtSq[2]; 
  edge_incr[9] = 0.7071067811865474*nuVtSq[1]*temp_edge[11]+0.7071067811865475*nuVtSq[0]*temp_edge[9]; 
  edge_incr[10] = 0.421637021355784*nuVtSq[2]*temp_edge[10]+0.7071067811865475*nuVtSq[0]*temp_edge[10]+0.4216370213557839*nuVtSq[3]*temp_edge[6]+0.6210590034081187*nuVtSq[1]*temp_edge[6]+0.6210590034081187*nuVtSq[2]*temp_edge[3]+0.7071067811865474*temp_edge[2]*nuVtSq[3]; 
  edge_incr[11] = 0.6324555320336759*nuVtSq[2]*temp_edge[11]+0.7071067811865475*nuVtSq[0]*temp_edge[11]+0.7071067811865474*nuVtSq[1]*temp_edge[9]; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
  out[6] += edge_incr[6]*rdvSq4+diff_incr[6]*rdvSq4+vol_incr[6]; 
  out[7] += edge_incr[7]*rdvSq4+diff_incr[7]*rdvSq4+vol_incr[7]; 
  out[8] += edge_incr[8]*rdvSq4+diff_incr[8]*rdvSq4+vol_incr[8]; 
  out[9] += edge_incr[9]*rdvSq4+diff_incr[9]*rdvSq4+vol_incr[9]; 
  out[10] += edge_incr[10]*rdvSq4+diff_incr[10]*rdvSq4+vol_incr[10]; 
  out[11] += edge_incr[11]*rdvSq4+diff_incr[11]*rdvSq4+vol_incr[11]; 
} 
