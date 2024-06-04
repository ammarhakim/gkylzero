#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_1x2v_ser_p1(const double *dv_min, const double *vmap, const double *vmapSq, 
    double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
    const double *bmag, double* GKYL_RESTRICT drag_rad_nu_surf, double* GKYL_RESTRICT drag_rad_nu) 
{ 
  // dv_min: minimum cell lengths in velocity space.
  // vmap: velocity space mapping, from computational to physical.
  // vmapSq: squared velocity space mapping.
  // charge: elementary charge.
  // mass: mass of electron.
  // a, alpha, beta, gamma, v0: Radiation fitting parameters.
  // bmag: magnetic field amplitude.
  // drag_rad_nu_surf: Radiation drag in direction dir at corresponding surfaces.
  //                   Note: Each cell owns their *lower* edge surface evaluation.
  // drag_rad_nu: Radiation drag direction dir volume expansion.

  double scaled_v0 = v0/sqrt(mass/(2.0*fabs(charge)));
  double c_const = 8.0*fabs(charge)/sqrt(M_PI)*pow(2*fabs(charge)/mass,gamma/2.0); 
  double const_mult = a*(alpha+beta)/c_const;
  double vmag = 0.0;
  double vmag_surf = 0.0;

  double vsqnu_surf_nodal[6] = {0.0};
  double vsqnu_nodal[12] = {0.0};
  vmag_surf = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[0] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  vmag_surf = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[1] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  vmag_surf = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[2] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  vmag_surf = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[3] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  vmag_surf = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[4] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  vmag_surf = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
  vsqnu_surf_nodal[5] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));

  drag_rad_nu_surf[0] = 0.2777777777777778*vsqnu_surf_nodal[5]+0.4444444444444444*vsqnu_surf_nodal[4]+0.2777777777777778*vsqnu_surf_nodal[3]+0.2777777777777778*vsqnu_surf_nodal[2]+0.4444444444444444*vsqnu_surf_nodal[1]+0.2777777777777778*vsqnu_surf_nodal[0]; 
  drag_rad_nu_surf[1] = 0.2777777777777778*vsqnu_surf_nodal[5]+0.4444444444444444*vsqnu_surf_nodal[4]+0.2777777777777778*vsqnu_surf_nodal[3]-0.2777777777777778*vsqnu_surf_nodal[2]-0.4444444444444444*vsqnu_surf_nodal[1]-0.2777777777777778*vsqnu_surf_nodal[0]; 
  drag_rad_nu_surf[2] = 0.3726779962499649*vsqnu_surf_nodal[5]-0.3726779962499649*vsqnu_surf_nodal[3]+0.3726779962499649*vsqnu_surf_nodal[2]-0.3726779962499649*vsqnu_surf_nodal[0]; 
  drag_rad_nu_surf[3] = 0.3726779962499649*vsqnu_surf_nodal[5]-0.3726779962499649*vsqnu_surf_nodal[3]-0.3726779962499649*vsqnu_surf_nodal[2]+0.3726779962499649*vsqnu_surf_nodal[0]; 
  drag_rad_nu_surf[4] = 0.2484519974999766*vsqnu_surf_nodal[5]-0.4969039949999532*vsqnu_surf_nodal[4]+0.2484519974999766*vsqnu_surf_nodal[3]+0.2484519974999766*vsqnu_surf_nodal[2]-0.4969039949999532*vsqnu_surf_nodal[1]+0.2484519974999766*vsqnu_surf_nodal[0]; 
  drag_rad_nu_surf[5] = 0.2484519974999767*vsqnu_surf_nodal[5]-0.4969039949999535*vsqnu_surf_nodal[4]+0.2484519974999767*vsqnu_surf_nodal[3]-0.2484519974999767*vsqnu_surf_nodal[2]+0.4969039949999535*vsqnu_surf_nodal[1]-0.2484519974999767*vsqnu_surf_nodal[0]; 

  vmag = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[0] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[1] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[2] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[3] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*(0.7071067811865475*(0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[4] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[5] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[6] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[7] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.7071067811865475*(bmag[1]+bmag[0])*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3]))/mass);
  vsqnu_nodal[8] = 2.0*(0.7071067811865475*vmap[2]-0.7071067811865475*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]-0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.4999999999999999*(bmag[1]+bmag[0])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[9] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*(0.4999999999999999*(bmag[1]+bmag[0])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[10] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  vmag = sqrt((0.6324555320336759*vmapSq[2]+0.9486832980505137*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*(0.4999999999999999*(bmag[1]+bmag[0])*(vmap[3]+vmap[2]))/mass);
  vsqnu_nodal[11] = 2.0*(0.7071067811865475*(vmap[3]+vmap[2]))*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));

  drag_rad_nu[0] = 0.1964185503295965*vsqnu_nodal[11]+0.3142696805273545*vsqnu_nodal[10]+0.1964185503295965*vsqnu_nodal[9]+0.1964185503295965*vsqnu_nodal[8]+0.3142696805273545*vsqnu_nodal[7]+0.1964185503295965*vsqnu_nodal[6]+0.1964185503295965*vsqnu_nodal[5]+0.3142696805273545*vsqnu_nodal[4]+0.1964185503295965*vsqnu_nodal[3]+0.1964185503295965*vsqnu_nodal[2]+0.3142696805273545*vsqnu_nodal[1]+0.1964185503295965*vsqnu_nodal[0]; 
  drag_rad_nu[1] = 0.1964185503295965*vsqnu_nodal[11]+0.3142696805273545*vsqnu_nodal[10]+0.1964185503295965*vsqnu_nodal[9]+0.1964185503295965*vsqnu_nodal[8]+0.3142696805273545*vsqnu_nodal[7]+0.1964185503295965*vsqnu_nodal[6]-0.1964185503295965*vsqnu_nodal[5]-0.3142696805273545*vsqnu_nodal[4]-0.1964185503295965*vsqnu_nodal[3]-0.1964185503295965*vsqnu_nodal[2]-0.3142696805273545*vsqnu_nodal[1]-0.1964185503295965*vsqnu_nodal[0]; 
  drag_rad_nu[2] = 0.2635231383473649*vsqnu_nodal[11]-0.2635231383473649*vsqnu_nodal[9]+0.2635231383473649*vsqnu_nodal[8]-0.2635231383473649*vsqnu_nodal[6]+0.2635231383473649*vsqnu_nodal[5]-0.2635231383473649*vsqnu_nodal[3]+0.2635231383473649*vsqnu_nodal[2]-0.2635231383473649*vsqnu_nodal[0]; 
  drag_rad_nu[3] = 0.1964185503295965*vsqnu_nodal[11]+0.3142696805273545*vsqnu_nodal[10]+0.1964185503295965*vsqnu_nodal[9]-0.1964185503295965*vsqnu_nodal[8]-0.3142696805273545*vsqnu_nodal[7]-0.1964185503295965*vsqnu_nodal[6]+0.1964185503295965*vsqnu_nodal[5]+0.3142696805273545*vsqnu_nodal[4]+0.1964185503295965*vsqnu_nodal[3]-0.1964185503295965*vsqnu_nodal[2]-0.3142696805273545*vsqnu_nodal[1]-0.1964185503295965*vsqnu_nodal[0]; 
  drag_rad_nu[4] = 0.2635231383473649*vsqnu_nodal[11]-0.2635231383473649*vsqnu_nodal[9]+0.2635231383473649*vsqnu_nodal[8]-0.2635231383473649*vsqnu_nodal[6]-0.2635231383473649*vsqnu_nodal[5]+0.2635231383473649*vsqnu_nodal[3]-0.2635231383473649*vsqnu_nodal[2]+0.2635231383473649*vsqnu_nodal[0]; 
  drag_rad_nu[5] = 0.1964185503295965*vsqnu_nodal[11]+0.3142696805273545*vsqnu_nodal[10]+0.1964185503295965*vsqnu_nodal[9]-0.1964185503295965*vsqnu_nodal[8]-0.3142696805273545*vsqnu_nodal[7]-0.1964185503295965*vsqnu_nodal[6]-0.1964185503295965*vsqnu_nodal[5]-0.3142696805273545*vsqnu_nodal[4]-0.1964185503295965*vsqnu_nodal[3]+0.1964185503295965*vsqnu_nodal[2]+0.3142696805273545*vsqnu_nodal[1]+0.1964185503295965*vsqnu_nodal[0]; 
  drag_rad_nu[6] = 0.2635231383473649*vsqnu_nodal[11]-0.2635231383473649*vsqnu_nodal[9]-0.2635231383473649*vsqnu_nodal[8]+0.2635231383473649*vsqnu_nodal[6]+0.2635231383473649*vsqnu_nodal[5]-0.2635231383473649*vsqnu_nodal[3]-0.2635231383473649*vsqnu_nodal[2]+0.2635231383473649*vsqnu_nodal[0]; 
  drag_rad_nu[7] = 0.2635231383473649*vsqnu_nodal[11]-0.2635231383473649*vsqnu_nodal[9]-0.2635231383473649*vsqnu_nodal[8]+0.2635231383473649*vsqnu_nodal[6]-0.2635231383473649*vsqnu_nodal[5]+0.2635231383473649*vsqnu_nodal[3]+0.2635231383473649*vsqnu_nodal[2]-0.2635231383473649*vsqnu_nodal[0]; 
  drag_rad_nu[8] = 0.1756820922315766*vsqnu_nodal[11]-0.3513641844631533*vsqnu_nodal[10]+0.1756820922315766*vsqnu_nodal[9]+0.1756820922315766*vsqnu_nodal[8]-0.3513641844631533*vsqnu_nodal[7]+0.1756820922315766*vsqnu_nodal[6]+0.1756820922315766*vsqnu_nodal[5]-0.3513641844631533*vsqnu_nodal[4]+0.1756820922315766*vsqnu_nodal[3]+0.1756820922315766*vsqnu_nodal[2]-0.3513641844631533*vsqnu_nodal[1]+0.1756820922315766*vsqnu_nodal[0]; 
  drag_rad_nu[9] = 0.1756820922315767*vsqnu_nodal[11]-0.3513641844631534*vsqnu_nodal[10]+0.1756820922315767*vsqnu_nodal[9]+0.1756820922315767*vsqnu_nodal[8]-0.3513641844631534*vsqnu_nodal[7]+0.1756820922315767*vsqnu_nodal[6]-0.1756820922315767*vsqnu_nodal[5]+0.3513641844631534*vsqnu_nodal[4]-0.1756820922315767*vsqnu_nodal[3]-0.1756820922315767*vsqnu_nodal[2]+0.3513641844631534*vsqnu_nodal[1]-0.1756820922315767*vsqnu_nodal[0]; 
  drag_rad_nu[10] = 0.1756820922315767*vsqnu_nodal[11]-0.3513641844631534*vsqnu_nodal[10]+0.1756820922315767*vsqnu_nodal[9]-0.1756820922315767*vsqnu_nodal[8]+0.3513641844631534*vsqnu_nodal[7]-0.1756820922315767*vsqnu_nodal[6]+0.1756820922315767*vsqnu_nodal[5]-0.3513641844631534*vsqnu_nodal[4]+0.1756820922315767*vsqnu_nodal[3]-0.1756820922315767*vsqnu_nodal[2]+0.3513641844631534*vsqnu_nodal[1]-0.1756820922315767*vsqnu_nodal[0]; 
  drag_rad_nu[11] = 0.1756820922315766*vsqnu_nodal[11]-0.3513641844631533*vsqnu_nodal[10]+0.1756820922315766*vsqnu_nodal[9]-0.1756820922315766*vsqnu_nodal[8]+0.3513641844631533*vsqnu_nodal[7]-0.1756820922315766*vsqnu_nodal[6]-0.1756820922315766*vsqnu_nodal[5]+0.3513641844631533*vsqnu_nodal[4]-0.1756820922315766*vsqnu_nodal[3]+0.1756820922315766*vsqnu_nodal[2]-0.3513641844631533*vsqnu_nodal[1]+0.1756820922315766*vsqnu_nodal[0]; 

  double vsqnu_corner_vol_nod[12];
  double vsqnu_corner_surf_nod[6];
  vsqnu_corner_vol_nod[0] = 0.02516526685334564*vsqnu_nodal[11]-0.08931639747704079*vsqnu_nodal[10]+0.1981257268392565*vsqnu_nodal[9]-0.09391805448271529*vsqnu_nodal[8]+0.3333333333333333*vsqnu_nodal[7]-0.7394152788506181*vsqnu_nodal[6]-0.09391805448271529*vsqnu_nodal[5]+0.3333333333333333*vsqnu_nodal[4]-0.7394152788506181*vsqnu_nodal[3]+0.3505069510775155*vsqnu_nodal[2]-1.244016935856292*vsqnu_nodal[1]+2.759535388563216*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[1] = (-0.09391805448271529*vsqnu_nodal[11])+0.3333333333333333*vsqnu_nodal[10]-0.7394152788506181*vsqnu_nodal[9]+0.02516526685334564*vsqnu_nodal[8]-0.08931639747704079*vsqnu_nodal[7]+0.1981257268392565*vsqnu_nodal[6]+0.3505069510775155*vsqnu_nodal[5]-1.244016935856292*vsqnu_nodal[4]+2.759535388563216*vsqnu_nodal[3]-0.09391805448271529*vsqnu_nodal[2]+0.3333333333333333*vsqnu_nodal[1]-0.7394152788506181*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[2] = 0.1339745962155614*vsqnu_nodal[10]-0.5*vsqnu_nodal[7]-0.5*vsqnu_nodal[4]+1.866025403784439*vsqnu_nodal[1]; 
  vsqnu_corner_vol_nod[3] = (-0.5*vsqnu_nodal[10])+0.1339745962155614*vsqnu_nodal[7]+1.866025403784439*vsqnu_nodal[4]-0.5*vsqnu_nodal[1]; 
  vsqnu_corner_vol_nod[4] = 0.1981257268392565*vsqnu_nodal[11]-0.08931639747704079*vsqnu_nodal[10]+0.02516526685334564*vsqnu_nodal[9]-0.7394152788506181*vsqnu_nodal[8]+0.3333333333333333*vsqnu_nodal[7]-0.09391805448271529*vsqnu_nodal[6]-0.7394152788506181*vsqnu_nodal[5]+0.3333333333333333*vsqnu_nodal[4]-0.09391805448271529*vsqnu_nodal[3]+2.759535388563216*vsqnu_nodal[2]-1.244016935856292*vsqnu_nodal[1]+0.3505069510775155*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[5] = (-0.7394152788506181*vsqnu_nodal[11])+0.3333333333333333*vsqnu_nodal[10]-0.09391805448271529*vsqnu_nodal[9]+0.1981257268392565*vsqnu_nodal[8]-0.08931639747704079*vsqnu_nodal[7]+0.02516526685334564*vsqnu_nodal[6]+2.759535388563216*vsqnu_nodal[5]-1.244016935856292*vsqnu_nodal[4]+0.3505069510775155*vsqnu_nodal[3]-0.7394152788506181*vsqnu_nodal[2]+0.3333333333333333*vsqnu_nodal[1]-0.09391805448271529*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[6] = (-0.09391805448271529*vsqnu_nodal[11])+0.3333333333333333*vsqnu_nodal[10]-0.7394152788506181*vsqnu_nodal[9]+0.3505069510775155*vsqnu_nodal[8]-1.244016935856292*vsqnu_nodal[7]+2.759535388563216*vsqnu_nodal[6]+0.02516526685334564*vsqnu_nodal[5]-0.08931639747704079*vsqnu_nodal[4]+0.1981257268392565*vsqnu_nodal[3]-0.09391805448271529*vsqnu_nodal[2]+0.3333333333333333*vsqnu_nodal[1]-0.7394152788506181*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[7] = 0.3505069510775155*vsqnu_nodal[11]-1.244016935856292*vsqnu_nodal[10]+2.759535388563216*vsqnu_nodal[9]-0.09391805448271529*vsqnu_nodal[8]+0.3333333333333333*vsqnu_nodal[7]-0.7394152788506181*vsqnu_nodal[6]-0.09391805448271529*vsqnu_nodal[5]+0.3333333333333333*vsqnu_nodal[4]-0.7394152788506181*vsqnu_nodal[3]+0.02516526685334564*vsqnu_nodal[2]-0.08931639747704079*vsqnu_nodal[1]+0.1981257268392565*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[8] = (-0.5*vsqnu_nodal[10])+1.866025403784439*vsqnu_nodal[7]+0.1339745962155614*vsqnu_nodal[4]-0.5*vsqnu_nodal[1]; 
  vsqnu_corner_vol_nod[9] = 1.866025403784439*vsqnu_nodal[10]-0.5*vsqnu_nodal[7]-0.5*vsqnu_nodal[4]+0.1339745962155614*vsqnu_nodal[1]; 
  vsqnu_corner_vol_nod[10] = (-0.7394152788506181*vsqnu_nodal[11])+0.3333333333333333*vsqnu_nodal[10]-0.09391805448271529*vsqnu_nodal[9]+2.759535388563216*vsqnu_nodal[8]-1.244016935856292*vsqnu_nodal[7]+0.3505069510775155*vsqnu_nodal[6]+0.1981257268392565*vsqnu_nodal[5]-0.08931639747704079*vsqnu_nodal[4]+0.02516526685334564*vsqnu_nodal[3]-0.7394152788506181*vsqnu_nodal[2]+0.3333333333333333*vsqnu_nodal[1]-0.09391805448271529*vsqnu_nodal[0]; 
  vsqnu_corner_vol_nod[11] = 2.759535388563216*vsqnu_nodal[11]-1.244016935856292*vsqnu_nodal[10]+0.3505069510775155*vsqnu_nodal[9]-0.7394152788506181*vsqnu_nodal[8]+0.3333333333333333*vsqnu_nodal[7]-0.09391805448271529*vsqnu_nodal[6]-0.7394152788506181*vsqnu_nodal[5]+0.3333333333333333*vsqnu_nodal[4]-0.09391805448271529*vsqnu_nodal[3]+0.1981257268392565*vsqnu_nodal[2]-0.08931639747704079*vsqnu_nodal[1]+0.02516526685334564*vsqnu_nodal[0]; 
  vsqnu_corner_surf_nod[0] = (-0.0687527876293696*vsqnu_surf_nodal[5])+0.2440169358562925*vsqnu_surf_nodal[4]-0.5412895520113615*vsqnu_surf_nodal[3]+0.2565888965948001*vsqnu_surf_nodal[2]-0.9106836025229592*vsqnu_surf_nodal[1]+2.020120109712598*vsqnu_surf_nodal[0]; 
  vsqnu_corner_surf_nod[1] = 1.366025403784439*vsqnu_surf_nodal[1]-0.3660254037844386*vsqnu_surf_nodal[4]; 
  vsqnu_corner_surf_nod[2] = (-0.5412895520113615*vsqnu_surf_nodal[5])+0.2440169358562925*vsqnu_surf_nodal[4]-0.0687527876293696*vsqnu_surf_nodal[3]+2.020120109712598*vsqnu_surf_nodal[2]-0.9106836025229592*vsqnu_surf_nodal[1]+0.2565888965948001*vsqnu_surf_nodal[0]; 
  vsqnu_corner_surf_nod[3] = 0.2565888965948001*vsqnu_surf_nodal[5]-0.9106836025229592*vsqnu_surf_nodal[4]+2.020120109712598*vsqnu_surf_nodal[3]-0.0687527876293696*vsqnu_surf_nodal[2]+0.2440169358562925*vsqnu_surf_nodal[1]-0.5412895520113615*vsqnu_surf_nodal[0]; 
  vsqnu_corner_surf_nod[4] = 1.366025403784439*vsqnu_surf_nodal[4]-0.3660254037844386*vsqnu_surf_nodal[1]; 
  vsqnu_corner_surf_nod[5] = 2.020120109712598*vsqnu_surf_nodal[5]-0.9106836025229592*vsqnu_surf_nodal[4]+0.2565888965948001*vsqnu_surf_nodal[3]-0.5412895520113615*vsqnu_surf_nodal[2]+0.2440169358562925*vsqnu_surf_nodal[1]-0.0687527876293696*vsqnu_surf_nodal[0]; 

  // If vsqnu < 0 at corner nodes, revaluate it at this node.
  if (vsqnu_corner_vol_nod[0] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[0] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[1] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[1] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[2] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[2] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[3] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[3] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[4] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[4] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[5] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[5] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[6] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[6] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[7] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[7] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[8] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[8] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[9] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[9] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[10] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_vol_nod[10] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }
  if (vsqnu_corner_vol_nod[11] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2]))/mass);
    vsqnu_corner_vol_nod[11] = fabs(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])<0.1*dv_min[1]? 0.0 : 2.0*(1.224744871391589*vmap[3]+0.7071067811865475*vmap[2])*((a/c_const*(alpha+beta)*pow(vmag,gamma))/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta)));
  }

  if (vsqnu_corner_surf_nod[0] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[0] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }
  if (vsqnu_corner_surf_nod[1] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[1] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }
  if (vsqnu_corner_surf_nod[2] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((0.7071067811865475*bmag[0]-1.224744871391589*bmag[1])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[2] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }
  if (vsqnu_corner_surf_nod[3] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]-1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[3] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }
  if (vsqnu_corner_surf_nod[4] < 0.0) {
    vmag = sqrt((0.7071067811865475*vmapSq[0]-0.7905694150420947*vmapSq[2])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[4] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }
  if (vsqnu_corner_surf_nod[5] < 0.0) {
    vmag = sqrt((1.58113883008419*vmapSq[2]+1.224744871391589*vmapSq[1]+0.7071067811865475*vmapSq[0])+2.0*((1.224744871391589*bmag[1]+0.7071067811865475*bmag[0])*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3]))/mass);
    vsqnu_corner_surf_nod[5] = fabs(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])<0.1*dv_min[1]? 0.0 : 2.0*(0.7071067811865475*vmap[2]-1.224744871391589*vmap[3])*((a/c_const*(alpha+beta)*pow(vmag_surf,gamma))/(beta*pow(vmag_surf/scaled_v0,-alpha)+alpha*pow(vmag_surf/scaled_v0,beta)));
  }

  drag_rad_nu[0] = 0.1178511301977579*vsqnu_corner_vol_nod[11]+0.1178511301977579*vsqnu_corner_vol_nod[10]+0.4714045207910317*vsqnu_corner_vol_nod[9]+0.4714045207910317*vsqnu_corner_vol_nod[8]+0.1178511301977579*vsqnu_corner_vol_nod[7]+0.1178511301977579*vsqnu_corner_vol_nod[6]+0.1178511301977579*vsqnu_corner_vol_nod[5]+0.1178511301977579*vsqnu_corner_vol_nod[4]+0.4714045207910317*vsqnu_corner_vol_nod[3]+0.4714045207910317*vsqnu_corner_vol_nod[2]+0.1178511301977579*vsqnu_corner_vol_nod[1]+0.1178511301977579*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[1] = 0.06804138174397717*vsqnu_corner_vol_nod[11]+0.06804138174397717*vsqnu_corner_vol_nod[10]+0.2721655269759087*vsqnu_corner_vol_nod[9]+0.2721655269759087*vsqnu_corner_vol_nod[8]+0.06804138174397717*vsqnu_corner_vol_nod[7]+0.06804138174397717*vsqnu_corner_vol_nod[6]-0.06804138174397717*vsqnu_corner_vol_nod[5]-0.06804138174397717*vsqnu_corner_vol_nod[4]-0.2721655269759087*vsqnu_corner_vol_nod[3]-0.2721655269759087*vsqnu_corner_vol_nod[2]-0.06804138174397717*vsqnu_corner_vol_nod[1]-0.06804138174397717*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[2] = 0.2041241452319315*vsqnu_corner_vol_nod[11]+0.2041241452319315*vsqnu_corner_vol_nod[10]-0.2041241452319315*vsqnu_corner_vol_nod[7]-0.2041241452319315*vsqnu_corner_vol_nod[6]+0.2041241452319315*vsqnu_corner_vol_nod[5]+0.2041241452319315*vsqnu_corner_vol_nod[4]-0.2041241452319315*vsqnu_corner_vol_nod[1]-0.2041241452319315*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[3] = 0.06804138174397717*vsqnu_corner_vol_nod[11]-0.06804138174397717*vsqnu_corner_vol_nod[10]+0.2721655269759087*vsqnu_corner_vol_nod[9]-0.2721655269759087*vsqnu_corner_vol_nod[8]+0.06804138174397717*vsqnu_corner_vol_nod[7]-0.06804138174397717*vsqnu_corner_vol_nod[6]+0.06804138174397717*vsqnu_corner_vol_nod[5]-0.06804138174397717*vsqnu_corner_vol_nod[4]+0.2721655269759087*vsqnu_corner_vol_nod[3]-0.2721655269759087*vsqnu_corner_vol_nod[2]+0.06804138174397717*vsqnu_corner_vol_nod[1]-0.06804138174397717*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[4] = 0.1178511301977579*vsqnu_corner_vol_nod[11]+0.1178511301977579*vsqnu_corner_vol_nod[10]-0.1178511301977579*vsqnu_corner_vol_nod[7]-0.1178511301977579*vsqnu_corner_vol_nod[6]-0.1178511301977579*vsqnu_corner_vol_nod[5]-0.1178511301977579*vsqnu_corner_vol_nod[4]+0.1178511301977579*vsqnu_corner_vol_nod[1]+0.1178511301977579*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[5] = 0.0392837100659193*vsqnu_corner_vol_nod[11]-0.0392837100659193*vsqnu_corner_vol_nod[10]+0.1571348402636772*vsqnu_corner_vol_nod[9]-0.1571348402636772*vsqnu_corner_vol_nod[8]+0.0392837100659193*vsqnu_corner_vol_nod[7]-0.0392837100659193*vsqnu_corner_vol_nod[6]-0.0392837100659193*vsqnu_corner_vol_nod[5]+0.0392837100659193*vsqnu_corner_vol_nod[4]-0.1571348402636772*vsqnu_corner_vol_nod[3]+0.1571348402636772*vsqnu_corner_vol_nod[2]-0.0392837100659193*vsqnu_corner_vol_nod[1]+0.0392837100659193*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[6] = 0.1178511301977579*vsqnu_corner_vol_nod[11]-0.1178511301977579*vsqnu_corner_vol_nod[10]-0.1178511301977579*vsqnu_corner_vol_nod[7]+0.1178511301977579*vsqnu_corner_vol_nod[6]+0.1178511301977579*vsqnu_corner_vol_nod[5]-0.1178511301977579*vsqnu_corner_vol_nod[4]-0.1178511301977579*vsqnu_corner_vol_nod[1]+0.1178511301977579*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[7] = 0.06804138174397717*vsqnu_corner_vol_nod[11]-0.06804138174397717*vsqnu_corner_vol_nod[10]-0.06804138174397717*vsqnu_corner_vol_nod[7]+0.06804138174397717*vsqnu_corner_vol_nod[6]-0.06804138174397717*vsqnu_corner_vol_nod[5]+0.06804138174397717*vsqnu_corner_vol_nod[4]+0.06804138174397717*vsqnu_corner_vol_nod[1]-0.06804138174397717*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[8] = 0.105409255338946*vsqnu_corner_vol_nod[11]+0.105409255338946*vsqnu_corner_vol_nod[10]-0.210818510677892*vsqnu_corner_vol_nod[9]-0.210818510677892*vsqnu_corner_vol_nod[8]+0.105409255338946*vsqnu_corner_vol_nod[7]+0.105409255338946*vsqnu_corner_vol_nod[6]+0.105409255338946*vsqnu_corner_vol_nod[5]+0.105409255338946*vsqnu_corner_vol_nod[4]-0.210818510677892*vsqnu_corner_vol_nod[3]-0.210818510677892*vsqnu_corner_vol_nod[2]+0.105409255338946*vsqnu_corner_vol_nod[1]+0.105409255338946*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[9] = 0.06085806194501844*vsqnu_corner_vol_nod[11]+0.06085806194501844*vsqnu_corner_vol_nod[10]-0.1217161238900369*vsqnu_corner_vol_nod[9]-0.1217161238900369*vsqnu_corner_vol_nod[8]+0.06085806194501844*vsqnu_corner_vol_nod[7]+0.06085806194501844*vsqnu_corner_vol_nod[6]-0.06085806194501844*vsqnu_corner_vol_nod[5]-0.06085806194501844*vsqnu_corner_vol_nod[4]+0.1217161238900369*vsqnu_corner_vol_nod[3]+0.1217161238900369*vsqnu_corner_vol_nod[2]-0.06085806194501844*vsqnu_corner_vol_nod[1]-0.06085806194501844*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[10] = 0.06085806194501844*vsqnu_corner_vol_nod[11]-0.06085806194501844*vsqnu_corner_vol_nod[10]-0.1217161238900369*vsqnu_corner_vol_nod[9]+0.1217161238900369*vsqnu_corner_vol_nod[8]+0.06085806194501844*vsqnu_corner_vol_nod[7]-0.06085806194501844*vsqnu_corner_vol_nod[6]+0.06085806194501844*vsqnu_corner_vol_nod[5]-0.06085806194501844*vsqnu_corner_vol_nod[4]-0.1217161238900369*vsqnu_corner_vol_nod[3]+0.1217161238900369*vsqnu_corner_vol_nod[2]+0.06085806194501844*vsqnu_corner_vol_nod[1]-0.06085806194501844*vsqnu_corner_vol_nod[0]; 
  drag_rad_nu[11] = 0.03513641844631532*vsqnu_corner_vol_nod[11]-0.03513641844631532*vsqnu_corner_vol_nod[10]-0.07027283689263064*vsqnu_corner_vol_nod[9]+0.07027283689263064*vsqnu_corner_vol_nod[8]+0.03513641844631532*vsqnu_corner_vol_nod[7]-0.03513641844631532*vsqnu_corner_vol_nod[6]-0.03513641844631532*vsqnu_corner_vol_nod[5]+0.03513641844631532*vsqnu_corner_vol_nod[4]+0.07027283689263064*vsqnu_corner_vol_nod[3]-0.07027283689263064*vsqnu_corner_vol_nod[2]-0.03513641844631532*vsqnu_corner_vol_nod[1]+0.03513641844631532*vsqnu_corner_vol_nod[0]; 

  drag_rad_nu_surf[0] = 0.1666666666666667*vsqnu_corner_surf_nod[5]+0.6666666666666666*vsqnu_corner_surf_nod[4]+0.1666666666666667*vsqnu_corner_surf_nod[3]+0.1666666666666667*vsqnu_corner_surf_nod[2]+0.6666666666666666*vsqnu_corner_surf_nod[1]+0.1666666666666667*vsqnu_corner_surf_nod[0]; 
  drag_rad_nu_surf[1] = 0.09622504486493764*vsqnu_corner_surf_nod[5]+0.3849001794597506*vsqnu_corner_surf_nod[4]+0.09622504486493764*vsqnu_corner_surf_nod[3]-0.09622504486493764*vsqnu_corner_surf_nod[2]-0.3849001794597506*vsqnu_corner_surf_nod[1]-0.09622504486493764*vsqnu_corner_surf_nod[0]; 
  drag_rad_nu_surf[2] = 0.2886751345948129*vsqnu_corner_surf_nod[5]-0.2886751345948129*vsqnu_corner_surf_nod[3]+0.2886751345948129*vsqnu_corner_surf_nod[2]-0.2886751345948129*vsqnu_corner_surf_nod[0]; 
  drag_rad_nu_surf[3] = 0.1666666666666667*vsqnu_corner_surf_nod[5]-0.1666666666666667*vsqnu_corner_surf_nod[3]-0.1666666666666667*vsqnu_corner_surf_nod[2]+0.1666666666666667*vsqnu_corner_surf_nod[0]; 
  drag_rad_nu_surf[4] = 0.149071198499986*vsqnu_corner_surf_nod[5]-0.2981423969999719*vsqnu_corner_surf_nod[4]+0.149071198499986*vsqnu_corner_surf_nod[3]+0.149071198499986*vsqnu_corner_surf_nod[2]-0.2981423969999719*vsqnu_corner_surf_nod[1]+0.149071198499986*vsqnu_corner_surf_nod[0]; 
  drag_rad_nu_surf[5] = 0.08606629658238703*vsqnu_corner_surf_nod[5]-0.1721325931647741*vsqnu_corner_surf_nod[4]+0.08606629658238703*vsqnu_corner_surf_nod[3]-0.08606629658238703*vsqnu_corner_surf_nod[2]+0.1721325931647741*vsqnu_corner_surf_nod[1]-0.08606629658238703*vsqnu_corner_surf_nod[0]; 

} 
