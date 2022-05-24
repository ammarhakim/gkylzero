// Gkyl ------------------------------------------------------------------------
//
// Header file for FEM Poisson solver.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 
 
int fem_poisson_num_nodes_global_1x_ser_p1_periodicx(const int *numCells);
int fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(const int *numCells);

void fem_poisson_local_to_global_1x_ser_p1_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_Uxperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p1_Uxnonperiodicx(const int *numCells, const int *idx, long *globalIdxs);


int fem_poisson_num_nodes_global_2x_ser_p1_periodicx_periodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells);

void fem_poisson_local_to_global_2x_ser_p1_periodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxnonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxnonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_periodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_periodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_nonperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxnonperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p1_Uxnonperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);


int fem_poisson_num_nodes_global_3x_ser_p1_periodicx_periodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_periodicx_periodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells);

void fem_poisson_local_to_global_3x_ser_p1_periodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_periodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_nonperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p1_Uxnonperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);



int fem_poisson_num_nodes_global_1x_ser_p2_periodicx(const int *numCells);
int fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(const int *numCells);

void fem_poisson_local_to_global_1x_ser_p2_periodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_Uxperiodicx(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_1x_ser_p2_Uxnonperiodicx(const int *numCells, const int *idx, long *globalIdxs);


int fem_poisson_num_nodes_global_2x_ser_p2_periodicx_periodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p2_periodicx_nonperiodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_periodicy(const int *numCells);
int fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_nonperiodicy(const int *numCells);

void fem_poisson_local_to_global_2x_ser_p2_periodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_periodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_nonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_nonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxnonperiodicx_periodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxnonperiodicx_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_periodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_periodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_nonperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_nonperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxnonperiodicx_Uyperiodicy(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_2x_ser_p2_Uxnonperiodicx_Uynonperiodicy(const int *numCells, const int *idx, long *globalIdxs);


int fem_poisson_num_nodes_global_3x_ser_p2_periodicx_periodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_periodicx_periodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_periodicx_nonperiodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_periodicx_nonperiodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_periodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_periodicy_nonperiodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_nonperiodicy_periodicz(const int *numCells);
int fem_poisson_num_nodes_global_3x_ser_p2_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells);

void fem_poisson_local_to_global_3x_ser_p2_periodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_periodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_periodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_nonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_nonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uyperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uyperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uynonperiodicy_periodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uynonperiodicy_nonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_periodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_periodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_nonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_nonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_periodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_nonperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uyperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uyperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uynonperiodicy_Uzperiodicz(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global_3x_ser_p2_Uxnonperiodicx_Uynonperiodicy_Uznonperiodicz(const int *numCells, const int *idx, long *globalIdxs);



EXTERN_C_END 
