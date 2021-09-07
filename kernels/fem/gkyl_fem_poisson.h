// Gkyl ------------------------------------------------------------------------
//
// Header file for FEM Poisson solver.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once
 
void fem_poisson_num_nodes_global1xser_periodicxP1(const int *numCells);
void fem_poisson_num_nodes_global1xser_nonPeriodicxP1(const int *numCells);

void fem_poisson_local_to_global1xser_periodicxP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_nonPeriodicxP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_UxperiodicxP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_UxnonPeriodicxP1(const int *numCells, const int *idx, long *globalIdxs);


void fem_poisson_num_nodes_global2xser_periodicxperiodicyP1(const int *numCells);
void fem_poisson_num_nodes_global2xser_periodicxnonPeriodicyP1(const int *numCells);
void fem_poisson_num_nodes_global2xser_nonPeriodicxperiodicyP1(const int *numCells);
void fem_poisson_num_nodes_global2xser_nonPeriodicxnonPeriodicyP1(const int *numCells);

void fem_poisson_local_to_global2xser_periodicxperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxnonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxnonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxnonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxnonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxUyperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxUynonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxUyperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxUynonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxUyperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxUynonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxUyperiodicyP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxUynonPeriodicyP1(const int *numCells, const int *idx, long *globalIdxs);


void fem_poisson_num_nodes_global3xser_periodicxperiodicyperiodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxperiodicynonPeriodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxnonPeriodicyperiodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxnonPeriodicynonPeriodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxperiodicyperiodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxperiodicynonPeriodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxnonPeriodicyperiodiczP1(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxnonPeriodicynonPeriodiczP1(const int *numCells);

void fem_poisson_local_to_global3xser_periodicxperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicynonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyUzperiodiczP1(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyUznonPeriodiczP1(const int *numCells, const int *idx, long *globalIdxs);



void fem_poisson_num_nodes_global1xser_periodicxP2(const int *numCells);
void fem_poisson_num_nodes_global1xser_nonPeriodicxP2(const int *numCells);

void fem_poisson_local_to_global1xser_periodicxP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_nonPeriodicxP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_UxperiodicxP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global1xser_UxnonPeriodicxP2(const int *numCells, const int *idx, long *globalIdxs);


void fem_poisson_num_nodes_global2xser_periodicxperiodicyP2(const int *numCells);
void fem_poisson_num_nodes_global2xser_periodicxnonPeriodicyP2(const int *numCells);
void fem_poisson_num_nodes_global2xser_nonPeriodicxperiodicyP2(const int *numCells);
void fem_poisson_num_nodes_global2xser_nonPeriodicxnonPeriodicyP2(const int *numCells);

void fem_poisson_local_to_global2xser_periodicxperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxnonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxnonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxnonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxnonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxUyperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_periodicxUynonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxUyperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_nonPeriodicxUynonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxUyperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxperiodicxUynonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxUyperiodicyP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global2xser_UxnonPeriodicxUynonPeriodicyP2(const int *numCells, const int *idx, long *globalIdxs);


void fem_poisson_num_nodes_global3xser_periodicxperiodicyperiodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxperiodicynonPeriodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxnonPeriodicyperiodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_periodicxnonPeriodicynonPeriodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxperiodicyperiodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxperiodicynonPeriodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxnonPeriodicyperiodiczP2(const int *numCells);
void fem_poisson_num_nodes_global3xser_nonPeriodicxnonPeriodicynonPeriodiczP2(const int *numCells);

void fem_poisson_local_to_global3xser_periodicxperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicynonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxnonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxnonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxnonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxnonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUyperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_periodicxUynonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUyperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_nonPeriodicxUynonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUyperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxperiodicxUynonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUyperiodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyUzperiodiczP2(const int *numCells, const int *idx, long *globalIdxs);
void fem_poisson_local_to_global3xser_UxnonPeriodicxUynonPeriodicyUznonPeriodiczP2(const int *numCells, const int *idx, long *globalIdxs);



